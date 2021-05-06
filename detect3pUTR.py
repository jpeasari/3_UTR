#!/usr/bin/env python

# detect3pUTR.py
# Script to detect 3' UTR extensions from gene models and aligned RNA-seq data

import sys
import subprocess

# Expected parameters are:
#   1. GFF wih gene annotations
#   2. BAM with aligned reads to genome
#   3. OUT file, a tsv file with the detected 3' UTRs, GFF format

gff  = sys.argv[1]
bam  = sys.argv[2]
out  = sys.argv[3]

var  = 0.1 # percentage of variation, used to define when a variation drop/raise
size = 1500 # default extension to detect UTR, for yeast, max size was 1400
samtools_bin = 'samtools' # path to samtools if not in PATH

def defineDrop(region, strand):
    # Function to define where the signal is going to 0 or has an increment > var
    # It uses samtools depth to define the signal (depth coverage over region)
    # Return: last position before the drop

    global bam
    # Samtools command to be run with subprocess
    cmd  = [ samtools_bin, "depth", "-a", "-r", region, bam ]
    # cmd_out is the resulting table, columns are:
    #   1. chr
    #   2. position
    #   3. depth
    cmd_out = subprocess.Popen(
                    cmd, 
                    stdout = subprocess.PIPE, 
                    stderr = subprocess.STDOUT, 
                    text=True
                )
    
    # parsing the output, it considers the strand direction
    # where "+" will sort from low to high, "-" is reversed
    deps = []
    for line in cmd_out.stdout.readlines():
        if strand == "-":
            deps.insert(0, line.rstrip())
        else:
            deps.append(line.rstrip())

    # process all positions, lim defines variation limit,
    # curr is the current position
    lim = None
    curr = 0
    for x in range(0, len(deps)):
        line = deps[x]
        ln = line.split("\t")
        pos = ln[1]
        dep = int(ln[2])
        if dep == 0:  # signal drops to 0
            return pos

        if lim:
            if dep >= lim:  # signal raise beyond variation
                return pos
        else:
            lim = dep * (1 + var)
            curr = pos
    return curr 


## Main ##

# load genes coordinates, stored in "genes"
genes = {}
count = 0
print ("loading gene coordinates from {}".format(gff))
with open(gff, "r") as gh:
    for line in gh:
        if line.startswith("#"):  # skip headers
            continue
        ln = line.split("\t")
        if len(ln) < 8: # skip incompleted lines
            continue
        if ln[2] == "gene": # gene, mRNA, exon and CDS have the same coordinates
            count += 1
            chrom  = ln[0]
            start  = ln[3]
            end    = ln[4]
            strand = ln[6]
            info   = ln[8]
            inf    = info.split(";")
            gene   = inf[0].replace("ID=gene-", "") # extract gene name
            genes[gene] = {}
            genes[gene]["chrom"]  = chrom
            genes[gene]["start"]  = start
            genes[gene]["end"]    = end
            genes[gene]["strand"] = strand

print ("loaded {} genes".format(count))
print ("analyzing BAM {}".format(bam))

# iterates over all genes, write output if 3' UTR is detected
oh = open(out, "w")
for gene in genes.keys():
    chrom  = genes[gene]["chrom"]
    start  = genes[gene]["start"]
    end    = genes[gene]["end"]
    strand = genes[gene]["strand"]

    # as format is GFF3, adding gene information
    utr_info = "ID=3pUTR-{};gene={};locus_tag={}".format(gene, gene, gene)

    # define utr_start and utr_end
    if strand == "-":
        c1 = int(start) - size
        if c1 < 0:
            c1 = 0 
        c2 = start
        region = "{}:{}-{}".format(chrom, str(c1), c2)
        utr_end   = start
        utr_start = defineDrop(region, strand)
    else:
        region = "{}:{}-{}".format(chrom, end, str(int(end) + size))
        utr_start = end
        utr_end   = defineDrop(region, strand)

    if utr_end == utr_start: # for genes without signal
        continue

    # write output
    oh.write("\t".join([chrom, "UTRpred", "3_prime_UTR", utr_start, utr_end, ".", strand, ".", utr_info]))
    oh.write("\n")

oh.close() 
