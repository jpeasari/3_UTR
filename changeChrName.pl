#!/usr/bin/perl

use strict;
use warnings;

$ARGV[2] or die "use: perl changeChrName.pl COLNUM INPUT OUTPUT\n";

my $col = shift @ARGV;
my $in  = shift @ARGV;
my $out = shift @ARGV;
my %chr = (
    "NC_001133.9" => "chrI",
    "NC_001134.8" => "chrII",
    "NC_001135.5" => "chrIII",
    "NC_001136.10" => "chrIV",
    "NC_001137.3" => "chrV",
    "NC_001138.5" => "chrVI",
    "NC_001139.9" => "chrVII",
    "NC_001140.6" => "chrVIII",
    "NC_001141.2" => "chrIX",
    "NC_001142.9" => "chrX",
    "NC_001143.9" => "chrXI",
    "NC_001144.5" => "chrXII",
    "NC_001145.3" => "chrXIII",
    "NC_001146.8" => "chrXIV",
    "NC_001147.6" => "chrXV",
    "NC_001148.4" => "chrXVI",
    "NC_001224.1" => "chrMT",
    "Reference"   => "Reference"
);

open (my $ih, "<",  $in) or die;
open (my $oh, ">", $out) or die;
while (<$ih>) {
    chomp;
    my @line = split(/\t/, $_);
    if (defined $chr{ $line[$col] }) {
        $line[$col] = $chr{ $line[$col] };
    }
    else {
        die "undefined chr $line[$col]\n";
    }
    print $oh join "\t", @line;
    print $oh "\n";
}
close $ih;
close $oh;