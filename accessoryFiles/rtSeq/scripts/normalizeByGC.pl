#!/usr/bin/perl
use strict;

my $infileFactors = $ARGV[0];
my $infileGCtab   = $ARGV[1];
my $maxSim        = $ARGV[2];

#######################################################
## Get CFs
open CF, $infileFactors;
my ($hiQ, $loQ, %cf);

while (<CF>){
  chomp;
  my ($GC,$median);
  ($GC,$median,$loQ,$hiQ) = split(/\t/,$_);
  next if ($_ =~ /hiQ/);
  $cf{$GC} = $median;
}

close CF;

#######################################################
## Normalize ALL
my $out      = $infileGCtab; $out      =~ s/\.tab/.GCcorrected.bedgraph/;

open IN      , $infileGCtab;
open OUT     , ">", $out;

while (<IN>){
  chomp;
  my ($cs, $from, $to, $sim, $real, $nonRep, $GC) = split(/\t/,$_);

  next if ($_ =~ /NA/ || $real <= 0 );
  next unless ($cf{$GC});

  ## CORRECT FOR GC AND MAPABILITY
  next unless ($sim);
  my $cover = ($real / $cf{$GC}) * ($maxSim / $sim) ;

  next if (($cover < $loQ || $cover > $hiQ) || $cs !~ /^chr[0123456789]+$/);
  next if ($sim / $maxSim < 0.75);

  ## USE GENOME VALS IF THERE ARE NO READS SELECTED FROM A VERY SHORT CHROMOSOME (i.e. chrM)

  print OUT      join("\t",$cs,$from,$to,$cover)."\n";
}

close IN;
close OUT;
