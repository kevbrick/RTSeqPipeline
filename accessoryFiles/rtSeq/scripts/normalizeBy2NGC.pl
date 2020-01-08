#!/usr/bin/perl
use strict;

my $infileFactors            = $ARGV[0];
my $infileGCtab              = $ARGV[1];
my $maxSim                   = $ARGV[2];
my $useLowCoverageCorrection = $ARGV[3];
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

open IN    , $infileGCtab;
open OUT   , ">", $out;
open GCOUT , ">", 'GCData.tab';

while (<IN>){
  chomp;
  my ($cs, $from, $to, $real, $sim, $GC) = split(/\t/,$_);

  next if ($_ =~ /NA/ || $real <= 0 );
  next unless ($cf{$GC});

  ## CORRECT FOR GC AND MAPABILITY
  next unless ($sim);
  my ($cover);

  if ($useLowCoverageCorrection){
    $cover = ($real / $cf{$GC}) * ($maxSim / $sim) ;
  }else{
    $cover = ($real / $cf{$GC}) ;
  }

  #next if (($cover < $loQ || $cover > $hiQ) || $cs !~ /^chr[0123456789]+$/);
  next if ($cs !~ /^chr[0123456789XYM]+$/);

  next if (($sim/$maxSim) < 0.60);
  next if (($sim/$maxSim) > 1.10);

  ## USE GENOME VALS IF THERE ARE NO READS SELECTED FROM A VERY SHORT CHROMOSOME (i.e. chrM)

  print OUT   join("\t",$cs,$from,$to,$cover)."\n";
  print GCOUT join("\t",$GC,$real,($real / $cf{$GC}),$cover)."\n";
}

close IN;
close OUT;
close GCOUT;
