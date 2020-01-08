#!/usr/bin/perl
use strict;

my $infileFactors = $ARGV[0];
my $infileMeans   = $ARGV[1];
my $infileGCtab   = $ARGV[2];
my $maxSim        = $ARGV[3];

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
## Get Chrom Means
open CM, $infileMeans;
my (%chrom);

while (<CM>){
  next if ($_ =~ /median/);
  chomp;
  my @F = split(/\t/,$_);
  $chrom{$F[0]}->{mean}   = $F[1];
  $chrom{$F[0]}->{median} = $F[2];
}

close CM;

#######################################################
## Normalize ALL
my $out      = $infileGCtab; $out      =~ s/\.tab/.GCcorrected.bedgraph/;
my $outCSMN  = $infileGCtab; $outCSMN  =~ s/\.tab/.normByChromMean.bedgraph/;
my $outCSMED = $infileGCtab; $outCSMED =~ s/\.tab/.normByChromMedian.bedgraph/;
my $outWGMN  = $infileGCtab; $outWGMN  =~ s/\.tab/.normByGenomeMean.bedgraph/;
my $outWGMED = $infileGCtab; $outWGMED =~ s/\.tab/.normByGenomeMedian.bedgraph/;

open IN      , $infileGCtab;
open OUT     , ">", $out;
open OUTCSMN , ">", $outCSMN;
open OUTCSMED, ">", $outCSMED;
open OUTWGMN , ">", $outWGMN;
open OUTWGMED, ">", $outWGMED;

while (<IN>){
  chomp;
  my ($cs, $from, $to, $sim, $real, $nonRep, $GC) = split(/\t/,$_);

  next if ($_ =~ /NA/ || $real <= 0 );
  next unless ($cf{$GC});

  ## CORRECT FOR GC AND MAPABILITY
  next unless ($sim);
  my $cover = ($real / $cf{$GC}) * ($maxSim / $sim) ;

  next if (($cover < $loQ || $cover > $hiQ ) && $cs !~ /^chr[0123456789]+$/);
  next if ($sim / $maxSim < 0.50);
  ## USE GENOME VALS IF THERE ARE NO READS SELECTED FROM A VERY SHORT CHROMOSOME (i.e. chrM)
  $chrom{$cs}->{mean}   = $chrom{"all"}->{mean} unless ($chrom{$cs}->{mean});
  $chrom{$cs}->{median} = $chrom{"all"}->{median} unless ($chrom{$cs}->{median});

  print OUT      join("\t",$cs,$from,$to,$cover)."\n";
  print OUTCSMN  join("\t",$cs,$from,$to,$cover/$chrom{$cs}->{mean})."\n";
  print OUTCSMED join("\t",$cs,$from,$to,$cover/$chrom{$cs}->{median})."\n";
  print OUTWGMN  join("\t",$cs,$from,$to,$cover/$chrom{"all"}->{median})."\n";
  print OUTWGMED join("\t",$cs,$from,$to,$cover/$chrom{"all"}->{median})."\n";
}

close IN;
close OUT;
close OUTCSMN;
close OUTCSMED;
close OUTWGMN;
close OUTWGMED;
