use strict;
use Statistics::Descriptive;

my ($bgIN, $name) = @ARGV;

## READ ALL INTERVAL VALUES
open BG, $bgIN;

my (%cov,@bg);
while (<BG>){
  chomp;
  my ($cs,$from,$to,$cover) = split(/\t/,$_);
  if ($cover){
    push @{$cov{'all'}}, $cover if ($cs !~ /chr[XYMZW]/);
    push @{$cov{$cs}}, $cover;
  }

  push @bg, $_;
}

close BG;

## GET STATS
my (%csStat);
for my $cs(keys(%cov)){
  my $stat = Statistics::Descriptive::Full->new();

  $stat->add_data(@{$cov{$cs}});
  $csStat{$cs}->{mean} = $stat->mean;
  $csStat{$cs}->{med}  = $stat->median;
  $csStat{$cs}->{sd}   = $stat->standard_deviation;
  $csStat{$cs}->{hiPC} = $stat->percentile(99.5);
  $csStat{$cs}->{loPC}   = $stat->percentile(0.5);
}

## NORMALIZE BGs
open ALLMN,  '>', $name.".normByGenomeMean.penultimateBG";
open ALLMED, '>', $name.".normByGenomeMedian.penultimateBG";
open CSMN,   '>', $name.".normByChromMean.penultimateBG";
open CSMED,  '>', $name.".normByChromMedian.penultimateBG";

for (@bg){
  chomp;
  my ($cs,$from,$to,$cover) = split(/\t/,$_);

  my ($allMean,$allMed,$csMean,$csMed) = (0,0,0,0);

  ## Check that val is within +- 2s.d.s
  #if ($cover > ($csStat{all}->{mean} - 2*$csStat{all}->{sd}) || ($cover < $csStat{all}->{mean} + 2*$csStat{all}->{sd})){
  if ($cover){
    if (($cover > $csStat{all}->{loPC}) && ($cover < $csStat{all}->{hiPC})) {
      my $allMeanI = ($cover && $csStat{all}->{mean})?log2($cover/$csStat{all}->{mean}):0;
      my $allMedI  = ($cover && $csStat{all}->{med})?log2($cover/$csStat{all}->{med}):0;
      my $csMeanI  = ($cover && $csStat{$cs}->{mean})?log2($cover/$csStat{$cs}->{mean}):0;
      my $csMedI   = ($cover && $csStat{$cs}->{med})?log2($cover/$csStat{$cs}->{med}):0;

      $allMean  = (abs($allMeanI) <= 1.5)?$allMeanI:0;
      $allMed   = (abs($allMedI)  <= 1.5)?$allMedI:0;
      $csMean   = (abs($csMeanI)  <= 1.5)?$csMeanI:0;
      $csMed    = (abs($csMedI)  <= 1.5)?$csMedI:0;
    }else{
      if (($cover < $csStat{all}->{loPC})) {
        my $allMeanI = ($cover && $csStat{all}->{mean})?log2($csStat{'all'}->{loPC}/$csStat{all}->{mean}):0;
        my $allMedI  = ($cover && $csStat{all}->{med})?log2($csStat{'all'}->{loPC}/$csStat{all}->{med}):0;
        my $csMeanI  = ($cover && $csStat{$cs}->{mean})?log2($csStat{all}->{loPC}/$csStat{$cs}->{mean}):0;
        my $csMedI   = ($cover && $csStat{$cs}->{med})?log2($csStat{all}->{loPC}/$csStat{$cs}->{med}):0;

        $allMean  = (abs($allMeanI) <= 1.5)?$allMeanI:0;
        $allMed   = (abs($allMedI)  <= 1.5)?$allMedI:0;
        $csMean   = (abs($csMeanI)  <= 1.5)?$csMeanI:0;
        $csMed    = (abs($csMedI)  <= 1.5)?$csMedI:0;
      }

      if (($cover > $csStat{all}->{hiPC})) {
        my $allMeanI = ($cover && $csStat{all}->{mean})?log2($csStat{'all'}->{hiPC}/$csStat{all}->{mean}):0;
        my $allMedI  = ($cover && $csStat{all}->{med})?log2($csStat{'all'}->{hiPC}/$csStat{all}->{med}):0;
        my $csMeanI  = ($cover && $csStat{$cs}->{mean})?log2($csStat{all}->{hiPC}/$csStat{$cs}->{mean}):0;
        my $csMedI   = ($cover && $csStat{$cs}->{med})?log2($csStat{all}->{hiPC}/$csStat{$cs}->{med}):0;

        $allMean  = (abs($allMeanI) <= 1.5)?$allMeanI:0;
        $allMed   = (abs($allMedI)  <= 1.5)?$allMedI:0;
        $csMean   = (abs($csMeanI)  <= 1.5)?$csMeanI:0;
        $csMed    = (abs($csMedI)  <= 1.5)?$csMedI:0;
      }
    }
  }

  print ALLMN  join("\t",$cs,$from,$to,$allMean)."\n";
  print ALLMED join("\t",$cs,$from,$to,$allMed)."\n";
  print CSMN   join("\t",$cs,$from,$to,$csMean)."\n";
  print CSMED  join("\t",$cs,$from,$to,$csMed)."\n";
}

close ALLMN;
close ALLMED;
close CSMN;
close CSMED;

######################################################
sub log2{
  my $v = shift;
  return(log($v)/log(2));
}
