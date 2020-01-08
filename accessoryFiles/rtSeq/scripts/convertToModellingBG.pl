use strict;

my %csSz;

open FAI, $ARGV[1];
while (<FAI>){
  chomp;
  my ($cs,$sz,$tot) = split(/\t/,$_);
  $csSz{$cs}=$sz;
}
close FAI;

open BG, $ARGV[0];
while (<BG>){
  chomp;
  my ($cs,$from,$to,$score) = split(/\t/,$_);
  for (my $x = $from; $x < $to; $x += 10000){
    print join("\t",$cs,$x,($x+9999),$score)."\n" unless ($x+9999 > $csSz{$cs});
  }
}
close BG;
