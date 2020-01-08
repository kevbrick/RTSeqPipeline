use strict;
use List::Util qw (min reduce); 

my $bam = $ARGV[0];

my @arr; 
if ($ARGV[1]){
	@arr = split(/,/,$ARGV[1]);
}else{
	@arr = (50,76,150);
}

open my $PIPE, '-|', "samtools view ".($bam =~ /\.sam$/?" -S ":"")." $bam |cut -f10 |head -n 50000 |tail -n 5000";

my %cnt; 

while (<$PIPE>){
	chomp;
	my $f = length($_);
        my @a2 = map { abs($_ - $f) } @arr;
	my $mina2 = min(@a2); 
	my $sz; 
	for my $n(0..$#arr){$sz = $arr[$n] if ($a2[$n] == $mina2)}; 
	$cnt{$sz}++; 
}

my $highest = List::Util::reduce { $cnt{$b} > $cnt{$a} ? $b : $a } keys %cnt;
 
print $highest."\n";
