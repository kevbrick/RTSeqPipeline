#!/bin/bash

chr=$1
g=$2

psr50="$NXF_PIPEDIR/accessoryFiles/rtSeq/$g/pseudoreads/50bpReads_1bpStep/"$chr".bam"
psr150="$NXF_PIPEDIR/accessoryFiles/rtSeq/$g/pseudoreads/150bpReads_1bpStep/"$chr".bam"

fai="$GENOMES/$g/genome.fa.fai"
grep -w $chr $fai >$chr.fai

bedtools makewindows -g $chr.fai -w 101 -s 1 | intersectBed -v -sorted -a - -b $NXF_PIPEDIR/accessoryFiles/rtSeq/$g/blacklist/$g.blacklist.bed | perl -lane 'print join("\t",@F) unless (($F[2]-$F[1]) != 101)' > $chr.win101.bed

bedtools nuc -fi $GENOMES/$g/BWAIndex/version0.7.10/genome.fa -bed $chr.win101.bed -C | grep ^chr | cut -f1-3,5 | perl -M"Math::Round" -pi -e 's/(0\.\d\d+)/round($1*100)/e' |cut -f4 >$chr.win101.GC.tab

intersectBed -a $chr.win101.bed  -b $psr50  -c -sorted -g $fai |perl -pi -e 's/\\./0/g' |cut -f4 >$chr.psr50.tab
intersectBed -a $chr.win101.bed  -b $psr150 -c -sorted -g $fai |perl -pi -e 's/\\./0/g' |cut -f4 >$chr.psr150.tab
