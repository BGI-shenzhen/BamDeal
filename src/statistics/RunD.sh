#!/bin/sh
#$ -S /bin/sh
#Version1.0	hewm@genomics.org.cn	2018-05-17
echo Start Time : 
date
source	/ifshk7/BC_PS/heweiming/06.Delelop/SourceMe2.sh	
../../BamDeal	statistics	DeteSV	-InFile	../../../../TestBamDeal/02.bam/out.sort.bam	-OutFile	outFix	-Insert	280	-SD	40	
echo End Time : 
date
