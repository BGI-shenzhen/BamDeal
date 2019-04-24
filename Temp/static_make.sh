#!/bin/sh
#$ -S /bin/sh
#Version1.0	hewm@genomics.org.cn	2018-05-4
echo Start Time : 
date
g++	--std=c++11	-g	-O3	-o	BamDeal_Linux	BamDeal.cpp	-lhts	-lncurses	-lm	-lpthread	-lboost_thread	-lboost_system	-lz	-I	/ifshk4/BC_PUB/biosoft/newblc/01.Usr/include/	-lrt	-L	/ifshk4/BC_PUB/biosoft/newblc/01.Usr/lib	-L./src/lib	-static	-llzma	-lbz2	
echo End Time : 
date
