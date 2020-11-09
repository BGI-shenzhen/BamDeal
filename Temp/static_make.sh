#!/bin/sh
#$ -S /bin/sh
#Version1.0	P_bc_rd@genomics.org.cn	2020-06-05
echo Start Time : 
date
/share/app/gcc-5.2.0/bin/g++	--std=c++11	-g	-O3	-I/hwfswh2/BC_PUB/Software/06.Develop/02.SpeUsr/include/	-L/hwfswh2/BC_PUB/Software/06.Develop/02.SpeUsr/lib/	-o	BamDeal_Linux	BamDeal.cpp	-lhts	-lncurses	-lm	-lpthread	-lboost_thread	-lz	-lrt	-static	-llzma	-lbz2	-I/hwfswh2/BC_PUB/Software/06.Develop/02.SpeUsr/include/	-L/hwfswh2/BC_PUB/Software/06.Develop/02.SpeUsr/lib/	
echo End Time : 
date
