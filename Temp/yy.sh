#!/bin/sh
#$ -S /bin/sh
#Version1.0	  hewm@genomics.org.cn	  2022-12-14
echo Start Time : 
date
/share/app/gcc-5.2.0/bin/g++	  --std=c++11	  -g	  -O3	  BamDeal.cpp	  -lcurl	  -lncurses	  -lhts	  -lm	  -lpthread	  -lboost_thread	  -lz	  -lrt	  -static	  -llzma	  -lbz2	  -I	/home/heweiming/PS/dede/include/  	  -L	 /home/heweiming/PS/dede/lib/ 	  -lcurl	  -lcrypto	  -ldl	  -L	  /home/heweiming/PS/dede/Lib/	  -o	  BamDeal_Linux	  
echo End Time : 
date
