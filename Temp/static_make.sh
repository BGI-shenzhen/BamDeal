#!/bin/sh
#$ -S /bin/sh
#Version1.0	P_bc_rd@genomics.org.cn	2020-06-05
echo Start Time : 
date
#/share/app/gcc-5.2.0/bin/g++	--std=c++11	-g	-O3	-I/hwfswh2/BC_PUB/Software/06.Develop/02.SpeUsr/include/	-L/hwfswh2/BC_PUB/Software/06.Develop/02.SpeUsr/lib/	-o	BamDeal_Linux	BamDeal.cpp	-lhts	-lncurses	-lm	-lpthread	-lboost_thread	-lz	-lrt	-static	-llzma	-lbz2	-I/hwfswh2/BC_PUB/Software/06.Develop/02.SpeUsr/include/	-L/hwfswh2/BC_PUB/Software/06.Develop/02.SpeUsr/lib/	
../gcc-6.4.0/bin/g++  --std=c++11     -g      -O3     -I  /home/BCadmin/03.Soft_ALL/Dele/samtools-1.12/include/  -L /home/BCadmin/03.Soft_ALL/Dele/samtools-1.12/lib/   -o      BamDeal_Linux   BamDeal.cpp     -lhts      -lm     -lpthread   -lz  -static  -L/data/storage05/BC_PUB/Software/03.Soft_ALL/curl-7.77.0/lib   -L /data/storage05/BC_PUB/Software/01.Usr/lib/  -L /usr/lib/  -L/data/storage05/BC_PUB/Software/03.Soft_ALL/curl-7.77.0/lib    -lcurl  -lcrypto  -lboost_thread    -lrt  -ldl   -llzma  -lbz2  -lcares
echo End Time : 
date
