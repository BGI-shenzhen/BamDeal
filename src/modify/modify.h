#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <zlib.h>
#include "bamRand.h"
#include "bamCat.h"
#include "bamFilter.h"
#include "bamShiftQ.h"
#include "bamSplit.h"
#include "bamSubChr.h"
#include "bamAssign.h"
#include "bamLimit.h"



using namespace std;

int bamRand_main(int argc, char **argv);
int bam_Filter_main(int argc, char **argv);
int bam_ShiftQ_main(int argc, char **argv);
int bamSplit_main(int argc, char *argv[]) ;
int bamCat_main(int argc, char *argv[]);
int bam_SubChr_main(int argc, char **argv);
int bamAssign_main(int argc, char *argv[]);
int bamLimit_main(int argc, char *argv[]);

static int  modify_usage ()
{
	cerr<<""
		"\n"
		"\t\tbamFilter       filter low quality read in bam\n"
		"\t\tbamSplit        split single/muti-Bam by chr\n"
		"\t\tbamAssign       split single/muti-Bam by assign chr\n"
		"\t\tbamCat          Merge/Cat muti (diff header) bam to one bam\n"
		"\t\tbamRand         random out partly of bam read\n"
		"\t\tbamSubChr       extract or remove chr(s) from SAM/BAM\n"
		"\t\tbamShiftQ       modify seq Phred quality in bam\n"
		"\t\tbamLimit        Limit big bam to muti subbam by fix line\n"
		"\n"
		"\t\tHelp            Show this help\n"
		"\n";
	return 1;
}

int modify_main(int argc, char *argv[])
{
	if (argc < 2) { return modify_usage(); }
	else if (strcmp(argv[1], "bamRand") == 0) { return bamRand_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "bamFilter") == 0) { return bam_Filter_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "bamSplit") == 0) { return bamSplit_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "bamAssign") == 0) { return bamAssign_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "bamSubChr") == 0) { return bam_SubChr_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "bamCat") == 0) { return bamCat_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "bamShiftQ") == 0) { return bam_ShiftQ_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "bamLimit") == 0) { return bamLimit_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "Help") == 0 || strcmp(argv[1], "help") == 0 || strcmp(argv[1], "?")== 0 || ( argv[1][0] == '-' &&( argv[1][1] =='h' || argv[1][1] =='H' || argv[1][1] =='?' ) )  || strcmp(argv[1], "less") == 0 )
	{
		return modify_usage();
	}
	else
	{
		cerr<<"convert [main] unrecognized command "<<argv[1]<<endl;
		return 1;
	}
	return 0;
}

///////// swimming in the sky and flying in the sea ////////////


