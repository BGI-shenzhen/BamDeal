#ifndef bamStatistics_H_
#define bamStatistics_H_


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <zlib.h>
#include "BamCoverage.h"
#include "BamLowDepthRegion.h"
#include "BamDeteCNV.h"
#include "BamDeteSV.h"
#include "BamBasesCount.h"
#include "BamSiteRead.h"



using namespace std;

int bamSV_main(int argc, char *argv[]);
int bamCov_main(int argc, char *argv[]) ;
int bamLowDepth_main(int argc, char *argv[]);
int bamCNV_main(int argc, char *argv[]);
int BaseCount_main(int argc, char *argv[]);

static int  statistics_usage ()
{
	cerr<<""
		"\n"
		"\t\tCoverage            Calculate Genome Coverage/Depth/GC Dis based Bam\n"
		"\t\tBasesCount          Calculate Genome every Site's four base Depth\n"
		"\t\tDeteCNV             Detect CNV/Deletion Region by merge Depth info based Bam\n"
		"\t\tDeteSV              Detect SV by Pair End Read insert size in Bam\n"
		"\t\tLowDepth            GiveOut bed file of low Depth Region(may BigDeletion)\n"
		"\t\tSiteRead            Extract the read (with base) covering the SNP site\n"
		"\n"
		"\t\tHelp                Show this help\n"
		"\n";
	return 1;
}

int statistics_main(int argc, char *argv[])
{
	if (argc < 2) { return statistics_usage(); }
	else if (strcmp(argv[1], "Coverage") == 0) { return bamCov_main(argc-1, argv+1);}
	else if (strcmp(argv[1], "LowDepth") == 0) { return bamLowDepth_main(argc-1, argv+1);}
	else if (strcmp(argv[1], "DeteSV") == 0) { return bamSV_main(argc-1, argv+1);}
	else if (strcmp(argv[1], "BasesCount") == 0) { return BaseCount_main(argc-1, argv+1);}
	else if (strcmp(argv[1], "DeteCNV") == 0) { return bamCNV_main(argc-1, argv+1);}
	else if (strcmp(argv[1], "SiteRead") == 0) { return SiteRead_main(argc-1, argv+1);}
	else if (strcmp(argv[1], "Help") == 0 || strcmp(argv[1], "help") == 0 || strcmp(argv[1], "?")== 0 || ( argv[1][0] == '-' &&( argv[1][1] =='h' || argv[1][1] =='H' || argv[1][1] =='?' ) )  || strcmp(argv[1], "less") == 0 )
	{
		return statistics_usage();
	}
	else
	{
		cerr<<"convert [main] unrecognized command "<<argv[1]<<endl;
		return 1;
	}
	return 0;
}

#endif 
///////// swimming in the sky and flying in the sea ////////////


