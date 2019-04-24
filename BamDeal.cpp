#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include "src/convert/covert.h"
#include "src/modify/modify.h"
#include "src/statistics/statistics.h"
#include "src/visualize/visualize.h"

using namespace std;

int convert_main(int argc, char *argv[]) ;
int modify_main(int argc, char *argv[]);
int statistics_main(int argc, char *argv[]);
int visualize_main(int argc, char *argv[]);

static int  AllTools_usage ()
{
	cerr <<"Program: BamDeal\nVersion: 0.22\thewm2008@gmail.com\t"<<__DATE__<<endl;

	cerr<<""
		"\n"
		"\tUsage:\n\n"
		"\t\tconvert        convert tools\n"
		"\t\tmodify         modify tools\n"
		"\t\tstatistics     statistics analysis tools\n"
		"\t\tvisualize      visualize tools for bam\n"
		"\n"
		"\t\tHelp           Show help in detail\n"
		"\n";
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) { return AllTools_usage(); }
	else if (strcmp(argv[1], "convert") == 0) { return convert_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "statistics") == 0) { return statistics_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "visualize") == 0) { return visualize_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "modify") == 0) { return modify_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "Help") == 0 || strcmp(argv[1], "help") == 0 || strcmp(argv[1], "?")== 0 || ( argv[1][0] == '-' &&( argv[1][1] =='h' || argv[1][1] =='H' || argv[1][1] =='?' ) )  || strcmp(argv[1], "less") == 0 )
	{
		return  AllTools_usage();
	}
	else
	{
		cerr<<"BamDeal [main] unrecognized command "<<argv[1]<<endl;
		return 1;
	}

	return 0;
}


