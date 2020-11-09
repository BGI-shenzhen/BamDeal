#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <zlib.h>
#include "soap2bam.h"
#include "bam2fq.h"
#include "bam2fa.h"
#include "bam2soap.h"



using namespace std;

int soap2bam_main(int argc, char **argv);
int Bam2fq_main(int argc, char **argv);
int Bam2fa_main(int argc, char **argv);
int Bam2Soap_main(int argc, char **argv);

static int  covert_usage ()
{
	cerr<<""
		"\n"
		"\t\tsoap2bam        soap    -->  bam/sam Format\n"
		"\t\tbam2soap        bam/sam -->  soap    Format\n"
		"\t\tbam2fq          bam/sam -->  Fastq   Format\n"
		"\t\tbam2fa          bam/sam -->  Fasta   Format\n"
		"\n"
		"\t\tHelp            Show this help\n"
		"\n";
	return 1;
}


int convert_main(int argc, char *argv[])
{
	if (argc < 2) { return covert_usage(); }
	else if (strcmp(argv[1], "soap2bam") == 0) { return soap2bam_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "bam2soap") == 0) { return Bam2Soap_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "bam2fq") == 0) { return Bam2fq_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "bam2fa") == 0) { return Bam2fa_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "Help") == 0 || strcmp(argv[1], "help") == 0 || strcmp(argv[1], "?")== 0 || ( argv[1][0] == '-' &&( argv[1][1] =='h' || argv[1][1] =='H' || argv[1][1] =='?' ) )  || strcmp(argv[1], "less") == 0 )
	{
		return covert_usage();
	}
	else
	{
		cerr<<"convert [main] unrecognized command "<<argv[1]<<endl;
		return 1;
	}
	return 0;
}

///////// swimming in the sky and flying in the sea ////////////


