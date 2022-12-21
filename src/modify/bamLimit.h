
#ifndef BamLimit_H_
#define BamLimit_H_


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include "../ALL/comm.h"
#include "../ALL/DataClass.h"
#include <htslib/sam.h>
#include <cstdlib>
using namespace std;
typedef long long llong ;

/*
   int  print_usage_Limit()
   {
   cout <<""
   "\n"
   "\tUsage: bamLimit -InPut <in.bam> -OutPut <outFix>\n"
   "\n"
   "\t\t-InPut        <str>   InPut sam/bam File\n"
   "\t\t-OutDir       <str>   OutPut bam OutFix\n"
   "\n"
   "\t\t-MaxNum       <int>   Max Read Number for each bam[1000000000]\n"
   "\n"
   "\t\t-help                 show this help\n" 
   "\n";
   return 1;
   }
   */

void print_usage_Limit()
{
	cout<<""
		"\n"
		"\tUsage: bamLimit -i  <in.bam> -o ./\n"
		"\n"
		"\t\t-i    <str>    input SAM/BAM file\n"
		"\t\t-o    <str>    output directory, default [PWD]\n"
		"\n"
		"\t\t-n    <int>    max read number for each bam[1000000000]\n"
		"\n"
		"\t\t-h             show more details for help\n"
		"\n"
		"\n";
}

void More_print_usage_Limit()
{
	cout<<""
		"\n"
		"\t\t1. bamLimit -i  <in.bam> -o ./ -n N\n"
		"\t\tThis will split the input BAM file and output the splitting files in current directory. Each output has at most N reads, default value of N is 1000000000.\n"
		"\n"
		"\n"
		"\n";
}

int parse_Acmd_Limit(int argc, char **argv, In3str1v * para_Xam04 )
{
	if (argc <2 ) {print_usage_Limit();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InPut" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_Xam04->InStr1=argv[i];
		}
		else if (flag  ==  "OutDir")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0; }
			i++;
			para_Xam04->InStr2=argv[i];
		}
		else if (flag  ==  "MaxNum")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0; }
			i++;
			para_Xam04->InInt=atoi(argv[i]);
		}
		else if (flag  == "help" ||  flag  == "h")
		{
			More_print_usage_Limit();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((para_Xam04->InStr1).empty() ||   (para_Xam04->InStr2).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;    
	}
	return 1 ;
}



//programme entry
///////// swimming in the sky and flying in the sea ////////////
//int main(int argc, char **argv)
int bamLimit_main(int argc, char **argv)
{
	In3str1v * para_Xam04 = new In3str1v;
	para_Xam04->InInt=1000000000;
	if (parse_Acmd_Limit(argc, argv , para_Xam04 )==0)
	{
		delete para_Xam04 ; 
		return 0 ;
	}

	bam_hdr_t *header;
	bam1_t *aln = bam_init1();
	samFile *in = sam_open((para_Xam04->InStr1).c_str(), "r");


	header = sam_hdr_read(in);
	string AAOUT=(para_Xam04->InStr2)+".1.bam";
	samFile *outR1 = sam_open(AAOUT.c_str(), "wb");
	if (sam_hdr_write(outR1, header) < 0) 
	{
		fprintf(stderr, "Error writing output.\n");
		exit(-1);
	}
	int X=0;
	int Flag=0;
	int Cout=1;

	while (sam_read1(in, header, aln) >= 0)
	{
		if  (Flag>=(para_Xam04->InInt))
		{
			Flag=0;
			sam_close(outR1);
			Cout++;
			string AAOUT=(para_Xam04->InStr2)+"."+Int2Str(Cout)+".bam";
			outR1 = sam_open(AAOUT.c_str(), "wb");
			if (sam_hdr_write(outR1, header) < 0)
			{
				fprintf(stderr, "Error writing output.\n");
				exit(-1);
			}
			X=sam_write1(outR1, header, aln);
			Flag++;
		}
		else
		{
			X=sam_write1(outR1, header, aln);
			Flag++;
		}
	}

	sam_close(in);
	sam_close(outR1);


	delete para_Xam04 ;
	return 0;
}
#endif 
///////// swimming in the sky and flying in the sea ////////////
