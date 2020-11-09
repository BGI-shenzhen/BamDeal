#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include "../ALL/comm.h"
#include <htslib/sam.h>
#include <cstdlib>
#include "../ALL/gzstream/gzstream.h"
using namespace std;
typedef long long llong ;


class Para_FF {
	public:
		string input_Sam ;
		string OutPut ;
		llong start ;
		llong end ;
		string chr;
		int minLeng ;
		int minMapQ ;
		int maxHit ;
		bool TF ;

		Para_FF()
		{
			start=0;
			end=1000000000;
			minLeng=30 ;
			minMapQ =15 ;
			maxHit=1 ;
			TF=false;
			chr="";
		}
};


void print_Ausage_FF()
{
	cout<<""
		"\n"
		"\tUsage: BamFilter -i <in.bam> -o <out.bam>\n"
		"\n"
		"\t\t-i    <str>    input  SAM/BAM file\n"
		"\t\t-o    <str>    output BAM file\n"
		"\n"
		"\t\t-q    <int>    the quality to filter reads, default [15]\n"
		"\t\t-l    <int>    the length to filter reads, default [30]\n"
		"\t\t-s    <int>    the beginning of interval containing the 1-based leftmost mapping position of first matching base, default [0]\n"
		"\t\t-e    <int>    the end of interval containing the 1-based leftmost mapping position of first matching base, default [1e9]\n"
		"\t\t-c    <str>    specify the chromosome to output, default [all chromosomes]\n"
		"\t\t-d             remove the duplicate read\n"
		"\n"
		"\t\t-h             show more details for help\n"
		"\n";
}

void More_Help_FF()
{
	cout<<""
		"\n"
		"\t\t1. BamFilter -i <in.bam> -o AAA\n"
		"\t\t   This will remove the aligned reads whose quality lower than 15 or length shorter than 30bp and output the reads left to the file named AAA in current directory.\n"
		"\n"
		"\t\t2. BamFilter -i <in.bam> -o AAA -q Q -l L -s S -e E -c ChrX\n"
		"\t\t   This will remove the aligned reads whose quality lower than Q or length shorter than L bp or the 1-based leftmost mapping position of first matching base does not locate within [S,E], and output the reads left of ChrX to the file named AAA in current directory.\n"
		"\n";
}


int parse_Acmd_FF(int argc, char **argv, Para_FF * para_FF)
{
	if (argc <2 ) {print_Ausage_FF();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InPut" ||  flag  == "i")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_FF->input_Sam=argv[i];
		}
		else if (flag  ==  "OutPut" ||  flag  == "o" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0; }
			i++;
			para_FF->OutPut=argv[i];
		}
		else if (flag  ==  "Chr" ||  flag  == "c" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0; }
			i++;
			para_FF->chr=argv[i];
		}
		else if (flag  ==  "MinMapQ" ||  flag  == "q" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_FF->minMapQ=atoi(argv[i]);
		}
		else if (flag  ==  "MinLeng" ||  flag  == "l")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_FF->minLeng=atoi(argv[i]);
		}
		else if (flag  ==  "End" ||  flag  == "e")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_FF->end=atol(argv[i]);
		}
		else if (flag  == "d")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_FF->TF=true;
		}
		else if (flag  ==  "Start" ||  flag  == "s")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0; }
			i++;
			para_FF->start=atol(argv[i]);
		}
		else if (flag  == "help" ||  flag  == "h" )
		{
			More_Help_FF();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((para_FF->input_Sam).empty() ||   (para_FF->OutPut).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;    
	}

	return 1 ;
}



//programme entry
///////// swimming in the sky and flying in the sea ////////////
//int main(int argc, char **argv)
int bam_Filter_main(int argc, char **argv)
{
	Para_FF * para_FF = new Para_FF;
	if (parse_Acmd_FF(argc, argv , para_FF )==0)
	{
		delete para_FF ; 
		return 0 ;
	}

	bam_hdr_t *header;
	bam1_t *aln = bam_init1();
	samFile *in = sam_open((para_FF->input_Sam).c_str(), "r");

	samFile *outR1 = sam_open(para_FF->OutPut.c_str(), "wb");

	header = sam_hdr_read(in);

	if (sam_hdr_write(outR1, header) < 0) 
	{
		fprintf(stderr, "Error writing output.\n");
		exit(-1);
	}
	int tmp=0;
	if (para_FF->TF)
	{
		while (sam_read1(in, header, aln) >= 0)
		{
			if ( (aln->core).qual < (para_FF->minMapQ) ) 
			{
				continue ;
			}
			if ( (aln->core).pos < ((para_FF->start))  ||    (aln->core).pos > (para_FF->end) ) 
			{
				continue ;
			}
			if  ((aln->core).l_qseq   < ( para_FF->minLeng ))
			{
				continue ;
			}
			if ( (aln->core).flag   & 0x400 )  { continue ; }

			string chrID="*";
			if ((aln->core).tid >= 0)
			{ // chr
				chrID=header->target_name[(aln->core).tid];
			}

			if (!(para_FF->chr).empty()   &&   chrID!=(para_FF->chr))
			{
				continue ;
			}

			tmp=sam_write1(outR1, header, aln);
		}
	}
	else
	{
		while (sam_read1(in, header, aln) >= 0)
		{
			if ( (aln->core).qual < (para_FF->minMapQ) ) 
			{
				continue ;
			}
			if ( (aln->core).pos < ((para_FF->start))  ||    (aln->core).pos > (para_FF->end) ) 
			{
				continue ;
			}
			if  ((aln->core).l_qseq   < ( para_FF->minLeng ))
			{
				continue ;
			}
			string chrID="*";
			if ((aln->core).tid >= 0)
			{ // chr
				chrID=header->target_name[(aln->core).tid];
			}

			if (!(para_FF->chr).empty()   &&   chrID!=(para_FF->chr))
			{
				continue ;
			}

			tmp=sam_write1(outR1, header, aln);
		}

	}

	sam_close(in);
	sam_close(outR1);


	delete para_FF ;
	return 0;
}

///////// swimming in the sky and flying in the sea ////////////
