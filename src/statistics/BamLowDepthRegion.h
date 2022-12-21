#ifndef bamLDR_H_
#define bamLDR_H_
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstring>
#include <list>
#include <map>
#include <iomanip>
#include "../ALL/comm.h"
#include "../ALL/DataClass.h"
#include <cstdlib>
#include <zlib.h>
#include <htslib/sam.h>
#include <htslib/kseq.h>
#include <stdio.h>

using namespace std;
typedef unsigned long long ubit64_t;
//KSEQ_INIT(gzFile, gzread)

/*
int  bamLDR_help()
{
	cout <<""
		"\n"
		"\tUsage: LowDepth  -InList     <bam.list>  -OutPut  <out.bed>\n"
		"\tUsage: LowDepth  -InDepthFa  <Depth.fa.gz>  -OutPut  <out.bed>\n"
		"\n"
		"\t\t-InList      <str>     Input Bam/Sam File List\n"
		"\t\t-InFile      <str>     Input Bam/Sam File File[repeat]\n"
		"\t\t-InDepthF    <str>     In DepthFA File,bamCoverage OutFile\n"
		"\t\t-OutPut      <str>     OutPut Bed Region File\n"
		"\n"
		"\t\t-LowDepth    <int>     Regard SiteDepth < X as LowDepth[2]\n"
		"\t\t-MinLength   <int>     Filter too short region [1000]\n"
		"\t\t-MinQ        <int>     Ignore too low mapQ read[10]\n"
		"\n"
		"\t\t-help                  Show this help [hewm2008 v1.01]\n"
		"\n";
	return 1;
}
*/


void bamLDR_help()
{
cout<<""
"\n"
"\tUsage: LowDepth -l  <bam.list>  -o  <out.bed>\n"
"\tUsage: LowDepth -d  <Depth.fa.gz> -o <out.bed>\n"
"\n"
"\t\t-i    <str>    input SAM/BAM files, delimited by space\n"
"\t\t-l    <str>    input list of SAM/BAM files\n"
"\t\t-d    <str>    depth along site in reference FASTA\n"
"\t\t-o    <str>    output bed region file\n"
"\n"
"\t\t-x    <int>    set the minimum value of low depth,default[2]\n"
"\t\t-s    <int>    the length to filter short region, default [1000]\n"
"\t\t-q    <int>    ignore too low mapQ read, default [10]\n"
"\n"
"\t\t-h             show more details for help\n"
"\n"
"\n";
}

void More_bamLDR_help()
{
cout<<""
"\n"
"\t\t1. LowDepth -i  < A.bam B.bam > -o  <out.bed> -q Q -s S\n"
"\t\tLowDepth -l  <bam.list>  -o  <out.bed> -q Q -s S\n"
"\t\tThis will generate one compressed file (bed file of low depth region) named out.bed.gz in current directory.\n"
"\t\t(1.1)      the reads with quality lower than Q will be removed from analysis, default value of Q is 10.\n"
"\t\t(1.2)      the bed region length lower than S will be filtered.\n"
"\t\t2. LowDepth -d  <Depth.fa.gz> -o <out.bed> -q Q\n"
"\t\tThis operation is the same with the example above but with different input format.\n"
"\t\t(2.1) the input file <Depth.fa.gz> shows the depth along the reference FASTA. User could get this file with the function Coverage in the Statistics module. Below isan example of the usage of this command. More details could be found with the -h in this command.\n"
"\t\tbamdeal statistics Coverage -i <in.bam> -o AAA -r <Ref.fa>\n"
"\n"
"\n";
}
int bamLDR_help01(int argc, char **argv , In3str1v * paraFA04 )
{
	if (argc <=1 ) {bamLDR_help();return 0;}
	int file_count=0;
	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InList" ||  flag  == "List"   ||  flag  == "l")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			string A=argv[i];
			file_count+=(ReadList(A ,(paraFA04->List)));
		}
		else if (flag  == "InFile" ||  flag  == "i")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			string A=argv[i];
			(paraFA04->List).push_back(A);
			file_count++;
			bool RunT=true;
			while(RunT)
			{
				if ( (i+1) < argc  && (argv[i+1][0]!='-'))
				{
					i++;
					A=argv[i];
					(paraFA04->List).push_back(A);
					file_count++;
				}
				else
				{
					RunT=false;
				}
			}
		}

		else if (flag  ==  "InDepthF"  ||  flag  == "d" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr3=argv[i];
		}
		else if (flag  ==  "OutPut"  ||  flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr2=argv[i];
		}
		else if (flag  ==  "LowDepth"   ||  flag  == "x" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InInt=atoi(argv[i]);
		}
		else if (flag  ==  "MinQ" ||  flag  == "q")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InInt2=atoi(argv[i]);
		}
		else if (flag  ==  "MinLength" ||  flag  == "s"  )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InF=atof(argv[i]);
		}
		else if (flag == "help" ||  flag  == "h")
		{
			More_bamLDR_help();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ( (paraFA04->InStr2).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	if  ( (file_count  < 1 )  &&  (paraFA04->InStr3).empty())
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;		
	}
	(paraFA04->InStr2)=add_Asuffix(paraFA04->InStr2);
	return 1 ;
}


int bamLowDepth_main(int argc, char *argv[])
//int main(int argc, char *argv[])
{
	In3str1v *paraFA04 = new In3str1v;
	paraFA04->InInt=2;
	paraFA04->InInt2=10;
	paraFA04->InF=1000;
	if ((bamLDR_help01(argc, argv, paraFA04)==0))
	{
		delete paraFA04 ;
		return 1 ;
	}

	int MinLength=int(paraFA04->InF);

	ogzstream  OUT ((paraFA04->InStr2).c_str());

	if((!OUT.good()))
	{
		cerr << "open OUT File error: "<<(paraFA04->InStr2)<<endl;
		delete  paraFA04 ; return  1;
	}


	if  ((paraFA04->InStr3).empty())
	{

	}
	else
	{

		igzstream FAIN ((paraFA04->InStr3).c_str(),ifstream::in);

		bool S_E=false;
		int Start=0;
		int End=0;
		int Count=0;
		string chrName="";
		int Depth=0;
		getline(FAIN,chrName);
		chrName=chrName.substr(1);

		while(!FAIN.eof())
		{
			string  line ;
			getline(FAIN,line);
			if (line.length()<=0)  { continue  ; }
			if (line[0] == '>')
			{
				if (S_E)
				{
					if ((End-Start)>MinLength)
					{
						OUT<<chrName<<"\t"<<Start<<"\t"<<End<<endl;
					}
				}
				chrName=line.substr(1);
				S_E=false;
				Count=0;
			}
			else
			{
				istringstream isone (line,istringstream::in);
				while(isone>>Depth)
				{				 		
					Count++;
					if ( (!S_E)   &&  ( Depth > (paraFA04->InInt)) )
					{
						continue ;
					}
					else if ( (S_E)   &&  ( Depth  > (paraFA04->InInt))  )
					{
						S_E=false ;
						if ((End-Start)>MinLength)
						{
							OUT<<chrName<<"\t"<<Start<<"\t"<<End<<"\n";
						}
					}
					else if   ( (!S_E)   &&  ( Depth  < (paraFA04->InInt)) )
					{
						S_E=true ;
						Start=Count ;  End=Count;
					}
					else 
					{
						End=Count;
					}
				}
			}


		}

		if (S_E)
		{
			if ((End-Start)>MinLength)
			{
				OUT<<chrName<<"\t"<<Start<<"\t"<<End<<endl;
			}
		}

		OUT.close();
		FAIN.close();
		delete  paraFA04 ; return  0;
	}





	string  BamPath=(paraFA04->List)[0];
	bam_hdr_t *header;
	samFile *BamIn = hts_open(BamPath.c_str(), "r");
	header = sam_hdr_read(BamIn);
	sam_close(BamIn);


	cout<<"begin new the memory ...\n";

	unsigned short int **depth = new unsigned short int*[(header->n_targets)]; //开辟行  
	for(int i = 0; i < (header->n_targets); i++)  
	{
		int CC=(header->target_len[i])+500;
		depth[i] = new unsigned short int [CC]; //开辟列  
		for (int32_t j =0 ; j< CC ; j++)
		{
			depth[i][j]=0;
		}
	}

	cout<<"new the memory done"<<endl;


	int FileNum=(paraFA04->List).size();
	for (int kij=0 ;  kij< FileNum ;  kij++ )
	{
		string  line = (paraFA04->List)[kij];
		if (line.length()<=0)  { continue  ; }

		cout <<"Begin Bam :"<<line<<endl;
		bam_hdr_t *headerA;
		bam1_t *aln = bam_init1();

		samFile *InBam = sam_open(line.c_str(), "r");
		headerA = sam_hdr_read(InBam);

		if ((header->n_targets)!=(headerA->n_targets))
		{
			cerr<<"Error1 :Bam header should be the same\n"<<BamPath<<"\n"<<line<<endl;
			for(int i = 0; i <(header->n_targets); i++)
			{
				delete[] depth[i];  
			}
			delete[] depth;
			delete paraFA04 ;
			return 1;
		}
		bool NoSameBam=false;
		for(int i = 0; i < (header->n_targets); i++) 
		{
			if (strcmp(header->target_name[i],headerA->target_name[i])!=0)
			{
				NoSameBam=true;
				cerr<<header->target_name[i]<<"\t"<<headerA->target_name[i]<<endl;
				break ;
			}
		}
		if (NoSameBam)
		{
			cerr<<"Error2 :Bam header should be the same\n"<<BamPath<<"\n"<<line<<endl;
			for(int i = 0; i <(header->n_targets); i++)
			{
				delete[] depth[i];  
			}
			delete[] depth;
			delete paraFA04 ;
			return 1;
		}


		uint32_t *cigar;
		while (sam_read1(InBam, header, aln) >= 0)
		{
			if ((aln->core).tid < 0) {continue ;}
			if ( (aln->core).qual < (paraFA04->InInt))
			{
				continue ;
			}

			cigar = bam_get_cigar(aln);
			int32_t This=((aln->core).pos);
			for(int i=0; i < aln->core.n_cigar;++i)
			{				
				int cig=bam_cigar_op(cigar[i]);
				int ncig = bam_cigar_oplen(cigar[i]);

				if  (cig==BAM_CMATCH || cig==BAM_CEQUAL  || cig==BAM_CDIFF )
				{
					int end=This+ncig;
					for (  ; This<end;This++)
					{
						depth[(aln->core).tid][This]++;
					}
					continue;
				}
				else if  (cig==BAM_CDEL ||   cig==BAM_CREF_SKIP )
				{
					This=This+ncig;
				}
			}
		}

		sam_close(InBam);
		bam_destroy1(aln);
		bam_hdr_destroy(headerA);
	}

	cout <<"ALL Bam Read done"<<endl;

	for(int i = 0; i <(header->n_targets); i++)
	{
		bool S_E=false ;
		int Start=0;
		int End=0;
		for (int j =0 ; j< (header->target_len[i]) ; j++)
		{
			if ( (!S_E)   &&  ( depth[i][j] > (paraFA04->InInt)) )
			{
				continue ;
			}
			else if ( (S_E)   &&  ( depth[i][j] > (paraFA04->InInt))  )
			{
				S_E=false ;
				if ((End-Start)>MinLength)
				{
					OUT<<(header->target_name[i])<<"\t"<<Start<<"\t"<<End<<"\n";
				}
			}
			else if   ( (!S_E)   &&  ( depth[i][j] < (paraFA04->InInt)) )
			{
				S_E=true ;
				Start=j ;  End=j;
			}
			else 
			{
				End=j;
			}
		}
		if (S_E)
		{
			if ((End-Start)>MinLength)
			{
				OUT<<(header->target_name[i])<<"\t"<<Start<<"\t"<<End<<endl;
			}
		}
	}

	OUT.close();

	//释放开辟的资源  
	for(int i = 0; i <(header->n_targets); i++)
	{
		delete[] depth[i];  
	}
	delete[] depth;  


	bam_hdr_destroy(header);
	delete paraFA04 ;
	return 0;
}
#endif // bamLDR_H_  //
///////// swimming in the sky and flying in the sea ////////////





