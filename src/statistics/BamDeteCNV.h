#ifndef bamCNV_H_
#define bamCNV_H_
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
#include "Bmath/pnorm.c"
#include <cstdlib>
#include <zlib.h>
#include <htslib/sam.h>
#include <htslib/kseq.h>
#include <stdio.h>
#include <algorithm>
#include <math.h>

using namespace std;
typedef unsigned long long ubit64_t;
//KSEQ_INIT(gzFile, gzread)
/*
int  bamCNV_help()
{
	cout <<""
		"\n"
		"\tUsage: DeteCNV  -List  <bam.list>  -OutPut  <outPrefix>\n"
		"\n"
		"\t\t-InList      <str>     Input Bam/Sam File List\n"
		"\t\t-InFile      <str>     Input Bam/Sam File File[repeat]\n"
		"\t\t-Ref         <str>     In Ref.fa If Want Out Ref N-base ratio info\n"
		"\t\t-OutPut      <str>     OutPut CNV File Prefix\n"
		"\n"
		"\t\t-Ratio       <float>   DepthRatio to judge breakpoint of merge adjacent[0.45]\n"
		"\t\t-ChrDepth              Use Chr Depth as mean depth,default whole genome\n"
		"\n"
		"\t\t-MinLength   <int>     Min Length of CNV length [1800]\n"
		"\t\t-PValue      <float>   PValue of CNV Depth bias[0.02]\n"
		"\t\t-MinQ        <int>     Filter the low mapQ read in bam[10]\n"
		"\n"
		"\t\t-help                  Show this help [hewm2008]\n"
		"\n";
	return 1;
}
*/
void bamCNV_help()
{
cout<<""
"\n"
"\tUsage: DeteCNV  -l  <bam.list>  -r  <Ref.fa>  -o  <outPrefix>\n"
"\tUsage: DeteCNV  -i  <A.bam B.bam>  -r  <Ref.fa>  -o  <outPrefix>\n"
"\n"
"\t\t-i    <str>      input SAM/BAM files, delimited by space\n"
"\t\t-l    <str>      input list of SAM/BAM files\n"
"\t\t-o    <str>      prefix of output file\n"
"\t\t-r    <str>      input reference FASTA to get N-base ratio \n"
"\n"
"\t\t-f    <float>    depthRatio to judge breakpoint of merge adjacent[0.45]\n"
"\t\t-c               for each chromosome, use its own mean of depth into calculation\n"
"\t\t                 default would use the mean of depth of the whole genome\n"
"\n"
"\t\t-m    <int>      set the minimum length of CNV, default [1800]\n"
"\t\t-p    <float>    p-value of CNV depth bias, default [0.02]\n"
"\t\t-q    <int>      the quality to filter reads, default [10]\n"
"\n"
"\t\t-h               show more details for help\n"
"\n"
"\n";
}



void More_bamCNV_help()
{
cout<<""
"\n"
"\t\t1. DeteCNV -i <A.bam B.bam> -r  <Ref.fa> -o AAA -q Q -m M\n"
"\t\tDeteCNV -l <bam.list> -r  <Ref.fa> -o AAA -q Q -m M\n"
"\t\tThis will generate two files (raw cnv file and filtered cnv file) and output the result to current directory. Output files are named with the prefix AAA.\n"
"\t\t(1.1)      the reads with quality lower than Q will be removed from analysis, default value of Q is 10.\n"
"\t\t(1.2)      it will not be considered a CNV with the minimum length less than M, default value of M is 1800.\n"
"\n"
"\n"
"\n";
}


int bamCNV_help01(int argc, char **argv , In3str1v * paraFA04 )
{
	if (argc <=1 ) {bamCNV_help();return 0;}
	(paraFA04->InStr1)="0.02";
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

		if (flag  == "InList" ||  flag  == "List"  ||  flag  == "l")
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
		else if (flag  ==  "OutPut"  ||  flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr2=argv[i];
		}
		else if (flag  ==  "Ref"   ||  flag  == "r")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr3=argv[i];
		}
		else if (flag  ==  "PValue" ||  flag  == "p")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			(paraFA04->InStr1)=argv[i];
		}
		else if (flag  ==  "Ratio"  ||  flag  == "f")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InF=atof(argv[i]);
		}
		else if (flag  ==  "MinLength"  ||  flag  == "m")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InInt2=atoi(argv[i]);
		}
		else if (flag  ==  "MinQ"  ||  flag  == "q")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InInt=atoi(argv[i]);
		}
		else if (flag  ==  "ChrDepth" ||  flag  == "c" )
		{
			paraFA04->TF=false ;
		}
		else if (flag == "help" ||  flag  == "h")
		{
			More_bamCNV_help();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((file_count<1) || (paraFA04->InStr2).empty() ||  (paraFA04->InStr3).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	(paraFA04->InStr2)=add_Asuffix(paraFA04->InStr2);
	return 1 ;
}


int bamCNV_main(int argc, char *argv[])
//int main(int argc, char *argv[])
{
	In3str1v *paraFA04 = new In3str1v;
	paraFA04->InInt=10;
	paraFA04->InInt2=1800;
	paraFA04->InF=0.45;
	if ((bamCNV_help01(argc, argv, paraFA04)==0))
	{
		delete paraFA04 ;
		return 0 ;
	}

	double FilterPValue=atof((paraFA04->InStr1).c_str()) ;
	if  (FilterPValue >0.5 )
	{
		FilterPValue=0.5;
		cout<<"warning P value should no too big, reset 0.5"<<endl;
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
		for (int32_t j =0 ; j< CC  ; j++)
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
			delete paraFA04;
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


	string PrefixO=(paraFA04->InStr2).substr(0,(paraFA04->InStr2).length()-3);
	string StatF=PrefixO+".tmp";
	ofstream  OUTSTAT (StatF.c_str());

	cout<<"#ID\tLength\tCoverageBase\tTotalDepth\tCoverage%\tMeanDepth"<<endl;

	map <unsigned short int,ubit64_t> MapDepthSum ;
	map <unsigned short int,ubit64_t> :: iterator Th_mapiT ;
	ubit64_t  Sum_chrLen=0;
	ubit64_t Sum_CNVerBase=0;
	ubit64_t Sum_Depth=0;


	struct Region
	{
		int start;
		int end ;
		double depth;
	};

	vector<Region> ThisChrRegion ;
	vector<Region> ArryRegion ;


	for(int i = 0; i <(header->n_targets); i++)
	{
		map <unsigned short int,ubit64_t> MapDepthChr;
		map <unsigned short int,ubit64_t> :: iterator mapIt_Depth;
		Region One ;
		One.start=0;One.end=0;One.depth=depth[i][0];
		ArryRegion.clear();


		for (int j =0 ; j< (header->target_len[i]) ; j++)
		{
			mapIt_Depth=MapDepthChr.find(depth[i][j]);
			if  (mapIt_Depth==MapDepthChr.end())
			{
				MapDepthChr.insert( map <unsigned short int,ubit64_t>  :: value_type (depth[i][j],1));
			}
			else
			{
				(mapIt_Depth->second)++;
			}

			if (One.depth==depth[i][j])
			{
				One.end=j;
			}
			else
			{
				ArryRegion.push_back(One);
				One.start=j;One.end=j;One.depth=depth[i][j];
			}
		}

		ArryRegion.push_back(One);


		ubit64_t CNVBase=0;
		ubit64_t TotalDepth=0;
		for (mapIt_Depth=MapDepthChr.begin(); mapIt_Depth!=MapDepthChr.end(); mapIt_Depth++)
		{
			if ((mapIt_Depth->first)!=0)
			{
				CNVBase+=mapIt_Depth->second;
			}
			Th_mapiT=MapDepthSum.find(mapIt_Depth->first);
			TotalDepth+=(mapIt_Depth->second)*(mapIt_Depth->first);
			if (Th_mapiT==MapDepthSum.end())
			{
				MapDepthSum.insert( map <unsigned short int,ubit64_t>  :: value_type (mapIt_Depth->first,mapIt_Depth->second));
			}
			else
			{
				(Th_mapiT->second)+=(mapIt_Depth->second);
			}
		}

		double CNVerage=CNVBase*100.0/(header->target_len[i]);
		double MeanDepth=TotalDepth*1.0/(header->target_len[i]);
		Sum_chrLen+=(header->target_len[i]);
		Sum_CNVerBase+=CNVBase ;
		Sum_Depth+=TotalDepth;
		cout<<(header->target_name[i])<<"\t"<<(header->target_len[i])<<"\t"<<CNVBase<<"\t"<<TotalDepth<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\t"<<CNVerage<<"\t"<<MeanDepth<<endl;


		ThisChrRegion.clear();

		bool RunFlag=false ;
		int ii=0;
		int ZV_size = 0 ;
		bool VecAB=true ;
		for ( ii=0 ; ii< 11 ; ii++ )
		{
			if (RunFlag) {break ;}
			RunFlag=true ;
			if (VecAB)
			{
				ThisChrRegion.clear();
				ZV_size =ArryRegion.size();
				One.start=ArryRegion[0].start;
				One.end=ArryRegion[0].end;
				One.depth=ArryRegion[0].depth;
				VecAB=false ;
				for (int j=1 ; j< ZV_size ; j++)
				{
					if  (fabs(One.depth-ArryRegion[j].depth)<1.5)
					{
						int LengthA=ArryRegion[j].end-ArryRegion[j].start+1;
						int LengthB=One.end-One.start+1;
						One.depth=(One.depth*LengthB+ArryRegion[j].depth*LengthA)/(LengthA+LengthB);
						One.end=ArryRegion[j].end;
						RunFlag=false ;
					}
					else 
					{
						ThisChrRegion.push_back(One);
						One.start=ArryRegion[j].start;
						One.end=ArryRegion[j].end;
						One.depth=ArryRegion[j].depth;
					}
				}
				ThisChrRegion.push_back(One);
			}
			else
			{
				ArryRegion.clear();
				VecAB=true;
				ZV_size =ThisChrRegion.size();
				One.start=ThisChrRegion[0].start;
				One.end=ThisChrRegion[0].end;
				One.depth=ThisChrRegion[0].depth;

				for (int j=1 ; j< ZV_size ; j++)
				{
					if  (fabs(One.depth-ThisChrRegion[j].depth)<1.5)
					{
						int LengthA=ThisChrRegion[j].end-ThisChrRegion[j].start+1;
						int LengthB=One.end-One.start+1;
						One.depth=(One.depth*LengthB+ThisChrRegion[j].depth*LengthA)/(LengthA+LengthB);
						One.end=ThisChrRegion[j].end;
						RunFlag=false ;
					}
					else 
					{
						ArryRegion.push_back(One);
						One.start=ThisChrRegion[j].start;
						One.end=ThisChrRegion[j].end;
						One.depth=ThisChrRegion[j].depth;
					}
				}
				ArryRegion.push_back(One);
			}
		}


		RunFlag=false ;
		double AA=0.8 ;
		double binDouble=0.1;
		if ( (paraFA04->InF) >  AA )
		{
			AA=paraFA04->InF;
		}
		else
		{
			binDouble=(AA-(paraFA04->InF))/10;
		}		
		double BB=1/AA;

		for (   ; ii< 20 ; ii++)
		{
			if  (AA > paraFA04->InF )
			{
				RunFlag=false ; 
				AA=AA-binDouble;
			}
			BB=1/AA;
			if (RunFlag) {break ;}			
			RunFlag=true;

			if (VecAB)
			{
				ThisChrRegion.clear();
				VecAB=false ;
				ZV_size =ArryRegion.size();
				One.start=ArryRegion[0].start;
				One.end=ArryRegion[0].end;
				One.depth=ArryRegion[0].depth;
				for (int j=1 ; j< ZV_size ; j++)
				{
					double ratio = 0;
					if (ArryRegion[j].depth>0)
					{
						ratio=One.depth/ArryRegion[j].depth;
					}
					if  (ratio > AA  && ratio < BB )
					{
						int LengthA=ArryRegion[j].end-ArryRegion[j].start+1;
						int LengthB=One.end-One.start+1;
						One.depth=(One.depth*LengthB+ArryRegion[j].depth*LengthA)/(LengthA+LengthB);
						One.end=ArryRegion[j].end;
						RunFlag=false ;
					}
					else 
					{
						ThisChrRegion.push_back(One);
						One.start=ArryRegion[j].start;
						One.end=ArryRegion[j].end;
						One.depth=ArryRegion[j].depth;
					}
				}
				ThisChrRegion.push_back(One);
			}
			else
			{
				ArryRegion.clear();
				VecAB=true ;
				ZV_size =ThisChrRegion.size();
				One.start=ThisChrRegion[0].start;
				One.end=ThisChrRegion[0].end;
				One.depth=ThisChrRegion[0].depth;

				for (int j=1 ; j< ZV_size ; j++)
				{
					double ratio = 0;
					if (ArryRegion[j].depth>0)
					{
						ratio=One.depth/ThisChrRegion[j].depth;
					}
					if  (ratio > AA  && ratio < BB )
					{
						int LengthA=ThisChrRegion[j].end-ThisChrRegion[j].start+1;
						int LengthB=One.end-One.start+1;
						One.depth=(One.depth*LengthB+ThisChrRegion[j].depth*LengthA)/(LengthA+LengthB);
						One.end=ThisChrRegion[j].end;
						RunFlag=false;
					}
					else 
					{
						ArryRegion.push_back(One);
						One.start=ThisChrRegion[j].start;
						One.end=ThisChrRegion[j].end;
						One.depth=ThisChrRegion[j].depth;
					}
				}
				ArryRegion.push_back(One);
			}
		}




		///                meger Low Depth  /// 
		AA=MeanDepth*0.1;
		if (VecAB)
		{
			ThisChrRegion.clear();
			ZV_size =ArryRegion.size();
			One.start=ArryRegion[0].start;
			One.end=ArryRegion[0].end;
			One.depth=ArryRegion[0].depth;
			VecAB=false ;
			for (int j=1 ; j< ZV_size ; j++)
			{
				if  ( One.depth < AA  && ArryRegion[j].depth < AA )
				{
					int LengthA=ArryRegion[j].end-ArryRegion[j].start+1;
					int LengthB=One.end-One.start+1;
					One.depth=(One.depth*LengthB+ArryRegion[j].depth*LengthA)/(LengthA+LengthB);
					One.end=ArryRegion[j].end;
				}
				else 
				{
					ThisChrRegion.push_back(One);
					One.start=ArryRegion[j].start;
					One.end=ArryRegion[j].end;
					One.depth=ArryRegion[j].depth;
				}
			}
			ThisChrRegion.push_back(One);
		}
		else
		{
			ArryRegion.clear();VecAB=true ;
			ZV_size =ThisChrRegion.size();
			One.start=ThisChrRegion[0].start;
			One.end=ThisChrRegion[0].end;
			One.depth=ThisChrRegion[0].depth;
			for (int j=1 ; j< ZV_size ; j++)
			{
				if  ( One.depth < AA  && ThisChrRegion[j].depth < AA )
				{
					int LengthA=ThisChrRegion[j].end-ThisChrRegion[j].start+1;
					int LengthB=One.end-One.start+1;
					One.depth=(One.depth*LengthB+ThisChrRegion[j].depth*LengthA)/(LengthA+LengthB);
					One.end=ThisChrRegion[j].end;
				}
				else 
				{
					ArryRegion.push_back(One);
					One.start=ThisChrRegion[j].start;
					One.end=ThisChrRegion[j].end;
					One.depth=ThisChrRegion[j].depth;
				}
			}
			ArryRegion.push_back(One);
		}



		RunFlag=false ;
		AA=paraFA04->InF;
		BB=1/AA;

		for (  ; ii< 100; ii++)
		{
			if (RunFlag) {break ;}			
			RunFlag=true;
			if (VecAB)
			{
				ThisChrRegion.clear();VecAB=false;
				ZV_size =ArryRegion.size();
				One.start=ArryRegion[0].start;
				One.end=ArryRegion[0].end;
				One.depth=ArryRegion[0].depth;
				for (int j=1 ; j< ZV_size; j++)
				{
					double ratio = 0;
					if (ArryRegion[j].depth>0)
					{
						ratio=One.depth/ArryRegion[j].depth;
					}
					if  (ratio > AA  && ratio < BB )
					{
						int LengthA=ArryRegion[j].end-ArryRegion[j].start+1;
						int LengthB=One.end-One.start+1;
						One.depth=(One.depth*LengthB+ArryRegion[j].depth*LengthA)/(LengthA+LengthB);
						One.end=ArryRegion[j].end;
						RunFlag=false ;
					}
					else 
					{
						ThisChrRegion.push_back(One);
						One.start=ArryRegion[j].start;
						One.end=ArryRegion[j].end;
						One.depth=ArryRegion[j].depth;
					}
				}
				ThisChrRegion.push_back(One);
				ArryRegion.clear();
			}
			else
			{
				ArryRegion.clear();VecAB=true;
				ZV_size =ThisChrRegion.size();
				One.start=ThisChrRegion[0].start;
				One.end=ThisChrRegion[0].end;
				One.depth=ThisChrRegion[0].depth;

				for (int j=1 ; j< ZV_size ; j++)
				{
					double ratio = 0;
					if (ArryRegion[j].depth>0)
					{
						ratio=One.depth/ThisChrRegion[j].depth;
					}
					if  (ratio > AA  && ratio < BB )
					{
						int LengthA=ThisChrRegion[j].end-ThisChrRegion[j].start+1;
						int LengthB=One.end-One.start+1;
						One.depth=(One.depth*LengthB+ThisChrRegion[j].depth*LengthA)/(LengthA+LengthB);
						One.end=ThisChrRegion[j].end;
						RunFlag=false;
					}
					else 
					{
						ArryRegion.push_back(One);
						One.start=ThisChrRegion[j].start;
						One.end=ThisChrRegion[j].end;
						One.depth=ThisChrRegion[j].depth;
					}
				}
				ArryRegion.push_back(One);
				ThisChrRegion.clear();
			}
		}




		double CNV_N;
		if (VecAB)
		{
			ZV_size=ArryRegion.size();
			for (int j=0 ; j< ZV_size ; j++)
			{
				CNV_N=ArryRegion[j].depth/MeanDepth;
				OUTSTAT<<(header->target_name[i])<<"\t"<<ArryRegion[j].start<<"\t"<<ArryRegion[j].end<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<ArryRegion[j].depth<<"\t"<<CNV_N<<"\n";
			}
		}
		else
		{
			ZV_size=ThisChrRegion.size();
			for (int j=0 ; j< ZV_size ; j++)
			{
				CNV_N=ThisChrRegion[j].depth/MeanDepth;
				OUTSTAT<<(header->target_name[i])<<"\t"<<ThisChrRegion[j].start<<"\t"<<ThisChrRegion[j].end<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<ThisChrRegion[j].depth<<"\t"<<CNV_N<<"\n";
			}
		}
	}
	double CNVerage=Sum_CNVerBase*100.0/Sum_chrLen ;
	double MeanDepth=Sum_Depth*1.0/Sum_chrLen ;
	cout<<"#Genome\t"<<Sum_chrLen<<"\t"<<Sum_CNVerBase<<"\t"<<Sum_Depth<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\t"<<CNVerage<<"\t"<<MeanDepth<<endl;


	OUTSTAT.close();
	for(int i = 0; i <(header->n_targets); i++)
	{
		delete[] depth[i];  
	}
	delete[] depth;
	bam_hdr_destroy(header);











	//                             //
	//                             //
	//                             //
	//                             //
	//                             //


	long  double  Lsd=0.0 ;  ubit64_t L_count=0;
	long  double  Rsd=0.0 ;  ubit64_t R_count=0;
	long  double diff=0.0 ;
	long  double this_count=0.0 ;
	string tempsss;
	int Start , End;
	double DepthRegion;


	igzstream INB ((StatF).c_str(),ifstream::in);
	if(!INB.good())
	{
		cerr << "open IN File error: "<<StatF<<endl;
		return 1;
	}

	if (paraFA04->TF)
	{
		while(!INB.eof())
		{
			string  line ;
			getline(INB,line);
			if ( (line.length()<=0)  ||  (line[0] == '#') )  { continue ;}
			istringstream  stream  (line,istringstream::in);
			stream >> tempsss>>Start>>End ;
			stream >>DepthRegion  ;
			this_count=DepthRegion/MeanDepth;
			diff=this_count-1.0 ;
			int NumC=(End-Start+1);
			if  (diff>0)
			{
				R_count+=NumC;
				Rsd+=(diff*diff)*NumC;
			}
			else if  (diff<0)
			{
				L_count+=NumC;
				Lsd+=(diff*diff)*NumC;
			}
		}
	}
	else
	{
		while(!INB.eof())
		{
			string  line ;
			getline(INB,line);
			if ( (line.length()<=0)  ||  (line[0] == '#') )  { continue ;}
			istringstream  stream  (line,istringstream::in);
			stream >> tempsss>>Start>>End ;
			stream >>DepthRegion>> this_count ;
			diff=this_count-1.0 ;
			int NumC=(End-Start+1);
			if  (diff>0)
			{
				R_count+=NumC;
				Rsd+=(diff*diff)*NumC;
			}
			else if  (diff<0)
			{
				L_count+=NumC;
				Lsd+=(diff*diff)*NumC;
			}
		}

	}

	INB.close();
	Rsd=sqrt(Rsd/R_count);
	Lsd=sqrt(Lsd/L_count);


	gzFile fp;
	kseq_t *seq;
	int l;
	fp = gzopen((paraFA04->InStr3).c_str(), "r");
	seq = kseq_init(fp);
	map <string ,string> ChrSeq ;

	while ((l = kseq_read(seq)) >= 0)
	{
		string RefSeq=(seq->seq.s);
		string chr=(seq->name.s);
		ChrSeq.insert(map <string,string>  ::value_type(chr,RefSeq));	
	}
	kseq_destroy(seq);
	gzclose(fp);

	map <string ,string > ::iterator Mapit;

	string StatRaw=PrefixO+".raw.gz";
	ogzstream OUT  (StatRaw.c_str());
	OUT<<"#Chr\tStart\tEnd\tLength\tMeanDepth\tCNV\tNratio\tZScore\tPValue\n";
	OUT<<"##MeanDepth: "<<MeanDepth<<" ; Rsd: " <<Rsd <<" ; Lsd:" <<Lsd<<endl;

	string StatCNV=PrefixO+".cnv.gz";
	ogzstream OUTCNV  (StatCNV.c_str());
	OUTCNV<<"#Chr\tStart\tEnd\tLength\tMeanDepth\tCNV\tNratio\tZScore\tPValue\n";


	igzstream INC ((StatF).c_str(),ifstream::in);

	if(!INC.good())
	{
		cerr << "open IN File error: "<<StatF<<endl;
		return 1;
	}

	double ZScore,PValue;

	if (paraFA04->TF)
	{
		while(!INC.eof())
		{
			string  line ;
			getline(INC,line);
			if ( (line.length()<=0)  ||  (line[0] == '#') )  { continue ;}
			istringstream  stream  (line,istringstream::in);
			stream >>tempsss>>Start>>End ;
			stream >>DepthRegion  ;
			this_count=DepthRegion/MeanDepth;
			diff=this_count-1.0 ;
			int NumC=(End-Start+1);
			Mapit=ChrSeq.find(tempsss);
			int Ascii[256] = {0};
			for(int ix=Start ; ix<=End ; ix++)
			{
				Ascii[(Mapit->second)[ix]]++;
			}
			int N_Num=Ascii['N']+Ascii['n'];
			double Nratio=N_Num*100.0/NumC;
			PValue=1.0;
			if  (diff>0)
			{
				PValue=(1-pnorm5(this_count, 1.0 , Rsd  , 1 , 0 ));
				ZScore=diff/Rsd;
			}
			else if  (diff<0)
			{
				PValue=(pnorm5(this_count, 1.0 , Lsd  , 1 , 0 ));
				ZScore=diff/Lsd;
			}
			else
			{
				PValue=0.99;
				ZScore=0.0;
			}
			OUT<<tempsss<<"\t"<<Start<<"\t"<<End<<"\t"<<NumC<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<DepthRegion<<"\t"<<this_count<<"\t"<<Nratio<<"\t"<<ZScore<<"\t"<<setprecision(6)<<PValue<<"\n";
			if (PValue< FilterPValue  && NumC > (paraFA04->InInt2)  &&  Nratio<10 )
			{
				OUTCNV<<tempsss<<"\t"<<Start<<"\t"<<End<<"\t"<<NumC<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<DepthRegion<<"\t"<<this_count<<"\t"<<Nratio<<"\t"<<ZScore<<"\t"<<setprecision(6)<<PValue<<"\n";
			}
		}
	}
	else
	{
		while(!INC.eof())
		{
			string  line ;
			getline(INC,line);
			if ( (line.length()<=0)  ||  (line[0] == '#') )  { continue ;}
			istringstream  stream  (line,istringstream::in);
			stream >>tempsss>>Start>>End ;
			stream >>DepthRegion >>this_count ;
			diff=this_count-1.0 ;
			int NumC=(End-Start+1);
			Mapit=ChrSeq.find(tempsss);
			int Ascii[256] = {0};
			for(int ix=Start ; ix<=End ; ix++)
			{
				Ascii[(Mapit->second)[ix]]++;
			}

			int N_Num=Ascii['N']+Ascii['n'];
			double Nratio=N_Num*100.0/NumC;

			if  (diff>0)
			{
				PValue=(1-pnorm5(this_count, 1.0 , Rsd  , 1 , 0 ));
				ZScore=diff/Rsd;
			}
			else if  (diff<0)
			{
				PValue=(pnorm5(this_count, 1.0 , Lsd  , 1 , 0 ));
				ZScore=diff/Lsd;
			}
			else
			{
				PValue=0.99;
				ZScore=0.0;
			}

			OUT<<tempsss<<"\t"<<Start<<"\t"<<End<<"\t"<<NumC<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<DepthRegion<<"\t"<<this_count<<"\t"<<Nratio<<"\t"<<ZScore<<"\t"<<setprecision(6)<<PValue<<"\n";
			if (PValue< FilterPValue  && NumC > (paraFA04->InInt2)  &&  Nratio<10 )
			{
				OUTCNV<<tempsss<<"\t"<<Start<<"\t"<<End<<"\t"<<NumC<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<DepthRegion<<"\t"<<this_count<<"\t"<<Nratio<<"\t"<<ZScore<<"\t"<<setprecision(6)<<PValue<<"\n";
			}
		}

	}

	INC.close();
	OUTCNV.close();	
	OUT.close();


	string Rmrm="rm -rf "+ StatF ;
	std::system(Rmrm.c_str()) ;

	delete paraFA04 ;
	return 0;
}
#endif // bamCNV_H_  //
///////// swimming in the sky and flying in the sea ////////////




