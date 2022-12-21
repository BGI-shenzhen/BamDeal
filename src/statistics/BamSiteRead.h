#ifndef bamSiteRead_H_
#define bamSiteRead_H_
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
void SiteRead_help()
{
	cout<<""
		"\n"
		"\tUsage:   BamSiteRead -l <bam.list> -b SNP.vcf -o  <outPrefix>\n"
		"\tUsage:   BamSiteRead -i < A.bam B.bam > -b SNP.vcf -o <outPrefix>\n"
		"\n"
		"\t\t-b     <str>     input SNP Site(chr site)\n"
		"\t\t-i     <str>     input SAM/BAM files, delimited by space\n"
		"\t\t-l     <str>     input list of SAM/BAM files\n"
		"\t\t-o     <str>     prefix of output file\n"
		"\n"
		"\t\t-q     <int>     the quality to filter reads, default [10]\n"
		"\t\t-d               Filter the PCR or optical duplicate read(DUP)\n"
		"\t\t-s               Filter the secondary alignment read(SECONDARY)\n"
		"\t\t-c               Filter the not passing quality controls read(QCFAIL)\n"
		"\n"
		"\t\t-h               show more details for help\n"
		"\n";
}

void More_SiteRead_help()
{
	cout<<""
		"\n"
		"\t\t1. SiteRead -i <A.bam B.bam> -b SNP.vcf -o AAA -q -1\n"
		"\t\tSiteRead -l <bam.list> -b SNP.vcf  -o AAA -q -1\n"
		"\t\tThis will generate Read(with read and base info covering the SNP site).\n"
		"\t\t...\n"
		"\n"
		"\n";
}

int SiteRead_help01(int argc, char **argv , In3str1v * paraFA04 )
{
	if (argc <=1 ) {SiteRead_help();return 0;}
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

		if (flag  == "InList" ||  flag  == "List" || flag  == "l")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			string A=argv[i];
			file_count+=(ReadList(A ,(paraFA04->List)));
		}
		else if (flag  == "InFile" || flag  == "i" )
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
		else if (flag  ==  "OutPut"  || flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr2=argv[i];
		}
		else if (flag  ==  "Bed"  || flag  == "b")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr1=argv[i];
		}
		else if (flag  ==  "MinQ"  || flag  == "q")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InInt=atoi(argv[i]);
		}
		else if (flag  ==  "RmDup"  || flag  == "d" )
		{
			paraFA04->TF=false ;
		}
		else if (flag  ==  "RmSec"   ||  flag  =="s")
		{
			paraFA04->TF2=false ;
		}
		else if (flag  ==  "RmQCFAIL"   ||  flag  =="c")
		{
			paraFA04->InInt2=1 ;
		}
		else if (flag == "help"  || flag  == "h")
		{
			More_SiteRead_help();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if ( (file_count<1) || (paraFA04->InStr2).empty()  ||  (paraFA04->InStr1).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	(paraFA04->InStr2)=add_Asuffix(paraFA04->InStr2);
	return 1 ;
}


int SiteRead_main(int argc, char *argv[])
//int main(int argc, char *argv[])
{
	In3str1v *paraFA04 = new In3str1v;
	paraFA04->InInt=10;
	if ((SiteRead_help01(argc, argv, paraFA04)==0))
	{
		delete paraFA04 ;
		return 0 ;
	}

	(paraFA04->InStr2)=(paraFA04->InStr2).substr(0,(paraFA04->InStr2).length()-3);
	string path=(paraFA04->InStr2);	
	string ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);

	string PrefixO=path;

	if (ext == "Read" || ext == "read")
	{
		PrefixO=path.substr(0,path.rfind('.'));
	}

	string  outFaDepth=PrefixO+".Read.gz";
	ogzstream  OUT (outFaDepth.c_str());

	if((!OUT.good()))
	{
		cerr << "open OUT File error: "<<outFaDepth<<endl;
		delete  paraFA04 ; return  1;
	}

	OUT<<"#CHROM\tPOS\tBase\tReadID\tMapPOS\tMapQ\tSeq\tBaseSeqQ\n";
	uint8_t Base[16] = {0,65,67,0,71,0,0,0,84,0,0,0,0,0,0,78};


	string  BamPath=(paraFA04->List)[0];
	bam_hdr_t *header;
	samFile *BamIn = hts_open(BamPath.c_str(), "r");
	header = sam_hdr_read(BamIn);
	sam_close(BamIn);


	cout<<"begin new the memory ...\n";
	bool QCFAIL=true;
	if ((paraFA04->InInt2)!=0) {QCFAIL=false;}

	bool **depth = new bool *[(header->n_targets)]; //开辟行  
	map <string,int> Name2Int ;
	for(int i = 0; i < (header->n_targets); i++)
	{
		int CC=(header->target_len[i])+500;
		depth[i] = new bool [CC]; //开辟列  
		for (int32_t j =0 ; j< CC ; j++)
		{
			depth[i][j]=false;
		}
		string Name=header->target_name[i];
		Name2Int[Name]=i;
	}

	cout<<"new the memory done"<<endl;

	igzstream VCFIN ((paraFA04->InStr1).c_str(),ifstream::in);
	if (VCFIN.fail())
	{
		cerr << "open VCFIN File error: "<<(paraFA04->InStr1)<<endl;
		delete  paraFA04 ; return 1;
	}
	int Position=0;
	while(!VCFIN.eof())
	{
		string  line ;
		getline(VCFIN,line);
		if (line.length()<=0 )  { continue  ; }
		else if ((line[0] =='#'))
		{
			continue ;
		}
		istringstream isone (line,istringstream::in);
		string chr_id;
		isone>>chr_id>>Position;
		Position--;
		depth[Name2Int[chr_id]][Position]=true;
	}
	VCFIN.close();
	int shiftQ=33;
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
		uint8_t *seq;
		int iref =0;
		int iseq=0;  int  arryFlag=0;


		//#define BAM_CMATCH      0
		//#define BAM_CINS        1
		//#define BAM_CDEL        2
		//#define BAM_CREF_SKIP   3
		//#define BAM_CSOFT_CLIP  4
		//#define BAM_CHARD_CLIP  5
		//#define BAM_CPAD        6
		//#define BAM_CEQUAL      7
		//#define BAM_CDIFF       8


		if ((paraFA04->TF)   &&  (paraFA04->TF2) &&  (QCFAIL) )
		{
			while (sam_read1(InBam, header, aln) >= 0)
			{
				if ((aln->core).tid < 0) {continue ;}
				if ( (aln->core).qual < (paraFA04->InInt)){continue ;}
				cigar = bam_get_cigar(aln);
				seq = bam_get_seq(aln);
				iref=((aln->core).pos);
				iseq=0;
				for(int i=0; i < aln->core.n_cigar;++i)
				{				
					int cig=bam_cigar_op(cigar[i]);
					int ncig = bam_cigar_oplen(cigar[i]);
					if (cig==BAM_CMATCH || cig==BAM_CEQUAL  || cig==BAM_CDIFF)
					{
						int end=iref+ncig;
						for (  ; iref<end;iref++)
						{
							if (depth[(aln->core).tid][iref])
							{
								string readID=bam_get_qname(aln);
								if  ( (((aln)->core.flag&64) != 0))
								{
									readID=readID + "/1" ;
								}
								else if ((((aln)->core.flag&128) != 0) )
								{
									readID=readID + "/2" ;
								}
								char BaseNow=Base[bam_seqi(seq, iseq)];

								uint8_t  *seqQ=bam_get_qual(aln);
								uint8_t  *seq=bam_get_seq(aln);
								string BaseSeq="";
								string BaseSeqQ="";
								for(int io=0; io < (aln->core).l_qseq ; io++)
								{
									BaseSeqQ=BaseSeqQ+char(seqQ[io]+shiftQ);
									BaseSeq=BaseSeq+char(Base[bam_seqi(seq, io)]);
								}
								Position=iref+1;uint32_t MQ=(aln->core).qual;
								OUT<<header->target_name[(aln->core).tid]<<"\t"<<Position<<"\t"<<BaseNow<<"\t"<<readID<<"\t"<<(aln->core).pos<<"\t"<<MQ<<"\t"<<BaseSeq<<"\t"<<BaseSeqQ<<"\n";
							}
							iseq++;
						}
						continue;
					}
					else if  (cig==BAM_CDEL ||   cig==BAM_CREF_SKIP )
					{
						iref=iref+ncig;
						continue;
					}
					else if(cig==BAM_CINS || cig==BAM_CSOFT_CLIP )
					{
						iseq += ncig;
						continue;
					}
					//        BAM_CHARD_CLIP   BAM_CPAD
				}
			}

		}
		else if ((!paraFA04->TF)   &&  (paraFA04->TF2)   &&  (QCFAIL) )
		{
			while (sam_read1(InBam, header, aln) >= 0)
			{
				if ((aln->core).tid < 0) {continue ;}
				if ( (aln->core).qual < (paraFA04->InInt)){	continue ;}
				if ( (aln->core).flag   & 0x400 )  { continue ; }
				cigar = bam_get_cigar(aln);
				seq = bam_get_seq(aln);
				iref=((aln->core).pos);
				iseq=0;
				for(int i=0; i < aln->core.n_cigar;++i)
				{				
					int cig=bam_cigar_op(cigar[i]);
					int ncig = bam_cigar_oplen(cigar[i]);
					if (cig==BAM_CMATCH || cig==BAM_CEQUAL  || cig==BAM_CDIFF)
					{
						int end=iref+ncig;
						for (  ; iref<end;iref++)
						{
							if (depth[(aln->core).tid][iref])
							{
								string readID=bam_get_qname(aln);
								if  ( (((aln)->core.flag&64) != 0))
								{
									readID=readID + "/1" ;
								}
								else if ((((aln)->core.flag&128) != 0) )
								{
									readID=readID + "/2" ;
								}
								char BaseNow=Base[bam_seqi(seq, iseq)];
								uint8_t  *seqQ=bam_get_qual(aln);
								uint8_t  *seq=bam_get_seq(aln);
								string BaseSeq="";
								string BaseSeqQ="";
								for(int io=0; io < (aln->core).l_qseq ; io++)
								{
									BaseSeqQ=BaseSeqQ+char(seqQ[io]+shiftQ);
									BaseSeq=BaseSeq+char(Base[bam_seqi(seq, io)]);
								}							
								Position=iref+1;uint32_t MQ=(aln->core).qual;
								OUT<<header->target_name[(aln->core).tid]<<"\t"<<Position<<"\t"<<BaseNow<<"\t"<<readID<<"\t"<<(aln->core).pos<<"\t"<<MQ<<"\t"<<BaseSeq<<"\t"<<BaseSeqQ<<"\n";
							}

							iseq++;
						}
						continue;
					}
					else if  (cig==BAM_CDEL ||   cig==BAM_CREF_SKIP )
					{
						iref=iref+ncig;
						continue;
					}
					else if(cig==BAM_CINS || cig==BAM_CSOFT_CLIP )
					{
						iseq += ncig;
						continue;
					}
				}
			}
		}
		else if ((!paraFA04->TF)   &&  (!paraFA04->TF2)    &&  (QCFAIL)  )
		{
			while (sam_read1(InBam, header, aln) >= 0)
			{
				if ((aln->core).tid < 0) {continue ;}
				if ( (aln->core).qual < (paraFA04->InInt)){	continue ;}
				if ( (aln->core).flag   & 0x400 )  { continue ; }
				if ( (aln->core).flag   & 0x100 )  { continue ; }
				cigar = bam_get_cigar(aln);
				seq = bam_get_seq(aln);
				iref=((aln->core).pos);
				iseq=0;
				for(int i=0; i < aln->core.n_cigar;++i)
				{				
					int cig=bam_cigar_op(cigar[i]);
					int ncig = bam_cigar_oplen(cigar[i]);
					if (cig==BAM_CMATCH || cig==BAM_CEQUAL  || cig==BAM_CDIFF)
					{
						int end=iref+ncig;
						for (  ; iref<end;iref++)
						{
							if (depth[(aln->core).tid][iref])
							{
								string readID=bam_get_qname(aln);
								if  ( (((aln)->core.flag&64) != 0))
								{
									readID=readID + "/1" ;
								}
								else if ((((aln)->core.flag&128) != 0) )
								{
									readID=readID + "/2" ;
								}
								char BaseNow=Base[bam_seqi(seq, iseq)];
								uint8_t  *seqQ=bam_get_qual(aln);
								uint8_t  *seq=bam_get_seq(aln);
								string BaseSeq="";
								string BaseSeqQ="";
								for(int io=0; io < (aln->core).l_qseq ; io++)
								{
									BaseSeqQ=BaseSeqQ+char(seqQ[io]+shiftQ);
									BaseSeq=BaseSeq+char(Base[bam_seqi(seq, io)]);
								}
								Position=iref+1;uint32_t MQ=(aln->core).qual;
								OUT<<header->target_name[(aln->core).tid]<<"\t"<<Position<<"\t"<<BaseNow<<"\t"<<readID<<"\t"<<(aln->core).pos<<"\t"<<MQ<<"\t"<<BaseSeq<<"\t"<<BaseSeqQ<<"\n";

							}


							iseq++;
						}
						continue;
					}
					else if  (cig==BAM_CDEL ||   cig==BAM_CREF_SKIP )
					{
						iref=iref+ncig;
						continue;
					}
					else if(cig==BAM_CINS || cig==BAM_CSOFT_CLIP )
					{
						iseq += ncig;
						continue;
					}
				}
			}
		}
		else if ((paraFA04->TF)   &&  (!paraFA04->TF2)    &&  (QCFAIL)  )
		{
			while (sam_read1(InBam, header, aln) >= 0)
			{
				if ((aln->core).tid < 0) {continue ;}
				if ( (aln->core).qual < (paraFA04->InInt)){	continue ;}
				if ( (aln->core).flag   & 0x100 )  { continue ; }
				cigar = bam_get_cigar(aln);
				seq = bam_get_seq(aln);
				iref=((aln->core).pos);
				iseq=0;
				for(int i=0; i < aln->core.n_cigar;++i)
				{				
					int cig=bam_cigar_op(cigar[i]);
					int ncig = bam_cigar_oplen(cigar[i]);
					if (cig==BAM_CMATCH || cig==BAM_CEQUAL  || cig==BAM_CDIFF)
					{
						int end=iref+ncig;
						for (  ; iref<end;iref++)
						{
							if (depth[(aln->core).tid][iref])
							{
								string readID=bam_get_qname(aln);
								if  ( (((aln)->core.flag&64) != 0))
								{
									readID=readID + "/1" ;
								}
								else if ((((aln)->core.flag&128) != 0) )
								{
									readID=readID + "/2" ;
								}
								char BaseNow=Base[bam_seqi(seq, iseq)];
								uint8_t  *seqQ=bam_get_qual(aln);
								uint8_t  *seq=bam_get_seq(aln);
								string BaseSeq="";
								string BaseSeqQ="";
								for(int io=0; io < (aln->core).l_qseq ; io++)
								{
									BaseSeqQ=BaseSeqQ+char(seqQ[io]+shiftQ);
									BaseSeq=BaseSeq+char(Base[bam_seqi(seq, io)]);
								}
								Position=iref+1;uint32_t MQ=(aln->core).qual;
								OUT<<header->target_name[(aln->core).tid]<<"\t"<<Position<<"\t"<<BaseNow<<"\t"<<readID<<"\t"<<(aln->core).pos<<"\t"<<MQ<<"\t"<<BaseSeq<<"\t"<<BaseSeqQ<<"\n";

							}



							iseq++;
						}
						continue;
					}
					else if  (cig==BAM_CDEL ||   cig==BAM_CREF_SKIP )
					{
						iref=iref+ncig;
						continue;
					}
					else if(cig==BAM_CINS || cig==BAM_CSOFT_CLIP )
					{
						iseq += ncig;
						continue;
					}
				}
			}
		}
		else if ((!paraFA04->TF)   &&  (!paraFA04->TF2)    &&  (!QCFAIL)  )
		{
			while (sam_read1(InBam, header, aln) >= 0)
			{
				if ((aln->core).tid < 0) {continue ;}
				if ( (aln->core).qual < (paraFA04->InInt)){	continue ;}
				if ( (aln->core).flag   & 0x400 )  { continue ; }
				if ( (aln->core).flag   & 0x100 )  { continue ; }
				if ( (aln->core).flag   & 0x200 )  { continue ; }
				cigar = bam_get_cigar(aln);
				seq = bam_get_seq(aln);
				iref=((aln->core).pos);
				iseq=0;
				for(int i=0; i < aln->core.n_cigar;++i)
				{				
					int cig=bam_cigar_op(cigar[i]);
					int ncig = bam_cigar_oplen(cigar[i]);
					if (cig==BAM_CMATCH || cig==BAM_CEQUAL  || cig==BAM_CDIFF)
					{
						int end=iref+ncig;
						for (  ; iref<end;iref++)
						{
							if (depth[(aln->core).tid][iref])
							{
								string readID=bam_get_qname(aln);
								if  ( (((aln)->core.flag&64) != 0))
								{
									readID=readID + "/1" ;
								}
								else if ((((aln)->core.flag&128) != 0) )
								{
									readID=readID + "/2" ;
								}
								char BaseNow=Base[bam_seqi(seq, iseq)];
								uint8_t  *seqQ=bam_get_qual(aln);
								uint8_t  *seq=bam_get_seq(aln);
								string BaseSeq="";
								string BaseSeqQ="";
								for(int io=0; io < (aln->core).l_qseq ; io++)
								{
									BaseSeqQ=BaseSeqQ+char(seqQ[io]+shiftQ);
									BaseSeq=BaseSeq+char(Base[bam_seqi(seq, io)]);
								}
								Position=iref+1;uint32_t MQ=(aln->core).qual;
								OUT<<header->target_name[(aln->core).tid]<<"\t"<<Position<<"\t"<<BaseNow<<"\t"<<readID<<"\t"<<(aln->core).pos<<"\t"<<MQ<<"\t"<<BaseSeq<<"\t"<<BaseSeqQ<<"\n";



							}


							iseq++;
						}
						continue;
					}
					else if  (cig==BAM_CDEL ||   cig==BAM_CREF_SKIP )
					{
						iref=iref+ncig;
						continue;
					}
					else if(cig==BAM_CINS || cig==BAM_CSOFT_CLIP )
					{
						iseq += ncig;
						continue;
					}
				}
			}
		}
		else if ((!paraFA04->TF)   &&  (paraFA04->TF2)    &&  (!QCFAIL)  )
		{
			while (sam_read1(InBam, header, aln) >= 0)
			{
				if ((aln->core).tid < 0) {continue ;}
				if ( (aln->core).qual < (paraFA04->InInt)){	continue ;}
				if ( (aln->core).flag   & 0x400 )  { continue ; }
				if ( (aln->core).flag   & 0x200 )  { continue ; }
				cigar = bam_get_cigar(aln);
				seq = bam_get_seq(aln);
				iref=((aln->core).pos);
				iseq=0;
				for(int i=0; i < aln->core.n_cigar;++i)
				{				
					int cig=bam_cigar_op(cigar[i]);
					int ncig = bam_cigar_oplen(cigar[i]);
					if (cig==BAM_CMATCH || cig==BAM_CEQUAL  || cig==BAM_CDIFF)
					{
						int end=iref+ncig;
						for (  ; iref<end;iref++)
						{
							if (depth[(aln->core).tid][iref])
							{
								string readID=bam_get_qname(aln);
								if  ( (((aln)->core.flag&64) != 0))
								{
									readID=readID + "/1" ;
								}
								else if ((((aln)->core.flag&128) != 0) )
								{
									readID=readID + "/2" ;
								}
								char BaseNow=Base[bam_seqi(seq, iseq)];
								uint8_t  *seqQ=bam_get_qual(aln);
								uint8_t  *seq=bam_get_seq(aln);
								string BaseSeq="";
								string BaseSeqQ="";
								for(int io=0; io < (aln->core).l_qseq ; io++)
								{
									BaseSeqQ=BaseSeqQ+char(seqQ[io]+shiftQ);
									BaseSeq=BaseSeq+char(Base[bam_seqi(seq, io)]);
								}
								Position=iref+1;uint32_t MQ=(aln->core).qual;
								OUT<<header->target_name[(aln->core).tid]<<"\t"<<Position<<"\t"<<BaseNow<<"\t"<<readID<<"\t"<<(aln->core).pos<<"\t"<<MQ<<"\t"<<BaseSeq<<"\t"<<BaseSeqQ<<"\n";



							}


							iseq++;
						}
						continue;
					}
					else if  (cig==BAM_CDEL ||   cig==BAM_CREF_SKIP )
					{
						iref=iref+ncig;
						continue;
					}
					else if(cig==BAM_CINS || cig==BAM_CSOFT_CLIP )
					{
						iseq += ncig;
						continue;
					}
				}
			}
		}
		else if ((paraFA04->TF)   &&  (!paraFA04->TF2)    &&  (!QCFAIL)  )
		{
			while (sam_read1(InBam, header, aln) >= 0)
			{
				if ((aln->core).tid < 0) {continue ;}
				if ( (aln->core).qual < (paraFA04->InInt)){	continue ;}
				if ( (aln->core).flag   & 0x100 )  { continue ; }
				if ( (aln->core).flag   & 0x200 )  { continue ; }
				cigar = bam_get_cigar(aln);
				seq = bam_get_seq(aln);
				iref=((aln->core).pos);
				iseq=0;
				for(int i=0; i < aln->core.n_cigar;++i)
				{				
					int cig=bam_cigar_op(cigar[i]);
					int ncig = bam_cigar_oplen(cigar[i]);
					if (cig==BAM_CMATCH || cig==BAM_CEQUAL  || cig==BAM_CDIFF)
					{
						int end=iref+ncig;
						for (  ; iref<end;iref++)
						{
							if (depth[(aln->core).tid][iref])
							{
								string readID=bam_get_qname(aln);
								if  ( (((aln)->core.flag&64) != 0))
								{
									readID=readID + "/1" ;
								}
								else if ((((aln)->core.flag&128) != 0) )
								{
									readID=readID + "/2" ;
								}
								char BaseNow=Base[bam_seqi(seq, iseq)];
								uint8_t  *seqQ=bam_get_qual(aln);
								uint8_t  *seq=bam_get_seq(aln);
								string BaseSeq="";
								string BaseSeqQ="";
								for(int io=0; io < (aln->core).l_qseq ; io++)
								{
									BaseSeqQ=BaseSeqQ+char(seqQ[io]+shiftQ);
									BaseSeq=BaseSeq+char(Base[bam_seqi(seq, io)]);
								}
								Position=iref+1;uint32_t MQ=(aln->core).qual;
								OUT<<header->target_name[(aln->core).tid]<<"\t"<<Position<<"\t"<<BaseNow<<"\t"<<readID<<"\t"<<(aln->core).pos<<"\t"<<MQ<<"\t"<<BaseSeq<<"\t"<<BaseSeqQ<<"\n";



							}


							iseq++;
						}
						continue;
					}
					else if  (cig==BAM_CDEL ||   cig==BAM_CREF_SKIP )
					{
						iref=iref+ncig;
						continue;
					}
					else if(cig==BAM_CINS || cig==BAM_CSOFT_CLIP )
					{
						iseq += ncig;
						continue;
					}
				}
			}
		}
		else if ((paraFA04->TF)   &&  (paraFA04->TF2)    &&  (!QCFAIL)  )
		{
			while (sam_read1(InBam, header, aln) >= 0)
			{
				if ((aln->core).tid < 0) {continue ;}
				if ( (aln->core).qual < (paraFA04->InInt)){	continue ;}
				if ( (aln->core).flag   & 0x200 )  { continue ; }
				cigar = bam_get_cigar(aln);
				seq = bam_get_seq(aln);
				iref=((aln->core).pos);
				iseq=0;
				for(int i=0; i < aln->core.n_cigar;++i)
				{				
					int cig=bam_cigar_op(cigar[i]);
					int ncig = bam_cigar_oplen(cigar[i]);
					if (cig==BAM_CMATCH || cig==BAM_CEQUAL  || cig==BAM_CDIFF)
					{
						int end=iref+ncig;
						for (  ; iref<end;iref++)
						{
							if (depth[(aln->core).tid][iref])
							{
								string readID=bam_get_qname(aln);
								if  ( (((aln)->core.flag&64) != 0))
								{
									readID=readID + "/1" ;
								}
								else if ((((aln)->core.flag&128) != 0) )
								{
									readID=readID + "/2" ;
								}
								char BaseNow=Base[bam_seqi(seq, iseq)];
								uint8_t  *seqQ=bam_get_qual(aln);
								uint8_t  *seq=bam_get_seq(aln);
								string BaseSeq="";
								string BaseSeqQ="";
								for(int io=0; io < (aln->core).l_qseq ; io++)
								{
									BaseSeqQ=BaseSeqQ+char(seqQ[io]+shiftQ);
									BaseSeq=BaseSeq+char(Base[bam_seqi(seq, io)]);
								}
								Position=iref+1;uint32_t MQ=(aln->core).qual;
								OUT<<header->target_name[(aln->core).tid]<<"\t"<<Position<<"\t"<<BaseNow<<"\t"<<readID<<"\t"<<(aln->core).pos<<"\t"<<MQ<<"\t"<<BaseSeq<<"\t"<<BaseSeqQ<<"\n";



							}


							iseq++;
						}
						continue;
					}
					else if  (cig==BAM_CDEL ||   cig==BAM_CREF_SKIP )
					{
						iref=iref+ncig;
						continue;
					}
					else if(cig==BAM_CINS || cig==BAM_CSOFT_CLIP )
					{
						iseq += ncig;
						continue;
					}
				}
			}
		}

		sam_close(InBam);
		bam_destroy1(aln);
		bam_hdr_destroy(headerA);
	}

	cout <<"ALL Bam Read done"<<endl;


	OUT.close();
	for(int i = 0; i <(header->n_targets); i++)
	{
		delete[] depth[i];  
	}
	delete[] depth;  

	bam_hdr_destroy(header);
	delete paraFA04 ;
	return 0;
}
#endif // bamSiteRead_H_  //
///////// swimming in the sky and flying in the sea ////////////





