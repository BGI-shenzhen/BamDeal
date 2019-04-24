#ifndef BamSplit_H_
#define BamSplit_H_
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
#include <htslib/kstring.h>
#include <stdio.h>

using namespace std;


void BamSplit_help()
{
	cout<<""
		"\n"
		"\tUsage: bamSplit  -l <bam.list>\n"
		"\tUsage: bamSplit  -i A.bam  B.bam\n"
		"\n"
		"\t\t-i    <str>     input SAM/BAM file(s)\n"
		"\t\t-l    <str>     Input list of SAM/BAM files\n"
		"\n"
		"\t\t-o    <str>     output directory, default [PWD]\n"
		"\t\t-s              set output files in SAM format, default in BAM format.\n"
		"\t\t-q    <int>     reads with quality lower than this would be classified to unmap.bam, default [10]\n"
		"\t\t-r              reset output files headers by remove the chromosomes not in the output files\n"
		"\n"
		"\t\t-h              show more details for help [hewm2008 v1.04]\n"
		"\n";
}


void More_HelpBamSplit()
{
	cout<<""
		"\n"
		"\t\t1. BamSplit -i A.bam B.bam -q Q\n"
		"\t\t   This will split A.bam and B.bam by chromosome and output the splitting BAM files to current directory. \n"
		"\t\t   (1.1)  A.bam and B.bam should have same header.\n"
		"\t\t   (1.2)  The splitting files would keep the same headers as input files.\n"
		"\t\t   (1.3)  Reads in A.bam and B.bam with quality lower than Q would be outputted to the file unmap.bam.\n"
		"\n"
		"\t\t2. BamSplit -l bam.list -r -s \n"
		"\t\t   This will split the files in bam.list by chromosome and output the splitting SAM files to current directory. \n"
		"\t\t   (2.1)  files in bam.list should have same header. \n"
		"\t\t          For example, if user has two BAM files A.bam and B.bam to split, bam.list should be formatted as:\n"
		"\t\t            ./A.bam\n"
		"\t\t            ./B.bam\n"
		"\t\t   (2.2)  with -s added, the output files would be compressed SAM files.\n"
		"\t\t   (2.3)  with -r added, the header of each output file would keep its own chromosome.\n"
		"\n";
}



int BamSplit_help01(int argc, char **argv , In3str1v * paraBamS )
{
	if (argc <2 ) {BamSplit_help();return 0;}
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

		if (flag  == "InList" || flag  == "l" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			string A=argv[i];
			file_count+=(ReadList(A ,(paraBamS->List)));
		}
		else if (flag  == "InFile"  || flag  == "i" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			string A=argv[i];
			(paraBamS->List).push_back(A);
			file_count++;
			bool RunT=true;
			while(RunT)
			{
				if ( (i+1) < argc  && (argv[i+1][0]!='-'))
				{
					i++;
					A=argv[i];
					(paraBamS->List).push_back(A);
					file_count++;
				}
				else
				{
					RunT=false;
				}
			}
		}
		else if (flag  ==  "OutDir"    || flag  == "o" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraBamS->InStr2=argv[i];			
		}
		else if (flag  ==  "MinQ"  || flag  == "q")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraBamS->InInt=atoi(argv[i]);
		}
		else if (flag  ==  "ReSetHead"  || flag  == "r")
		{
			paraBamS->TF=true ;
		}
		else if (flag  ==  "OutSam"  || flag  == "s" )
		{
			paraBamS->TF2=false ;
		}
		else if (flag == "help"  || flag  == "h")
		{
			More_HelpBamSplit();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((file_count<1) )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	return 1 ;
}


int Split2Sam( In3str1v * paraBamS )
{

	string  BamPath=(paraBamS->List)[0];

	bam_hdr_t *header;
	samFile *BamIn = hts_open(BamPath.c_str(), "r");
	header = sam_hdr_read(BamIn);
	sam_close(BamIn);

	ogzstream *Soap2Chr = new ogzstream[(header->n_targets)];
	ogzstream UnMapOut ;
	string UUUtmp=(paraBamS->InStr2)+"/Unmap.sam.gz";
	UnMapOut.open(UUUtmp.c_str());
	for (int j=0; j<(header->n_targets) ; j++)
	{
		UnMapOut<<"@SQ\tSN:"<<(header->target_name[j])<<"\tLN:"<<(header->target_len[j])<<endl;
	}


	if  (!(paraBamS->TF))
	{

		for (int i=0; i<(header->n_targets) ; i++)
		{
			string aaa=(paraBamS->InStr2)+"/"+(header->target_name[i])+".sam.gz" ;
			Soap2Chr[i].open(aaa.c_str());
			if  (!Soap2Chr[i].good())
			{
				cerr<<"Can't open follow output:\n"<<aaa<<endl;
			}
			for (int j=0; j<(header->n_targets) ; j++)
			{
				Soap2Chr[i]<<"@SQ\tSN:"<<(header->target_name[j])<<"\tLN:"<<(header->target_len[j])<<endl;
			}
		}
	}
	else
	{
		for (int i=0; i<(header->n_targets) ; i++)
		{
			string aaa=(paraBamS->InStr2)+"/"+(header->target_name[i])+".sam.gz" ;
			Soap2Chr[i].open(aaa.c_str());
			if  (!Soap2Chr[i].good())
			{
				cerr<<"Can't open follow output:\n"<<aaa<<endl;
			}
			Soap2Chr[i]<<"@SQ\tSN:"<<(header->target_name[i])<<"\tLN:"<<(header->target_len[i])<<endl;
		}
	}

	int FileNum=(paraBamS->List).size();
	for (int ii=0 ; ii< FileNum ; ii++)
	{
		string  line =(paraBamS->List)[ii];
		if (line.length()<=0)  { continue  ; }		
		cout <<"Begin Bam :"<<line<<endl;
		bam_hdr_t *headerA;
		bam1_t *aln = bam_init1();

		samFile *InBam = sam_open(line.c_str(), "r");
		headerA = sam_hdr_read(InBam);

		if ((header->n_targets)!=(headerA->n_targets))
		{
			cerr<<"Error1 :Bam header should be the same\n"<<BamPath<<"\n"<<line<<endl;
			delete paraBamS ;
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
			delete paraBamS ;
			return 1;
		}

		while (sam_read1(InBam, header, aln) >= 0)
		{
			kstring_t Kstr = { 0, 0, NULL };
			int A=sam_format1(header, aln, &Kstr);

			if ( ((aln->core).tid < 0)  || (aln->core).qual < (paraBamS->InInt) )
			{
				UnMapOut<<Kstr.s<<endl;
			}
			else
			{
				Soap2Chr[(aln->core).tid]<<Kstr.s<<endl;
			}
		}
		sam_close(InBam);
		bam_destroy1(aln);
		bam_hdr_destroy(headerA);
	}



	UnMapOut.close();
	for (int i=0; i<(header->n_targets) ; i++)
	{
		Soap2Chr[i].close() ;
	}

	bam_hdr_destroy(header);
	delete [] Soap2Chr ;

	return 1 ;
}




int Split2Bam( In3str1v * paraBamS )
{

	string  BamPath=(paraBamS->List)[0];

	bam_hdr_t *header;
	samFile *BamIn = hts_open(BamPath.c_str(), "r");
	header = sam_hdr_read(BamIn);
	sam_close(BamIn);

	htsFile **OUTBam = new  htsFile *[(header->n_targets)];
	bam_hdr_t **OUTheader=  new  bam_hdr_t *[(header->n_targets)];

	htsFile *UUM;
	string aaa=(paraBamS->InStr2)+"/UnMap.bam" ;
	UUM=hts_open(aaa.c_str(), "wb");
	int tmp=sam_hdr_write(UUM, header);


	for (int i=0; i<(header->n_targets) ; i++)
	{
		string aaa=(paraBamS->InStr2)+"/"+(header->target_name[i])+".bam" ;
		OUTBam[i]=hts_open(aaa.c_str(), "wb");
		OUTheader[i]=NULL;

		kstring_t str = { 0, 0, NULL };
		kputsn("@SQ\tSN:", 7, &str);
		kputs((header->target_name[i]), &str);
		kputsn("\tLN:", 4, &str);
		kputw(header->target_len[i], &str);
		kputc('\n', &str);

		OUTheader[i]= sam_hdr_parse(str.l, str.s);
		OUTheader[i]->l_text = str.l; OUTheader[i]->text = str.s;
		OUTheader[i]=sam_hdr_sanitise(OUTheader[i]);
		if (sam_hdr_write(OUTBam[i], OUTheader[i]) < 0)
		{
			cerr<<"Can't Write bam Header:\n"<<aaa<<endl;
			exit(1);
		}
	}
	//sam_close(OUTBam[0]);  return 1 ;

	int FileNum=(paraBamS->List).size();
	for (int ii=0 ; ii< FileNum ; ii++)
	{
		string  line =(paraBamS->List)[ii];
		if (line.length()<=0)  { continue  ; }

		cout <<"Begin Bam :"<<line<<endl;
		bam_hdr_t *headerA;
		bam1_t *aln = bam_init1();

		samFile *InBam = sam_open(line.c_str(), "r");
		headerA = sam_hdr_read(InBam);

		if ((header->n_targets)!=(headerA->n_targets))
		{
			cerr<<"Error1 :Bam header should be the same\n"<<BamPath<<"\n"<<line<<endl;
			delete paraBamS ;
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
			delete paraBamS ;
			return 1;
		}

		while (sam_read1(InBam, header, aln) >= 0)
		{
			int32_t i=(aln->core).tid ;
			if ( ((aln->core).qual < (paraBamS->InInt))  || (i < 0 ) )
			{
				int AA=bam_write1(UUM->fp.bgzf, aln);
			}			
			else
			{
				
				if ((aln->core).tid==(aln->core).mtid)
				{
					(aln->core).mtid=0;
				}
				else if  ( (aln->core).mtid >0 )
				{
					(aln->core).mtid=-1;
				}
				(aln->core).tid=0;

				int AA=bam_write1(OUTBam[i]->fp.bgzf, aln);
			}
		}
		sam_close(InBam);
		bam_destroy1(aln);

		bam_hdr_destroy(headerA);
	}


	sam_close(UUM);
	for (int i=0; i<(header->n_targets) ; i++)
	{
		sam_close(OUTBam[i]);
		bam_hdr_destroy(OUTheader[i]);
	}

	bam_hdr_destroy(header);
	delete [] OUTheader ;
	delete [] OUTBam ;
	return 1 ;
}





int Split2BamSameHeader( In3str1v * paraBamS )
{

	string  BamPath=(paraBamS->List)[0];

	bam_hdr_t *header;
	samFile *BamIn = hts_open(BamPath.c_str(), "r");
	header = sam_hdr_read(BamIn);
	sam_close(BamIn);

	htsFile **OUTBam = new  htsFile *[(header->n_targets)];

	htsFile *UUM;
	string aaa=(paraBamS->InStr2)+"/UnMap.bam" ;
	UUM=hts_open(aaa.c_str(), "wb");
	int tmp=sam_hdr_write(UUM, header);

	for (int i=0; i<(header->n_targets) ; i++)
	{
		aaa=(paraBamS->InStr2)+"/"+(header->target_name[i])+".bam" ;
		OUTBam[i]=hts_open(aaa.c_str(), "wb");		
		if (sam_hdr_write(OUTBam[i], header) < 0)
		{
			cerr<<"Can't Write bam Header:\n"<<aaa<<endl;
			exit(1);
		}
	}
	//sam_close(OUTBam[0]);  return 1 ;

	int FileNum=(paraBamS->List).size();
	for (int ii=0 ; ii< FileNum ; ii++)
	{
		string  line =(paraBamS->List)[ii];
		if (line.length()<=0)  { continue  ; }

		cout <<"Begin Bam :"<<line<<endl;
		bam_hdr_t *headerA;
		bam1_t *aln = bam_init1();

		samFile *InBam = sam_open(line.c_str(), "r");
		headerA = sam_hdr_read(InBam);

		if ((header->n_targets)!=(headerA->n_targets))
		{
			cerr<<"Error1 :Bam header should be the same\n"<<BamPath<<"\n"<<line<<endl;
			delete paraBamS ;
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
			delete paraBamS ;
			return 1;
		}

		while (sam_read1(InBam, header, aln) >= 0)
		{
			int32_t i=(aln->core).tid ;
			if (( (aln->core).qual < (paraBamS->InInt)  )   ||  (i < 0 ))
			{
				int A=bam_write1(UUM->fp.bgzf, aln);
			}
			else
			{
				int A=bam_write1(OUTBam[i]->fp.bgzf, aln);
			}
		}
		sam_close(InBam);
		bam_destroy1(aln);

		bam_hdr_destroy(headerA);
	}

	sam_close(UUM);

	for (int i=0; i<(header->n_targets) ; i++)
	{
		sam_close(OUTBam[i]);
	}

	bam_hdr_destroy(header);
	delete [] OUTBam ;
	return 1 ;
}








int bamSplit_main(int argc, char *argv[])
{
	In3str1v *paraBamS = new In3str1v;
	paraBamS->InInt=10;
	paraBamS->TF=false ;
	paraBamS->InStr2="./" ;
	if ((BamSplit_help01(argc, argv, paraBamS)==0))
	{
		delete paraBamS ;
		return 0 ;
	}

	if (paraBamS->TF2)
	{
		if  (paraBamS->TF)
		{
			Split2Bam( paraBamS );
		}
		else
		{
			Split2BamSameHeader( paraBamS );
		}
	}
	else
	{
		Split2Sam( paraBamS );
	}

	delete paraBamS ;
	return 0;
}
#endif // BamSplit_H_  //
///////// swimming in the sky and flying in the sea ////////////





