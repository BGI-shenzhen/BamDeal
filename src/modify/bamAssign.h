#ifndef BamAssign_H_
#define BamAssign_H_
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

void BamAssign_help()
{
	cout<<""
		"\n"
		"\tUsage: bamAssign  -i  A.bam  B.bam  -a <assign.list> -r\n"
		"\tUsage: bamAssign  -l  <bam.list> -a <assign.list>\n"
		"\n"
		"\t\t-i     <str>     input SAM/BAM/CRAM file(s)\n"
		"\t\t-l     <str>     Input list of SAM/BAM/CRAM files\n"
		"\t\t-a     <str>     list indicating how to assign chromosomes to outputs\n"
		"\n"
		"\t\t-o     <str>     output directory, default [PWD]\n"
		"\t\t-q     <int>     reads with quality lower than this would be classified to unmap.bam, default [10]\n"
		"\t\t-r               reset output files headers by remove the chromosomes not in the output files\n"
		"\n"
		"\t\t-h               show more details for help [hewm2008 v1.02]\n"
		"\n"
		"\n";
}

void More_HelpBamAssign()
{
	cout<<""
		"\n"
		"\t\t1. bamAssign -i A.bam B.bam -a <assign.list> -q Q -r\n"
		"\t\t   bamAssign -l  <bam.list> -a <assign.list> -q Q -r\n"
		"\t\t   This will assign the chromosomes in the input SAM/BAM to different files according to the <assign.list> and output the results to current directory. \n"
		"\t\t    (1.1)  -l lists the input files. For example, if user has two input files A.bam and B.bam, bam.list should be formatted as:\n"
		"\t\t            ./A.bam\n"
		"\t\t            ./B.bam\n"
		"\t\t    (1.2)  Headers of the input SAM/BAMs should be the same.\n"
		"\t\t    (1.3)  The <assign.list> indicates how to assign chromosomes to outputs. \n"
		"\t\t           For example, if there are 5 chromosomes (Chr1, Chr2, Chr3, Chr4, Chr5) in the input, and the <assign.list> is formatted as:\n"
		"\t\t            Chr1  AAA\n"
		"\t\t            Chr2  AAA\n"
		"\t\t            Chr3  BBB\n"
		"\t\t            Chr4  BBB\n"
		"\t\t          bamAssign would output the reads of Chr1 and Chr2 to the file named AAA and output the reads of Chr3 and Chr4 to the file BBB. \n"
		"\t\t    (1.4) Chromosomes not in the <assign.list> (like Chr5 in this example) would be outputted to a file named NaAss.bam.\n"
		"\t\t    (1.5) Both unmapped reads and reads with quality lower than Q (default value of Q is 10) would be outputted to file named UnMap.bam.\n"
		"\t\t    (1.6) With -r added, the header of the output SAM/BAM will keep its own chromosomes. As the example above, header of AAA would contain only Chr1 and Chr2. Without -r, the header of output would be the same as the header of input.\n"
		"\n"
		"\n";
}

int BamAssign_help01(int argc, char **argv , In3str1v * paraBamS )
{
	if (argc <2 ) {BamAssign_help();return 0;}
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
		else if (flag  == "InFile" || flag  == "i")
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
		else if (flag  ==  "OutDir"  || flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraBamS->InStr2=argv[i];			
		}
		else if (flag  ==  "MinQ" || flag  == "q")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraBamS->InInt=atoi(argv[i]);
		}
		else if (flag  ==  "ReSetHead" || flag  == "r")
		{
			paraBamS->TF=false ;
		}
		else if (flag  ==  "chr2File"  || flag  == "a")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraBamS->InStr3=argv[i];
		}
		else if (flag == "help" || flag  == "h")
		{
			More_HelpBamAssign();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((file_count<1) ||  (paraBamS->InStr3).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}

	return 1 ;

}





int Assign2Bam( In3str1v * paraBamS )
{



	map <string , int > chr2Int;
	map <string , int > outFile2Int;
	map <string , int > :: iterator  Mapit ;

	igzstream LIST ((paraBamS->InStr3).c_str(),ifstream::in); // igzstream
	if (!LIST.good())
	{
		cerr << "open chr2outFile error: "<<(paraBamS->InStr3)<<endl;
		return  1 ;
	}
	int OutFileNum=0;
	int tmp=0;
	while(!LIST.eof())
	{
		string  line ;
		getline(LIST,line);
		if (line.length()<=0)  { continue  ; }
		if (line[0] == '#')  { continue  ; }
		vector<string> inf;
		split(line,inf," \t");
		Mapit=outFile2Int.find(inf[1]);
		if (Mapit!=outFile2Int.end())
		{
			tmp=Mapit->second;
		}
		else
		{
			tmp=OutFileNum;
			outFile2Int.insert( map <string ,int > :: value_type(inf[1],tmp));
			OutFileNum++;
		}

		Mapit=chr2Int.find(inf[0]);
		if (Mapit!=chr2Int.end())
		{
			if ((Mapit->second)!=tmp)
			{

				cerr<<"chr ID: "<<inf[0]<<"\tcan't to two outFile, please check it."<<endl;
				return 1;
			}
		}
		else
		{
			chr2Int.insert( map <string,int > :: value_type (inf[0],tmp) );
		}


	}
	LIST.close();





	string  BamPath=(paraBamS->List)[0];

	bam_hdr_t *header;
	samFile *BamIn = hts_open(BamPath.c_str(), "r");
	header = sam_hdr_read(BamIn);
	sam_close(BamIn);



	htsFile *UUM;
	string aaa=(paraBamS->InStr2)+"/UnMap.bam" ;
	UUM=hts_open(aaa.c_str(), "wb");
	tmp=sam_hdr_write(UUM, header);



	map <int,int > Int2IntMap ;
	map <int, int > :: iterator  MapitInt ;

	tmp=OutFileNum;
	for (int i=0; i<(header->n_targets) ; i++)
	{
		aaa=(header->target_name[i]);
		Mapit=chr2Int.find(aaa);
		if (Mapit!=chr2Int.end())
		{
			Int2IntMap.insert(map <int ,int > :: value_type(i,Mapit->second));
		}
		else
		{
			paraBamS->TF2=false;			
			Int2IntMap.insert( map <int ,int > :: value_type(i,tmp));
			cerr<<"warning :chr ID :"<<aaa<<" can't found in the chr2outfile List, will be assigned to the NaAss.bam"<<endl;
			chr2Int.insert( map <string ,int > :: value_type(aaa,tmp));
		}
	}

	if  (paraBamS->TF2)
	{
	}
	else
	{
		OutFileNum++;
		outFile2Int.insert(map <string ,int > :: value_type("NaAss",tmp));
	}


	htsFile **OUTBam = new  htsFile *[OutFileNum];


	if (paraBamS->TF)
	{
		for ( Mapit=outFile2Int.begin() ;Mapit!=outFile2Int.end()  ; Mapit++)
		{

			aaa=(paraBamS->InStr2)+"/"+(Mapit->first)+".bam" ;			
			OUTBam[(Mapit->second)]=hts_open(aaa.c_str(), "wb");
			tmp=sam_hdr_write(OUTBam[(Mapit->second)], header);
		}


		int FileNum=(paraBamS->List).size();
		for (int ii=0 ; ii< FileNum ; ii++)
		{
			string  line =(paraBamS->List)[ii];
			if (line.length()<=0)  { continue; }

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
				delete paraBamS;
				return 1;
			}

			while (sam_read1(InBam, header, aln) >= 0)
			{
				int  i=(aln->core).tid ;
				if ( ((aln->core).qual < (paraBamS->InInt))  || (i < 0 ) )
				{
					int AA=bam_write1(UUM->fp.bgzf, aln);
				}			
				else
				{
					MapitInt=Int2IntMap.find(i);
					int AA=bam_write1(OUTBam[MapitInt->second]->fp.bgzf, aln);
				}
			}
			sam_close(InBam);
			bam_destroy1(aln);

			bam_hdr_destroy(headerA);
		}


	}
	else
	{


		map <int,int > ChangIDMap ;

		for ( Mapit=outFile2Int.begin() ;Mapit!=outFile2Int.end()  ; Mapit++)
		{
			int Flag=0;
			kstring_t str = { 0, 0, NULL };
			for ( MapitInt=Int2IntMap.begin() ; MapitInt!=Int2IntMap.end(); MapitInt++)
			{
				if ((MapitInt->second)==(Mapit->second))
				{
					string chr=(header->target_name[MapitInt->first]);
					int TR=(header->target_len[MapitInt->first]);
					kputsn("@SQ\tSN:", 7, &str);
					kputs(chr.c_str(), &str);
					kputsn("\tLN:", 4, &str);
					kputw(TR, &str);
					kputc('\n', &str);
					ChangIDMap.insert(map <int,int > :: value_type(MapitInt->first,Flag));
					Flag++;
				}
			}

			aaa=(paraBamS->InStr2)+"/"+(Mapit->first)+".bam" ;			
			OUTBam[(Mapit->second)]=hts_open(aaa.c_str(), "wb");

			bam_hdr_t *NewHeader;
			NewHeader= sam_hdr_parse(str.l, str.s);
			NewHeader->l_text = str.l; NewHeader->text = str.s;
			NewHeader=sam_hdr_sanitise(NewHeader);			
			tmp=sam_hdr_write(OUTBam[(Mapit->second)], NewHeader);
			bam_hdr_destroy(NewHeader);

		}






		int FileNum=(paraBamS->List).size();
		for (int ii=0 ; ii< FileNum ; ii++)
		{
			string  line =(paraBamS->List)[ii];
			if (line.length()<=0)  { continue; }

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
				delete paraBamS;
				return 1;
			}

			while (sam_read1(InBam, header, aln) >= 0)
			{
				int  i=(aln->core).tid ;
				if ( ((aln->core).qual < (paraBamS->InInt))  || (i < 0 ) )
				{
					tmp=bam_write1(UUM->fp.bgzf, aln);
				}			
				else
				{

					if ((aln->core).tid>=0)
					{
						MapitInt=ChangIDMap.find((aln->core).tid);
						if (MapitInt!=ChangIDMap.end())
						{
							(aln->core).tid=MapitInt->second;
						}
						else
						{
							(aln->core).tid=-1;
						}
					}

					if ((aln->core).mtid>=0)
					{
						MapitInt=ChangIDMap.find((aln->core).mtid);
						if (MapitInt!=ChangIDMap.end())
						{
							(aln->core).mtid=MapitInt->second;
						}
						else
						{
							(aln->core).mtid=-1;
						}
					}

					MapitInt=Int2IntMap.find(i);
					tmp=bam_write1(OUTBam[MapitInt->second]->fp.bgzf, aln);

				}
			}
			sam_close(InBam);
			bam_destroy1(aln);

			bam_hdr_destroy(headerA);
		}















	}

	sam_close(UUM);
	for (int i=0; i<OutFileNum ; i++)
	{
		sam_close(OUTBam[i]);
	}

	delete [] OUTBam ;

	bam_hdr_destroy(header);


	return 1 ;
}






int bamAssign_main(int argc, char *argv[])
{
	In3str1v *paraBamS = new In3str1v;
	paraBamS->InInt=10;
	paraBamS->InStr2="./" ;
	if ((BamAssign_help01(argc, argv, paraBamS)==0))
	{
		delete paraBamS ;
		return 0 ;
	}

	Assign2Bam( paraBamS );

	delete paraBamS ;
	return 0;
}
#endif // BamAssign_H_  //
///////// swimming in the sky and flying in the sea ////////////





