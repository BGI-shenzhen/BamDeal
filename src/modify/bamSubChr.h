#ifndef BamSubChr_H_
#define BamSubChr_H_

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
#include <htslib/kstring.h>
#include <cstdlib>
#include <sys/select.h>

using namespace std;
typedef long long llong ;


void print_SubChr_Usage()
{
	cout<<""
		"\n"
		"\tUsage: SubChr -i <in.bam>  -d <list>  -o <out.bam>\n"
		"\tUsage: SubChr -i <in.bam>  -k <list>  -o <out.bam>\n"
		"\n"
		"\t\t-i      <str>     input SAM/BAM\n"
		"\t\t-o      <str>     output BAM\n"
		"\n"
		"\t\t-k      <str>     list of chromosomes to be kept\n"
		"\t\t-d      <str>     list of chromosomes to be deleted\n"
		"\n"
		"\t\t-u                remove unmapped reads\n"
		"\t\t-r                reset output headers by remove the chr(s) not in the out files\n"
		"\n"
		"\t\t-h                show more details for help\n"
		"\n";
}


void More_HelpSubChr()
{
	cout<<""
		"\n"
		"\t\t1. SubChr -i <in.bam> -d delete.list -o AAA -r\n"
		"\t\t   This will remove the chromosome(s) in the delete.list from input SAM/BAM and output the reads left to a file named AAA in current directory. For example, if user would like to remove Chr1 and Chr2 from input BAM file, the delete.list would be formatted as:\n"
		"\t\t     Chr1\n"
		"\t\t     Chr2\n"
		"\t\t  The header of the output SAM/BAM will also remove the name(s) of the chromosome(s) in the delete.list. As the example above, Chr1 and Chr2 would not show up in the header of output file.\n"
		"\n"
		"\t\t2. SubChr -i <in.bam> -k keep.list -o AAA -r\n"
		"\t\t   This will extract the chromosome(s) in the keep.list from input SAM/BAM and output them to a file named AAA in current directory. keep.list is formatted in the same way as delete.list shown above.\n"
		"\n"
		"\t\t3. SubChr -i <in.bam> -o AAA -u\n"
		"\t\t   This will remove the unmapped reads from input SAM/BAM and output the reads left to a file named AAA in current directory. \n"
		"\n";
}


int parse_Acmd_SubChr(int argc, char **argv, In3str1v * para_Xam04)
{
	if (argc <2 ) {print_SubChr_Usage();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","") ;

		if (flag  == "InPut"  || flag  == "i" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_Xam04->InStr1=argv[i];
		}
		else if (flag  ==  "OutPut" || flag  == "o" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0; }
			i++;
			para_Xam04->InStr2=argv[i];
		}
		else if (flag  ==  "KeepList"   || flag  == "k")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_Xam04->InStr3=(argv[i]);
			para_Xam04->InInt=(para_Xam04->InInt) | 0x1;
		}
		else if (flag  ==  "RmList"  || flag  == "d")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_Xam04->InStr3=(argv[i]);
			para_Xam04->InInt=(para_Xam04->InInt) | 0x2;
		}
		else if (flag  == "ReSetHead"  || flag  == "r")
		{
			para_Xam04->TF=false;
		}
		else if (flag  == "RmUnmap"  || flag  == "u")
		{
			para_Xam04->TF2=false;
		}
		else if (flag  == "help"  || flag  == "h")
		{
			More_HelpSubChr();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}

	if  ((para_Xam04->InStr1).empty() ||   (para_Xam04->InStr2).empty()  ||  (para_Xam04->InStr3).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;    
	}

	if (para_Xam04->InInt==3)
	{
		cerr<< "Para [ -RmList ]  and [ -KeepList ] can not use together"<<endl;
		return 0;
	}
	else if   (para_Xam04->InInt==0)
	{
		cerr<< "One of Para [ -RmList ]  and [ -KeepList ] must be seted"<<endl;
		return 0;
	}

	string path=(para_Xam04->InStr2);

	string ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);

	if (ext != "bam")
	{
		(para_Xam04->InStr2)=path+".bam";
	}

	return 1 ;


}


//programme entry
///////// swimming in the sky and flying in the sea ////////////
//int main(int argc, char **argv)
int bam_SubChr_main(int argc, char **argv)
{
	In3str1v * para_Xam04 = new In3str1v;
	if (parse_Acmd_SubChr(argc, argv , para_Xam04 )==0)
	{
		delete para_Xam04; 
		return 0;
	}

	samFile *outR1 = sam_open(para_Xam04->InStr2.c_str(), "wb");



	bam_hdr_t *header_SS;
	bam1_t *aln_SS = bam_init1();
	samFile *in_SS = sam_open((para_Xam04->InStr1).c_str(), "r");


	header_SS = sam_hdr_read(in_SS);


	igzstream LIST (para_Xam04->InStr3.c_str(),ifstream::in); // igzstream

	int soapfilecout=0;
	if (!LIST.good())
	{
		cerr << "open List error: "<<para_Xam04->InStr3<<endl;
	}

	map <string ,int > ChrInfo;


	while(!LIST.eof())
	{
		string  line ;
		getline(LIST,line);
		if (line.length()<=0)  { continue  ; }
		ChrInfo.insert( map <string ,int >  :: value_type (line,-1));	
		soapfilecout++;
	}
	LIST.close();


	string chrID ;
	int tmp=1;

	if ( para_Xam04->TF )
	{
		if (sam_hdr_write(outR1, header_SS) < 0) 
		{
			fprintf(stderr, "Error writing output.\n");
			exit(-1);
		}


		if ( (para_Xam04->InInt)==1 )
		{
			while((sam_read1(in_SS, header_SS, aln_SS) >= 0))
			{
				if ((aln_SS->core).tid >= 0)
				{ // chr
					chrID=header_SS->target_name[(aln_SS->core).tid];
					if (ChrInfo.find(chrID)!=ChrInfo.end())
					{
						tmp=bam_write1(outR1->fp.bgzf, aln_SS);
					}
				}
			}
		}
		else
		{
			while((sam_read1(in_SS, header_SS, aln_SS) >= 0))
			{
				if ((aln_SS->core).tid >= 0)
				{ // chr
					chrID=header_SS->target_name[(aln_SS->core).tid];
					if (ChrInfo.find(chrID)==ChrInfo.end())
					{
						tmp=bam_write1(outR1->fp.bgzf, aln_SS);
					}
				}
				else
				{
					if (para_Xam04->TF2)
					{
						tmp=bam_write1(outR1->fp.bgzf, aln_SS);
					}
				}
			}
		}

	}


	else
	{


		kstring_t str = { 0, 0, NULL };
		map <string , int > :: iterator MapIt ;
		int Flag=0;

		if ( (para_Xam04->InInt)==1)
		{
			for(int i = 0; i < (header_SS->n_targets); i++)
			{
				string chr=header_SS->target_name[i];
				int TR=header_SS->target_len[i];
				MapIt=ChrInfo.find(chr);
				if (MapIt!=ChrInfo.end())
				{
					kputsn("@SQ\tSN:", 7, &str);
					kputs((MapIt->first).c_str(), &str);
					kputsn("\tLN:", 4, &str);
					kputw(TR, &str);
					kputc('\n', &str);
					MapIt->second= Flag ;
					Flag++;
				}
			}



			bam_hdr_t *NewHeader;
			NewHeader= sam_hdr_parse(str.l, str.s);
			NewHeader->l_text = str.l; NewHeader->text = str.s;
			NewHeader=sam_hdr_sanitise(NewHeader);
			tmp=sam_hdr_write(outR1, NewHeader);
			bam_hdr_destroy(NewHeader);

			while((sam_read1(in_SS, header_SS, aln_SS) >= 0))
			{
				if ((aln_SS->core).tid >= 0)
				{ // chr
					chrID=header_SS->target_name[(aln_SS->core).tid];
					MapIt=ChrInfo.find(chrID);
					if (MapIt!=ChrInfo.end())
					{
						(aln_SS->core).tid=MapIt->second;
						if ((aln_SS->core).mtid>=0)
						{				
							chrID=header_SS->target_name[(aln_SS->core).mtid];
							MapIt=ChrInfo.find(chrID);
							if (MapIt!=ChrInfo.end())
							{
								(aln_SS->core).mtid=MapIt->second;
							}
							else
							{
								(aln_SS->core).mtid=-1;
							}
						}
						tmp=bam_write1(outR1->fp.bgzf, aln_SS);
					}
				}
			}



		}
		else
		{


			map <string , int >  NewChangID ;
			for(int i = 0; i < (header_SS->n_targets); i++)
			{
				string chr=header_SS->target_name[i];
				int TR=header_SS->target_len[i];
				MapIt=ChrInfo.find(chr);
				if (MapIt==ChrInfo.end())
				{
					kputsn("@SQ\tSN:", 7, &str);
					kputs(chr.c_str(), &str);
					kputsn("\tLN:", 4, &str);
					kputw(TR, &str);
					kputc('\n', &str);
					NewChangID.insert( map <string ,int >  :: value_type (chr,Flag));
					Flag++;
				}
			}


			bam_hdr_t *NewHeader; 

			NewHeader= sam_hdr_parse(str.l, str.s);
			NewHeader->l_text = str.l; NewHeader->text = str.s;
			NewHeader=sam_hdr_sanitise(NewHeader);
			tmp=sam_hdr_write(outR1, NewHeader);
			bam_hdr_destroy(NewHeader);

			while((sam_read1(in_SS, header_SS, aln_SS) >= 0))
			{
				if ((aln_SS->core).tid >= 0)
				{ // chr
					chrID=header_SS->target_name[(aln_SS->core).tid];

					MapIt=NewChangID.find(chrID);
					if (MapIt!=NewChangID.end())
					{
						(aln_SS->core).tid=MapIt->second;
						if ((aln_SS->core).mtid>=0)
						{				
							chrID=header_SS->target_name[(aln_SS->core).mtid];
							MapIt=NewChangID.find(chrID);
							if (MapIt!=NewChangID.end())
							{
								(aln_SS->core).mtid=MapIt->second;
							}
							else
							{
								(aln_SS->core).mtid=-1;
							}
						}
						tmp=bam_write1(outR1->fp.bgzf, aln_SS);
					}
				}
				else
				{
					if (para_Xam04->TF2)
					{
						if ((aln_SS->core).mtid>=0)
						{				
							chrID=header_SS->target_name[(aln_SS->core).mtid];
							MapIt=NewChangID.find(chrID);
							if (MapIt!=NewChangID.end())
							{
								(aln_SS->core).mtid=MapIt->second;
							}
							else
							{
								(aln_SS->core).mtid=-1;
							}
						}
						tmp=bam_write1(outR1->fp.bgzf, aln_SS);
					}
				}
			}
		}







	}

	sam_close(in_SS);

	sam_close(outR1);
	bam_destroy1(aln_SS);
	bam_hdr_destroy(header_SS);

	delete para_Xam04 ;
	return 0;
}
#endif
///////// swimming in the sky and flying in the sea ////////////
