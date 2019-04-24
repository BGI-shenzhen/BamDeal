#ifndef BamCat_H_
#define BamCat_H_
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

void BamCat_help()
{
	cout<<""
		"\n"
		"\tUsage: bamCat  -l  <sort.bam.list>  -o <out.sort.bam> -s\n"
		"\tUsage: bamCat  -i  A.bam B.bam  -o C.bam\n"
		"\n"
		"\t\t-i    <str>     input SAM/BAM/CRAM files(s)\n"
		"\t\t-l    <str>     input list of SAM/BAM/CRAM files\n"
		"\t\t-o    <str>     output BAM file\n"
		"\n"
		"\t\t-s              output sort bam file when all inputs were sorted\n"
		"\n"
		"\t\t-h              show more details for help\n"		
		"\n";
}

void More_Help_BamCat()
{
	cout<<""
		"\n"
		"\t\t1. bamCat -i A.bam B.bam -o AAA\n"
		"\t\t   bamCat -i A.sort.bam B.sort.sam -o AAA -s\n"
		"\t\t   bamCat -l bam.list -o AAA\n"
		"\t\t   bamCat -l bam.sort.list -o AAA -s\n"
		"\t\t   This will capture the input files into one file and output the result to a sorted BAM file named AAA in current directory. \n"
		"\t\t    (1.1)   -i input files could be all BAM, all SAM or the mix of two.\n"
		"\t\t    (1.2)   -l lists the input files. For example, if user has two input files A.bam and B.bam, bam.list should be formatted as:\n"
		"\t\t              ./A.bam\n"
		"\t\t              ./B.bam\n"
		"\t\t    (1.3)   Input files could have different headers. If the input files have the same chromosome with different length, warning would be given and the longest length would be kept in the header for that chromosome.\n"
		"\t\t    (1.4)  -s only works when the input files are sorted. Otherwise, the input files would be simply captured together.\n"
		"\n";
}

int BamCat_help01(int argc, char **argv , In3str1v * paraBamS )
{
	if (argc <2 ) {BamCat_help();return 0;}
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

		if (flag  == "InList"  ||  flag  == "l" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			string A=argv[i];
			file_count+=(ReadList(A ,(paraBamS->List)));
		}
		else if (flag  == "InFile"  ||  flag  == "i")
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
		else if (flag  ==  "OutFile" ||  flag  == "o" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraBamS->InStr2=argv[i];			
		}
		else if (flag  ==  "Merge"  ||  flag  == "Sorted" ||  flag  == "s")
		{
			paraBamS->TF=false;
		}
		else if (flag == "help" ||  flag  == "h")
		{
			More_Help_BamCat();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  (  (file_count<1) ||  (paraBamS->InStr2).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}

	string path=(paraBamS->InStr2);
	string ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);
	if (ext != "bam")
	{
		(paraBamS->InStr2)=path+".bam" ;
	}

	return 1 ;
}


void  run_Aread (  samFile *InBam , int & position , int & chrID,  htsFile *OUT ,  bam1_t *aln , bam_hdr_t *headerA  )
{
	int tmp=bam_write1(OUT->fp.bgzf, aln);

	for (int i =1 ; i>0 ; i++)
	{
		if (sam_read1(InBam, headerA, aln) >= 0)
		{
			if ((aln->core).pos<=position && (aln->core).tid == chrID )
			{
				tmp=bam_write1(OUT->fp.bgzf, aln);
			}
			else
			{				
				position=(aln->core).pos ;
				chrID=(aln->core).tid;
				break ;
			}
		}
		else
		{
			position=-8;
			chrID=-888;
			break ;
		}
	}
}



void  run_Aread (  samFile *InBam , int & position , int & chrID,  htsFile *OUT ,  bam1_t *aln , bam_hdr_t *headerA , map <string , int > & ChrInfo )
{
	string chrTmp ;
	map <string , int > :: iterator MapIt ;

	int tmp=bam_write1(OUT->fp.bgzf, aln);

	for (int i =1 ; i>0 ; i++)
	{
		if (sam_read1(InBam, headerA, aln) >= 0)
		{
			if ((aln->core).tid>=0)
			{
				chrTmp=(headerA->target_name[(aln->core).tid]);
				MapIt=ChrInfo.find(chrTmp);
				(aln->core).tid=MapIt->second;
			}

			if ( (aln->core).mtid >=0 )
			{
				chrTmp=((headerA)->target_name[(aln->core).mtid]);
				MapIt=ChrInfo.find(chrTmp);
				if (MapIt!=ChrInfo.end())
				{
					(aln->core).mtid=MapIt->second;
				}
				else
				{
					(aln->core).mtid=-1;
				}
			}

			if ((aln->core).pos<=position && (aln->core).tid == chrID )
			{
				tmp=bam_write1(OUT->fp.bgzf, aln);
			}
			else
			{				
				position=(aln->core).pos ;
				chrID=(aln->core).tid;
				break ;
			}
		}
		else
		{
			position=-8;
			chrID=-888;
			break ;
		}
	}
}


int chose_Afile ( const int Posi[] , const int ChrID[] ,  const int count )
{

	int IDD =-8888;
	int i=0;

	for (   ; i <count  ; i++)
	{
		if (ChrID[i]!=-888)
		{
			IDD=ChrID[i];
			break;
		}
	}

	for (   ; i <count  ; i++)
	{
		if (ChrID[i]!=-888)
		{
			if (IDD> ChrID[i])
			{
				IDD=ChrID[i];
			}
		}
	}

	int min_Aflag=-1;
	int start =-2 ;
	for ( i=0 ; i <count  ; i++)
	{
		if  ( ChrID[i]!=IDD )
		{
			continue ;
		}
		if (  start <0 )
		{
			min_Aflag=i ;
			start=Posi[i];
		}
		else if ( Posi[i]>-1 &&  Posi[i]<start)
		{
			min_Aflag=i ;
			start=Posi[i];
		}
	}

	if  (Posi[min_Aflag]<-7)
	{
		min_Aflag=-2;
	}
	return  min_Aflag ;
}




bool ALL_same_header( vector <string> List  )
{
	string  BamPath=List[0];
	bam_hdr_t *header;
	samFile *BamIn = hts_open(BamPath.c_str(), "r");
	header = sam_hdr_read(BamIn);
	sam_close(BamIn);

	int FileNum=List.size();
	for (int ii=1 ; ii< FileNum ; ii++)
	{
		string  line =List[ii];
		bam_hdr_t *headerA;
		samFile *InBam = sam_open(line.c_str(), "r");
		headerA = sam_hdr_read(InBam);
		sam_close(InBam);

		if ((header->n_targets)!=(headerA->n_targets))
		{
			bam_hdr_destroy(header);
			bam_hdr_destroy(headerA);
			return false;
		}
		bool NoSameBam=false;
		for(int i = 0; i < (header->n_targets); i++)
		{
			if ((strcmp(header->target_name[i],headerA->target_name[i])!=0)  ||  ((header->target_len[i])!=(headerA->target_len[i])))
			{
				NoSameBam=true;
				break ;
			}
		}
		if (NoSameBam)
		{
			bam_hdr_destroy(header);
			bam_hdr_destroy(headerA);
			return false;
		}
		bam_hdr_destroy(headerA);
	}

	bam_hdr_destroy(header);
	return true ;

}




void bamCat_WithDiffHead(In3str1v *paraBamS )
{

	htsFile *OUT;
	OUT=hts_open((paraBamS->InStr2).c_str(), "wb");

	map <string , int > ChrInfo ;
	map <string , int > :: iterator MapIt ;


	int FileNum=(paraBamS->List).size();


	for (int ii=0 ; ii< FileNum ; ii++)
	{
		string  line =(paraBamS->List)[ii];
		if (line.length()<=0)  { continue  ; }
		bam_hdr_t *headerA;
		samFile *InBam = sam_open(line.c_str(), "r");
		headerA = sam_hdr_read(InBam);
		bool NoSameBam=false;
		for(int i = 0; i < (headerA->n_targets); i++)
		{
			string chr=headerA->target_name[i];
			int TR=headerA->target_len[i];
			MapIt=ChrInfo.find(chr);
			if (MapIt==ChrInfo.end())
			{
				ChrInfo.insert( map <string ,int >  :: value_type (chr,TR));
			}
			else
			{
				if ((MapIt->second)==TR)
				{
					continue ;
				}
				cerr << "Warning! chr name : "<<chr<<" repeat muti-time with diff length: "<<(MapIt->second)<<" and "<<TR<<endl;
				if ((MapIt->second)<TR)
				{
					MapIt->second=TR;
				}
			}
		}
		bam_hdr_destroy(headerA);
		sam_close(InBam);
	}


	kstring_t str = { 0, 0, NULL };
	//	map <string , int > Chr2Int ;
	int Flag=0;

	for ( MapIt=ChrInfo.begin(); MapIt!=ChrInfo.end(); MapIt++)
	{
		kputsn("@SQ\tSN:", 7, &str);
		kputs((MapIt->first).c_str(), &str);
		kputsn("\tLN:", 4, &str);
		kputw(MapIt->second, &str);
		kputc('\n', &str);
		MapIt->second= Flag ;
		Flag++;
	}

	bam_hdr_t *NewHeader;
	NewHeader= sam_hdr_parse(str.l, str.s);
	NewHeader->l_text = str.l; NewHeader->text = str.s;
	NewHeader=sam_hdr_sanitise(NewHeader);

	int tmp=sam_hdr_write(OUT, NewHeader);
	string chrTmp ;
	if (paraBamS->TF)
	{
		for (int ii=0 ; ii< FileNum ; ii++)
		{
			string  line =(paraBamS->List)[ii];
			if (line.length()<=0)  { continue ;}
			cout <<"Begin Bam :"<<line<<endl;
			bam_hdr_t *headerA;
			bam1_t *aln = bam_init1();

			samFile *InBam = sam_open(line.c_str(), "r");
			headerA = sam_hdr_read(InBam);

			while (sam_read1(InBam, headerA, aln) >= 0)
			{
				if (  (aln->core).tid>=0 )
				{
					chrTmp=(headerA->target_name[(aln->core).tid]);
					MapIt=ChrInfo.find(chrTmp);
					(aln->core).tid=MapIt->second;
				}
				if ( (aln->core).mtid >=0 )
				{
					chrTmp=((headerA)->target_name[(aln->core).mtid]);
					MapIt=ChrInfo.find(chrTmp);
					if (MapIt!=ChrInfo.end())
					{
						(aln->core).mtid=MapIt->second;
					}
					else
					{
						(aln->core).mtid=-1;
					}

				}

				tmp=bam_write1(OUT->fp.bgzf, aln);
			}
			sam_close(InBam);
			bam_destroy1(aln);
			bam_hdr_destroy(headerA);
		}
	}
	else
	{

		samFile  **InBam = new samFile *[FileNum];
		bam_hdr_t **headerA =new bam_hdr_t *[FileNum];
		bam1_t  **aln= new bam1_t *[FileNum];
		int *Posi=new int[FileNum] ;
		int *ChrID=new int[FileNum] ;

		for (int ii=0 ; ii< FileNum ; ii++)
		{
			aln[ii]=bam_init1();
			InBam[ii]=sam_open((paraBamS->List)[ii].c_str(), "r");
			headerA[ii]=sam_hdr_read(InBam[ii]);
			if (sam_read1(InBam[ii], headerA[ii], aln[ii]) >= 0)
			{
				if ( (aln[ii]->core).tid >=0 )
				{
					chrTmp=((headerA[ii])->target_name[(aln[ii]->core).tid]);
					MapIt=ChrInfo.find(chrTmp);
					(aln[ii]->core).tid=MapIt->second;
				}

				if (  (aln[ii]->core).mtid >=0 )
				{
					chrTmp=((headerA[ii])->target_name[(aln[ii]->core).mtid]);
					MapIt=ChrInfo.find(chrTmp);

					if (MapIt!=ChrInfo.end())
					{
						(aln[ii]->core).mtid=MapIt->second;
					}
					else
					{
						(aln[ii]->core).mtid=-1;
					}

				}

				ChrID[ii]=(aln[ii]->core).tid;
				Posi[ii]=(aln[ii]->core).pos;
			}
			else
			{
				cerr<<"can Not read bam file : "<<(paraBamS->List)[ii]<<endl;
				Posi[ii]=-2;
				ChrID[ii]=-888;
			}
		}

		int file_Arun_Acout=chose_Afile(Posi ,ChrID , FileNum) ;

		while( file_Arun_Acout > -1 )
		{
			run_Aread (InBam[file_Arun_Acout] , Posi[file_Arun_Acout] , ChrID[file_Arun_Acout], OUT , aln[file_Arun_Acout] ,headerA[file_Arun_Acout] , ChrInfo );
			file_Arun_Acout=chose_Afile(Posi ,ChrID , FileNum) ;
		}

		for (int ii=0 ; ii< FileNum ; ii++)
		{
			sam_close(InBam[ii]);
			bam_destroy1(aln[ii]);
			bam_hdr_destroy(headerA[ii]);

		}


		delete [] Posi;
		delete [] ChrID;
		delete [] InBam ;
		delete [] aln ;
		delete [] headerA ;

	}


	sam_close(OUT);
	bam_hdr_destroy(NewHeader);

}








void bamCat_WithSameHead(In3str1v *paraBamS )
{

	htsFile *OUT;
	OUT=hts_open((paraBamS->InStr2).c_str(), "wb");
	bam_hdr_t *NewHeader;


	string  BamPath=(paraBamS->List)[0];
	samFile *BamIn = hts_open(BamPath.c_str(), "r");
	NewHeader = sam_hdr_read(BamIn);
	sam_close(BamIn);

	int tmp=sam_hdr_write(OUT, NewHeader);
	string chrTmp ;
	int FileNum=(paraBamS->List).size();

	if (paraBamS->TF)
	{
		for (int ii=0 ; ii< FileNum ; ii++)
		{
			string  line =(paraBamS->List)[ii];
			if (line.length()<=0)  { continue ;}
			cout <<"Begin Bam :"<<line<<endl;
			bam_hdr_t *headerA;
			bam1_t *aln = bam_init1();

			samFile *InBam = sam_open(line.c_str(), "r");
			headerA = sam_hdr_read(InBam);
			while (sam_read1(InBam, headerA, aln) >= 0)
			{
				tmp=bam_write1(OUT->fp.bgzf, aln);
			}
			sam_close(InBam);
			bam_destroy1(aln);
			bam_hdr_destroy(headerA);
		}
	}
	else
	{

		samFile  **InBam = new samFile *[FileNum];
		bam_hdr_t **headerA =new bam_hdr_t *[FileNum];
		bam1_t  **aln= new bam1_t *[FileNum];
		int *Posi=new int[FileNum];
		int *ChrID=new int[FileNum];

		for (int ii=0 ; ii< FileNum ; ii++)
		{
			aln[ii]=bam_init1();
			InBam[ii]=sam_open((paraBamS->List)[ii].c_str(), "r");
			headerA[ii]=sam_hdr_read(InBam[ii]);
			if (sam_read1(InBam[ii], headerA[ii], aln[ii]) >= 0)
			{
				ChrID[ii]=(aln[ii]->core).tid;
				Posi[ii]=(aln[ii]->core).pos;
			}
			else
			{
				cerr<<"can Not read bam file : "<<(paraBamS->List)[ii]<<endl;
				Posi[ii]=-2;
				ChrID[ii]=-888;
			}
		}

		int file_Arun_Acout=chose_Afile(Posi ,ChrID , FileNum) ;

		while( file_Arun_Acout > -1 )
		{
			run_Aread (InBam[file_Arun_Acout] , Posi[file_Arun_Acout] , ChrID[file_Arun_Acout], OUT , aln[file_Arun_Acout] ,headerA[file_Arun_Acout]);
			file_Arun_Acout=chose_Afile(Posi ,ChrID , FileNum) ;
		}

		for (int ii=0 ; ii< FileNum ; ii++)
		{
			sam_close(InBam[ii]);
			bam_destroy1(aln[ii]);
			bam_hdr_destroy(headerA[ii]);
		}


		delete [] Posi;
		delete [] ChrID;
		delete [] InBam ;
		delete [] aln ;
		delete [] headerA ;

	}



	sam_close(OUT);
	bam_hdr_destroy(NewHeader);



}


int bamCat_main(int argc, char *argv[])
{
	In3str1v *paraBamS = new In3str1v;

	if ((BamCat_help01(argc, argv, paraBamS)==0))
	{
		delete paraBamS ;
		return 0;
	}

	bool Samhead=ALL_same_header(paraBamS->List);


	if  (Samhead)
	{
		bamCat_WithSameHead( paraBamS ) ;
	}
	else
	{
		bamCat_WithDiffHead( paraBamS ) ;
	}



	delete paraBamS ;
	return 0;
}
#endif // BamCat_H_  //
///////// swimming in the sky and flying in the sea ////////////





