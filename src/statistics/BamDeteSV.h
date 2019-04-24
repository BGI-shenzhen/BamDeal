#ifndef bamSV_H_
#define bamSV_H_
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



int  bamSV_help()
{
	cout <<""
		"\n"
		"\tUsage: bamSV(beta)  -InFile A.bam B.bam  -OutFile outFix\n"
		"\n"
		"\t\t-InFile     <str>     In Bam/Sam File to Dete SV[repeat]\n"
		"\t\t-OutFile    <str>     OutFile Prefix\n"
		"\n"
		"\t\t-InList     <str>     In Bam File List\n"
		"\t\t-Insert     <int>     set insert[auto]\n"
		"\t\t-SD         <int>     set insert SD[auto]\n"
		"\t\t-Ref        <str>     In Ref.fa to get gap info to filter SV\n"
		"\t\t-MinDepth   <int>     minimum count of PE in one cluster[6]\n"
		"\n"
		"\t\t-help                 Show more help [hewm2008]\n"
		"\n";
	return 1;
}


int bamSV_help01(int argc, char **argv , In3str1v * paraBamS )
{
	if (argc <=2 ) {bamSV_help();return 0;}
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

		if (flag  == "InList")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			string A=argv[i];
			file_count+=(ReadList(A ,(paraBamS->List)));
		}
		else if (flag  == "InFile")
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
		else if (flag  ==  "OutFile")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraBamS->InStr2=argv[i];			
		}
		else if (flag  ==  "Ref")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraBamS->InStr1=argv[i];			
		}
		else if (flag  ==  "Insert")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraBamS->InInt=atoi(argv[i]);
		}
		else if (flag  ==  "SD")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraBamS->InInt2=atoi(argv[i]);
		}
		else if (flag  ==  "MinDepth")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraBamS->InF=atof(argv[i]);
		}




		else if (flag == "help")
		{
			bamSV_help();return 0;
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


	return 1 ;
}



int bamSV_main(int argc, char *argv[])
//int main(int argc, char *argv[])
{
	In3str1v *paraBamS = new In3str1v;
	paraBamS->InF=6.0;
	if ((bamSV_help01(argc, argv, paraBamS)==0))
	{
		delete paraBamS ;
		return 0;
	}

	string strPE=(paraBamS->InStr2)+".OutInsert.pe.gz";
	ogzstream  OUTPE (strPE.c_str());

	string strSV=(paraBamS->InStr2)+".sv";
	ofstream  OUTSV (strSV.c_str());


	OUTSV<<"#Chr\tSV_type\tLength\tClusterStartA\tClusterStartB\tClusterEndA\tClusterEndB\tFR\tSupRead\n";

	struct clusterPE 
	{
		int Start_A;
		int End_A;
		int Start_B;
		int End_B;
		int chrID ;
		char Flag[2];
		int Number ;
	};

	vector <clusterPE> Vcluster_A;
	int  ZV_sizeA=0;

	bam_hdr_t *headerALL;
	samFile *InBamALL = sam_open((paraBamS->List)[0].c_str(), "r");
	headerALL = sam_hdr_read(InBamALL);
	sam_close(InBamALL);

	int FileNum=(paraBamS->List).size();

	int Final_Insert=0; 

	for (int ii=0 ; ii<FileNum ; ii++)
	{
		string  line ;
		line =(paraBamS->List)[ii];
		if (line.length()<=0)  { continue  ; }
		vector <string> TmpVV ;
		split (line , TmpVV ," \t");
		line=TmpVV[0];

		int InsertThis = (paraBamS->InInt) ;
		int SD=(paraBamS->InInt2);
		if (TmpVV.size()>1)
		{
			InsertThis=atoi(TmpVV[1].c_str());
		}
		if (TmpVV.size()>2)
		{
			SD=atoi(TmpVV[2].c_str());			
		}

		bam1_t *aln = bam_init1();

		if (InsertThis<10  ||  SD<5 )
		{
			cout <<"Begin Bam :"<<line<<endl;
			bam_hdr_t *header;

			samFile *InBam = sam_open(line.c_str(), "r");
			header = sam_hdr_read(InBam);


			int tmp=0;
			map <int ,int > map_insert ;
			map <int ,int > :: iterator mapIt ;
			int Count=0;
			while (sam_read1(InBam, header, aln) >= 0 &&  Count <88888)
			{
				if ( ((aln->core).tid < 0) || ( (aln->core).qual < 10)   ) {continue ;}
				if  ( (((aln)->core.flag & 0x40) ==0 ))
				{
					continue ;
				}
				if  ( (((aln)->core.flag & 0x1) == 0))
				{
					continue ;
				}
				if  ( (((aln)->core.flag & 0x8) != 0))
				{
					continue ;
				}
				if  ((aln->core).tid != (aln->core).mtid)
				{
					continue ;
				}
				if ( (aln->core).isize<0)
				{
					(aln->core).isize=0-(aln->core).isize;
				}
				mapIt=map_insert.find((aln->core).isize);
				if (mapIt!=map_insert.end())
				{
					(mapIt->second)++;
				}
				else
				{
					map_insert.insert( map <int,int>  :: value_type ((aln->core).isize,1));
				}
				Count++;
			}
			sam_close(InBam);
			bam_hdr_destroy(header);

			int  max_y = 0;
			int  max_x = 0;
			mapIt=map_insert.begin();

			while(mapIt !=map_insert.end())
			{
				if (max_y <= (mapIt->second))
				{
					max_y= mapIt->second;
					max_x= mapIt->first;
				}
				mapIt++;
			}

			int Num=0;
			int SumDiff=0;
			int cutoffX=(max_x*100); 
			if (cutoffX>500000) {cutoffX=500000;}

			for ( mapIt=map_insert.begin();  mapIt !=map_insert.end() ; mapIt++ )
			{
				if ( (mapIt->first)>cutoffX)
				{
					continue;
				}
				int Diff=(mapIt->first)-max_x;
				Num+=(mapIt->second);
				SumDiff=(mapIt->second)*Diff*Diff;
			}
			SD=int(sqrt(SumDiff/Num));
			if (InsertThis<20)
			{
				InsertThis=max_x;
			}
		}

		cout<<"Insert Size Peaks is Estimated: "<<InsertThis<<" , and SD: "<<SD<<endl;
		cout<<"PE Insert Size out off Peaks +- 4*SD will be detected"<<endl;


		int MinXX=InsertThis-4*SD;
		int MaxXX=InsertThis+4*SD;
		if (Final_Insert < InsertThis)
		{
			Final_Insert=InsertThis ;
		}


		bam_hdr_t *headerC;

		samFile *InBamC = sam_open(line.c_str(), "r");
		headerC = sam_hdr_read(InBamC);

		while (sam_read1(InBamC, headerC, aln) >= 0 )
		{
			if ( ((aln->core).tid < 0) || ( (aln->core).qual < 10)   ) {continue ;}
			if  ( (((aln)->core.flag & 0x40) ==0))
			{
				continue ;
			}
			if  ( (((aln)->core.flag & 0x1) == 0))
			{
				continue ;
			}
			if  ( (((aln)->core.flag & 0x8) != 0))
			{
				continue ;
			}
			if  ((aln->core).tid != (aln->core).mtid)
			{
				continue ;
			}
			int Insert=(aln->core).isize;
			if ( Insert<0)
			{
				Insert=0-Insert;
			}
			if  (Insert > MaxXX  ||  Insert < MinXX)
			{
				clusterPE  ThisPE;
				ThisPE.Number=1;
				ThisPE.Flag[0]='+';
				ThisPE.Flag[1]='+';
				ThisPE.Flag[2]='\0';
				ThisPE.chrID=(aln->core).tid;

				if  (((aln)->core.flag & 0x10) ) {ThisPE.Flag[0]='-';}
				if  (((aln)->core.flag & 0x20) ) {ThisPE.Flag[1]='-';}
				string readID=bam_get_qname(aln);

				if  ((aln->core).pos < (aln->core).mpos)
				{
					ThisPE.Start_B=(aln->core).mpos;
					ThisPE.End_B=(aln->core).mpos;
					ThisPE.Start_A=(aln->core).pos;
					ThisPE.End_A=(aln->core).pos;
					OUTPE<<(headerC->target_name[(aln->core).tid])<<"\t"<<(aln->core).pos<<"\t"<<(headerC->target_name[(aln->core).mtid])<<"\t"<<(aln->core).mpos<<"\t"<<Insert<<"\t"<<ThisPE.Flag<<"\t"<<readID<<"\n";
				}
				else
				{
					char ttt=ThisPE.Flag[0]; ThisPE.Flag[0]=ThisPE.Flag[1]; ThisPE.Flag[1]=ttt ;
					ThisPE.Start_A=(aln->core).mpos;
					ThisPE.End_A=(aln->core).mpos;
					ThisPE.Start_B=(aln->core).pos;
					ThisPE.End_B=(aln->core).pos;
					OUTPE<<(headerC->target_name[(aln->core).mtid])<<"\t"<<(aln->core).mpos<<"\t"<<(headerC->target_name[(aln->core).tid])<<"\t"<<(aln->core).pos<<"\t"<<Insert<<"\t"<<ThisPE.Flag<<"\t"<<readID<<"\n";
				}


				ZV_sizeA=Vcluster_A.size();
				bool TF_End=true;
				for (int ii=0 ; ii< ZV_sizeA ; ii++)
				{
					if (  (ThisPE.chrID==Vcluster_A[ii].chrID)  &&  (strcmp(Vcluster_A[ii].Flag ,ThisPE.Flag)==0 )   &&  abs(ThisPE.Start_A-Vcluster_A[ii].Start_A)<(Final_Insert) &&  abs(ThisPE.Start_B-Vcluster_A[ii].Start_B)<(Final_Insert) )
					{
						(Vcluster_A[ii].Number)+=(ThisPE.Number);
						if (ThisPE.Start_A<Vcluster_A[ii].Start_A) {Vcluster_A[ii].Start_A=ThisPE.Start_A;}
						if (ThisPE.Start_B<Vcluster_A[ii].Start_B) {Vcluster_A[ii].Start_B=ThisPE.Start_B;}

						if (ThisPE.End_A>Vcluster_A[ii].End_A) {Vcluster_A[ii].End_A=ThisPE.End_A;}
						if (ThisPE.End_B>Vcluster_A[ii].End_B) {Vcluster_A[ii].End_B=ThisPE.End_B;}
						TF_End=false ;
						break;
					}
				}
				if  (TF_End)
				{
					Vcluster_A.push_back(ThisPE);
				}
			}
		}

		bam_destroy1(aln);
		sam_close(InBamC);
		bam_hdr_destroy(headerC);

	}



	OUTPE.close();


	int MinDepthF = int (paraBamS->InF);
	for (int ii=0 ; ii< ZV_sizeA ; ii++)
	{
		if ( (Vcluster_A[ii].Number) < MinDepthF )
		{
			continue ;
		}
		int Len=Vcluster_A[ii].Start_B-Vcluster_A[ii].End_A;

		if (Vcluster_A[ii].Flag[0] == '+' &&   Vcluster_A[ii].Flag[1] == '-' )
		{
			if  (Len<0) { continue ;}
			else
			{
				if  (Len>Final_Insert)
				{
					OUTSV<<(headerALL->target_name[Vcluster_A[ii].chrID])<<"\tDeletion\t"<<Len<<"\t"<<Vcluster_A[ii].Start_A<<"\t"<<Vcluster_A[ii].End_A<<"\t"<<Vcluster_A[ii].Start_B<<"\t"<<Vcluster_A[ii].End_B<<"\t"<<Vcluster_A[ii].Flag<<"\t"<<(Vcluster_A[ii].Number)<<endl;
				}
				else
				{
					OUTSV<<(headerALL->target_name[Vcluster_A[ii].chrID])<<"\tInsertion\t"<<Len<<"\t"<<Vcluster_A[ii].Start_A<<"\t"<<Vcluster_A[ii].End_A<<"\t"<<Vcluster_A[ii].Start_B<<"\t"<<Vcluster_A[ii].End_B<<"\t"<<Vcluster_A[ii].Flag<<"\t"<<(Vcluster_A[ii].Number)<<endl;
				}
			}
		}
		else
		{
			OUTSV<<(headerALL->target_name[Vcluster_A[ii].chrID])<<"\tNANAtmp\t"<<Len<<"\t"<<Vcluster_A[ii].Start_A<<"\t"<<Vcluster_A[ii].End_A<<"\t"<<Vcluster_A[ii].Start_B<<"\t"<<Vcluster_A[ii].End_B<<"\t"<<Vcluster_A[ii].Flag<<"\t"<<(Vcluster_A[ii].Number)<<endl;
		}
	}


	/*
	   vector <clusterPE> Vcluster_B;

	   int  ZV_sizeB=0;
	   bool VecAB=true ;
	   bool TF_End=true;

	   if ( (paraBamS->InInt) < 250 )  { (paraBamS->InInt)=250;}

	   for ( int cc=3 ; cc< (paraBamS->InInt) ; cc++)
	   {
	   if  (VecAB)
	   {
	   Vcluster_B.clear();
	   ZV_sizeA=Vcluster_A.size();
	   VecAB=false;
	   for (int ii=0 ; ii< ZV_sizeA  ; ii++)
	   {
	   ZV_sizeB=Vcluster_B.size();
	   TF_End=true;
	   for (int jj=0 ; jj< ZV_sizeB  ; jj++)
	   {
	   if (( Vcluster_B[jj].chrID==Vcluster_A[ii].chrID)  &&  (strcmp(Vcluster_A[ii].Flag ,Vcluster_B[jj].Flag)==0 )   &&  (abs(Vcluster_B[jj].Start_A-Vcluster_A[ii].Start_A)<cc ) && ( abs(Vcluster_B[jj].Start_B-Vcluster_A[ii].Start_B)<cc ))
	   {
	   (Vcluster_B[jj].Number)+=(Vcluster_A[ii].Number);
	   if (Vcluster_B[jj].Start_A>Vcluster_A[ii].Start_A) {Vcluster_B[jj].Start_A=Vcluster_A[ii].Start_A;}
	   if (Vcluster_B[jj].Start_B>Vcluster_A[ii].Start_B) {Vcluster_B[jj].Start_B=Vcluster_A[ii].Start_B;}

	   if (Vcluster_B[jj].End_A<Vcluster_A[ii].End_A) {Vcluster_B[jj].End_A=Vcluster_A[ii].End_A;}
	   if (Vcluster_B[jj].End_B<Vcluster_A[ii].End_B) {Vcluster_B[jj].End_B=Vcluster_A[ii].End_B;}
	   TF_End=false;
	   break;
	   }
	   }
	   if  (TF_End)
	   {
	   Vcluster_B.push_back(Vcluster_A[ii]);
	   }
	   }
	   }
	   else
	   {



	   Vcluster_A.clear();
	   ZV_sizeB=Vcluster_B.size();
	   VecAB=true;
	   for (int jj=0 ; jj< ZV_sizeB; jj++)
	   {
	   TF_End=true;
	   ZV_sizeA=Vcluster_A.size();
	   for (int ii=0 ; ii< ZV_sizeA ; ii++)
	   {
	   if ((Vcluster_B[jj].chrID==Vcluster_A[ii].chrID)  &&  (strcmp(Vcluster_A[ii].Flag ,Vcluster_B[jj].Flag)==0 )   &&  abs(Vcluster_B[jj].Start_A-Vcluster_A[ii].Start_A)<cc &&  abs(Vcluster_B[jj].Start_B-Vcluster_A[ii].Start_B)<cc)
	   {
	   (Vcluster_A[ii].Number)+=(Vcluster_B[jj].Number);
	   if (Vcluster_B[jj].Start_A<Vcluster_A[ii].Start_A) {Vcluster_A[ii].Start_A=Vcluster_B[jj].Start_A;}
	   if (Vcluster_B[jj].Start_B<Vcluster_A[ii].Start_B) {Vcluster_A[ii].Start_B=Vcluster_B[jj].Start_B;}

	   if (Vcluster_B[jj].End_A>Vcluster_A[ii].End_A) {Vcluster_A[ii].End_A=Vcluster_B[jj].End_A;}
	   if (Vcluster_B[jj].End_B>Vcluster_A[ii].End_B) {Vcluster_A[ii].End_B=Vcluster_B[jj].End_B;}
	   TF_End=false;
	   break;
	   }
	   }
	   if  (TF_End)
	   {
	   Vcluster_A.push_back(Vcluster_B[jj]);
	   }
	   }

}
}




int MinDepthF = int (paraBamS->InF);
if  (VecAB)
{
	ZV_sizeA=Vcluster_A.size();
	for (int ii=0 ; ii< ZV_sizeA ; ii++)
	{
		if ( (Vcluster_A[ii].Number) < MinDepthF )
		{
			continue ;
		}
		int Len=Vcluster_A[ii].Start_B-Vcluster_A[ii].End_A;

		if (Vcluster_A[ii].Flag[0] == '+' &&   Vcluster_A[ii].Flag[1] == '-' )
		{
			if  (Len<0) { continue ;}
			else
			{
				if  (Len>Final_Insert)
				{
					OUTSV<<(headerALL->target_name[Vcluster_A[ii].chrID])<<"\tDeletion\t"<<Len<<"\t"<<Vcluster_A[ii].Start_A<<"\t"<<Vcluster_A[ii].End_A<<"\t"<<Vcluster_A[ii].Start_B<<"\t"<<Vcluster_A[ii].End_B<<"\t"<<Vcluster_A[ii].Flag<<"\t"<<(Vcluster_A[ii].Number)<<endl;
				}
				else
				{
					OUTSV<<(headerALL->target_name[Vcluster_A[ii].chrID])<<"\tInsertion\t"<<Len<<"\t"<<Vcluster_A[ii].Start_A<<"\t"<<Vcluster_A[ii].End_A<<"\t"<<Vcluster_A[ii].Start_B<<"\t"<<Vcluster_A[ii].End_B<<"\t"<<Vcluster_A[ii].Flag<<"\t"<<(Vcluster_A[ii].Number)<<endl;
				}
			}
		}
		else
		{
			OUTSV<<(headerALL->target_name[Vcluster_A[ii].chrID])<<"\tNANAtmp\t"<<Len<<"\t"<<Vcluster_A[ii].Start_A<<"\t"<<Vcluster_A[ii].End_A<<"\t"<<Vcluster_A[ii].Start_B<<"\t"<<Vcluster_A[ii].End_B<<"\t"<<Vcluster_A[ii].Flag<<"\t"<<(Vcluster_A[ii].Number)<<endl;
		}
	}
}
else
{
	ZV_sizeB=Vcluster_B.size();
	for (int jj=0 ; jj< ZV_sizeB ; jj++)
	{
		if ( (Vcluster_B[jj].Number) < MinDepthF )
		{
			continue ;
		}
		int Len=Vcluster_B[jj].Start_B-Vcluster_B[jj].End_A;
		if (Vcluster_B[jj].Flag[0] == '+' &&   Vcluster_B[jj].Flag[1] == '-' )
		{
			if  (Len<0) { continue ;}
			else
			{
				if  (Len>Final_Insert)
				{
					OUTSV<<(headerALL->target_name[Vcluster_B[jj].chrID])<<"\tDeletion\t"<<Len<<"\t"<<Vcluster_B[jj].Start_A<<"\t"<<Vcluster_B[jj].End_A<<"\t"<<Vcluster_B[jj].Start_B<<"\t"<<Vcluster_B[jj].End_B<<"\t"<<Vcluster_B[jj].Flag<<"\t"<<(Vcluster_B[jj].Number)<<endl;
				}
				else
				{
					OUTSV<<(headerALL->target_name[Vcluster_B[jj].chrID])<<"\tInsertion\t"<<Len<<"\t"<<Vcluster_B[jj].Start_A<<"\t"<<Vcluster_B[jj].End_A<<"\t"<<Vcluster_B[jj].Start_B<<"\t"<<Vcluster_B[jj].End_B<<"\t"<<Vcluster_B[jj].Flag<<"\t"<<(Vcluster_B[jj].Number)<<endl;
				}
			}
		}
		else
		{
			OUTSV<<(headerALL->target_name[Vcluster_B[jj].chrID])<<"\tNAtmp\t"<<Len<<"\t"<<Vcluster_B[jj].Start_A<<"\t"<<Vcluster_B[jj].End_A<<"\t"<<Vcluster_B[jj].Start_B<<"\t"<<Vcluster_B[jj].End_B<<"\t"<<Vcluster_B[jj].Flag<<"\t"<<(Vcluster_B[jj].Number)<<endl;
		}
	}
}


*/



OUTSV.close();
bam_hdr_destroy(headerALL);
delete paraBamS ;
return 0;


}

#endif // bamSV_H_  //


///////// swimming in the sky and flying in the sea ////////////





