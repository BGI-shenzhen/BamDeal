#ifndef bamDepGC_H_
#define bamDepGC_H_
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
int  bamDeGCshow_help()
{
	cout <<""
		"\n"
		"\tUsage: DepthGC  -InList  <bam.list>  -Ref  <Ref.fa> -OutPut  <outPrefix>\n"
		"\tUsage: DepthGC  -InWigF  <In.DepthGC.gz>   -OutPut  <outPrefix>\n"
		"\n"
		"\t\t-InList    <str>     Input Bam/Sam File List\n"
		"\t\t-InFile    <str>     Input Bam/Sam File File[repeat]\n"
		"\t\t-InWigF    <str>     Input Depth-GC wig File,BamCoverage OutPut\n"
		"\t\t-OutPut    <str>     OutPut Prefix File for Depth-GC Plot pdf\n"
		"\n"
		"\t\t-Ref       <str>     In Ref.fa Want Out Depth-GC wig info\n"
		"\t\t-Windows   <int>     Windows size for Depth-GC wig[10000]\n"
		"\t\t-MinQ      <int>     Ignore too low mapQ read[10]\n"
		"\t\t-MaxYaxis  <int>     Max Y axis to plot the Result [3*meanD]\n"
		"\t\t-keepR               Keep the Rscript to RePlot the Figure\n"
		"\n"
		"\t\t-help                Show this help [hewm2008 v1.20]\n"
		"\n";
	return 1;
}
*/

void bamDeGCshow_help()
{
cout<<""
"\n"
"\tUsage: DepthGC  -l  <bam.list>  -r  <Ref.fa> -o  <outPrefix>\n"
"\tUsage: DepthGC  -f  <DepthGC.wig.gz>   -o  <outPrefix>\n"
"\n"
"\t\t-i      <str>     input SAM/BAM files, delimited by space\n"
"\t\t-l      <str>     input list of SAM/BAM files\n"
"\t\t-o      <str>     prefix of output file\n"
"\n"
"\t\t-r      <str>     input reference FASTA\n"
"\t\t-f      <str>     file containing depth and GC content in each window. This file is one of the output files of Bamdeal statistics Coverage.\n"
"\t\t-w      <int>     window size to calculate base frequency, default [10000]\n"
"\t\t-q      <int>     reads with quality lower than this will be filtered, default [10]\n"
"\t\t-y      <int>     maximum of y axis of the plot, default [3*mean of depth]\n"
"\t\t-k                output the Rscript used to generate plots\n"
"\n"
"\t\t-h                show more details for help\n"
"\n"
"\n";
}


void More_bamDeGCshow_help()
{
cout<<""
"\n"
"\t\t1. DepthGC -i <A.bam B.bam> -r <Ref.fa> -o AAA -k -q Q\n"
"\t\tDepthGC -l <bam.list> -r <Ref.fa> -o AAA -k -q Q\n"
"\t\tThis will generate the plot of depth v.s GC content by window. The output file will be named with the prefix AAA and output to current directory.\n"
"\t\t(1.1) -l lists the input files. For example, if user has two input files A.bam and B.bam, bam.list should be formatted as:\n"
"\t\t./A.bam\n"
"\t\t./B.bam\n"
"\t\t(1.2) the reads with quality lower than Q will be removed, default value of Q is 10.\n"
"\n"
"\t\t2. DepthGC -f <DepthGC.wig.gz> -o AAA -k -q Q\n"
"\t\tThis operation is the same with the example above but with different input format.\n"
"\t\t(2.1) the input file <DepthGC.wig.gz> shows depth and GC content by window. User could get this file with the function Coverage in the Statistics module. Below is an example of the usage of this command. More details could be found with the -h in this command.\n"
"\t\tbamdeal statistics Coverage -i <in.bam> -o AAA -r <Ref.fa>\n"
"\n"
"\n";
}


int bamDeGCshow_help01(int argc, char **argv , In3str1v * paraFA04 )
{
	if (argc <=1 ) {bamDeGCshow_help();return 0;}
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

		if (flag  == "InList" ||  flag  == "List"   || flag == "l")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			string A=argv[i];
			file_count+=(ReadList(A ,(paraFA04->List)));
		}
		else if (flag  == "InFile"   || flag == "i")
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
		else if (flag  ==  "InWigF"   || flag == "f")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			(paraFA04->InStr1)=argv[i];
		}	
		else if (flag  ==  "OutPut"   || flag == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr2=argv[i];
		}
		else if (flag  ==  "Ref"   || flag == "r")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr3=argv[i];
		}
		else if (flag  ==  "Windows"  || flag == "w")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InF=atof(argv[i]);
		}
		else if (flag  ==  "MinQ"   || flag == "q")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InInt=atoi(argv[i]);
		}
		else if (flag  ==  "MaxYaxis"   || flag == "y")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InInt2=atoi(argv[i]);
		}
		else if (flag  ==  "keepR"   || flag == "k" )
		{
			paraFA04->TF=false ;
		}
		else if (flag == "help"  || flag == "h" )
		{
		   More_bamDeGCshow_help();return 0;
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
	
	if (!(paraFA04->InStr1).empty())
	{
		return 1 ;
	}

	if ( (file_count<1)   ||  (paraFA04->InStr3).empty())
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}

	//	(paraFA04->InStr2)=add_Asuffix(paraFA04->InStr2);
	return 1 ;
}


void DrawDepthGC( string File , In3str1v * paraFA04 )
{

	if  ((paraFA04->InInt2)<2)
	{
		igzstream LLT (File.c_str(),ifstream::in); // igzstream

		if (!LLT.good())
		{
			cerr << "open List error: "<<File<<endl;
		}

		string name ;
		double Depth ;
		double GC ;
		double SumDepth=0;
		double Count =0;
		while(!LLT.eof())
		{
			string  line ;
			getline(LLT,line);
			if (line.length()<=0)  { continue  ; }
			if (line[0] ==  '#')  { continue  ; }
			if (line[line.length()-1] ==  'A')  { continue  ; }
			istringstream isone (line,istringstream::in);
			isone>>name>>Depth;
			SumDepth+=Depth;
			Count++;
		}
		LLT.close();
		(paraFA04->InInt2)=(int(SumDepth/Count)+1)*3;
	}



	igzstream LIST (File.c_str(),ifstream::in); // igzstream
	if (!LIST.good())
	{
		cerr << "open List error: "<<File<<endl;
	}



	string name ;
	double Depth ;
	double GC ;
	map <double ,map <double , int> > hash_map;

	while(!LIST.eof())
	{
		string  line ;
		getline(LIST,line);
		if (line.length()<=0)  { continue  ; }
		if (line[0] ==  '#')  { continue  ; }
		if (line[line.length()-1] ==  'A')  { continue  ; }
		istringstream isone (line,istringstream::in);
		isone>>name>>Depth>>GC;
		if (Depth>(paraFA04->InInt2))
		{
			Depth=(paraFA04->InInt2);
		}
		hash_map[GC][Depth]++;
	}
	LIST.close();

	string OutPlotTmp=(paraFA04->InStr2)+".tmp";
	ofstream  OUTPLOT (OutPlotTmp.c_str());

	OUTPLOT<<"GC_ratio\tDepth\tNumbers\n";
	map <double,map<double, int> > ::const_iterator outerit=hash_map.begin();
	map<double,int > ::const_iterator innerit ;


	for(outerit=hash_map.begin() ; outerit!=hash_map.end(); outerit++)
	{
		for(innerit=outerit->second.begin() ; innerit!=outerit->second.end();  innerit++)
		{
			OUTPLOT<<outerit->first<<"\t"<<innerit->first<<"\t"<<innerit->second<<endl;
		}
	}
	OUTPLOT.close();

	string OutPlotr=(paraFA04->InStr2)+".r";
	ofstream  OUTR (OutPlotr.c_str());

	OUTR<<""
		"\n"
		"library(ggplot2);\nlibrary(gridExtra);\nlibrary(ggExtra)\n" 
		"read.table(\""<<OutPlotTmp<<"\",header=T)->tab;\n"
		"GC_ratio<-tab[,1];\n"
		"Depth<-tab[,2];\n"
		"Numbers<-tab[,3];\n"
		"df <- data.frame(GC_ratio , Depth );\n"
		"pdf(\""<<(paraFA04->InStr2)<<".pdf\");\n"
		"p <- ggplot(df, aes(GC_ratio, Depth,colour=Numbers)) + geom_point(size=1.5) +theme(legend.direction=\"horizontal\",legend.position= \"bottom\",legend.title = element_text(face = \"bold\",size=7));\n"
		"ggExtra::ggMarginal(p, type = \"histogram\");\n"
		"dev.off();\n"
		"png(\""<<(paraFA04->InStr2)<<".png\");\n"
		"ggExtra::ggMarginal(p, type = \"histogram\");\n"
		"\ndev.off();"<<endl;
	OUTR.close();


	char   buf[2048]={'\0'};
	string cc="which  Rscript  2> /dev/null ";
	memset( buf, '\0', sizeof(buf) );
	FILE   *stream ;
	stream=popen(cc.c_str(),"r") ;
	fread( buf, sizeof(char), sizeof(buf), stream);
	string binPath=buf;
	binPath=binPath.substr(0,binPath.length()-1);
	if (binPath == "" )
	{
		cout <<"\twarning: can't find the [Rscript] in your $PATH ; no png Figure Out"<<endl;
		cout <<"\t\tRscript "<<OutPlotr<<endl;
	}
	else
	{
		if (paraFA04->TF)
		{
			cc=binPath+"\t"+OutPlotr+"  ; rm  -rf  "+ OutPlotr +"  "+ OutPlotTmp ;
			std::system(cc.c_str()) ;
		}
		else
		{
			cc=binPath+"\t"+OutPlotr ;
			std::system(cc.c_str()) ;
			cout <<"\t\tRePlot : Rscript  "<<OutPlotr<<endl;
		}
	}

}


int bamDepthGC_main(int argc, char *argv[])
//int main(int argc, char *argv[])
{
	In3str1v *paraFA04 = new In3str1v;
	paraFA04->InInt=10;
	paraFA04->InF=10000;
	if ((bamDeGCshow_help01(argc, argv, paraFA04)==0))
	{
		delete paraFA04 ;
		return 0 ;
	}

	if (!(paraFA04->InStr1).empty())
	{
		DrawDepthGC( paraFA04->InStr1 ,  paraFA04 ) ;
		delete paraFA04 ;
		return 1 ;
	}

	string GCDepth=(paraFA04->InStr2);
	GCDepth=GCDepth+".DepthGC.wig.gz";
	ogzstream  OUTGC (GCDepth.c_str());

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
		while (sam_read1(InBam, header, aln)>= 0)
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
				if  (cig==BAM_CMATCH || cig==BAM_CEQUAL || cig==BAM_CDIFF)
				{
					int end=This+ncig;
					for (  ; This<end;This++ )
					{
						depth[(aln->core).tid][This]++;
					}
				}
				else if  (cig==BAM_CDEL ||   cig==BAM_CREF_SKIP )
				{
					This=This+ncig;
				}
			}
			/*///
			  for (int32_t j=0 ; j< ((aln->core).l_qseq) ; j++)
			  {
			  int32_t  This=((aln->core).pos)+j;
			  depth[(aln->core).tid][This]++;
			  }
			  *////
		}
		sam_close(InBam);
		bam_destroy1(aln);
		bam_hdr_destroy(headerA);
	}

	cout <<"ALL Bam Read done"<<endl;

	OUTGC<<"##Bigin\tDepth(X)\tGC(%)\n";

	gzFile fp;
	kseq_t *seq;
	int l;
	fp = gzopen((paraFA04->InStr3).c_str(), "r");
	seq = kseq_init(fp);
	int Windows=int(paraFA04->InF);
	int cutOffLen=int(Windows*0.4);

	while ((l = kseq_read(seq)) >= 0)
	{
		int  NumIntSca=0;
		for(NumIntSca=0; NumIntSca <(header->n_targets); NumIntSca++)
		{
			if ( strcmp(header->target_name[NumIntSca],(seq->name.s)) ==0)
			{
				break;
			}
		}
		int  LastDis=(header->target_len[NumIntSca])-(Windows);
		if (LastDis<0) {continue;}
		OUTGC<<"#>"<<(seq->name.s)<<"\tWindowsSite:\t"<<(Windows)<<endl;

		for (int j =0 ; j< LastDis ; j+=(Windows))
		{
			int  Ascii[256] = {0};
			int TotalDepth=0;
			for (int kk=0 ; kk<(Windows) ; kk++)
			{
				int site=j+kk;
				Ascii[seq->seq.s[site]]++;
				TotalDepth+=depth[NumIntSca][site];
			}
			int  N_Aleng= Ascii['N']+Ascii['n'] ;
			int  eff_Acout=(Windows)-N_Aleng ;
			if  (eff_Acout<cutOffLen)
			{
				OUTGC<<(j+1)<<"\tNA\tNA"<<endl;
				continue ;
			}
			int  GC_Aleng=Ascii['C']+Ascii['G'] + Ascii['c']+Ascii['g'] ;
			double MeanDepth=TotalDepth*1.0/eff_Acout;
			double GC_rio=GC_Aleng*100.0/eff_Acout;
			OUTGC<<(j+1)<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\t"<<MeanDepth<<"\t"<<GC_rio<<endl;

		}
	}

	OUTGC.close();
	kseq_destroy(seq);
	gzclose(fp);


	DrawDepthGC(  GCDepth ,  paraFA04 ) ;



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
#endif // bamDepGC_H_  //
///////// swimming in the sky and flying in the sea ////////////





