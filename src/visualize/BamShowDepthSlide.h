#ifndef bamDepSlid_H_
#define bamDepSlid_H_
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


int  bamDepSliding_help()
{
	cout <<""
		"\n"
		"\tUsage: DepthSlide  -InList  <bam.list>  -Ref  <Ref.fa> -OutPut  <outPrefix>\n"
		"\tUsage: DepthSlide  -InWigF  <DepthGC.wig.gz>   -OutPut  <outPrefix> -DrawChr chr1,chr2\n"
		"\n"
		"\t\t-InList      <str>     Input Bam/Sam File List\n"
		"\t\t-InFile      <str>     Input Bam/Sam File File[repeat]\n"
		"\t\t-InWigF      <str>     Input Depth-GC wig File,BamCoverage OutPut\n"
		"\t\t-OutPut      <str>     OutPut Prefix File for Depth-Slide pdf\n"
		"\n"
		"\t\t-Ref         <str>     In Ref.fa Want Out Depth-GC wig info\n"
		"\t\t-Windows     <int>     Windows size for Depth-GC wig[10000]\n"
		"\t\t-SlideRatio  <float>   Windows slidint ratio(0,1] default 1\n"
		"\t\t-MinQ        <int>     Ignore too low mapQ read[10]\n"
		"\n"
		"\t\t-DrawChr     <str>     only draw some chr Depth sliding [ALL_chr]\n"
		"\t\t-MaxY        <int>     MaxDepth Yaxis to plot Result [4*meanDepth]\n"
		"\t\t-keepR                 Keep the Rscript to RePlot the Figure\n"
		"\n"
		"\t\t-help                  Show this help [hewm2008 v1.20]\n"
		"\n";
	return 1;
}

int bamDepSliding_help01(int argc, char **argv , In3str1v * paraFA04 ,  string & OnlyChr , string & MaxY )
{
	if (argc <=2 ) {bamDepSliding_help();return 0;}
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

		if (flag  == "InList" ||  flag  == "List")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			string A=argv[i];
			file_count+=(ReadList(A ,(paraFA04->List)));
		}
		else if (flag  == "InFile")
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

		else if (flag  ==  "InWigF")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr1=argv[i];
			file_count++;
		}
		else if (flag  ==  "MaxY")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			MaxY=argv[i];
		}
		else if (flag  ==  "DrawChr")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			OnlyChr=argv[i];
		}	
		else if (flag  ==  "OutPut")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr2=argv[i];
		}
		else if (flag  ==  "Ref")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr3=argv[i];
		}
		else if (flag  ==  "Windows")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InInt2=atoi(argv[i]);
		}
		else if (flag  ==  "MinQ")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InInt=atoi(argv[i]);
		}
		else if (flag  ==  "SlideRatio")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InF=atof(argv[i]);
		}
		else if (flag  ==  "keepR")
		{
			paraFA04->TF=false ;
		}
		else if (flag == "help")
		{
			bamDepSliding_help();return 0;
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


	if  ( (paraFA04->InInt2)<10 )
	{
		cerr<< " windows size too small"<<endl;
		return 0;
	}
	
	if (!((paraFA04->InStr1).empty()))
	{
		return 1 ;
	}

	if (  ( file_count<1 )  ||  (paraFA04->InStr3).empty())
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	
	if ( (paraFA04->InF) <=0  ||   (paraFA04->InF) > 1 )
	{
		cerr <<"SlideRatio  Must be the in  (0,1] "<<endl;
		return 0;
	}
	//	(paraFA04->InStr2)=add_Asuffix(paraFA04->InStr2);
	return 1 ;
}


void DrawDepthSliding( string File , In3str1v * paraFA04 , string OnlyChr , string MaxYtemp )
{

	igzstream LIST (File.c_str(),ifstream::in); // igzstream
	if (!LIST.good())
	{
		cerr << "open List error: "<<File<<endl;
	}

	string OutPlotTmp=(paraFA04->InStr2)+".tmp";
	ofstream  OUTPLOT (OutPlotTmp.c_str());
	OUTPLOT<<"Chr\tSite\tDepth\n";

	string chr="";
	getline(LIST,chr);
	string site="";
	int  Depth=0;
	//	string OnlyChr=(paraFA04->List)[0];
	//	string MaxYtemp=(paraFA04->List)[1];

	int Count=0;
	ubit64_t SumDepth=0;
	vector<string> ChrVec;
	int VecSize=0;
	if  (OnlyChr!="NA")
	{
		split(OnlyChr,ChrVec,",-");
		VecSize=ChrVec.size();
	}

	bool YYY=false;

	while(!LIST.eof())
	{
		string  line ;
		getline(LIST,line);
		if (line.length()<=0)  { continue  ; }
		istringstream isone (line,istringstream::in);
		if (line[0] ==  '#')
		{
			isone>>chr;
			chr=chr.substr(2);
			if  (OnlyChr!="NA")
			{
				YYY=true ;
				for (int iii=0 ; iii<VecSize ; iii++)
				{
					if (ChrVec[iii] == chr )
					{
						YYY=false;
						break;
					}
				}
			}
			continue ; 
		}
		if (YYY)
		{
			continue ;
		}
		isone>>site>>Depth;
		Count++;
		SumDepth+=Depth;
		OUTPLOT<<chr<<"\t"<<site<<"\t"<<Depth<<"\n";
	}
	LIST.close();
	OUTPLOT.close();

	int MaxYdepth=int(SumDepth/Count)*4;
	if  (MaxYtemp!="NA")
	{
		MaxYdepth=atoi(MaxYtemp.c_str());
	}

	string OutPlotr=(paraFA04->InStr2)+".r";
	ofstream  OUTR (OutPlotr.c_str());

	OUTR<<""
		"\n"
		"library(ggplot2);\n" 
		"read.table(\""<<OutPlotTmp<<"\",header=T)->data;\ndata = data[order(data$Chr),];\nchrlist = sort(levels(data$Chr))\n"
		"data$Site = data$Site/1000\n"
		"if(length(chrlist)>1){\n"
		"for(i in 2:length(chrlist))\n"
		"{\n"
		"\tadd = max(data$Site[which(data$Chr==chrlist[i-1])])\n"
		"\tnew_site = data$Site[which(data$Chr==chrlist[i])]+add\n"
		"\tdata$Site[which(data$Chr==chrlist[i])] = new_site\n"
		"}\n"
		"}\n"
		"data$Site = data$Site/1000\n"
		"x_axis = scale_x_continuous(name=\"Site by chromosome(Mb)\",breaks=as.vector(by(data$Site,data$Chr,median)),labels = sort(levels(data$Chr)))\n"
		"plot_theme = theme(axis.title=element_text(face='bold',size=10),axis.text = element_text(size=8),legend.position='none')\n"
		"p <- ggplot(data,aes(Site,Depth,col=Chr))+geom_point(size=1.25)+ylim(1,"<<MaxYdepth<<") + x_axis +plot_theme\n"
		"pdf(\""<<(paraFA04->InStr2)<<".pdf\",h=6,w=25);\n"
		"p\n"
		"dev.off()\n"
		"png(\""<<(paraFA04->InStr2)<<".png\",h=6,w=25);\n"
		"p\n"
		"dev.off()\n"<<endl;

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


int bamDepthSlide_main(int argc, char *argv[])
//int main(int argc, char *argv[])
{
	In3str1v *paraFA04 = new In3str1v;
	paraFA04->InInt=10;
	paraFA04->InF=1;
	paraFA04->InInt2=10000;
	string OnlyChr ="NA";
	string MaxYtemp="NA";
	if ((bamDepSliding_help01(argc, argv, paraFA04,OnlyChr,MaxYtemp)==0))
	{
		delete paraFA04;
		return 0;
	}
	if	(!(paraFA04->InStr1).empty())
	{
		DrawDepthSliding(paraFA04->InStr1,  paraFA04 ,OnlyChr,MaxYtemp );
		delete paraFA04;
		return 1;
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
				}
				else if  (cig==BAM_CDEL ||   cig==BAM_CREF_SKIP )
				{
					This=This+ncig;
				}
			}

			/*			
						for (int32_t j=0 ; j< ((aln->core).l_qseq) ; j++)
						{
						int32_t  This=((aln->core).pos)+j;
						depth[(aln->core).tid][This]++;
						}
						*/			
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
	int Windows=int(paraFA04->InInt2);
	int SlidingLen=int(Windows*(paraFA04->InF));
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
		for (int j =0 ; j< LastDis ; j+=SlidingLen)
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

	DrawDepthSliding( GCDepth , paraFA04 ,OnlyChr,MaxYtemp);

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
#endif // bamDepSlid_H_  //
///////// swimming in the sky and flying in the sea ////////////





