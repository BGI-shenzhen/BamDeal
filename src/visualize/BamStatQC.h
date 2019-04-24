#ifndef BamStatQC_H_
#define BamStatQC_H_
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
#include "../ALL/dealfun.cpp"
#include <cstdlib>
#include <zlib.h>
#include <htslib/sam.h>
#include <htslib/kseq.h>
#include <stdio.h>
#include <sys/stat.h>
#define MAX_Q 128

using namespace std;
typedef unsigned long long ubit64_t;
//KSEQ_INIT(gzFile, gzread)


void  BamStatQC_help()
{
	cout<<""
		"\n"
		"\tUsage: StatQC  -i  <in.bam>  -o ./\n"
		"\tUsage: StatQC  -l  <bam.list>\n"
		"\n"
		"\t\t-i    <str>     input SAM/BAM file(s)\n"
		"\t\t-l    <str>     input list of SAM/BAM files\n"
		"\n"
		"\t\t-o    <str>     output directory, default [PWD]\n"
		"\t\t-k              output the Rscript used to generate plots\n"
		"\n"
		"\t\t-h              show more details for help\n"
		"\n";
}

void More_HelpBamStatQC()
{
	cout<<""
		"\n"
		"\t\t1. StatQC -i A.bam B.bam -k\n"
		"\t\t   StatQC -l <bam.list> -k \n"
		"\t\t   This will generate four plots (GC with depth, insert size, base distribution and quality distribution) and output the result to current directory. Output files are named with the prefix bamQC.\n"
		"\t\t    (1.1) -l lists the input files. For example, if user has two input files A.bam and B.bam, bam.list should be formatted as:\n"
		"\t\t            ./A.bam\n"
		"\t\t            ./B.bam\n"
		"\n";
}


int BamStatQC_help01(int argc, char **argv , In3str1v * paraFA04)
{
	if (argc <2 ) {BamStatQC_help();return 0;}
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

		if (flag  == "InList" ||  flag  == "l" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			string A=argv[i];
			file_count+=(ReadList(A ,(paraFA04->List)));
		}
		else if (flag  ==  "InFile" ||  flag  == "i")
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
		else if (flag  ==  "OutDir" ||  flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr2=argv[i];
		}
		else if (flag  ==  "keepR" ||  flag  == "k")
		{
			paraFA04->TF=false ;
		}
		else if (flag == "help" ||  flag  == "h")
		{
			More_HelpBamStatQC();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}

	if  ( file_count ==0 )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}

	return 1 ;
}



int BamStatQC_main(int argc, char *argv[])
//int main(int argc, char *argv[])
{
	In3str1v *paraFA04 = new In3str1v;
	(paraFA04->InStr2)="./";
	paraFA04->InF=10000;
	if ((BamStatQC_help01(argc, argv, paraFA04)==0))
	{
		delete paraFA04;
		return 0;
	}


	int Raw_ASICC=Get_qual_Data ( (paraFA04->List)[0] );
	int ReadMaxL=0;
	Get_BamInfo((paraFA04->List)[0], ReadMaxL);


	mkdir((paraFA04->InStr2).c_str() , 0755 ) ;
	(paraFA04->InStr2)=(paraFA04->InStr2)+"/bamQC";

	//static	long int cycle_qual[MAX_L][MAX_Q];//√ø∏ˆcycleµƒ÷ ¡ø∑÷≤º
	ubit64_t  **cycle_qual = new ubit64_t *[ReadMaxL]; //¬ø¬™¬±⁄êŒ†
	for(int i = 0; i < ReadMaxL; i++)
	{
		cycle_qual[i] = new ubit64_t [MAX_Q]; //¬ø¬™¬±⁄ÅŒ†
		for (int32_t j =0 ; j< MAX_Q ; j++)
		{
			cycle_qual[i][j]=0;
		}
	}


	ubit64_t  **cycle_base = new ubit64_t *[ReadMaxL]; //¬ø¬™¬±⁄
	for(int i = 0; i < ReadMaxL; i++)
	{
		cycle_base[i] = new ubit64_t [5]; //¬ø¬™¬±⁄ÅŒ†
		for (int32_t j =0 ; j<5 ; j++)
		{
			cycle_base[i][j]=0;
		}
	}


	int FileNum=(paraFA04->List).size();
	ubit64_t  TotalBase=0;
	ubit64_t  TotalRead=0;
	ubit64_t  UnmapRead=0, UnmapBase=0;
	ubit64_t  LowQRead=0;
	ubit64_t  duplicate=0;
	ubit64_t  not_primary=0;
	map <int , int > InsertSize ;
	map <int , int > GCDepth;
	map <int , int > :: iterator MapIt  ;

	uint8_t UC[16] = {0,0,1,0,2,0,0,0,3,0,0,0,0,0,0,4};

	for (int ii=0  ; ii< FileNum ; ii++ )
	{
		string  line =(paraFA04->List)[ii];
		if (line.length()<=0)  { continue;  }

		cout <<"Begin Bam :"<<line<<endl;
		bam_hdr_t *header;
		bam1_t *aln = bam_init1();
		samFile *InBam = sam_open(line.c_str(), "r");
		header = sam_hdr_read(InBam);

		while (sam_read1(InBam, header, aln) >= 0)
		{
			TotalRead++; TotalBase+=((aln->core).l_qseq);
			if ((aln->core).tid < 0)
			{
				UnmapRead++ ;  UnmapBase+=((aln->core).l_qseq) ; LowQRead++; continue ;
			}

			if ( (aln->core).flag  & 0x2 )
			{
				MapIt=InsertSize.find((aln->core).isize);
				if (MapIt==InsertSize.end())
				{
					InsertSize.insert( map <int ,int >  :: value_type ((aln->core).isize,1));	
				}
				else
				{
					(MapIt->second)++;
				}
			}

			if ( (aln->core).flag   & 0x400 )
			{
				duplicate++;
			}
			if ( (aln->core).flag   & 0x100 )
			{
				not_primary++;
			}
			if (ReadMaxL==(aln->core).l_qseq)
			{
				uint8_t  *seq=bam_get_seq(aln);
				uint8_t  *seqQ=bam_get_qual(aln);
				int DepthGC[5]={0};
				int Tmp=0;
				for (int c=0;c<ReadMaxL;c++)
				{
					Tmp=UC[bam_seqi(seq, c)];
					++cycle_base[c][Tmp];
					DepthGC[Tmp]++;
					int q = seqQ[c]-Raw_ASICC;
					if(q<0) {q=0;}
					++cycle_qual[c][q];
				}				
				Tmp=DepthGC[1]+DepthGC[2];
				MapIt=GCDepth.find(Tmp);
				if (MapIt==GCDepth.end())
				{
					GCDepth.insert( map <int ,int >  :: value_type (Tmp,1));
				}
				else
				{
					(MapIt->second)++;
				}
			}

			if ( (aln->core).qual ==0 )
			{
				LowQRead++;
			}
		}
		sam_close(InBam);
		bam_destroy1(aln);
		bam_hdr_destroy(header);
	}

	cout <<"ALL Bam Read done"<<endl;





	ubit64_t proper_pair=0;

	string OutIS=(paraFA04->InStr2)+".InsertSize";
	ofstream  OUTIS (OutIS.c_str());

	MapIt=InsertSize.begin();
	int X_InMax=0;
	int Y_InMax=0;

	OUTIS<<"#Insert\tNumber\n";
	for ( ; MapIt!=InsertSize.end(); MapIt++)
	{
		OUTIS<<MapIt->first<<"\t"<<MapIt->second<<"\n";
		proper_pair+=(MapIt->second);
		if ( (MapIt->second)>Y_InMax)
		{
			Y_InMax=(MapIt->second);
			X_InMax=MapIt->first;
		}
	}
	OUTIS.close();


	string OutGCD=(paraFA04->InStr2)+".GC-depth";
	ofstream  OUTGCD (OutGCD.c_str());
	int X_max=0;
	int Y_max=0;
	MapIt=GCDepth.begin();
	OUTGCD<<"#GC\tNumber\n";

	for ( ; MapIt!=GCDepth.end(); MapIt++)
	{			
		double PMap=(MapIt->first)*100.0/ReadMaxL;
		OUTGCD<<setiosflags(ios::fixed)<<setprecision(2)<<PMap<<"\t"<<MapIt->second<<"\n";
		if (Y_max<=(MapIt->second))
		{
			Y_max=(MapIt->second);
			X_max=PMap;
		}
	}
	OUTGCD.close();






	string OutStat=(paraFA04->InStr2)+".stat";
	cerr<<OutStat<<endl;
	ofstream  OUTS (OutStat.c_str());
	ubit64_t base[5]={0};

	for(int c=0;c<ReadMaxL;c++)
	{
		base[0]+=cycle_base[c][0];
		base[1]+=cycle_base[c][1];
		base[2]+=cycle_base[c][2];
		base[3]+=cycle_base[c][3];
		base[4]+=cycle_base[c][4];
	}

	ubit64_t baseACTG=base[0]+base[1]+base[2]+base[3];



	ubit64_t qual[MAX_Q]={0};

	for(int c=0;c<ReadMaxL;c++)
	{
		for(int q=0;q<MAX_Q;q++)
		{
			qual[q]+=cycle_qual[c][q];
		}
	}


	int qMAX=0;
	ubit64_t q30=0 ;
	ubit64_t q20=0;
	for (int q=0;q<MAX_Q;q++)
	{
		if (qual[q])
		{
			qMAX=q;
			if (q>=30) q30+=qual[q];
			if (q>=20) q20+=qual[q];
		}
	}


	int Low=33+Raw_ASICC;

	ubit64_t MappedR=TotalRead-UnmapRead;
	double PMap=MappedR*100.0/TotalRead ;
	double PD=duplicate*100.0/TotalRead ;
	double Pproper_pair=proper_pair*100.0/TotalRead ;

	double Pnot_prim=not_primary*100.0/TotalRead ;	
	int MeanRL=TotalBase/TotalRead;

	ubit64_t MappedB=TotalBase-UnmapBase;
	double PMapBase=MappedB*100.0/TotalBase;
	OUTS<<"#Read\n\t##Total:\t"<<TotalRead<<"\n"<<setiosflags(ios::fixed)<<setprecision(2)
		<<"       \t##Mapped:\t"<<MappedR<<"\t"<<PMap<<"%"<<"\n"
		<<"       \t##Proper pair:\t"<<proper_pair<<"\t"<<Pproper_pair<<"%"<<"\n"
		<<"       \t##Duplicate:\t"<<duplicate<<"\t"<<PD<<"%"<<"\n"
		<<"       \t##Not_primary:\t"<<not_primary<<"\t"<<Pnot_prim<<"%"<<"\n"
		<<"       \t##Zero MQ:\t"<<LowQRead<<"\n"
		<<"       \t##Mean Lenth:\t"<<MeanRL<<"\n"
		<<"#Base\n\t##Total:\t"<<TotalBase<<"\n"
		<<"       \t##Mapped:\t"<<MappedB<<"\t"<<PMapBase<<"%"<<"\n"
		<<"       \t##Qshift:\t"<<Low<<"\n"
		<<"       \t##GC%:\t"<<(100.0*(base[1]+base[2]))/baseACTG<<"%\n"
		<<"       \t##Q20:\t"<<(100.0*q20)/TotalBase<<"%\n"
		<<"       \t##Q30:\t"<<(100.0*q30)/TotalBase<<"%\n";
	OUTS.close();






	string OutQC=(paraFA04->InStr2)+".QC";
	ofstream OUTFile (OutQC.c_str());

	baseACTG+=base[4];

	OUTFile<<"#            A     C     G     T     N ";

	for(int q=0;q<=qMAX;q++)
	{
		OUTFile<<setw(4)<<q<<' ';
	}

	OUTFile<<endl;
	OUTFile<<"#Total   ";

	for(int b=0;b<=4;b++)  
	{
		OUTFile<<setw(5)<<setiosflags(ios::fixed)<<setprecision(1)<<((100.0*base[b])/baseACTG)<<' ';
	}
	for(int q=0;q<=qMAX;q++) 
	{
		OUTFile<<setw(4)<<(int)(1000*((double)qual[q]/TotalBase))<<" ";
	}
	OUTFile<<endl;

	for(int c=0;c<ReadMaxL;c++)	
	{
		OUTFile<<"base "<<setw(3)<<c+1<<' ';
		for(int b=0;b<=4;b++)  {OUTFile<<setw(5)<<setiosflags(ios::fixed)<<setprecision(1)<<100*(double)cycle_base[c][b]/TotalRead<<' ';}
		for(int q=0;q<=qMAX;q++)
		{
			OUTFile<<setw(4)<<(int)(1000*((double)cycle_qual[c][q]/TotalRead))<<' ';
		}
		OUTFile<<endl;
	}

	OUTFile.close();



	char   buf[2048]={'\0'};
	string cc="which  Rscript  2> /dev/null ";
	memset( buf, '\0', sizeof(buf) );
	FILE   *stream ;
	stream=popen(cc.c_str(),"r") ;
	fread( buf, sizeof(char), sizeof(buf), stream);
	string binPath=buf;
	binPath=binPath.substr(0,binPath.length()-1);


	string OutPlotr=(paraFA04->InStr2)+".r";
	ofstream  OUTR (OutPlotr.c_str());
	OUTR<<""
		"\n"
		"read.table(\""<<OutGCD<<"\",header=T)->r;\n"
		"pdf(\""<<OutGCD<<".pdf\");\n"
		"plot(r[,1],r[,2],col=\"blue\",type=\"l\",ylab=\"Depth\",xlab=\"GC content(%)\")\n"
		"text(x ="<<X_max<<", y ="<<Y_max<<", col=\"red\",labels = \" Peak GC% :"<<X_max<<"\",font=2);\n"
		"dev.off();\n"

		"read.table(\""<<OutIS<<"\",header=T)->r;\n"
		"pdf(\""<<OutIS<<".pdf\");\n"
		"plot(r[,1],r[,2],col=\"blue\",type=\"l\",ylab=\"Number\",xlab=\"Insert Size\",xlim=c(0,"<<X_InMax*4<<"))\n"
		"text(x ="<<X_InMax<<", y ="<<Y_InMax<<", col=\"red\",labels = \" Peak Insert Size% :"<<X_InMax<<"\",font=2);\n"
		"dev.off();\n"

		"read.table(\""<<OutQC<<"\")->r;\n"
		"pdf(\""<<OutQC<<".base.pdf\",w=12,h=6);\n"
		"plot(r[,2],r[,3],col=\"red\",type=\"l\",ylim=c(0,50),xlab=\"Position along reads\",ylab=\"Percentage (%)\")\n"
		"lines(r[,2],r[,4],col=\"green\",lty=2)\n"
		"lines(r[,2],r[,5],col=\"blue\",lty=3)\n"
		"lines(r[,2],r[,6],col=\"Purple\",lty=4)\n"
		"lines(r[,2],r[,7],col=\"Cyan\",lty=5)\n"
		"abline(v=max(r[,2])/2,lty=2)\n"
		"legend(\"topright\",c(\"A\",\"C\",\"G\",\"T\",\"N\"),col=c(\"red\",\"green\",\"blue\",\"Purple\",\"Cyan\"),cex=1,lty=c(1,2,3,4,5),bty=\"n\");\n"
		"dev.off();\n"

		"r2 = r[,c(2,8:ncol(r))]\n"
		"colnames(r2)=c(\"read\",0:(ncol(r2)-2))\n"

		"library(reshape)\n"
		"library(ggplot2)\n"
		"r3 = melt(r2,id=\"read\")\n"
		"colnames(r3)=c(\"read\",\"qualities\",\"depth\")\n"
		"r4 = r3\n"
		"r4$depth[which(r4$depth==0)] = NA\n"
		"pdf(\""<<OutQC<<".qual.pdf\",w=12,h=6);\n"
		"p <- ggplot(r4,aes(read,qualities,colour=depth))+geom_point(shape=15) + scale_color_gradient(low='green',high='red',na.value='white')\n"
		"text_set = theme(axis.title = element_text(size=10),axis.text = element_text(size=8),legend.text = element_text(size=8),legend.title = element_text(size=8))\n"
		"axis_set_x = scale_x_continuous(expand = c(0,0))\n"
		"theme_set = theme(panel.background = element_rect(fill = \"white\", colour = \"black\"))\n"
		"p = p + theme_set + text_set + axis_set_x\n"
		"p\n"
		"dev.off()\n"


		<<endl;

	OUTR.close();




	if (binPath == "" )
	{
		cout <<"\twarning: can't find the [Rscript] in your $PATH ; no png Figure Out"<<endl;
		cout <<"\t\tRscript "<<OutPlotr<<endl;
	}

	else
	{

		cc=binPath+"\t"+OutPlotr ;
		if (paraFA04->TF)
		{
			cc=binPath+"\t"+OutPlotr+"  ; rm  -rf  "+ OutPlotr  ;
		}
		std::system(cc.c_str()) ;
	}




	for(int i = 0; i <MAX_Q; i++)
	{
		delete[] cycle_qual[i];
	}

	for(int i = 0; i <5; i++)
	{
		delete[] cycle_base[i];
	}

	delete[] cycle_qual;
	delete[] cycle_base;








	delete paraFA04 ;
	return 0;
}
#endif // BamStatQC_H_  //
///////// swimming in the sky and flying in the sea ////////////



