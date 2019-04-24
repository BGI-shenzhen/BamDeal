#ifndef bamDepth_H_
#define bamDepth_H_
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


void bamDepth_help()
{
	cout<<""
		"\n"
		"\tUsage: DepthCov -i A.bam B.bam -o <OutPrefix>\n"
		"\tUsage: DepthCov -d <ref.depth.fa>  -o <OutPrefix>\n"
		"\n"
		"\t\t-i     <str>      input SAM/BAM/CRAM file(s)\n"
		"\t\t-l     <str>      input list of SAM/BAM/CRAM files\n"
		"\t\t-o     <str>      prefix of output file\n"
		"\n"
		"\t\t-d     <str>      depth along site in reference FASTA\n"
		"\t\t-m     <int>      x-axis of the plot, default [4*meanDepth]\n"
		"\t\t-q     <int>      the quality to filter reads, default [10]\n"
		"\t\t-k                output the Rscript used to generate plots\n"
		"\n"
		"\t\t-h                show more details for help [hewm2008 v1.40]\n"
		"\n";
}

void More_Help_bamDepth()
{
	cout<<""
		"\n"
		"\t\t1. DepthCov -i A.bam B.bam -o AAA -k -q Q\n"
		"\t\t   DepthCov -l <bam.list> -o AAA -k -q Q\n"
		"\t\t   This will generate a plot with sequencing depth in x-axis and two y-axes on different sides (left and right): base proportion and accumulative coverage. The output file will be named with the prefix AAA and output to current directory. \n"
		"\t\t     (1.1) -l lists the input files. For example, if user has two input files A.bam and B.bam, bam.list should be formatted as:\n"
		"\t\t           ./A.bam\n"
		"\t\t           ./B.bam\n"
		"\t\t     (1.2) the reads with quality lower than Q will be removed, default value of Q is 10.\n"
		"\n"
		"\t\t2. DepthCov -d <depth.fa> -o AAA -k -q Q\n"
		"\t\t   This will generate a plot with sequencing depth in x-axis and two y-axes on different sides (left and right): base proportion and accumulative coverage. The output file will be named with the prefix AAA and output to current directory. \n"
		"\t\t    (2.1) the input file <depth.fa> shows the depth along the reference FASTA. User could get this file with the function Coverage in the Statistics module. Below is an example of the usage of this command. More details could be found with the -h in this command.\n"
		"\n";
}


int bamDepth_help01(int argc, char **argv , In3str1v * paraFA04 )
{
	if (argc <2 ) {bamDepth_help();return 0;}
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


		if (flag  == "InList" ||  flag  == "List" ||  flag  == "l" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			string A=argv[i];
			file_count+=(ReadList(A ,(paraFA04->List)));
		}
		else if (flag  == "InFile"||  flag  == "i")
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
		else if (flag  ==  "InDepthF" ||  flag  == "d")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr3=argv[i];
		}
		else if (flag  ==  "OutPut" ||  flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr2=argv[i];
		}
		else if (flag  ==  "MaxXaxis" ||  flag  == "m")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InInt2=atoi(argv[i]);
		}
		else if (flag  ==  "MinQ" ||  flag  == "q")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InInt=atoi(argv[i]);
		}
		else if (flag  ==  "keepR" ||  flag  == "k")
		{
			paraFA04->TF=false ;
		}
		else if (flag == "help" ||  flag  == "h")
		{
			More_Help_bamDepth();return 0;
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
	if ( (file_count<1) &&  (paraFA04->InStr3).empty())
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}

	return 1 ;
}


void Plot_Depth_Cov( string File ,In3str1v * paraFA04 )
{
	igzstream FAIN (File.c_str(),ifstream::in);

	string OStatFile=(paraFA04->InStr2)+".Depth_Cov";
	ofstream  OUT (OStatFile.c_str());

	ubit64_t Count=0;	
	map <unsigned short int,ubit64_t> MapDepth ;
	map <unsigned short int,ubit64_t> :: iterator mapIt_Depth ;

	string chrName="";
	int Depth=0;
	getline(FAIN,chrName);
	chrName=chrName.substr(1);
	ubit64_t TotolBase =0;
	while(!FAIN.eof())
	{
		string  line ;
		getline(FAIN,line);
		if (line.length()<=0)  { continue  ; }
		if (line[0] == '>')
		{
			chrName=line.substr(1);
			//			Count=0;
		}
		else
		{
			istringstream isone (line,istringstream::in);
			while(isone>>Depth)
			{
				Count++;
				TotolBase+=Depth;
				mapIt_Depth=MapDepth.find(Depth);
				if  (mapIt_Depth==MapDepth.end())
				{
					MapDepth.insert( map <unsigned short int,ubit64_t>  :: value_type (Depth,1));
				}
				else
				{
					(mapIt_Depth->second)++;
				}
			}
		}
	}

	FAIN.close();


	if  (( paraFA04->InInt2)<2)
	{
		( paraFA04->InInt2) =(int(TotolBase/Count)+1)*4 ;
	}

	OUT<<"##Depth\tDepthNumber\tCoveraged_Base\tDepthNumber(%)\tCoverage(%)\n";
	ubit64_t SumDtmp=0;
	ubit64_t DThis=0;
	double DA_ratio=0.0;
	double DB_ratio=0.0;
	ubit64_t EndCount=0;
	for (mapIt_Depth=MapDepth.begin(); mapIt_Depth!=MapDepth.end();mapIt_Depth++)
	{
		DA_ratio= (mapIt_Depth->second)*100.0/(Count);
		SumDtmp+=(mapIt_Depth->second);
		DB_ratio= SumDtmp*100.0/(Count);
		if  ((mapIt_Depth->first)==0)
		{
			SumDtmp=SumDtmp-(mapIt_Depth->second);
			OUT<<"#0\t"<<(mapIt_Depth->second)<<"\t"<<SumDtmp<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<DA_ratio<<"\t"<<DB_ratio<<"\n";
		}
		else if  ((mapIt_Depth->first)<(paraFA04->InInt2))
		{
			OUT<<(mapIt_Depth->first)<<"\t"<<(mapIt_Depth->second)<<"\t"<<SumDtmp<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<DA_ratio<<"\t"<<DB_ratio<<"\n";
		}
		else
		{
			EndCount+=(mapIt_Depth->second);
		}
	}

	if (EndCount!=0)
	{
		DA_ratio= (EndCount)*100.0/(Count);
		DB_ratio= SumDtmp*100.0/(Count);
		OUT<<(paraFA04->InInt2)<<"\t"<<EndCount<<"\t"<<SumDtmp<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<DA_ratio<<"\t"<<DB_ratio<<endl;
	}
	OUT.close();



	string OutPlotr=(paraFA04->InStr2)+".r";
	ofstream  OUTR (OutPlotr.c_str());

	OUTR<<""
		"\n"
		"read.table(\""<<OStatFile<<"\",header=T)->r;\n"
		"pdf(\""<<(paraFA04->InStr2)<<".pdf\");\n"
		"plot(r[,1],r[,4],col=\"red\",type=\"l\",lwd=4,xaxt=\"n\",ylab=\"\",xlab=\"Sequencing Depth(X)\",col.axis=\"red\")\n"
		"axis(2,col=\"red\",col.ticks=\"red\",col.axis=\"red\")\n"
		"mtext(\"Base Proportion (%)\",side=2,line=2.2,col=\"red\")\n"
		"par(new=T)\n"
		"plot(r[,1],r[,5],col=\"blue\",type=\"l\",lwd=4,ylab=\"\",xlab=\"\",yaxt=\"n\",ylim=c(0,100))\n"
		"axis(4,col=\"blue\",col.ticks=\"blue\",col.axis=\"blue\")\n"
		"mtext(\"Accumulative coverage (%)\",side=4,line=-1,col=\"blue\")\n"
		"dev.off();\n"

		"pdf(\""<<(paraFA04->InStr2)<<".png\");\n"
		"plot(r[,1],r[,4],col=\"red\",type=\"l\",lwd=4,xaxt=\"n\",ylab=\"\",xlab=\"Sequencing Depth(X)\",col.axis=\"red\")\n"
		"axis(2,col=\"red\",col.ticks=\"red\",col.axis=\"red\")\n"
		"mtext(\"Base Proportion (%)\",side=2,line=2.2,col=\"red\")\n"
		"par(new=T)\n"
		"plot(r[,1],r[,5],col=\"blue\",type=\"l\",lwd=4,ylab=\"\",xlab=\"\",yaxt=\"n\",ylim=c(0,100))\n"
		"axis(4,col=\"blue\",col.ticks=\"blue\",col.axis=\"blue\")\n"
		"mtext(\"Accumulative coverage (%)\",side=4,line=-1,col=\"blue\")\n"
		"dev.off();\n"
		<<endl;

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
			cc=binPath+"\t"+OutPlotr+"  ; rm  -rf  "+ OutPlotr  ;
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

int bamDepthCov_main(int argc, char *argv[])
//int main(int argc, char *argv[])
{
	In3str1v *paraFA04 = new In3str1v;
	paraFA04->InInt=10;

	if ((bamDepth_help01(argc, argv, paraFA04)==0))
	{
		delete paraFA04 ;
		return 0 ;
	}

	if  (!(paraFA04->InStr3).empty())
	{
		Plot_Depth_Cov( paraFA04->InStr3 , paraFA04 );
		delete paraFA04 ;
		return 1;
	}


	string DepthOFile=(paraFA04->InStr2)+".depth.fa.gz";
	ogzstream  OUT (DepthOFile.c_str());
	if((!OUT.good()))
	{
		cerr << "open OUT File error: "<<(paraFA04->InStr2)<<endl;
		delete  paraFA04 ; return  0;
	}

	string  BamPath= (paraFA04->List)[0];
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




	for(int i = 0; i <(header->n_targets); i++)
	{
		OUT<<">"<<(header->target_name[i]);		
		for (int j =0 ; j< (header->target_len[i]) ; j++)
		{
			if (j%50000==0)
			{
				OUT<<endl<<depth[i][j];
			}
			else
			{
				OUT<<" "<<depth[i][j];
			}
		}
		OUT<<endl ;
	}
	OUT.close();



	//释放开辟的资源  
	for(int i = 0; i <(header->n_targets); i++)
	{
		delete[] depth[i];  
	}
	delete[] depth;  


	Plot_Depth_Cov( DepthOFile  , paraFA04 );

	bam_hdr_destroy(header);
	delete paraFA04 ;
	return 0;
}
#endif // bamDepth_H_  //
///////// swimming in the sky and flying in the sea ////////////





