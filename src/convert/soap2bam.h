#ifndef soap2bam_H_
#define soap2bam_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <cstdlib>

#include <zlib.h>
#include <stdio.h>


#include "../ALL/comm.h"
#include "../ALL/DataClass.h"
#include "../ALL/dealfun.cpp"


//KSEQ_INIT(gzFile, gzread)

using namespace std;
typedef long long  llong ;
int msort_main(int argc, char **argv) ;


void print_Ausage_Formt01()
{
	cout <<""
		"\n"
		"\tUsage: soap2bam  -i <in.soap> -s <out.sam>\n"
		"\tUsage: soap2bam  -i <in.soap> -b <out.bam> -d Ref.dict\n"
		"\n"
		"\t\t-i   <str>   input  Soap file\n"
		"\t\t-b   <str>   output BAM file\n"
		"\t\t-s   <str>   output SAM file\n"
		"\n"
		"\t\t-d   <str>   input ref FASTA to get header for BAM\n"
		"\t\t-p           if soap is PairOut,for flag\n"
		"\t\t-Q   <int>   shift sequence quality score by [+31] or [-31] or [0], default [0]\n"
		"\t\t-g           all read in memory to search mate information\n"
		"\n"
		"\t\t-h           show more help [hewm2008 v1.04]\n"
		"\n";
}

void More_Help1()
{
	cout<<""
		"\n"
		"\t\t1. Soap2bam -i <in.soap> -s AAA\n"
		"\t\t   This will convert the SOAP to SAM and output to a compressed file named AAA in current directory\n"
		"\n"
		"\t\t2. Soap2bam -i <in.soap> -b AAA -d Ref.fa\n"
		"\t\t   This will convert the SOAP to BAM with the header from reference FASTA and output to a compressed file named AAA in current directory.\n"
		"\n";
}



int parse_Acmd_Formt01(int argc, char **argv , Para_Formt01 * para_Formt01)
{
	if (argc <=1 ) {print_Ausage_Formt01();return 0;}

	for(int i = 1; i < argc; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InSoap" || flag  == "i")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_Formt01->input_file=argv[i];
		}           
		else if (flag  ==  "OutSam" || flag  == "s")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_Formt01->OutSamFile=argv[i];
		}
		else if (flag  ==  "OutBam"  || flag  == "b")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_Formt01->OutBamFile=argv[i];
		}
		else if (flag  ==  "Dict" ||  flag  == "d")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_Formt01->Dict=argv[i];
		}
		else if (flag  ==  "ShiftQ"  ||  flag  == "Q")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_Formt01->shiftQ=atoi(argv[i]);
		}
		else if (flag  ==  "Pair"  ||  flag  == "p"  )
		{
			(para_Formt01->PE)=1;
		}
		else if (flag  ==  "NoOri"||  flag  == "g" )
		{
			(para_Formt01->sort)=1;
		}       
		else if (flag  == "help"  ||  flag  == "h")
		{
			More_Help1();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}


	if  ((para_Formt01->input_file).empty()  )
	{
		cerr<< "lack argument for the must [-i] "<<endl ;
		return 0;
	}
	if ( (para_Formt01->OutSamFile).empty()    &&    (para_Formt01->OutBamFile).empty() )
	{
		cerr<< "lack argument for the must [-s]  or [-b] "<<endl ;
		return 0;
	}
	else if (   (! (para_Formt01->OutSamFile).empty() ))
	{
		(para_Formt01->OutSamFile)=add_Asuffix( para_Formt01->OutSamFile );
	}
	else
	{
		if ((para_Formt01->Dict).empty())
		{
			cerr<< "[-b] must together with [-d] "<<endl ;
			return 0;
		}
		string path=(para_Formt01->OutBamFile);
		string ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);		
		if (ext != "bam")
		{
			(para_Formt01->OutBamFile)=path+".bam" ;
		}
	}
	return 1 ;
}





int soap2sam_mainfun ( Para_Formt01 * para_Formt01 )
{
	ogzstream OUT ((para_Formt01->OutSamFile).c_str());
	if(!OUT.good())
	{
		cerr << "open InputFile error: "<<(para_Formt01->OutSamFile)<<endl;
		return 0;
	}

	if (!(para_Formt01->Dict).empty())
	{
		string path=(para_Formt01->Dict);
		string ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);

		if (ext=="dict")
		{
			Write_Sam_head (para_Formt01->Dict,OUT) ;
		}
		else
		{
			int AA=path.rfind('/') ==string::npos ? path.length() : path.rfind('/')+1 ;
			string File_Pre =path.substr(0,AA );
			ext =path.substr(AA);
			ext=replace_all(ext,".fasta.gz","");
			ext=replace_all(ext,".fasta","");
			ext=replace_all(ext,".fa.gz","");
			ext=replace_all(ext,".fa","");
			string newDictPath=File_Pre+ext+".dict";

			if ( access(newDictPath.c_str(), 0) == 0 )
			{
				Write_Sam_head (newDictPath , OUT) ;
			}
			else
			{
				string newDictPathV2=(para_Formt01->Dict)+".dict";
				if ( access(newDictPathV2.c_str(), 0) == 0 )
				{
					Write_Sam_head (newDictPathV2 , OUT) ;
				}
				else
				{
					string newDictPathV3=(para_Formt01->Dict)+".chrlist";
					if ( access(newDictPathV3.c_str(), 0) == 0 )
					{
						ifstream IN (newDictPathV3.c_str(),ifstream::in);
						string  line ;
						getline(IN,line);
						while(!IN.eof())
						{
							getline(IN,line);
							if (line.length()<=0)  { continue ; }
							vector<string> inf;
							split(line,inf," \t");
							OUT<<"@SQ\tSN:"<<inf[0]<<"\tLN:"<<inf[1]<<endl;
						}
						IN.close();
					}
					else
					{
						gzFile fp;
						kseq_t *seq;
						int l;
						fp = gzopen((para_Formt01->Dict).c_str(), "r");
						seq = kseq_init(fp);
						while ((l = kseq_read(seq)) >= 0)
						{
							OUT<<"@SQ\tSN:"<<(seq->name.s)<<"\tLN:"<<(seq->seq.l)<<endl;
						}
						kseq_destroy(seq);
						gzclose(fp);
					}
				}
			}
		}
	}




	if ( (para_Formt01->shiftQ)==0)
	{



		if (para_Formt01->sort)
		{
			char * A = const_cast<char*>((para_Formt01->input_file).c_str());
			//char * B = const_cast<char*>(tmpo.c_str());
			char * ssTmp[3]={(char *)"msort", (char *)"-k1", A };
			file_t *db;
			string  outPut ;
			db=ReadParaAnd ( 3 , ssTmp , outPut) ;
			sort_file_lines(db);   

			SamLine  *sam1 = new SamLine ;
			SamLine  *sam2 = new SamLine ;
			SamLine  *samtmp = new SamLine ;

			for(int i=0;i<db->line_size;i++)
			{
				string line=db->lines[i].str ;
				if ( soap2sam (line ,  para_Formt01 , sam1 ))
				{
					if ( (sam2->IF)  && ( (sam1->RID) == (sam2->RID))) 
					{
						mating( sam1 ,  sam2  );
						sam2->Print(OUT); sam2->rm();              
						sam1->Print(OUT); sam1->rm();
					}
					else
					{
						if ((sam2->IF))
						{
							sam2->Flag |= 0x8 ;
							sam2->Print(OUT);
							sam2->Flag |= 0x8 ;
						}
						samtmp->copy(sam2);
						sam2->copy(sam1); 
						sam1->copy(samtmp);
					}
				}
			}

			if ((sam2->IF))
			{
				sam2->Print(OUT); sam2->rm();
			}
			free_file_lines(db); 
			delete sam1 ;
			delete sam2 ;
			delete samtmp ;
		}
		else
		{
			igzstream IN ((para_Formt01->input_file).c_str(),ifstream::in);

			if(!IN.good())
			{
				cerr << "open InputFile error: "<<(para_Formt01->input_file)<<endl;
				return 1;
			}

			SamLine  *sam1 = new SamLine ;
			SamLine  *sam2 = new SamLine ;
			SamLine  *samtmp = new SamLine ;

			while(!IN.eof())
			{
				string  line ;
				getline(IN,line);
				if ( soap2sam (line ,  para_Formt01 , sam1 ))
				{
					if ( (sam2->IF)  && ( (sam1->RID) == (sam2->RID) )) 
					{
						mating( sam1 ,  sam2  );
						sam2->Print(OUT); sam2->rm();              
						sam1->Print(OUT); sam1->rm();
					}
					else
					{
						if ((sam2->IF))
						{
							sam2->Flag |= 0x8 ;
							sam2->Print(OUT);
							sam2->Flag |= 0x8 ;
						}
						samtmp->copy(sam2);
						sam2->copy(sam1); 
						sam1->copy(samtmp);
					}
				}
			}

			if ((sam2->IF))
			{
				sam2->Print(OUT); sam2->rm();
			}

			delete sam1 ;
			delete sam2 ;
			delete samtmp ;
			IN.close();    
		}




	}

	else
	{



		if (para_Formt01->sort)
		{
			char * A = const_cast<char*>((para_Formt01->input_file).c_str());
			char * ssTmp[3]={(char *)"msort", (char *)"-k1", A };
			file_t *db;
			string  outPut ;
			db=ReadParaAnd ( 3 , ssTmp , outPut) ;
			sort_file_lines(db);   

			SamLine  *sam1 = new SamLine ;
			SamLine  *sam2 = new SamLine ;
			SamLine  *samtmp = new SamLine ;

			for(int i=0;i<db->line_size;i++)
			{
				string line=db->lines[i].str ;
				if ( soap2samQ (line ,  para_Formt01 , sam1 ))
				{
					if ( (sam2->IF)  && ( (sam1->RID) == (sam2->RID))) 
					{
						mating( sam1 ,  sam2  );
						sam2->Print(OUT); sam2->rm();              
						sam1->Print(OUT); sam1->rm();
					}
					else
					{
						if ((sam2->IF))
						{
							sam2->Flag |= 0x8 ;
							sam2->Print(OUT);
							sam2->Flag |= 0x8 ;
						}
						samtmp->copy(sam2);
						sam2->copy(sam1); 
						sam1->copy(samtmp);
					}
				}
			}

			if ((sam2->IF))
			{
				sam2->Print(OUT); sam2->rm();
			}
			free_file_lines(db); 
			delete sam1 ;
			delete sam2 ;
			delete samtmp ;
		}
		else
		{
			igzstream IN ((para_Formt01->input_file).c_str(),ifstream::in);

			if(!IN.good())
			{
				cerr << "open InputFile error: "<<(para_Formt01->input_file)<<endl;
				return 1;
			}

			SamLine  *sam1 = new SamLine ;
			SamLine  *sam2 = new SamLine ;
			SamLine  *samtmp = new SamLine ;

			while(!IN.eof())
			{
				string  line ;
				getline(IN,line);
				if ( soap2samQ (line ,  para_Formt01 , sam1 ))
				{
					if ( (sam2->IF)  && ( (sam1->RID) == (sam2->RID) )) 
					{
						mating( sam1 ,  sam2  );
						sam2->Print(OUT); sam2->rm();              
						sam1->Print(OUT); sam1->rm();
					}
					else
					{
						if ((sam2->IF))
						{
							sam2->Flag |= 0x8 ;
							sam2->Print(OUT);
							sam2->Flag |= 0x8 ;
						}
						samtmp->copy(sam2);
						sam2->copy(sam1); 
						sam1->copy(samtmp);
					}
				}
			}

			if ((sam2->IF))
			{
				sam2->Print(OUT); sam2->rm();
			}

			delete sam1 ;
			delete sam2 ;
			delete samtmp ;
			IN.close();    
		}



	}

	OUT.close();
	return 1;
}




int soap2bam_mainfun ( Para_Formt01 * para_Formt01 )
{

	htsFile *OUTBam = hts_open((para_Formt01->OutBamFile).c_str(), "wb");
	bam_hdr_t *header=NULL;

	string path=(para_Formt01->Dict);

	string ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);

	if (ext=="dict")
	{
		ifstream IN (path.c_str(),ifstream::in);
		kstring_t str = { 0, 0, NULL };
		while(!IN.eof())
		{
			string  line ;
			getline(IN,line);
			if (line.length()<=0)  { continue ; }
			kputs(line.c_str(), &str);
			kputc('\n', &str);
		}
		IN.close();
		//		if (str.l == 0) kputsn("", 0, &str);
		header= sam_hdr_parse(str.l, str.s);
		header->l_text = str.l; header->text = str.s;
		header=sam_hdr_sanitise(header);
	}
	else
	{
		int AA=path.rfind('/') ==string::npos ? path.length() : path.rfind('/')+1 ;
		string File_Pre =path.substr(0,AA);
		ext =path.substr(AA);
		ext=replace_all(ext,".fasta.gz","");
		ext=replace_all(ext,".fasta","");
		ext=replace_all(ext,".fa.gz","");
		ext=replace_all(ext,".fa","");
		ext=replace_all(ext,".dict","");
		string newDictPath=File_Pre+ext+".dict";
		if ( access(newDictPath.c_str(), 0) == 0 )
		{

			ifstream IN (newDictPath.c_str(),ifstream::in);
			kstring_t str = { 0, 0, NULL };
			while(!IN.eof())
			{
				string  line ;
				getline(IN,line);
				if (line.length()<=0)  { continue ; }
				kputs(line.c_str(), &str);
				kputc('\n', &str);

			}
			IN.close();
			header= sam_hdr_parse(str.l, str.s);
			header->l_text = str.l; header->text = str.s;
			header=sam_hdr_sanitise(header);
		}

		string newDictPathV2=(para_Formt01->Dict)+".dict";
		if ( access(newDictPathV2.c_str(), 0) == 0 )
		{
			ifstream IN (newDictPathV2.c_str(),ifstream::in);
			kstring_t str = { 0, 0, NULL };
			while(!IN.eof())
			{
				string  line ;
				getline(IN,line);
				if (line.length()<=0)  { continue ; }
				kputs(line.c_str(), &str);
				kputc('\n', &str);
			}
			IN.close();
			header= sam_hdr_parse(str.l, str.s);
			header->l_text = str.l; header->text = str.s;
			header=sam_hdr_sanitise(header);
		}
		else
		{
			string newDictPathV3=(para_Formt01->Dict)+".chrlist";
			if ( access(newDictPathV3.c_str(), 0) == 0 )
			{

				ifstream IN (newDictPathV3.c_str(),ifstream::in);
				kstring_t str = { 0, 0, NULL };
				string  line ;
				getline(IN,line);
				while(!IN.eof())
				{
					getline(IN,line);
					if (line.length()<=0)  { continue ; }
					vector<string> inf;
					split(line,inf," \t");
					string Nline="@SQ\tSN:"+inf[0]+"\tLN:"+inf[1]+"\n";
					kputs(Nline.c_str(), &str);
				}
				IN.close();

				header= sam_hdr_parse(str.l, str.s);
				header->l_text = str.l; header->text = str.s;
				header=sam_hdr_sanitise(header);

			}
			else
			{
				header=Fa_hdr_read(para_Formt01->Dict);
			}
		}
	}


	if (sam_hdr_write(OUTBam, header) < 0)
	{
		fprintf(stderr, "Error writing output.\n");
		exit(-1);
	}

	bam1_t *aln = bam_init1();

	int tmp=0;



	if ( (para_Formt01->shiftQ)==0)
	{


		if (para_Formt01->sort)
		{
			char * A = const_cast<char*>((para_Formt01->input_file).c_str());
			//char * B = const_cast<char*>(tmpo.c_str());
			char * ssTmp[3]={(char *)"msort", (char *)"-k1", A };
			file_t *db;
			string  outPut ;
			db=ReadParaAnd ( 3 , ssTmp , outPut) ;
			sort_file_lines(db);   

			SamLine  *sam1 = new SamLine ;
			SamLine  *sam2 = new SamLine ;
			SamLine  *samtmp = new SamLine ;
			string OutStrA ,OutStrB;

			for(int i=0;i<db->line_size;i++)
			{
				string line=db->lines[i].str ;
				if ( soap2sam (line ,  para_Formt01 , sam1 ))
				{
					if ( (sam2->IF)  && ( (sam1->RID) == (sam2->RID))) 
					{
						mating( sam1 ,  sam2  );


						kstring_t str2 = { 0, 0, NULL };
						sam2->OUT2str(OutStrB);
						kputs(OutStrB.c_str(), &str2);
						sam2bam2( str2 , header , aln);
						tmp=sam_write1(OUTBam, header, aln);
						sam2->rm();     

						kstring_t str1 = { 0, 0, NULL };
						sam1->OUT2str(OutStrA);
						kputs(OutStrA.c_str(), &str1);
						sam2bam2( str1 , header , aln);
						tmp=sam_write1(OUTBam, header, aln);
						sam1->rm();     

					}
					else
					{
						if ((sam2->IF))
						{
							sam2->Flag |= 0x8 ;
							kstring_t str2 = { 0, 0, NULL };
							sam2->OUT2str(OutStrB);
							kputs(OutStrB.c_str(), &str2);
							sam2bam2( str2 , header , aln);
							tmp=sam_write1(OUTBam, header, aln);
							sam2->Flag |= 0x8 ;
						}
						samtmp->copy(sam2);
						sam2->copy(sam1); 
						sam1->copy(samtmp);
					}
				}
			}

			if ((sam2->IF))
			{
				kstring_t str2 = { 0, 0, NULL };
				sam2->OUT2str(OutStrB);
				kputs(OutStrB.c_str(), &str2);
				sam2bam2( str2 , header , aln);
				tmp=sam_write1(OUTBam, header, aln);
				sam2->rm();  
			}
			free_file_lines(db); 
			delete sam1 ;
			delete sam2 ;
			delete samtmp ;
		}

		else
		{
			igzstream IN ((para_Formt01->input_file).c_str(),ifstream::in);

			if(!IN.good())
			{
				cerr << "open InputFile error: "<<(para_Formt01->input_file)<<endl;
				return 1;
			}

			SamLine  *sam1 = new SamLine ;
			SamLine  *sam2 = new SamLine ;
			SamLine  *samtmp = new SamLine ;
			string OutStrA,OutStrB;

			while(!IN.eof())
			{
				string  line ;
				getline(IN,line);
				if ( soap2sam (line ,  para_Formt01 , sam1 ))
				{
					if ( (sam2->IF)  && ( (sam1->RID) == (sam2->RID) )) 
					{
						mating( sam1 ,  sam2  );


						kstring_t str2 = { 0, 0, NULL };
						sam2->OUT2str(OutStrB);
						kputs(OutStrB.c_str(), &str2);
						sam2bam2( str2 , header , aln);
						tmp=sam_write1(OUTBam, header, aln);
						sam2->rm();     

						kstring_t str1 = { 0, 0, NULL };
						sam1->OUT2str(OutStrA);
						kputs(OutStrA.c_str(), &str1);
						sam2bam2( str1 , header , aln);
						tmp=sam_write1(OUTBam, header, aln);
						sam1->rm();     

					}
					else
					{
						if ((sam2->IF))
						{
							sam2->Flag |= 0x8 ;

							kstring_t str2 = { 0, 0, NULL };
							sam2->OUT2str(OutStrB);
							kputs(OutStrB.c_str(), &str2);
							sam2bam2( str2 , header , aln);
							tmp=sam_write1(OUTBam, header, aln);
							sam2->rm();  

							sam2->Flag |= 0x8 ;
						}
						samtmp->copy(sam2);
						sam2->copy(sam1); 
						sam1->copy(samtmp);
					}
				}
			}

			if ((sam2->IF))
			{
				kstring_t str2 = { 0, 0, NULL };
				sam2->OUT2str(OutStrB);
				kputs(OutStrB.c_str(), &str2);
				sam2bam2( str2 , header , aln);
				tmp=sam_write1(OUTBam, header, aln);
				sam2->rm();     
			}

			delete sam1 ;
			delete sam2 ;
			delete samtmp ;
			IN.close();    
		}


	}
	else
	{




		if (para_Formt01->sort)
		{
			char * A = const_cast<char*>((para_Formt01->input_file).c_str());
			//char * B = const_cast<char*>(tmpo.c_str());
			char * ssTmp[3]={(char *)"msort", (char *)"-k1", A };
			file_t *db;
			string  outPut ;
			db=ReadParaAnd ( 3 , ssTmp , outPut) ;
			sort_file_lines(db);   

			SamLine  *sam1 = new SamLine ;
			SamLine  *sam2 = new SamLine ;
			SamLine  *samtmp = new SamLine ;
			string OutStrA ,OutStrB;

			for(int i=0;i<db->line_size;i++)
			{
				string line=db->lines[i].str ;
				if ( soap2samQ (line ,  para_Formt01 , sam1 ))
				{
					if ( (sam2->IF)  && ( (sam1->RID) == (sam2->RID))) 
					{
						mating( sam1 ,  sam2  );


						kstring_t str2 = { 0, 0, NULL };
						sam2->OUT2str(OutStrB);
						kputs(OutStrB.c_str(), &str2);
						sam2bam2( str2 , header , aln);
						tmp=sam_write1(OUTBam, header, aln);
						sam2->rm();     

						kstring_t str1 = { 0, 0, NULL };
						sam1->OUT2str(OutStrA);
						kputs(OutStrA.c_str(), &str1);
						sam2bam2( str1 , header , aln);
						tmp=sam_write1(OUTBam, header, aln);
						sam1->rm();     

					}
					else
					{
						if ((sam2->IF))
						{
							sam2->Flag |= 0x8 ;
							kstring_t str2 = { 0, 0, NULL };
							sam2->OUT2str(OutStrB);
							kputs(OutStrB.c_str(), &str2);
							sam2bam2( str2 , header , aln);
							tmp=sam_write1(OUTBam, header, aln);
							sam2->Flag |= 0x8 ;
						}
						samtmp->copy(sam2);
						sam2->copy(sam1); 
						sam1->copy(samtmp);
					}
				}
			}

			if ((sam2->IF))
			{
				kstring_t str2 = { 0, 0, NULL };
				sam2->OUT2str(OutStrB);
				kputs(OutStrB.c_str(), &str2);
				sam2bam2( str2 , header , aln);
				tmp=sam_write1(OUTBam, header, aln);
				sam2->rm();  
			}
			free_file_lines(db); 
			delete sam1 ;
			delete sam2 ;
			delete samtmp ;
		}

		else
		{
			igzstream IN ((para_Formt01->input_file).c_str(),ifstream::in);

			if(!IN.good())
			{
				cerr << "open InputFile error: "<<(para_Formt01->input_file)<<endl;
				return 1;
			}

			SamLine  *sam1 = new SamLine ;
			SamLine  *sam2 = new SamLine ;
			SamLine  *samtmp = new SamLine ;
			string OutStrA,OutStrB;

			while(!IN.eof())
			{
				string  line ;
				getline(IN,line);
				if ( soap2samQ (line ,  para_Formt01 , sam1 ))
				{
					if ( (sam2->IF)  && ( (sam1->RID) == (sam2->RID) )) 
					{
						mating( sam1 ,  sam2  );



						kstring_t str2 = { 0, 0, NULL };
						sam2->OUT2str(OutStrB);
						kputs(OutStrB.c_str(), &str2);
						sam2bam2( str2 , header , aln);
						tmp=sam_write1(OUTBam, header, aln);
						sam2->rm();

						kstring_t str1 = { 0, 0, NULL };
						sam1->OUT2str(OutStrA);
						kputs(OutStrA.c_str(), &str1);
						sam2bam2( str1 , header , aln);
						tmp=sam_write1(OUTBam, header, aln);
						sam1->rm();

					}
					else
					{
						if ((sam2->IF))
						{
							sam2->Flag |= 0x8 ;

							kstring_t str2 = { 0, 0, NULL };
							sam2->OUT2str(OutStrB);
							kputs(OutStrB.c_str(), &str2);
							sam2bam2( str2 , header , aln);
							tmp=sam_write1(OUTBam, header, aln);
							sam2->rm();  

							sam2->Flag |= 0x8 ;
						}
						samtmp->copy(sam2);
						sam2->copy(sam1); 
						sam1->copy(samtmp);
					}
				}
			}

			if ((sam2->IF))
			{
				kstring_t str2 = { 0, 0, NULL };
				sam2->OUT2str(OutStrB);
				kputs(OutStrB.c_str(), &str2);
				sam2bam2( str2 , header , aln);
				tmp=sam_write1(OUTBam, header, aln);
				sam2->rm();
			}

			delete sam1 ;
			delete sam2 ;
			delete samtmp ;
			IN.close();
		}







	}

	bam_hdr_destroy(header);
	bam_destroy1(aln);
	sam_close(OUTBam);
	return 1;
}





int soap2bam_main(int argc, char *argv[])
{
	Para_Formt01 * para_Formt01 = new Para_Formt01;
	if (parse_Acmd_Formt01(argc, argv , para_Formt01 )==0)
	{
		delete  para_Formt01 ;
		return 1;
	}


	if (!((para_Formt01->OutBamFile).empty()))
	{
		soap2bam_mainfun (  para_Formt01 );
	}
	else
	{
		soap2sam_mainfun (  para_Formt01 );
	}

	delete para_Formt01 ;
	return 0;

}


///////// swimming in the sky and flying in the sea ////////////
//


#endif // soap2bam_H_


////////////////////////swimming in the sea & flying in the sky //////////////////



