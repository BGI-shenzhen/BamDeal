
////////////////////////swimming in the sea & flying in the sky //////////////////


/*
 * DataClass.h
 *
 *  Created on: 2011-11-21
 *      Author: hewm@genomics.org.cn
 */

#ifndef DataClass_H_
#define DataClass_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "comm.h"

using namespace std;


///////////////// q_seq  for site ////////


class In3str1v {
	public:
		bool  TF ;
		int   InInt ;
		int   InInt2 ;
		bool TF2 ;
		double InF ;
		string InStr1 ;
		string InStr2 ;
		string InStr3 ;
		vector <string> List ;
		In3str1v()
		{
			InStr1="";
			InStr2="";
			InStr3="";
			TF=true ;
			TF2=true ;
			InInt=0 ;
			InInt2=0;
			InF=0.0;
		}
};



class Para_Formt01 {
	public:
		string input_file ;
		string OutSamFile ;
		string OutBamFile ;
		int PE ;
		int sort ;
		int shiftQ ;
		string Dict ;
		Para_Formt01()
		{
			input_file="" ;
			OutSamFile="" ;
			OutBamFile="" ;
			Dict="";
			PE=0 ;
			sort=0;
			shiftQ=0;
		}
};

class SamLine
{
	public:
		string RID;
		int Flag ;
		string seq;
		string Qseq;
		string cigar;
		string chr;
		llong position ;
		int mapQ ;
		string NM_i ;
		string MD ;
		bool IF ;
		int isize ;
		llong coor ;
		string XorD ;
		void copy( SamLine *  B )
		{
			RID=B->RID;   Flag=B->Flag;    seq=B->seq ;
			Qseq=B->Qseq; cigar=B->cigar; chr=B->chr;
			position=B->position ; mapQ=B->mapQ;  IF=B->IF;
			NM_i=B->NM_i ;MD=B->MD ; isize=B->isize;
			coor=B->coor; XorD=B->XorD;
		}
		SamLine()
		{
			RID="";  seq=""; Qseq=""; cigar=""; chr ="";
			Flag=0 ; position=0; mapQ =0 ; IF=false;
			NM_i=""; MD="";
			isize=0; coor=0; XorD="*";
		}
		void rm ()
		{
			RID="";  seq=""; Qseq=""; cigar=""; chr ="";
			Flag=0 ; position=0; mapQ =0 ; IF=false;
			NM_i=""; MD="";
			isize=0; coor=0; XorD="*";
			//return 1 ;
		}
		void Print(ogzstream & OUT )
		{
			OUT<<RID<<"\t"<<Flag<<"\t"<<chr<<"\t"<<position<<"\t"<<mapQ<<"\t"<<cigar<<"\t"<<XorD<<"\t"<<coor<<"\t"<<isize<<"\t"<<seq<<"\t"<<Qseq<<"\t"<<NM_i<<"\t"<<MD<<"\n";
		}
		void OUT2str(string & OutStr )
		{
			OutStr=RID+"\t"+Int2Str(Flag)+"\t"+chr+"\t"+Int2Str(position)+"\t"+Int2Str(mapQ)+"\t"+cigar+"\t"+XorD+"\t"+Int2Str(coor)+"\t"+Int2Str(isize)+"\t"+seq+"\t"+Qseq+"\t"+NM_i+"\t"+MD;
		}
};





////////swimming in the sky and flying in the sea *///////////





#endif /* DataClass_H_ */

//////////////// swimming in the sky and flying in the sea ////////////////



//////////////// swimming in the sky and flying in the sea ////////////////
