
#ifndef sameFun_H_
#define sameFun_H_


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <ctime>
#include <cstdlib>
#include <math.h>
#include <stdio.h>
#include <zlib.h>

#include "../ALL/msort/msort.h"
#include "../ALL/msort/sort_funs.c"
#include "../ALL/msort/stdhashc.cc"
#include "../ALL/msort/msort_h.h"
#include "../ALL/comm.h"

#include <htslib/sam.h>
#include <htslib/khash.h>
#include <htslib/kseq.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/kstring.h>
#include <htslib/hts_endian.h>

#ifdef UINT32_MAX

#else
    #define UINT32_MAX  (4294967295U)
#endif


//KSEQ_INIT(gzFile, gzread)


void Write_Sam_head (string InHead , ogzstream & OUT )
{

	igzstream DD (InHead.c_str(),ifstream::in);

	if(!DD.good())
	{
		cerr << "open InputFile error: "<<InHead<<endl;
		return  ;
	}
	while(!DD.eof())
	{
		string  line ;
		getline(DD,line);
		if (line.length()<2)
		{
			continue ;
		}
		OUT<<line<<endl;
	}
	DD.close();
}




///*///
bam_hdr_t * sam_hdr_sanitise(bam_hdr_t *h)
	//static bam_hdr_t * sam_hdr_sanitise(bam_hdr_t *h)
{
	if (!h)
		return NULL;

	if (h->l_text == 0)
		return h;

	uint32_t i, lnum = 0;
	char *cp = h->text, last = '\n';
	for (i = 0; i < h->l_text; i++) {
		if (cp[i] == 0)
			break;

		if (last == '\n') {
			lnum++;
			if (cp[i] != '@') {
				hts_log_error("Malformed SAM header at line %u", lnum);
				bam_hdr_destroy(h);
				return NULL;
			}
		}

		last = cp[i];
	}

	if (i < h->l_text) { // Early nul found.  Complain if not just padding.
		uint32_t j = i;
		while (j < h->l_text && cp[j] == '\0') j++;
		if (j < h->l_text)
			hts_log_warning("Unexpected NUL character in header. Possibly truncated");
	}
	if (last != '\n') {
		hts_log_warning("Missing trailing newline on SAM header. Possibly truncated");

		if (h->l_text == UINT32_MAX) {
			hts_log_error("No room for extra newline");
			bam_hdr_destroy(h);
			return NULL;
		}

		if (i >= h->l_text - 1) {
			cp = (char*) (realloc(h->text, (size_t) h->l_text+2));
			if (!cp) {
				bam_hdr_destroy(h);
				return NULL;
			}
			h->text = cp;
		}
		cp[i++] = '\n';

		if (h->l_text < i)
			h->l_text = i;
		cp[h->l_text] = '\0';
	}

	return h;
}
////*/

void  Get_BamInfo  ( string InFile , int & MaxRL)
{
	bam_hdr_t *header;
	bam1_t *aln = bam_init1();
	samFile *in = sam_open(InFile.c_str(), "r");
	header = sam_hdr_read(in);
	int  Count=0;
	uint8_t MaxPhredQ=0;
	uint8_t MinPhredQ=250;
	MaxRL=0;
	while  ( (sam_read1(in, header, aln) >= 0)   &&  (Count<16888) )
	{
		if ( (aln->core).qual < 10 )
		{
			continue ;
		}
		if ((aln->core).l_qseq> MaxRL )
		{
			MaxRL=(aln->core).l_qseq ;
		}
	}

	sam_close(in);
	bam_hdr_destroy(header);
	bam_destroy1(aln);

}




int Get_qual_Data ( string InFile)
{

	bam_hdr_t *header;
	bam1_t *aln = bam_init1();
	samFile *in = sam_open(InFile.c_str(), "r");
	header = sam_hdr_read(in);

	int  Count=0;
	uint8_t MaxPhredQ=0;
	uint8_t MinPhredQ=250;
	while  ( (sam_read1(in, header, aln) >= 0)   &&  (Count<16888) )
	{
		if ( (aln->core).qual < 10 )
		{
			continue ;
		}
		if  ((aln->core).l_qseq   < 10 )
		{
			continue ;
		}
		uint8_t  * seqQ=bam_get_qual(aln);
		Count++;
		for(int i=0; i < (aln->core).l_qseq ; i++)
		{
			if (MinPhredQ > seqQ[i])
			{
				MinPhredQ= seqQ[i] ;
			}
			if  (MaxPhredQ < seqQ[i])
			{
				MaxPhredQ = seqQ[i] ;
			}
		}
	}

	sam_close(in);
	bam_hdr_destroy(header);
	bam_destroy1(aln);

	int ASCII_raw=0;

	if(MinPhredQ >= 0 && MinPhredQ <= 42 && MaxPhredQ >= 0 && MaxPhredQ <= 42)
	{
	}
	else if  (MinPhredQ >= 31 && MinPhredQ <= 73 && MaxPhredQ >= 31 && MaxPhredQ <=73)
	{
		ASCII_raw=31 ;
	}
	return ASCII_raw ;
}



bam_hdr_t *Fa_hdr_read( string RefPath)
{
	gzFile fp;
	kseq_t *seq;
	bam_hdr_t *h = NULL;

	kstring_t str = { 0, 0, NULL };

	int l;
	fp = gzopen(RefPath.c_str(), "r");
	seq = kseq_init(fp);

	while ((l = kseq_read(seq)) >= 0)
	{
		stringstream   sstrm;
		sstrm  <<  (seq->seq.l);

		kputs("@SQ\tSN:", &str);
		kputs(seq->name.s, &str);
		kputs("\tLN:", &str);
		kputl(seq->seq.l, &str);
		kputc('\n', &str);
	}

	if (str.l == 0) kputsn("", 0, &str);
	h = sam_hdr_parse(str.l, str.s);
	h->l_text = str.l; h->text = str.s;

	kseq_destroy(seq);
	gzclose(fp);

	return sam_hdr_sanitise(h);
}



int sam2bam2( kstring_t str , bam_hdr_t *h, bam1_t *b)
{
	int ret;
	ret = sam_parse1(&str, h, b);
	if (ret < 0) {
		cerr<<"same the wrong at this "<<str.s ;
	}
	return 1;
}





int matmis( string line ,string & Base ,int & data )
{    
	vector<string> inf;
	split(line,inf,"->");
	if (inf.size()<2)
	{
		return 0 ;
	}
	else
	{
		Base=inf[0];
		int L=inf[1].length();
		int i=0 ;
		for ( i=0; i<L ; i++ )
		{
			char N= (inf[1][i]) ;
			if ( N  == 'A' ||  N == 'C' ||  N == 'T'  ||  N == 'G' )
			{
				break ;
			}
		}
		string B=inf[1].substr(0,i);
		data=atoi(B.c_str());
		return 1 ;
	}
	//    if ($ARGV[0] =~ /^([ACGT])->(\d+)/i);
}

//*//
int mating ( SamLine * A ,  SamLine * B )
{
	int Insert=0;
	if (  (A->chr) != "*" && ( (A->chr) ==(B->chr))  )  //# then calculate $isize
	{
		llong  x1=( (A->Flag) & 0x10 ) ? (A->position)+((A->seq).length()) : (A->position);
		llong  x2=( (B->Flag) & 0x10 ) ? (B->position)+((B->seq).length()) : (B->position);
		Insert=x2-x1;
	}
	//# update mate coordinate
	if  ( (B->chr) != "*" )
	{
		A->coor=B->position;
		A->isize=Insert;
		A->XorD= ((B->chr) == (A->chr)) ? "=" : (B->chr) ; 
		if ((B->Flag) & 0x10 )
		{
			A->Flag |= 0x20 ;
		}
	}
	else
	{
		A->Flag |= 0x8;
	}

	if  ( (A->chr) != "*" )
	{
		B->coor=A->position;
		B->isize=0-Insert;
		B->XorD= ((A->chr) == (B->chr)) ? "=" : (A->chr) ; 
		if ((A->Flag) & 0x10 )
		{
			B->Flag |= 0x20 ;
		}
	}
	else
	{
		B->Flag |= 0x8;
	}    
	return 1;
}



int soap2sam(string line , Para_Formt01 * para_Formt01 , SamLine *sam )
{
	vector<string> inf;
	split(line,inf," \t");
	int Length=inf.size();
	if ( Length < 9 || (inf[0].empty()) )
	{
		return 0 ;
	}
	sam->rm();
	int A=inf[3][0];
	if ( A>57 || A<48 ) // fix SOAP-2.1.x bugs
	{
		vector<string> tmp (Length-1);
		tmp[0]=inf[0]; tmp[1]=inf[1];   tmp[2]=inf[2]; 
		for(int ii=4 ; ii<Length ; ii++ )
		{
			tmp[ii-1]=tmp[ii];
		}
		inf.clear();
		inf=tmp ;
	}
	// readID //
	int RID_length=inf[0].length() ;
	//    sam->RID=getID(inf[0]);
	//*
	sam->RID=(inf[0]); 
	if ( RID_length >2 )
	{
		if ( inf[0][RID_length-2]== '/'  &&  ( (inf[0][RID_length-1]== '1')  || (inf[0][RID_length-1]== '2') ))
		{
			sam->RID=inf[0].substr(0, RID_length-2);                   
		}
	}
	///*///

	// initial flag (will be updated later)
	sam->Flag =0 ; 

	(sam->Flag) |= 1 | 1<<(inf[4] == "a" ? 6 : 7);
	if  (para_Formt01->PE)
	{
		(sam->Flag)  |= 2;
	}
	if (inf[6] == "-" )
	{           
		(sam->Flag) |= 0x10 ;
	}

	// # read & quality
	sam->seq=inf[1];
	int RLength=inf[1].length();
	int QLength=inf[2].length();
	sam->Qseq = (QLength> RLength ) ? inf[2].substr(0,RLength): inf[2] ;
	//*/////
	QLength=(sam->Qseq).length();
	for(string::size_type ix=0; ix<QLength; ix++)
	{
		(sam->Qseq)[ix]-=(para_Formt01->shiftQ) ; ////// chang the "#"  Quli to "B" Quli
	}
	///*////
	// cigar 
	sam->cigar= Int2Str(RLength)+"M";

	// # coor
	sam->chr = inf[7]; sam->position = atoi(inf[8].c_str());

	// mapQ
	sam->mapQ  = (inf[3] == "1" ) ? 30 : 0;

	// # mate coordinate
	//$s->[6] = '*'; $s->[7] = $s->[8] = 0;

	// # aux 
	sam->NM_i="NM:i:"+inf[9];
	sam->MD="";
	A=atoi(inf[9].c_str());
	if (A)
	{
		string Base ;
		int data ;
		map <int ,string  > MAP ;
		for (int ii=10 ; ii<Length ; ii++ )
		{
			if (matmis( inf[ii], Base , data))
			{
				MAP.insert(map < int,string > :: value_type (data,Base));
			}
		}
		int a=0;
		map <int , string> :: iterator it=MAP.begin();
		for(it=MAP.begin() ; it!=MAP.end(); it++)
		{
			int c =(it->first) - a ;
			(sam->MD) += Int2Str(c) + (it->second);
			a += (c  + 1);
		}
		(sam->MD) += Int2Str(RLength-a);
	}
	else
	{
		(sam->MD)=Int2Str(RLength);
	}
	(sam->MD)="MD:Z:"+(sam->MD);
	sam->IF=true ;
	return  1;
}
///////////////////




//*
int soap2samQ(string line , Para_Formt01 * para_Formt01 , SamLine *sam )
{
	vector<string> inf;
	split(line,inf," \t");
	int Length=inf.size();
	if ( Length < 9 || (inf[0].empty()) )
	{
		return 0 ;
	}
	sam->rm();
	int A=inf[3][0];
	if ( A>57 || A<48 ) // fix SOAP-2.1.x bugs
	{
		vector<string> tmp (Length-1);
		tmp[0]=inf[0]; tmp[1]=inf[1];   tmp[2]=inf[2]; 
		for(int ii=4 ; ii<Length ; ii++ )
		{
			tmp[ii-1]=tmp[ii];
		}
		inf.clear();
		inf=tmp ;
	}
	// readID //
	int RID_length=inf[0].length() ;
	//    sam->RID=getID(inf[0]);
	//*
	sam->RID=(inf[0]); 
	if ( RID_length >2 )
	{
		if ( inf[0][RID_length-2]== '/'  &&  ( (inf[0][RID_length-1]== '1')  || (inf[0][RID_length-1]== '2') ))
		{
			sam->RID=inf[0].substr(0, RID_length-2);                   
		}
	}
	///*///

	// initial flag (will be updated later)
	sam->Flag =0 ; 

	(sam->Flag) |= 1 | 1<<(inf[4] == "a" ? 6 : 7);
	if  (para_Formt01->PE)
	{
		(sam->Flag)  |= 2;
	}
	if (inf[6] == "-" )
	{           
		(sam->Flag) |= 0x10 ;
	}

	// # read & quality
	sam->seq=inf[1];
	int RLength=inf[1].length();
	int QLength=inf[2].length();
	sam->Qseq = (QLength> RLength ) ? inf[2].substr(0,RLength): inf[2] ;


	// cigar 
	sam->cigar= Int2Str(RLength)+"M";

	// # coor
	sam->chr = inf[7]; sam->position = atoi(inf[8].c_str());

	// mapQ
	sam->mapQ  = (inf[3] == "1" ) ? 30 : 0;

	// # mate coordinate
	//$s->[6] = '*'; $s->[7] = $s->[8] = 0;

	// # aux 
	sam->NM_i="NM:i:"+inf[9];
	sam->MD="";
	A=atoi(inf[9].c_str());
	if (A)
	{
		string Base ;
		int data ;
		map <int ,string  > MAP ;
		for (int ii=10 ; ii<Length ; ii++ )
		{
			if (matmis( inf[ii], Base , data))
			{
				MAP.insert(map < int,string > :: value_type (data,Base));
			}
		}
		int a=0;
		map <int , string> :: iterator it=MAP.begin();
		for(it=MAP.begin() ; it!=MAP.end(); it++)
		{
			int c =(it->first) - a ;
			(sam->MD) += Int2Str(c) + (it->second);
			a += (c  + 1);
		}
		(sam->MD) += Int2Str(RLength-a);
	}
	else
	{
		(sam->MD)=Int2Str(RLength);
	}
	(sam->MD)="MD:Z:"+(sam->MD);
	sam->IF=true ;
	return  1;
}
///////////////////





#endif 
///////// swimming in the sky and flying in the sea ////////////


