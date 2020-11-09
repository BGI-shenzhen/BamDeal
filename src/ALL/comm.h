#ifndef comm_H_
#define comm_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <cstdlib>
#include <math.h>
#include <stdio.h>
#include  <zlib.h>
#include "gzstream/gzstream.c"
#include "kseq.h"


/*///
#ifndef DBL_AMANT_ADIG
#define DBL_AMANT_ADIG 35
#endif

#define MAX_AFILE_ALIST_ALEN 512
#define MAX_ALIST_ANAME_ALEN 1024
#define MAX_ACHR_ANAME_ALEN 128
///*////
/*///
  const char decode[16] = {'N','A','C','N','G','N','N','N','T','N','N','N','N','N','N','N'};
  const int code[10] = {0,5,15,10,1,3,2,7,6,11};
  const int rev_Acode[16] = {0,4,6,5,4,1,8,7,6,8,3,9,5,7,9,2};
  const char abbv[17] = {'A','M','W','R','M','C','Y','S','W','Y','T','K','R','S','K','G','N'};

//*/

using namespace std;

typedef long long llong ;
typedef unsigned long long ubit64_t;

KSEQ_INIT(gzFile, gzread)

	/*
	   template < class TT > 
	   string  Int2Str ( TT  A )
	   {
	   stringstream   sstrm ;
	   sstrm  <<  A ;
	   return  sstrm.str();
	   }
	   */

int ReadList (string soaplist  ,   vector <string> & Soap_AStat  )
{
	igzstream LIST (soaplist.c_str(),ifstream::in); // igzstream
	int soapfilecout=0 ;
	if (!LIST.good())
	{
		cerr << "open List error: "<<soaplist<<endl;
		return  soapfilecout ;
	}
	while(!LIST.eof())
	{
		string  line ;
		getline(LIST,line);
		if (line.length()<=0)  { continue  ; }
		Soap_AStat.push_back(line);
		soapfilecout++;
	}
	LIST.close();
	return  soapfilecout ;
}

inline void  LogLackArg( string  flag )
{
	cerr << "\t\tLack Argument for [ -"<<flag<<" ]"<<endl;    
}

inline string add_Asuffix ( string path )
{
	string ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);
	if (ext != "gz")
	{
		path=path+".gz" ; 
	}
	return path ;
}


string &  replace_all(string &  str,const  string &  old_Avalue,const string &  new_Avalue)
{
	while(true)   {
		string::size_type  pos(0);
		if(   (pos=str.find(old_Avalue))!=string::npos   )
			str.replace(pos,old_Avalue.length(),new_Avalue);
		else   break;
	}
	return   str;
}

string &   replace_all_distinct(string&   str,const   string &   old_Avalue,const   string &   new_Avalue)
{
	for(string::size_type   pos(0);   pos!=string::npos;   pos+=new_Avalue.length())   {
		if(   (pos=str.find(old_Avalue,pos))!=string::npos   )
			str.replace(pos,old_Avalue.length(),new_Avalue);
		else   break;
	}
	return   str;
}

inline string getID (string ID )
{
	string ext =ID.substr(0,ID.rfind('#')==string::npos ? ID.length() : ID.rfind('#')) ;
	return ext ;
}

void split(const string& str,vector<string>& tokens,  const string& delimiters = " ")
{
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	string::size_type pos     = str.find_first_of(delimiters, lastPos);
	while (string::npos != pos || string::npos != lastPos)
	{
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		lastPos = str.find_first_not_of(delimiters, pos);
		pos = str.find_first_of(delimiters, lastPos);
	}
}


string Int2Str (size_t A )
{
	stringstream   sstrm ;
	sstrm  <<  A ;
	return  sstrm.str();
}

string Int2Str (int A )
{
	stringstream   sstrm ;
	sstrm  <<  A ;
	return  sstrm.str();
}

string Int2Str (llong  A )
{
	stringstream   sstrm ;
	sstrm  <<  A ;
	return  sstrm.str();
}

inline void Swap ( int & x ,int & y)
{
	int tmp=y;
	y=x;
	x=tmp;
}


int stat_str_base(string str , llong * Map , llong Leng )
{
	for(llong ix=0 ; ix<Leng ; ix++)
	{
		Map[str[ix]]++;
	}
	return 1 ;
}


int   stat_AMisMatch ( string str , string chr , map <string , llong > &  MisMatch )
{
	string mis1 = str.substr(str.rfind('\t')==string::npos ?str.length():str.rfind('\t') + 1);
	string::size_type s_Aend =mis1.size();
	llong ATGC=0 ;
	for(string::size_type ix=0 ; ix<s_Aend ; ix++)
	{
		if (mis1[ix] == 'A' || mis1[ix] == 'T' || mis1[ix] == 'C' || mis1[ix] == 'G')
		{
			ATGC++;
		}
	}
	MisMatch[chr]+=ATGC ;
	return ATGC ;
}

int ReadFaSeq (string RefFile , map <string,ubit64_t> & ChrLeng , map <string,string> & Seq)
{
	gzFile fp;
	kseq_t *seq;
	fp = gzopen(RefFile.c_str(), "r");
	seq = kseq_init(fp);
	int l=0, File_Acount=0 ;
	while ( (l = kseq_read(seq)) >= 0 )
	{
		string chr=seq->name.s;
		ubit64_t length_chr=ubit64_t(seq->seq.l);
		ChrLeng[chr]=length_chr;
		Seq[chr]=(seq->seq.s);
		File_Acount++;
	}
	kseq_destroy(seq);
	gzclose(fp);
	return  File_Acount ;
}

void Display( string seq , string ID , int linecut=100  )
{
	long sumleng=seq.length();
	int Endline=int(sumleng/linecut);
	cout<<">"<<ID<<endl ;
	for (int i=0 ; i< Endline ; i++ )
	{
		cout<<seq.substr(i*linecut,linecut)<<"\n";
	}
	long  AA=Endline*linecut ;
	if (sumleng> AA )
	{
		cout<<seq.substr(AA)<<endl;
	}
}



void Display( string seq , string ID , ogzstream  &  OUT , int linecut=100  )
{
	long sumleng=seq.length();
	int Endline=int(sumleng/linecut);
	OUT<<">"<<ID<<endl ;
	for (int i=0 ; i< Endline ; i++ )
	{
		OUT<<seq.substr(i*linecut,linecut)<<"\n";
	}
	long  AA=Endline*linecut ;
	if (sumleng> AA )
	{
		OUT<<seq.substr(AA)<<endl;
	}
}


void Display( string seq , ogzstream  &  OUT , int linecut=100  )
{
	long sumleng=seq.length();
	int Endline=int(sumleng/linecut);
	for (int i=0 ; i< Endline ; i++ )
	{
		OUT<<seq.substr(i*linecut,linecut)<<"\n";
	}
	long  AA=Endline*linecut ;
	if (sumleng> AA )
	{
		OUT<<seq.substr(AA)<<endl;
	}
}


int FaCutLine ( string FaPath)
{
	int linecut=100 ;
	igzstream INRef ((FaPath).c_str(),ifstream::in);
	if (INRef.fail())
	{
		cerr << "open File error: "<<(FaPath)<<endl;
		return  0;
	}

	string tmp ;
	int tmplinecut=0;
	for (int A=1 ; A<168 && (!INRef.eof())  ; A++ )
	{
		getline(INRef,tmp);
		if (tmp.length()<2)
		{
			continue ;
		}
		if (tmplinecut < tmp.length())
		{
			tmplinecut=tmp.length();
		}
	}
	INRef.close();


	if (tmplinecut<30)
	{
		linecut=100;
	}
	else if (tmplinecut>1000)
	{
		linecut=100 ;
	}
	else
	{
		linecut=tmplinecut;
	}

	return linecut ;

}


/*
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
*/

int GetShiftQ ( string FQPath )
{
	igzstream INFQ ((FQPath).c_str(),ifstream::in);
	if (INFQ.fail())
	{
		cerr << "open File error: "<<(FQPath)<<endl;
		return 64 ;
	}
	string tmp ;
	int minQ=50000;
	int maxQ=0;
	for (int A=1 ; A<46888 && (!INFQ.eof())  ; A++ )
	{
		getline(INFQ,tmp);
		if (tmp.length()<=0)  { continue  ; }

		if(A%4!=0)
		{
			continue;
		}
		string::size_type SeqQLength =tmp.size();
		for(int i=0 ; i<SeqQLength ; i++)
		{
			if(minQ>tmp[i])
			{
				minQ=tmp[i];
			}
			if(maxQ<tmp[i])
			{
				maxQ=tmp[i];
			}
		}
	}
	INFQ.close();

	if(minQ >= 33 &&  minQ <= 78  &&  maxQ >= 33 && maxQ <=78 )
	{
		return 33;
	}
	else if (minQ >= 64  &&  minQ <= 108  &&  maxQ >= 64 && maxQ <= 108)
	{
		return 64;
	}
	else
	{
		return 64 ;
	}
}


int GetShiftQSoap ( string SOAPPath ,int row )
{
	igzstream INSOAP ((SOAPPath).c_str(),ifstream::in);
	if (INSOAP.fail())
	{
		cerr << "open File error: "<<(SOAPPath)<<endl;
		return 64 ;
	}
	int minQ=50000;
	int maxQ=0;
	for ( int A=1 ; ( A<16888 && (!INSOAP.eof()))  ; A++ )
	{
		string line ;
		getline(INSOAP,line);
		if (line.length()<=0)  { continue  ; }

		vector<string> inf;
		split(line,inf," \t");

		string::size_type SeqQLength =(inf[row]).size();
		for(int i=0 ; i<SeqQLength ; i++)
		{
			if(minQ>(inf[row])[i])
			{
				minQ=(inf[row])[i];
			}
			if(maxQ<(inf[row])[i])
			{
				maxQ=(inf[row])[i];
			}
		}
	}
	INSOAP.close();
	if(minQ >= 33 &&  minQ <= 76  &&  maxQ >= 33 && maxQ <=76 )
	{
		return 33;
	}
	else if (minQ >= 64  &&  minQ <= 107  &&  maxQ >= 64 && maxQ <= 107)
	{
		return 64;
	}
	else
	{
		return 64 ;
	}
}

#endif // comm_H_  ;

