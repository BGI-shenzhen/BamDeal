#ifndef Msort_H_
#define Msort_H_

#include <iostream>
#include <fstream>
#include "msort.h"
#include <regex.h>
#include <string.h>
#include <string>
#include "../gzstream/gzstream.h"
#include <zlib.h>
using namespace std;
static file_t *file_db = NULL;

file_t* load_file_lines(FILE *in, int block_rows, int *line_breaker_table, int *col_breaker_table,
		int *col_idxs, int *col_starts, int *col_ends, int *col_types, hashci_t **enums,
		int *orders, int *str_cases, int col_size){
	file_t *db;
	line_t *line;
	sortkey_t *key;
	int default_order, default_case, default_type;
	long long i, j, m, n, k, b;
	int offset, flag;
#define BUFFER_SIZE (32 * 1024)
	char buffer[BUFFER_SIZE], *ptr, str[1024];
	db = (file_t*)malloc(sizeof(file_t));
	memset(db, 0, sizeof(file_t));
	db->file_size = 0;
	db->file_cap  = BUFFER_SIZE;
	db->text = (char*)malloc(db->file_cap + 1);
	while((n = fread(buffer, 1, BUFFER_SIZE, in)) > 0){
		if(db->file_size + n >= db->file_cap){
			db->file_cap = db->file_size + 64 * 1024 * 1024LLU + n;
			db->text = (char*)realloc(db->text, db->file_cap + 1LLU);
		}
		memcpy(db->text + db->file_size, buffer, n);
		db->file_size += n;
	}
	if(db->file_size < db->file_cap) db->text = (char*)realloc(db->text, db->file_size + 1);
	db->file_cap = db->file_size;
	db->text[db->file_size] = 0;
	db->line_size = 0;
	db->line_cap  = 1024;
	db->lines     = (line_t*)malloc(sizeof(line_t) * db->line_cap);
	i = 0;
	while(i<db->file_size){
		j = i;
		b = block_rows;
		while(j < db->file_size && b > 0){
			k = j;
			while(j < db->file_size && !line_breaker_table[(int)db->text[j]]){
				if(db->text[j] <= 0){
					db->text[j] = '*';
				}
				j ++;
			}
			if(j == k && b == block_rows){
			} else {
				b --;
			}
			j ++;
		}
		if(b == block_rows) break;
		j --;
		db->text[j] = 0;
		if(db->line_size == db->line_cap){
			db->line_cap = (size_t)(db->line_cap * 1.2) + 32;
			db->lines = (line_t*)realloc(db->lines, sizeof(line_t) * db->line_cap);
		}
		db->lines[db->line_size].line_num = db->line_size;
		db->lines[db->line_size].str      = db->text + i;
		db->lines[db->line_size].str_len  = j - i;
		db->line_size ++;
		i = j + 1;
	}
	if(db->line_size < db->line_cap) db->lines = (line_t*)realloc(db->lines, sizeof(line_t) * db->line_size);
	db->line_cap = db->line_size;
	if(db->line_size < 1){
		db->key_size = 0;
		db->key_cap  = 0;
		db->keys = NULL;
		db->param.size = 0;
		return db;
	}
	if(col_size <= 0){
		default_order = orders[0];
		default_case  = str_cases[0];
		default_type  = col_types[0];
		col_size = 0;
		j = 0;
		for(i=0;i<db->lines[0].str_len+1;i++){
			if(col_breaker_table[(int)db->lines[0].str[i]]){
				if(i > j){ 
					col_idxs[col_size]  = col_size;
					col_types[col_size] = default_type;
					str_cases[col_size] = default_case;
					orders[col_size]    = default_order;
					col_size ++;
				}
				j = i + 1;
			}
		}
	}
	for(i=0;i<MAX_COLS;i++) db->col_map[i] = -1;
	for(i=0;i<col_size;i++){
		if(col_idxs[i] >= MAX_COLS){
			fprintf(stderr, " -- Sort columns bigger than %d in %s -- %s:%d --\n", MAX_COLS, __FUNCTION__, __FILE__, __LINE__);
			exit(0);
		}
		db->col_map[col_idxs[i]] = i;
	}
	db->param.size  = col_size;
	db->param.cap   = col_size;
	db->param.col_idxs = (int*)malloc(sizeof(int) * col_size);
	memcpy(db->param.col_idxs, col_idxs, sizeof(int) * col_size);
	db->param.col_orders = (int*)malloc(sizeof(int) * col_size);
	memcpy(db->param.col_orders, orders, sizeof(int) * col_size);
	db->param.cmp_funs = (cmp_fun*)malloc(sizeof(cmp_fun) * col_size);
	db->param.cmp_ops = (int*)malloc(sizeof(int) * col_size);
	for(i=0;i<col_size;i++){
		switch(col_types[i]){
			case COL_TYPE_STR:
				if(orders[i] == ORDER_ASC){
					if(str_cases[i] == STR_CASE_YES){
						db->param.cmp_funs[i] = (cmp_fun)cmp_str_asc;
						db->param.cmp_ops[i] = CMP_STR_ASC;
					} else {
						db->param.cmp_funs[i] = (cmp_fun)cmp_str_case_asc;
						db->param.cmp_ops[i] = CMP_STR_CASE_ASC;
					}
				} else {
					if(str_cases[i] == STR_CASE_YES){
						db->param.cmp_funs[i] = (cmp_fun)cmp_str_dsc;
						db->param.cmp_ops[i] = CMP_STR_DSC;
					} else {
						db->param.cmp_funs[i] = (cmp_fun)cmp_str_case_dsc;
						db->param.cmp_ops[i] = CMP_STR_CASE_DSC;
					}
				}
				break;
			case COL_TYPE_ENUM:
			case COL_TYPE_NUM:
				if(orders[i] == ORDER_ASC){
					db->param.cmp_funs[i] = (cmp_fun)cmp_num_asc;
					db->param.cmp_ops[i] = CMP_NUM_ASC;
				} else {
					db->param.cmp_funs[i] = (cmp_fun)cmp_num_dsc;
					db->param.cmp_ops[i] = CMP_NUM_DSC;
				}
				break;
			case COL_TYPE_MIX:
				if(orders[i] == ORDER_ASC){
					if(str_cases[i] == STR_CASE_YES){
						db->param.cmp_funs[i] = (cmp_fun)cmp_mix_asc;
						db->param.cmp_ops[i] = CMP_MIX_ASC;
					} else {
						db->param.cmp_funs[i] = (cmp_fun)cmp_mix_case_asc;
						db->param.cmp_ops[i] = CMP_MIX_CASE_ASC;
					}
				} else {
					if(str_cases[i] == STR_CASE_YES){
						db->param.cmp_funs[i] = (cmp_fun)cmp_mix_dsc;
						db->param.cmp_ops[i] = CMP_MIX_DSC;
					} else {
						db->param.cmp_funs[i] = (cmp_fun)cmp_mix_case_dsc;
						db->param.cmp_ops[i] = CMP_MIX_CASE_DSC;
					}
				}
		}
	}
	db->key_size = 0;
	db->key_cap  = db->line_size * col_size;
	db->keys     = (sortkey_t*)malloc(sizeof(sortkey_t) * db->key_cap);
	memset(db->keys, 0, sizeof(sortkey_t) * db->key_cap);
	for(i=0;i<db->line_size;i++){
		line = db->lines + i;
		n = 0;
		m = 0;
		k = 0;
		for(j=0;j<line->str_len+1;j++){
			if(col_breaker_table[(int)line->str[j]]){
				if(j > n){
					k = db->col_map[m];
					if(k >= 0){
						key = db->keys + i * col_size + k;
						offset = (col_starts && col_starts[k])? col_starts[k]:0;
						key->str_val = line->str + n + offset;
						if(col_ends && col_ends[k]){
							if(j - n < col_ends[k]){
								key->num_val = j - n - offset;
							} else {
								key->num_val = col_ends[k] - offset;
							}
						} else {
							key->num_val = j - n - offset;
						}
						switch(col_types[k]){
							case COL_TYPE_STR:
								break;
							case COL_TYPE_NUM:
								memcpy(str, key->str_val, (size_t)key->num_val);
								str[(int)key->num_val] = 0;
								key->num_val = atof(str);
								break;
							case COL_TYPE_MIX:
								ptr = key->str_val;
								flag = 1;
								while(ptr < key->str_val + (int)key->num_val){
									if(*ptr >= '0' && *ptr <= '9'){
									} else {
										flag = 0;
										break;
									}
									ptr ++;
								}
								if(flag){
									memcpy(str, key->str_val, (size_t)key->num_val);
									str[(int)key->num_val] = 0;
									key->num_val = atof(str);
									key->str_val = NULL;
								}
								break;
							case COL_TYPE_ENUM:
								memcpy(str, key->str_val, (size_t)key->num_val);
								str[(int)key->num_val] = 0;
								if(str_cases[k] == STR_CASE_NO){
									ptr = str;
									while(*ptr){
										if(*ptr >= 'a' && *ptr <= 'z'){
											*ptr = *ptr + 'A' - 'a';
										}
										ptr ++;
									}
								}
								if(hci_get(enums[k], str, &offset)){
									key->num_val = offset;
								} else {
									key->num_val = 0;
								}
								break;
						}
					}
					m ++;
				}
				n = j + 1;
			}
		}
	}
	return db;
}

void free_file_lines(file_t *db){
	free(db->text);
	free(db->lines);
	free(db->keys);
	if(db->param.size){
		free(db->param.col_idxs);
		free(db->param.cmp_funs);
		free(db->param.col_orders);
	}
	free(db);
}

int _old_cmp_lines(const void *e1, const void *e2){
	line_t *line1, *line2;
	sortkey_t *key1, *key2;
	int i, ret;
	line1 = (line_t*)e1;
	line2 = (line_t*)e2;
	ret = 0;
	for(i=0;i<file_db->param.size;i++){
		key1 = file_db->keys + line1->line_num * file_db->param.size + i;
		key2 = file_db->keys + line2->line_num * file_db->param.size + i;
		if(key1->num_val == 0){
			if(key2->num_val == 0){
				return 0;
			} else {
				return (file_db->param.col_orders[i] == ORDER_ASC)? -1:1;
			}
		} else if(key2->num_val == 0){
			return (file_db->param.col_orders[i] == ORDER_ASC)? 1:-1;
		}
		ret = ((cmp_fun)(file_db->param.cmp_funs[i]))(key1, key2);
		if(ret == 0) continue;
		return ret;
	}
	return ret;
}

int cmp_lines(const void *e1, const void *e2){
	line_t *line1, *line2;
	sortkey_t *key1, *key2;
	int i, ret;
	line1 = (line_t*)e1;
	line2 = (line_t*)e2;
	ret = 0;
	for(i=0;i<file_db->param.size;i++){
		key1 = file_db->keys + line1->line_num * file_db->param.size + i;
		key2 = file_db->keys + line2->line_num * file_db->param.size + i;
		/*
		   if(key1->num_val == 0){
		   if(key2->num_val == 0){
		   return 0;
		   } else {
		   return (file_db->param.col_orders[i] == ORDER_ASC)? -1:1;
		   }
		   } else if(key2->num_val == 0){
		   return (file_db->param.col_orders[i] == ORDER_ASC)? 1:-1;
		   }
		   */
		switch(file_db->param.cmp_ops[i]){
			case CMP_STR_ASC:
				ret = my_strncmp(key1->str_val, key2->str_val, key1->num_val, key2->num_val);
				break;
			case CMP_STR_DSC:
				ret = my_strncmp(key2->str_val, key1->str_val, key2->num_val, key1->num_val);
				break;
			case CMP_STR_CASE_ASC:
				ret = my_strncasecmp(key1->str_val, key2->str_val, key1->num_val, key2->num_val);
				break;
			case CMP_STR_CASE_DSC:
				ret = my_strncasecmp(key2->str_val, key1->str_val, key2->num_val, key2->num_val);
				break;
			case CMP_NUM_ASC:
				if(key1->num_val == key2->num_val){
					ret = 0;
				} else if(key1->num_val < key2->num_val){
					return -1;;
				} else {
					return 1;;
				}
				break;
			case CMP_NUM_DSC:
				if(key1->num_val == key2->num_val){
					ret = 0;
				} else if(key1->num_val < key2->num_val){
					return 1;;
				} else {
					return -1;;
				}
				break;
			case CMP_MIX_ASC:
				if(key1->str_val){
					if(key2->str_val){
						ret = my_strncmp(key1->str_val, key2->str_val, key1->num_val, key2->num_val);
					} else {
						return 1;;
					}
				} else {
					if(key2->str_val){
						return -1;;
					} else {
						if(key1->num_val == key2->num_val){
							ret = 0;
						} else if(key1->num_val < key2->num_val){
							return -1;;
						} else {
							return 1;;
						}
					}
				}
				break;
			case CMP_MIX_DSC:
				if(key1->str_val){
					if(key2->str_val){
						ret = my_strncmp(key2->str_val, key1->str_val, key1->num_val, key2->num_val);
					} else {
						return -1;;
					}
				} else {
					if(key2->str_val){
						return 1;;
					} else {
						if(key1->num_val == key2->num_val){
							ret = 0;
						} else if(key1->num_val < key2->num_val){
							return 1;;
						} else {
							return -1;;
						}
					}
				}
				break;
			case CMP_MIX_CASE_ASC:
				if(key1->str_val){
					if(key2->str_val){
						ret = my_strncasecmp(key1->str_val, key2->str_val, key1->num_val, key2->num_val);
					} else {
						return 1;;
					}
				} else {
					if(key2->str_val){
						return -1;;
					} else {
						if(key1->num_val == key2->num_val){
							ret = 0;
						} else if(key1->num_val < key2->num_val){
							return -1;;
						} else {
							return 1;;
						}
					}
				}
				break;
			case CMP_MIX_CASE_DSC:
				if(key1->str_val){
					if(key2->str_val){
						ret = my_strncasecmp(key2->str_val, key1->str_val, key1->num_val, key2->num_val);
					} else {
						return -1;;
					}
				} else {
					if(key2->str_val){
						return 1;;
					} else {
						if(key1->num_val == key2->num_val){
							ret = 0;
						} else if(key1->num_val < key2->num_val){
							return 1;;
						} else {
							return -1;;
						}
					}
				}
				break;
			default:
				ret = ((cmp_fun)(file_db->param.cmp_funs[i]))(key1, key2);
		}
		if(ret == 0) continue;
		return ret;
	}
	return ret;
}

void sort_file_lines(file_t *db){
	file_db = db;
	qsort(db->lines, db->line_size, sizeof(line_t), cmp_lines);
}

/*
   void print_file_lines(file_t *db){
   int i;
   for(i=0;i<db->line_size;i++){
   printf("%s\n", db->lines[i].str);
   }
   }
   */

int PrintFileLines ( file_t *db , string outPut )
{
	if (outPut.empty() )    
	{
		for(int i=0;i<db->line_size;i++)
		{
			cout<<db->lines[i].str<<endl;
		}
	}
	else
	{
		string ext =outPut.substr(outPut.rfind('.') ==string::npos ? outPut.length() : outPut.rfind('.')+ 1) ;
		if (ext!="gz")
		{
			ofstream OUT (outPut.c_str());
			if(!OUT.good())
			{
				cerr << "open OUT File error: "<<outPut<<endl;
				exit(1);
				return 0 ;
			}

			for(int i=0;i<db->line_size;i++)
			{
				OUT<<db->lines[i].str<<endl;
			}
			OUT.close();
		}
		else
		{
			ogzstream OUT (outPut.c_str());
			if(!OUT.good())
			{
				cerr << "open OUT File error: "<<outPut<<endl;
				exit(1);
				return 0 ;
			}
			for(int i=0;i<db->line_size;i++)
			{
				OUT<<db->lines[i].str<<endl;
			}
			OUT.close();
		}
	}
	return 1;
}



hashci_t* init_builtin_enums(){
	hashci_t *builtin;
	char name[120];
	int i, v;
	v = 1;
	builtin = hci_init();
	for(i=1;i<=80;i++){
		sprintf(name, "CHR%d", i);
		hci_put(builtin, name, v);
		v ++;
	}
	sprintf(name, "CHRX");
	hci_put(builtin, name, v);
	v ++;
	sprintf(name, "CHRY");
	hci_put(builtin, name, v);
	v ++;
	sprintf(name, "CHRM");
	hci_put(builtin, name, v);
	v ++;
	return builtin;
}

void usage(char *prog){
	printf(
			//	"msort 2, sort file rows by multiple field, wriiten by  Ruan Jue <ruanjue@gmail.com>\n"
			" msort  -k8 -kn9 In.soap  -o Out.soap.gz\n"
			"Usage: %s [-hrfn] [-t <filed_separators>]  [-k [rfmneb]<col>[start-end]{enum1,...},...] [<file>]\n"
			"Options:\n"
			"   -h  display this document\n"
			"   -l  specify line brokers, you can define more than one characters\n"
			"        defaultly, they are <newline>\n"
			"   -L  treat N line as a block, defaultly N = 1\n"
			"   -t  specify field separators, you can define more than one characters\n"
			"        defaultly, they are <space>, <tab>\n"
			"   -r  reverse sort, it will be overwritten if fileds are specified\n"
			"   -f  ignore character`s case, ...\n"
			"   -n  treat fileds as number, ...\n"
			"   -m  treat fileds as number or string, ...\n"
			"   -o  OutPut the files[#.gz or not],Otherwise to STDOUT\n"
			"   -k  specify fileds to be sorted, eg.\n"
			"         rn10[2-6]: reverse sort 2-6 (include 6) of column 10, treat it as number\n"
			"         n11: sort column 11, treat it as number\n"
			"         f2[6-]: sort 6-end of column 2, treat it as string, ignore case\n"
			"         rfm3: reverse sort column 3, treat it as string or number, ignore string case\n"
			"         f3{red green blue}: sort column 3 , treat it as enum, ignore string case\n"
			"         b3: sort column 3, treat it as build-in enum, such as chromosomes\n"
			"         \t\tRuan Jue <ruanjue@gmail.com>\tmodify:hewm\n",prog
			);
	exit(0);
}

file_t * ReadParaAnd (int argc, char **argv , string & outPut ) {
	//int main(int argc, char **argv){
	if (argc<2) {usage(argv[0]) ;} //hewem
	file_t *db;
	int col_idxs[MAX_COLS];
	int col_starts[MAX_COLS];
	int col_ends[MAX_COLS];
	int col_types[MAX_COLS];
	int orders[MAX_COLS];
	int str_cases[MAX_COLS];
	hashci_t *enums[MAX_COLS];
	int i, j, m, n, c, col_size, col, start, end, block_rows, default_order, default_type, default_case;
	int memlimit;
	int col_breaker_table[256];
	int line_breaker_table[256];
	char *ptr, msg[1024];
	FILE *in;
	regex_t regex;
	regmatch_t matches[120];
	default_order = ORDER_ASC;
	default_type  = COL_TYPE_STR;
	default_case  = STR_CASE_YES;
	memset(col_starts, 0, MAX_COLS * sizeof(int));
	memset(col_ends, 0, MAX_COLS * sizeof(int));
	memset(enums, 0, sizeof(hashci_t*) * MAX_COLS);
	col_size = 0;
	memset(col_breaker_table, 0, sizeof(int) * 256);
	memset(line_breaker_table, 0, sizeof(int) * 256);
	col_breaker_table[' ']    = 1;
	col_breaker_table['\t']   = 1;
	line_breaker_table['\n'] = 1;
	line_breaker_table['\r'] = 1;
	block_rows = 1;
	memlimit = 2040;
	regcomp(&regex, "([rnmfeb]*)([0-9]+)(\\[([0-9]+)(-([0-9]+)?)?\\])?(\\{([^}]+)\\})?", REG_EXTENDED);
	while((c = getopt(argc, argv, "hrmnfdM:t:l:L:k:o:")) != -1){
		switch(c){
			case 'h':
				usage(argv[0]);
				break;
			case 'o':
				outPut =optarg ;
				break;
			case 'r':
				default_order = ORDER_DSC;
				break;
			case 'm':
				default_type  = COL_TYPE_MIX;
				break;
			case 'f':
				default_case  = STR_CASE_NO;
				break;
			case 'n':
				default_type  = COL_TYPE_NUM;
				break;
			case 'M':
				memlimit = atoi(optarg);
			case 't':
				j = 0;
				memset(col_breaker_table, 0, sizeof(int) * 256);
				for(i=0;i<(int)strlen(optarg);i++){
					if(optarg[i] == '\\'){
						if(j){
							col_breaker_table['\\'] = 1;
							j = 0;
						} else {
							j = 1;
						}
					} else {
						if(j){
							switch(optarg[i]){
								case 't':
									col_breaker_table['\t'] = 1;
									break;
								case 'n':
									col_breaker_table['\n'] = 1;
									break;
								case 'r':
									col_breaker_table['\r'] = 1;
									break;
								default:
									col_breaker_table[(int)optarg[i]] = 1;
							}
						} else {
							col_breaker_table[(int)optarg[i]] = 1;
						}
					}
				}
				break;
			case 'l':
				j = 0;
				memset(line_breaker_table, 0, sizeof(int) * 256);
				for(i=0;i<(int)strlen(optarg);i++){
					if(optarg[i] == '\\'){
						if(j){
							line_breaker_table['\\'] = 1;
							j = 0;
						} else {
							j = 1;
						}
					} else {
						if(j){
							switch(optarg[i]){
								case 't':
									line_breaker_table['\t'] = 1;
									break;
								case 'n':
									line_breaker_table['\n'] = 1;
									break;
								case 'r':
									line_breaker_table['\r'] = 1;
									break;
								default:
									line_breaker_table[(int)optarg[i]] = 1;
							}
						} else {
							line_breaker_table[(int)optarg[i]] = 1;
						}
					}
				}
				break;
			case 'L':
				block_rows = atoi(optarg);
				assert(block_rows >= 1);
				break;
			case 'k':
				j = 0;
				for(i=0;i<(int)strlen(optarg)+1;i++){
					if(optarg[i] == ',' || i == (int)strlen(optarg)){
						start = end = 0;
						int errcode;
						if((errcode = regexec(&regex, optarg + j, 120, matches, 0))){
							regerror(errcode, &regex, msg, 1023);
							fprintf(stderr, " -- bad field format '%s' %s\n", optarg + j, msg);
							j = i + 1;
							continue;
						}
						ptr = optarg + j + matches[2].rm_eo;
						col = strtol(optarg + j + matches[2].rm_so, &ptr, 10);
						if(col < 1){
							fprintf(stderr, " -- Error column number is 1 based, but find %d in %s -- %s:%d --\n", col,  __FUNCTION__, __FILE__, __LINE__);
							usage(argv[0]);
						}
						if(matches[4].rm_eo > matches[4].rm_so){
							ptr = optarg + j + matches[4].rm_eo;
							start = strtol(optarg + j + matches[4].rm_so, &ptr, 10);
						}
						if(matches[6].rm_eo > matches[6].rm_so){
							ptr = optarg + j + matches[6].rm_eo;
							end = strtol(optarg + j + matches[6].rm_so, &ptr, 10);
						}
						if(start > 0) start --;
						else if(start < 0) start = 0;
						if(end && end < start) end = start;
						col_idxs[col_size]   = col - 1;
						col_starts[col_size] = start;
						col_ends[col_size]   = end;
						col_types[col_size]  = COL_TYPE_STR;
						orders[col_size]     = ORDER_ASC;
						str_cases[col_size]  = STR_CASE_NO;
						for(col=matches[1].rm_so;col<matches[1].rm_eo;col++){
							switch(optarg[j + col]){
								case 'r':
									orders[col_size] = ORDER_DSC;
									break;
								case 'f':
									str_cases[col_size] = STR_CASE_YES;
									break;
								case 'n':
									col_types[col_size] = COL_TYPE_NUM;
									break;
								case 'm':
									col_types[col_size] = COL_TYPE_MIX;
									break;
								case 'b':
									col_types[col_size] = COL_TYPE_ENUM;
									enums[col_size] = init_builtin_enums();
									break;
								case 'e':
									col_types[col_size] = COL_TYPE_ENUM;
									break;
								default:
									fprintf(stderr, "Unknown modifier '%c' for field\n", optarg[j + col]);
							}
						}
						if(matches[8].rm_so < matches[8].rm_eo){
							col_types[col_size] = COL_TYPE_ENUM;
							n = matches[8].rm_so;
							col = 0;
							enums[col_size] = hci_init();
							for(m=matches[8].rm_so;m<=matches[8].rm_eo;m++){
								if(m == matches[8].rm_eo || optarg[j + m] == ' '){
									while(n < m && optarg[j + n] == ' ') n ++;
									if(m <= n) continue;
									memcpy(msg, optarg + j + n, m - n);
									msg[m-n] = 0;
									if(str_cases[col_size] == STR_CASE_NO){
										ptr = msg;
										while(*ptr){
											if(*ptr >= 'a' && *ptr <= 'z'){
												*ptr = *ptr + 'A' - 'a';
											}
											ptr ++;
										}
									}
									col ++;
									hci_put(enums[col_size], msg, col);
									n = m + 1;
								}
							}
						}
						col_size ++;
						j = i + 1;
					}
				}
				break;
		}
	}
	regfree(&regex);
	if(col_size == 0){
		col_types[0] = default_type;
		orders[0]    = default_order;
		str_cases[0] = default_case;
	}
	col_breaker_table[0] = 1;
	if(argc > optind){
		if(strcmp(argv[optind], "-") == 0){
			in = stdin;
		} else {
			if ( strlen(argv[optind]) > 3 && strcmp(argv[optind] + strlen(argv[optind]) - 3, ".gz") == 0)
			{
				char *cmd;
				cmd = (char*)malloc(strlen(argv[optind]) + 20);
				sprintf(cmd, "gzip -dc %s", argv[optind]);
				in = popen(cmd, "r");
				free(cmd);
			}
			else
			{
				in = fopen(argv[optind], "r");
			}
			if(in == NULL){
				perror("open file");
				return NULL ;
			}
		}
	} else {
		in = stdin;
	}
	db = load_file_lines(in, block_rows, line_breaker_table, col_breaker_table, col_idxs, col_starts, col_ends, col_types, enums, orders, str_cases, col_size);
	if(db == NULL){
		fprintf(stderr, "Please report below message to me <ruanjue@gmail.com>\n");
		fprintf(stderr, "-----------------------------------------\n");
		fprintf(stderr, " -- Error in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		fprintf(stderr, "-----------------------------------------\n");
		exit(1);
	}
	for(i=0;i<col_size;i++){
		if(enums[i]) hci_destroy(enums[i]);
	}
	if(in != stdin) fclose(in);
	return db ;
	/*
	   sort_file_lines(db);
	   print_file_lines(db);
	   free_file_lines(db);
	   return 0;
	   */
}

//int main(int argc, char **argv){
int msort_main(int argc, char **argv){
	file_t *db;
	string outPut ;
	db=ReadParaAnd ( argc, argv , outPut ) ;
	sort_file_lines(db);
	PrintFileLines(db ,outPut) ;
	free_file_lines(db);
}
#endif //Msort_H_

