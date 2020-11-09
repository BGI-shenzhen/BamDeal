#ifndef __MSORT_H_RJ
#define __MSORT_H_RJ

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include "stdhashc.h"

#define COL_TYPE_STR 0
#define COL_TYPE_NUM 1
#define COL_TYPE_MIX 2
#define COL_TYPE_ENUM 3

#define ORDER_ASC 0
#define ORDER_DSC 1

#define STR_CASE_YES 1
#define STR_CASE_NO  0

#define MAX_COLS 1024

typedef int (*cmp_fun)(const void *a, const void *b);

#define CMP_STR_ASC        1
#define CMP_STR_DSC        2
#define CMP_STR_CASE_ASC   3
#define CMP_STR_CASE_DSC   4
#define CMP_NUM_ASC        5
#define CMP_NUM_DSC        6
#define CMP_MIX_ASC        7
#define CMP_MIX_DSC        8
#define CMP_MIX_CASE_ASC   9
#define CMP_MIX_CASE_DSC   10
#define CMP_ENUM_ASC       11
#define CMP_ENUM_DSC       12
#define CMP_ENUM_CASE_ASC  13
#define CMP_ENUM_CASE_DSC  14

static inline int my_strncmp(char *s1, char *s2, double l1, double l2){
	int ret;
	size_t AL1=int(l1);
	size_t AL2=int(l2);
	if(l1 == l2){
		return strncmp(s1, s2, AL1);
	} else if(l1 < l2){
		if((ret = strncmp(s1, s2, AL1)) == 0){
			return -1;
		} else return ret;
	} else {
		if((ret = strncmp(s1, s2, AL2)) == 0){
			return 1;
		} else return ret;
	}
}


static inline int my_strncasecmp(char *s1, char *s2, double l1, double l2){
	int ret;
	size_t AL1=int(l1);
	size_t AL2=int(l2);
	if(l1 == l2){
		return strncasecmp(s1, s2, AL1);
	} else if(l1 < l2){
		if((ret = strncasecmp(s1, s2, AL1)) == 0){
			return -1;
		} else return ret;
	} else {
		if((ret = strncasecmp(s1, s2, AL2)) == 0){
			return 1;
		} else return ret;
	}
}

typedef struct {
	int *col_idxs;
	int *col_orders;
	cmp_fun *cmp_funs;
	int *cmp_ops;
	int size;
	int cap;
} sort_param;

typedef struct {
	double num_val;
	char *str_val;
} sortkey_t;

typedef struct {
	int line_num;
	char *str;
	int str_len;
} line_t;

typedef struct {
	char   *text;
	long long file_size;
	long long file_cap;
	line_t *lines;
	int line_size;
	int line_cap;
	sortkey_t  *keys;
	int key_size;
	int key_cap;
	sort_param param;
	int col_map[MAX_COLS];
} file_t;

int cmp_str_asc(const void *key1, const void *key2);
int cmp_str_dsc(const void *key1, const void *key2);
int cmp_str_case_asc(const void *key1, const void *key2);
int cmp_str_case_dsc(const void *key1, const void *key2);
int cmp_num_asc(const void *key1, const void *key2);
int cmp_num_dsc(const void *key1, const void *key2);
int cmp_mix_asc(const void *key1, const void *key2);
int cmp_mix_dsc(const void *key1, const void *key2);
int cmp_mix_case_asc(const void *key1, const void *key2);
int cmp_mix_case_dsc(const void *key1, const void *key2);

// static file_t *file_db = NULL;   //hewm

file_t* load_file_lines(FILE *in, int block_rows, int *line_breaker_table, int *col_breaker_table,
		int *col_idxs, int *col_starts, int *col_lends, int *col_types, hashci_t **enums, int *orders, int *str_cases, int col_size);

void sort_file_lines(file_t *db);

void print_file_lines(file_t *db);

void free_file_lines(file_t *db);

#endif
