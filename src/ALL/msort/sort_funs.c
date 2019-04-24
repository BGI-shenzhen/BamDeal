#include "msort.h"

int cmp_str_asc(const void *e1, const void *e2){
	sortkey_t *key1, *key2;
	key1  = (sortkey_t*)e1;
	key2  = (sortkey_t*)e2;
	return my_strncmp(key1->str_val, key2->str_val, key1->num_val, key2->num_val);
}

int cmp_str_case_asc(const void *e1, const void *e2){
	sortkey_t *key1, *key2;
	key1  = (sortkey_t*)e1;
	key2  = (sortkey_t*)e2;
	return my_strncasecmp(key1->str_val, key2->str_val, key1->num_val, key2->num_val);
}

int cmp_str_dsc(const void *e1, const void *e2){
	sortkey_t *key1, *key2;
	key1  = (sortkey_t*)e1;
	key2  = (sortkey_t*)e2;
	return my_strncmp(key2->str_val, key1->str_val, key2->num_val, key1->num_val);
}

int cmp_str_case_dsc(const void *e1, const void *e2){
	sortkey_t *key1, *key2;
	key1  = (sortkey_t*)e1;
	key2  = (sortkey_t*)e2;
	return my_strncasecmp(key2->str_val, key1->str_val, key2->num_val, key1->num_val);
}

int cmp_num_asc(const void *e1, const void *e2){
	sortkey_t *key1, *key2;
	key1  = (sortkey_t*)e1;
	key2  = (sortkey_t*)e2;
	if(key1->num_val == key2->num_val){
		return 0;
	} else if(key1->num_val < key2->num_val){
		return -1;
	} else {
		return 1;
	}
}

int cmp_num_dsc(const void *e1, const void *e2){
	sortkey_t *key1, *key2;
	key1  = (sortkey_t*)e1;
	key2  = (sortkey_t*)e2;
	if(key1->num_val == key2->num_val){
		return 0;
	} else if(key1->num_val < key2->num_val){
		return 1;
	} else {
		return -1;
	}
}

int cmp_mix_asc(const void *e1, const void *e2){
	sortkey_t *key1, *key2;
	key1  = (sortkey_t*)e1;
	key2  = (sortkey_t*)e2;
	if(key1->str_val){
		if(key2->str_val){
			return my_strncmp(key1->str_val, key2->str_val, key1->num_val, key2->num_val);
		} else {
			return 1;
		}
	} else {
		if(key2->str_val){
			return -1;
		} else {
			if(key1->num_val == key2->num_val){
				return 0;
			} else if(key1->num_val < key2->num_val){
				return -1;
			} else {
				return 1;
			}
		}
	}
}

int cmp_mix_dsc(const void *e1, const void *e2){
	sortkey_t *key1, *key2;
	key1  = (sortkey_t*)e1;
	key2  = (sortkey_t*)e2;
	if(key1->str_val){
		if(key2->str_val){
			return my_strncmp(key2->str_val, key1->str_val, key2->num_val, key1->num_val);
		} else {
			return -1;
		}
	} else {
		if(key2->str_val){
			return 1;
		} else {
			if(key1->num_val == key2->num_val){
				return 0;
			} else if(key1->num_val < key2->num_val){
				return 1;
			} else {
				return -1;
			}
		}
	}
}

int cmp_mix_case_asc(const void *e1, const void *e2){
	sortkey_t *key1, *key2;
	key1  = (sortkey_t*)e1;
	key2  = (sortkey_t*)e2;
	if(key1->str_val){
		if(key2->str_val){
			return my_strncasecmp(key1->str_val, key2->str_val, key1->num_val, key2->num_val);
		} else {
			return 1;
		}
	} else {
		if(key2->str_val){
			return -1;
		} else {
			if(key1->num_val == key2->num_val){
				return 0;
			} else if(key1->num_val < key2->num_val){
				return -1;
			} else {
				return 1;
			}
		}
	}
}

int cmp_mix_case_dsc(const void *e1, const void *e2){
	sortkey_t *key1, *key2;
	key1  = (sortkey_t*)e1;
	key2  = (sortkey_t*)e2;
	if(key1->str_val){
		if(key2->str_val){
			return my_strncasecmp(key2->str_val, key1->str_val, key2->num_val, key1->num_val);
		} else {
			return -1;
		}
	} else {
		if(key2->str_val){
			return 1;
		} else {
			if(key1->num_val == key2->num_val){
				return 0;
			} else if(key1->num_val < key2->num_val){
				return 1;
			} else {
				return -1;
			}
		}
	}
}

