/*
   stdhash -- standalone hash library

   Copyright (c) 2006, Heng Li <lh3lh3@gmail.com>

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/
/**
   Standalone hash library is a free C++ template library for hash tables.
   It implements open-address hashing with the "double hashing" technique,
   which makes this library different from most of other implementations,
   such as STL and TR1 hash based classes, which are all based on chained
   hashing. However, fewer implementations do not mean open-address hashing
   is less efficient. As a matter of fact, this hash library is at least as
   efficient, in terms of both time and space, as STL and TR1 libraries, and
   in some cases outperforms others.
 */

/*
  2007-07-18, 1.9.2: fixed a bug in hash_*_char when an element is erased.
 */

#ifndef LH3_STDHASH_H_
#define LH3_STDHASH_H_

#include <string.h>
#include <stdlib.h>

#define LH3_STDHASH_VERSION "1.9.2"

typedef unsigned int bit32_t;
typedef unsigned long long bit64_t;
typedef unsigned short bit16_t;
// even on 64-bit systems, hashint_t should be bit32_t
typedef bit32_t hashint_t;

/*
 * Hash functions
 */
// Do a web search "g_str_hash X31_HASH" for more information.
/** hash function for strings (char*) */
inline bit32_t __lh3_X31_hash_string(const char *s)
{
	bit32_t h = 0;
	for ( ; *s; s++)
		h = (h << 5) - h + *s;
	return h;
}
/** Jenkins' hash function for 32-bit integers. Not used in this library. */
inline bit32_t __lh3_Jenkins_hash_int(bit32_t key)
{
	key += (key << 12);
	key ^= (key >> 22);
	key += (key << 4);
	key ^= (key >> 9);
	key += (key << 10);
	key ^= (key >> 2);
	key += (key << 7);
	key ^= (key >> 12);
	return key;
}
/** Jenkins' hash function for 64-bit integers. Used when LH3_HASH_INT macro is set. */
inline bit64_t __lh3_Jenkins_hash_64(bit64_t key)
{
	key += ~(key << 32);
	key ^= (key >> 22);
	key += ~(key << 13);
	key ^= (key >> 8);
	key += (key << 3);
	key ^= (key >> 15);
	key += ~(key << 27);
	key ^= (key >> 31);
	return key;
}
/** Wang's hash function for 32-bit inegers. Used when LH3_HASH_INT macro is set. */
inline bit32_t __lh3_Wang_hash_int(bit32_t key)
{
	key += ~(key << 15);
	key ^=  (key >> 10);
	key +=  (key << 3);
	key ^=  (key >> 6);
	key += ~(key << 11);
	key ^=  (key >> 16);
	return key;
}
/** hash function for 16-bit integers */
inline bit32_t __lh3_hash_fun(bit16_t key)
{
#ifndef LH3_HASH_INT
	return __lh3_Wang_hash_int(bit32_t(key));
#else
	return bit32_t(key);
#endif
}
/** hash function for 32-bit integers */
inline bit32_t __lh3_hash_fun(bit32_t key)
{
#ifndef LH3_HASH_INT
	return __lh3_Wang_hash_int(key);
#else
	return key;
#endif
}
/** hash function for 64-bit integers */
inline bit32_t __lh3_hash_fun(bit64_t key)
{
#ifdef LH3_HASH_INT
	return bit32_t(__lh3_Jenkins_hash_64(key));
#else
	return bit32_t(key>>16) ^ bit32_t(key);
#endif
}
/** hash function for strings (char*) */
inline bit32_t __lh3_hash_fun(const char *key)
{
	return __lh3_X31_hash_string(key);
}

/*
 * Equality operator "=="
 */
/** "==" for 16-bit integers */
inline bool __lh3_key_equal(bit16_t a, bit16_t b)
{
	return a == b;
}
/** "==" for 32-bit integers */
inline bool __lh3_key_equal(bit32_t a, bit32_t b)
{
	return a == b;
}
/** "==" for 64-bit integers */
inline bool __lh3_key_equal(bit64_t a, bit64_t b)
{
	return a == b;
}
/** "==" for strings (char*) */
inline bool __lh3_key_equal(const char *a, const char *b)
{
	return (a && b && strcmp(a, b) == 0);
}

/*
 * Table for primes
 */
const int __lh3_HASH_PRIME_SIZE = 32;
static const bit32_t __lh3_prime_list[__lh3_HASH_PRIME_SIZE] =
{
  0ul,          3ul,          11ul,         23ul,         53ul,
  97ul,         193ul,        389ul,        769ul,        1543ul,
  3079ul,       6151ul,       12289ul,      24593ul,      49157ul,
  98317ul,      196613ul,     393241ul,     786433ul,     1572869ul,
  3145739ul,    6291469ul,    12582917ul,   25165843ul,   50331653ul,
  100663319ul,  201326611ul,  402653189ul,  805306457ul,  1610612741ul,
  3221225473ul, 4294967291ul
};

/** Threshold for rehashing */
const double __lh3_HASH_UPPER = 0.77;

/*
 * Constants and macros for retrieve/set flags "isempty" and "isdel".
 */
typedef bit32_t __lh3_flag_t;
const int __lh3_FLAG_SHIFT = 4;
const int __lh3_FLAG_MASK = 0xful;
const __lh3_flag_t __lh3_FLAG_DEFAULT = 0xaaaaaaaaul;

#define __lh3_isempty(flag, i) ((flag[i>>__lh3_FLAG_SHIFT]>>((i&__lh3_FLAG_MASK)<<1))&2)
#define __lh3_isdel(flag, i) ((flag[i>>__lh3_FLAG_SHIFT]>>((i&__lh3_FLAG_MASK)<<1))&1)
#define __lh3_isboth(flag, i) ((flag[i>>__lh3_FLAG_SHIFT]>>((i&__lh3_FLAG_MASK)<<1))&3)
#define __lh3_set_isdel_false(flag, i) (flag[i>>__lh3_FLAG_SHIFT]&=~(1ul<<((i&__lh3_FLAG_MASK)<<1)))
#define __lh3_set_isempty_false(flag, i) (flag[i>>__lh3_FLAG_SHIFT]&=~(2ul<<((i&__lh3_FLAG_MASK)<<1)))
#define __lh3_set_isboth_false(flag, i) (flag[i>>__lh3_FLAG_SHIFT]&=~(3ul<<((i&__lh3_FLAG_MASK)<<1)))
#define __lh3_set_isdel_true(flag, i) (flag[i>>__lh3_FLAG_SHIFT]|=1ul<<((i&__lh3_FLAG_MASK)<<1))

/*
 * Auxiliary functions for search/insert/erase.
 */
template <class keytype_t>
inline hashint_t __lh3_hash_search_aux(const keytype_t &key, hashint_t m, const keytype_t *keys, const __lh3_flag_t *flag)
{
	if (!m) return 0;
	hashint_t inc, k, i;
	k = __lh3_hash_fun(key);
	i = k % m;
	inc = 1 + k % (m - 1);
	hashint_t last = i;
	while (!__lh3_isempty(flag, i) && !__lh3_key_equal(keys[i], key)) {
		if (i + inc >= m) i = i + inc - m; // inc < m, and so never write this line as: "i += inc - m;"
		else i += inc;
		if (i == last) return m; // fail to find
	}
	return i;
}
template <class keytype_t>
inline hashint_t __lh3_hash_insert_aux(const keytype_t &key, hashint_t m, const keytype_t *keys, const __lh3_flag_t *flag)
{
	hashint_t inc, k, i, site;
	site = m;
	k = __lh3_hash_fun(key);
	i = k % m;
	inc = 1 + k % (m - 1);
	
	hashint_t last = i;
	while (!__lh3_isempty(flag, i) && !__lh3_key_equal(keys[i], key)) {
		if (__lh3_isdel(flag, i)) site = i;
		if (i + inc >= m) i = i + inc - m;
		else i += inc;
		if (i == last) return site;
	}
	if (__lh3_isempty(flag, i) && site != m) return site;
	else return i;
}
template <class keytype_t>
inline hashint_t __lh3_hash_erase_aux(const keytype_t &key, hashint_t m, const keytype_t *keys, __lh3_flag_t *flag)
{
	if (!m) return 0;
	hashint_t i;
	i = __lh3_hash_search_aux(key, m, keys, flag);
	if (i != m && !__lh3_isempty(flag, i)) {
		if (__lh3_isdel(flag, i)) return m; // has been deleted
		__lh3_set_isdel_true(flag, i); // set "isdel" flag as "true"
		return i;
	} else return m;
}

/** "iterator" class for "hash_set_char" and "hash_set_misc" */
template <class keytype_t>
class __lh3_hash_base_iterator
{
protected:
	hashint_t i;
	keytype_t *keys;
	__lh3_flag_t *flags;
public:
	__lh3_hash_base_iterator() {} // No initialization. This is unsafe, but reasonable use will not cause any problems.
	__lh3_hash_base_iterator(hashint_t _i, keytype_t *_keys, __lh3_flag_t *_flags) {
		i = _i; keys = _keys; flags = _flags;
	}
	inline const keytype_t &operator & () { return keys[i]; } // Keys should never be changed by an iterator.
	inline const keytype_t &key() { return keys[i]; } // an alias of the operator "&"
	inline bool operator != (const __lh3_hash_base_iterator &iter) { return i != iter.i; }
	inline bool operator == (const __lh3_hash_base_iterator &iter) { return i == iter.i; }
	inline bool operator < (const __lh3_hash_base_iterator &iter) { return i < iter.i; }
	inline bool operator > (const __lh3_hash_base_iterator &iter) { return i > iter.i; }
	inline void operator ++ () { ++i; }
	inline void operator ++ (int) { ++i; }
	inline void operator -- () { --i; }
	inline void operator -- (int) { --i; }
	inline bool isfilled() { return !__lh3_isboth(flags, i); }
	inline bool operator + () { return isfilled(); } // an alias of "isfilled()"
};

/** "iterator" class for "hash_map_char" and "hash_map_misc" */
template <class keytype_t, class valtype_t>
class __lh3_hash_val_iterator : public __lh3_hash_base_iterator<keytype_t>
{
protected:
	valtype_t *vals;
public:
	__lh3_hash_val_iterator() {}
	__lh3_hash_val_iterator(hashint_t _i, keytype_t *_keys, __lh3_flag_t *_flags, valtype_t *_vals) {
		this->i = _i; this->keys = _keys; this->flags = _flags; vals = _vals;
	}
	inline valtype_t &operator * () { return vals[this->i]; } // Values can be changed here.
	inline const valtype_t &value() { return vals[this->i]; } // the following two functions are alternatives to the operator "*".
	inline void value(const valtype_t &v) { vals[this->i] = v; }
};

/** Base class of all hash classes */
template <class keytype_t>
class __lh3_hash_base_class
{
protected:
	hashint_t n_capacity; /**< maximum size of the hash table */
	hashint_t n_size; /**< number of elements in hash table */
	hashint_t n_occupied; /**< number of cells that have not been flaged as "isempty" (n_capacity >= n_occupied >= n_size) */
	hashint_t upper_bound; /**< The upper bound. When n_occupied exceeds this, rehashing will be performed. */
	__lh3_flag_t *flags; /**< flag array which stores the status "isempty" or "isdel" of each hash cell. */
	keytype_t *keys; /**< array that stores hash keys */
	// return 0 for unchanged, 1 for empty, 2 for deleted
	inline int direct_insert_aux(const keytype_t &key, hashint_t m, keytype_t *K, __lh3_flag_t *F, hashint_t *i) {
		*i = __lh3_hash_insert_aux(key, m, K, F);
		if (__lh3_isempty(F, *i)) {
			K[*i] = key;
			__lh3_set_isboth_false(F, *i);
			return 1;
		} else if (__lh3_isdel(F, *i)) {
			K[*i] = key;
			__lh3_set_isboth_false(F, *i);
			return 2;
		} else return 0;
	}
	inline bool resize_aux1(hashint_t *new_capacity, __lh3_flag_t **new_flags) {
		hashint_t t;
		t = __lh3_HASH_PRIME_SIZE - 1;
		while (__lh3_prime_list[t] > *new_capacity) --t;
		*new_capacity = __lh3_prime_list[t+1];
		if (n_size >= hashint_t(*new_capacity * __lh3_HASH_UPPER + 0.5)) return false; // do not rehash
		keys = (keytype_t*)realloc(keys, *new_capacity * sizeof(keytype_t));
		if (keys == 0) return false; // insufficient memory?
		*new_flags = (__lh3_flag_t*)malloc(*new_capacity * sizeof(__lh3_flag_t));
		if (*new_flags == 0) { // insufficient memory?
			::free(*new_flags); return false;
		}
		for (t = 0; t < ((*new_capacity>>__lh3_FLAG_SHIFT) + 1); ++t)
			(*new_flags)[t] = __lh3_FLAG_DEFAULT;
		return true;
	}
	inline void resize_aux2(hashint_t new_capacity, __lh3_flag_t *new_flags) {
		::free(flags);
		flags = new_flags;
		n_capacity = new_capacity;
		n_occupied = n_size;
		upper_bound = hashint_t(n_capacity * __lh3_HASH_UPPER + 0.5);
	}
	/** Test whether rehashing is needed and perform rehashing if this is the fact. */
	inline void rehash() {
		if (n_occupied >= upper_bound) {
			if (n_capacity > (n_size<<1)) resize(n_capacity - 1); // do not enlarge
			else resize(n_capacity + 1); // enlarge the capacity
		}
	}
public:
	typedef __lh3_hash_base_iterator<keytype_t> iterator;
	__lh3_hash_base_class(void) {
		keys = 0; flags = 0;
		n_capacity = n_size = n_occupied = upper_bound = 0;;
	}
	~__lh3_hash_base_class(void) { ::free(keys); ::free(flags); }
	/** resize the hash table and perform rehashing */
	inline bool resize(hashint_t new_capacity) {
		__lh3_flag_t *new_flags;
		if (!resize_aux1(&new_capacity, &new_flags)) return false;
		for (hashint_t j = 0; j != n_capacity; ++j) {
			if (__lh3_isboth(flags, j) == 0) {
				keytype_t key = keys[j]; // take out the key
				__lh3_set_isdel_true(flags, j); // mark "deleted"
				while (1) {
					hashint_t inc, k, i;
					k = __lh3_hash_fun(key);
					i = k % new_capacity; // calculate the new position
					inc = 1 + k % (new_capacity - 1);
					while (!__lh3_isempty(new_flags, i)) {
						if (i + inc >= new_capacity) i = i + inc - new_capacity;
						else i += inc;
					}
					__lh3_set_isempty_false(new_flags, i);
					if (i < this->n_capacity && __lh3_isboth(flags, i) == 0) { // something is here
						{ keytype_t tmp = keys[i]; keys[i] = key; key = tmp; } // take it out
						__lh3_set_isdel_true(flags, i);
					} else { // put key and quit the loop
						keys[i] = key;
						break;
					}
				}
			}
		}
		resize_aux2(new_capacity, new_flags);
		return true;
	}
	/** get n_size */
	inline hashint_t size(void) const { return n_size; };
	/** get n_capacity */
	inline hashint_t capacity(void) const { return n_capacity; };
	/** the first iterator */
	inline iterator begin() { return iterator(0, keys, flags); }
	/** the last iterator */
	inline iterator end() { return iterator(n_capacity, keys, flags); }
	/** clear the hash table, but do not free the memory */
	inline void clear(void)	{
		if (flags) {
			for (hashint_t t = 0; t < ((n_capacity>>__lh3_FLAG_SHIFT) + 1); ++t)
				flags[t] = __lh3_FLAG_DEFAULT;
		}
		n_size = 0;
	}
	/** clear the hash table and free the memory */
	inline void free() {
		::free(keys); ::free(flags);
		keys = 0; flags = 0;
		n_capacity = n_size = n_occupied = upper_bound = 0;;
	}
};

/** hash_set_misc class */
template <class keytype_t>
class hash_set_misc : public __lh3_hash_base_class<keytype_t>
{
public:
	hash_set_misc(void) {};
	~hash_set_misc(void) {};
	/** search a key */
	inline bool find(const keytype_t &key) const {
		hashint_t i = __lh3_hash_search_aux(key, this->n_capacity, this->keys, this->flags);
		return (i != this->n_capacity && __lh3_isboth(this->flags, i) == 0)? true : false;
	}
	/** insert a key */
	inline bool insert(const keytype_t &key) {
		__lh3_hash_base_class<keytype_t>::rehash();
		hashint_t i;
		int ret = direct_insert_aux(key, this->n_capacity, this->keys, this->flags, &i);
		if (ret == 0) return true;
		if (ret == 1) { ++(this->n_size); ++(this->n_occupied); }
		else ++(this->n_size); // then ret == 2
		return false;
	}
	/** delete a key */
	inline bool erase(const keytype_t &key) {
		hashint_t i = __lh3_hash_erase_aux(key, this->n_capacity, this->keys, this->flags);
		if (i != this->n_capacity) {
			--(this->n_size);
			return true;
		} else return false;
	}
};

/** hash_map_misc class */
template <class keytype_t, class valtype_t>
class hash_map_misc : public hash_set_misc<keytype_t>
{
	valtype_t *vals;
	/** a copy of __lh3_hash_base_class<keytype_t>::rehash() */
	inline void rehash() {
		if (this->n_occupied >= this->upper_bound) {
			if (this->n_capacity > (this->n_size<<1)) resize(this->n_capacity - 1);
			else resize(this->n_capacity + 1);
		}
	}
public:
	hash_map_misc(void) { vals = 0; };
	~hash_map_misc(void) { ::free(vals); };
	typedef __lh3_hash_val_iterator<keytype_t, valtype_t> iterator;
	/** analogy of __lh3_hash_base_class<keytype_t>::resize(hashint_t) */
	inline bool resize(hashint_t new_capacity) {
		__lh3_flag_t *new_flags;
		if (!__lh3_hash_base_class<keytype_t>::resize_aux1(&new_capacity, &new_flags)) return false;
		vals = (valtype_t*)realloc(vals, sizeof(valtype_t) * new_capacity);
		if (vals == 0) { // insufficient enough memory?
			::free(new_flags);
			return false;
		}
		for (hashint_t j = 0; j != this->n_capacity; ++j) {
			if (__lh3_isboth(this->flags, j) == 0) {
				keytype_t key = this->keys[j]; // take out the key
				valtype_t val = vals[j];
				__lh3_set_isdel_true(this->flags, j); // mark "deleted"
				while (1) {
					hashint_t inc, k, i;
					k = __lh3_hash_fun(key);
					i = k % new_capacity; // calculate the new position
					inc = 1 + k % (new_capacity - 1);
					while (!__lh3_isempty(new_flags, i)) {
						if (i + inc >= new_capacity) i = i + inc - new_capacity;
						else i += inc;
					}
					__lh3_set_isempty_false(new_flags, i);
					if (i < this->n_capacity && __lh3_isboth(this->flags, i) == 0) { // something is here
						{ keytype_t tmp = this->keys[i]; this->keys[i] = key; key = tmp; } // take it out
						{ valtype_t tmp = vals[i]; vals[i] = val; val = tmp; } // take it out
						__lh3_set_isdel_true(this->flags, i);
					} else { // clear
						this->keys[i] = key;
						vals[i] = val;
						break;
					}
				}
			}
		}
		__lh3_hash_base_class<keytype_t>::resize_aux2(new_capacity, new_flags);
		return true;
	}
	inline bool find(const keytype_t &key, valtype_t *q) const {
		hashint_t i = __lh3_hash_search_aux(key, this->n_capacity, this->keys, this->flags);
		if (i != this->n_capacity && __lh3_isboth(this->flags, i) == 0) {
			*q = vals[i];
			return true;
		} else return false;
	}
	inline bool insert(const keytype_t &key, const valtype_t &val) {
		rehash();
		hashint_t i;
		int ret =this->direct_insert_aux(key, this->n_capacity, this->keys, this->flags, &i);
//int ret =direct_insert_aux(key, this->n_capacity, this->keys, this->flags, &i);
		vals[i] = val;
		if (ret == 0) return true;
		if (ret == 1) { ++(this->n_size); ++(this->n_occupied); }
		else ++(this->n_size); // then ret == 2
		return false;
	}
	inline bool insert(const keytype_t &key, valtype_t **q) {
		rehash();
		hashint_t i;
		int ret = direct_insert_aux(key, this->n_capacity, this->keys, this->flags, &i);
		*q = vals + i;
		if (ret == 0) return true;
		if (ret == 1) { ++(this->n_size); ++(this->n_occupied); }
		else ++(this->n_size); // then ret == 2
		return false;
	}
	inline bool erase(const keytype_t &key) { return hash_set_misc<keytype_t>::erase(key); }
	inline bool erase(const keytype_t &key, valtype_t **q) {
		hashint_t i = __lh3_hash_erase_aux(key, this->n_capacity, this->keys, this->flags);
		if (i != this->n_capacity) {
			--(this->n_size);
			*q = vals + i;
			return true;
		} else return false;
	}
	inline iterator begin() { return iterator(0, this->keys, this->flags, vals); }
	inline iterator end() { return iterator(this->n_capacity, this->keys, this->flags, vals); }
	inline void free() {
		hash_set_misc<keytype_t>::free(); ::free(vals); vals = 0;
	}
};

typedef char *__lh3_char_t;

/** hash_set_char class */
class hash_set_char : public __lh3_hash_base_class<char*>
{
protected:
	inline int insert_aux(const char *key, hashint_t m, char **K, __lh3_flag_t *F, hashint_t *i) {
		*i = __lh3_hash_insert_aux((const __lh3_char_t)key, m, K, F);
		if (__lh3_isempty(F, *i)) {
			K[*i] = strdup(key);
			__lh3_set_isboth_false(F, *i);
			return 1;
		} else if (__lh3_isdel(F, *i)) {
			K[*i] = strdup(key);
			__lh3_set_isboth_false(F, *i);
			return 2;
		} else return 0;
	}
public:
	hash_set_char(void) {};
	~hash_set_char(void) { clear(); };
	inline bool find(const char *key) const {
		hashint_t i = __lh3_hash_search_aux((const __lh3_char_t)key, this->n_capacity, this->keys, this->flags);
		return (i != this->n_capacity && __lh3_isboth(this->flags, i) == 0)? true : false;
	}
	inline bool insert(const char *key) {
		rehash();
		hashint_t i;
		int ret = insert_aux(key, this->n_capacity, this->keys, this->flags, &i);
		if (ret == 0) return true;
		if (ret == 1) { ++(this->n_size); ++(this->n_occupied); }
		else ++(this->n_size); // then ret == 2
		return false;
	}
	inline bool erase(const char *key) {
		hashint_t i = __lh3_hash_erase_aux((const __lh3_char_t)key, this->n_capacity, this->keys, this->flags);
		if (i != this->n_capacity) {
			::free(this->keys[i]); this->keys[i] = 0;
			--(this->n_size);
			return true;
		} else return false;
	}
	inline void clear(void) {
		for (hashint_t i = 0; i != this->n_capacity; ++i)
			if (!__lh3_isboth(this->flags, i)) {
				::free(this->keys[i]); this->keys[i] = 0;
			}
		__lh3_hash_base_class<char*>::clear();	
	}
	inline void free() {
		clear();
		__lh3_hash_base_class<char*>::free();
	}
};

/** hash_map_char class */
template <class valtype_t>
class hash_map_char : public hash_set_char
{
	valtype_t *vals;
	inline void rehash() {
		if (this->n_occupied >= this->upper_bound) {
			if (this->n_capacity > (this->n_size<<1)) resize(this->n_capacity - 1);
			else resize(this->n_capacity + 1);
		}
	}
public:
	typedef __lh3_hash_val_iterator<char*, valtype_t> iterator;
	hash_map_char(void) { vals = 0; };
	~hash_map_char(void) { clear(); ::free(vals); };
	inline valtype_t &val(hashint_t i) { return vals[i]; }
	inline void val(hashint_t i, const valtype_t &v) { return vals[i] = v; }
	inline bool resize(hashint_t new_capacity) {
		__lh3_flag_t *new_flags;
		if (!resize_aux1(&new_capacity, &new_flags)) return false;
		vals = (valtype_t*)realloc(vals, sizeof(valtype_t) * new_capacity);
		if (vals == 0) { // insufficient enough memory?
			::free(new_flags);
			return false;
		}
		for (hashint_t j = 0; j != this->n_capacity; ++j) {
			if (__lh3_isboth(this->flags, j) == 0) {
				char *key = this->keys[j]; // take out the key
				valtype_t val = vals[j];
				__lh3_set_isdel_true(this->flags, j); // mark "deleted"
				while (1) {
					hashint_t inc, k, i;
					k = __lh3_hash_fun(key);
					i = k % new_capacity; // calculate the new position
					inc = 1 + k % (new_capacity - 1);
					while (!__lh3_isempty(new_flags, i)) {
						if (i + inc >= new_capacity) i = i + inc - new_capacity;
						else i += inc;
					}
					__lh3_set_isempty_false(new_flags, i);
					if (i < this->n_capacity && __lh3_isboth(this->flags, i) == 0) { // something is here
						{ char* tmp = this->keys[i]; this->keys[i] = key; key = tmp; } // take it out
						{ valtype_t tmp = vals[i]; vals[i] = val; val = tmp; } // take it out
						__lh3_set_isdel_true(this->flags, i);
					} else { // clear
						this->keys[i] = key;
						vals[i] = val;
						break;
					}
				}
			}
		}
		resize_aux2(new_capacity, new_flags);
		return true;
	}
	inline bool find(const char *key, valtype_t *q) const {
		hashint_t i = __lh3_hash_search_aux((const __lh3_char_t)key, this->n_capacity, this->keys, this->flags);
		if (i != this->n_capacity && __lh3_isboth(this->flags, i) == 0) {
			*q = vals[i];
			return true;
		} else return false;
	}
	inline bool insert(const char *key, const valtype_t &val) {
		rehash();
		hashint_t i;
		int ret = insert_aux(key, this->n_capacity, this->keys, this->flags, &i);
		vals[i] = val;
		if (ret == 0) return true;
		if (ret == 1) { ++(this->n_size); ++(this->n_occupied); }
		else ++(this->n_size); // then ret == 2
		return false;
	}
	inline bool insert(const char *key, valtype_t **q) {
		rehash();
		hashint_t i;
		int ret = insert_aux(key, this->n_capacity, this->keys, this->flags, &i);
		*q = vals + i;
		if (ret == 0) return true;
		if (ret == 1) { ++(this->n_size); ++(this->n_occupied); }
		else ++(this->n_size); // then ret == 2
		return false;
	}
	inline bool erase(const char *key) { return hash_set_char::erase(key); }
	inline bool erase(const char *key, valtype_t **q) {
		hashint_t i = __lh3_hash_erase_aux((const __lh3_char_t)key, this->n_capacity, this->keys, this->flags);
		if (i != this->n_capacity) {
			::free(this->keys[i]); this->keys[i] = 0;
			--(this->n_size);
			*q = vals + i;
			return true;
		} else return false;
	}
	inline iterator begin() { return iterator(0, this->keys, this->flags, vals); }
	inline iterator end() { return iterator(this->n_capacity, this->keys, this->flags, vals); }
	inline void free() {
		hash_set_char::free(); ::free(vals); vals = 0;
	}
};

#endif // LH3_STDHASH_H_
