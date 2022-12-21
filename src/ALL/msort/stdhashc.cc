#include <stdlib.h>
#include "stdhash.hh"
#include "stdhashc.h"

// hashci

typedef hash_map_char<int> hashci_cpp_t;

hashci_t *hci_init()
{
	hashci_t *h = (hashci_t*)malloc(sizeof(hashci_t));
	h->ptr = new hashci_cpp_t;
	return h;
}
void hci_destroy(hashci_t *h)
{
	delete (hashci_cpp_t*)h->ptr;
	free(h);
}
int hci_put(hashci_t *h, const char *key, int value)
{
	return ((hashci_cpp_t*)h->ptr)->insert(key, value);
}
int hci_del(hashci_t *h, const char *key)
{
	return ((hashci_cpp_t*)h->ptr)->erase(key);
}
int hci_get(const hashci_t *h, const char *key, int *value)
{
	return ((hashci_cpp_t*)h->ptr)->find(key, value);
}
u_int32_t hci_size(hashci_t *h)
{
	return ((hashci_cpp_t*)h->ptr)->size();
}
u_int32_t hci_capacity(hashci_t *h)
{
	return ((hashci_cpp_t*)h->ptr)->capacity();
}
int hci_resize(hashci_t *h, u_int32_t new_size)
{
	return ((hashci_cpp_t*)h->ptr)->resize(new_size);
}
void hci_traverse(hashci_t *h, hashci_f func)
{
	hashci_cpp_t *hash = (hashci_cpp_t*)h->ptr;
	hashci_cpp_t::iterator iter;
	for (iter = hash->begin(); iter != hash->end(); ++iter) {
		if (iter.isfilled()) {
			func(iter.key(), iter.value());
		}
	}
}

// hashii

typedef hash_map_misc<bit32_t, int> hashii_cpp_t;

hashii_t *hii_init()
{
	hashii_t *h = (hashii_t*)malloc(sizeof(hashii_t));
	h->ptr = new hashii_cpp_t;
	return h;
}
void hii_destroy(hashii_t *h)
{
	delete (hashii_cpp_t*)h->ptr;
	free(h);
}
int hii_put(hashii_t *h, u_int32_t key, int value)
{
	return ((hashii_cpp_t*)h->ptr)->insert(key, value);
}
int hii_del(hashii_t *h, u_int32_t key)
{
	return ((hashii_cpp_t*)h->ptr)->erase(key);
}
int hii_get(const hashii_t *h, u_int32_t key, int *value)
{
	return ((hashii_cpp_t*)h->ptr)->find(key, value);
}
u_int32_t hii_size(hashii_t *h)
{
	return ((hashii_cpp_t*)h->ptr)->size();
}
u_int32_t hii_capacity(hashii_t *h)
{
	return ((hashii_cpp_t*)h->ptr)->capacity();
}
int hii_resize(hashii_t *h, u_int32_t new_size)
{
	return ((hashii_cpp_t*)h->ptr)->resize(new_size);
}
void hii_traverse(hashii_t *h, hashii_f func)
{
	hashii_cpp_t *hash = (hashii_cpp_t*)h->ptr;
	hashii_cpp_t::iterator iter;
	for (iter = hash->begin(); iter != hash->end(); ++iter) {
		if (iter.isfilled()) {
			func(iter.key(), iter.value());
		}
	}
}
