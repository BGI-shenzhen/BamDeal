#ifndef STDHASHC_H
#define STDHASHC_H

#include <sys/types.h>

/* hashci */

typedef struct
{
	void *ptr;
} hashci_t;

typedef int (*hashci_f)(const char *key, int value);

#ifdef __cplusplus
extern "C" {
#endif

	hashci_t *hci_init();
	void hci_destroy(hashci_t *h);
	int hci_put(hashci_t *h, const char *key, int value);
	int hci_del(hashci_t *h, const char *key);
    int hci_get(const hashci_t *h, const char *key, int *value);
	u_int32_t hci_size(hashci_t *h);
	u_int32_t hci_capacity(hashci_t *h);
	int hci_resize(hashci_t *h, u_int32_t new_size);
	void hci_traverse(hashci_t *h, hashci_f func);

#ifdef __cplusplus
}
#endif

/* hashii */

typedef struct
{
	void *ptr;
} hashii_t;

typedef int (*hashii_f)(u_int32_t key, int value);

#ifdef __cplusplus
extern "C" {
#endif

	hashii_t *hii_init();
	void hii_destroy(hashii_t *h);
	int hii_put(hashii_t *h, u_int32_t key, int value);
	int hii_del(hashii_t *h, u_int32_t key);
    int hii_get(const hashii_t *h, u_int32_t key, int *value);
	u_int32_t hii_size(hashii_t *h);
	u_int32_t hii_capacity(hashii_t *h);
	int hii_resize(hashii_t *h, u_int32_t new_size);
	void hii_traverse(hashii_t *h, hashii_f func);

#ifdef __cplusplus
}
#endif

#endif
