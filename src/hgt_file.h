#ifndef hgt_hgt_file_h
#define hgt_hgt_file_h
#include <stdio.h>

typedef struct _file_container hgt_file_container;
struct _file_container {
	FILE *p2;
	FILE *p3;
	FILE *p4;
	FILE *cov;
	FILE *ks;
	FILE *t2;
};

hgt_file_container *hgt_file_container_create(char *prefix);
void hgt_file_container_close(hgt_file_container *fc);
void hgt_file_container_flush(hgt_file_container *fc);
void hgt_file_container_destroy(hgt_file_container *fc);
void hgt_file_container_write_headers(hgt_file_container *fc);

#endif