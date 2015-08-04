#include "hgt_file.h"
#include <stdio.h>
#include <stdlib.h>

FILE *create_file(char *prefix, char * appdix, char *fmt);

FILE* create_file(char* prefix, char* appdix, char *file_format)
{
	char fn[100];
	int cx = snprintf(fn, 100, "%s.%s.%s", prefix, appdix, file_format);
	if (!(cx >= 0 && cx < 100)) {
		printf("could not create file name!\n");
		exit(EXIT_FAILURE);
	}
	FILE* f;
	f = fopen(fn, "w");
	return f;
}

hgt_file_container *hgt_file_container_create(char *prefix) {
	hgt_file_container *fc;
	fc = (hgt_file_container *)malloc(sizeof(hgt_file_container));
	fc->p2 = create_file(prefix, "p2", "txt");
	fc->p3 = create_file(prefix, "p3", "txt");
	fc->p4 = create_file(prefix, "p4", "txt");
	fc->cov = create_file(prefix, "cov", "txt");
	fc->ks = create_file(prefix, "ks", "txt");
	fc->t2 = create_file(prefix, "t2", "txt");
    fc->t3 = create_file(prefix, "t3", "txt");
    fc->t4 = create_file(prefix, "t4", "txt");
    fc->q2 = create_file(prefix, "q2", "txt");
    fc->q3 = create_file(prefix, "q3", "txt");
    fc->q4 = create_file(prefix, "q4", "txt");
	return fc;
}

void hgt_file_container_close(hgt_file_container *fc) {
	fclose(fc->p2);
	fclose(fc->p3);
	fclose(fc->p4);
	fclose(fc->cov);
	fclose(fc->ks);
	fclose(fc->t2);
    fclose(fc->t3);
    fclose(fc->t4);
    fclose(fc->q2);
    fclose(fc->q3);
    fclose(fc->q4);
}

void hgt_file_container_destroy(hgt_file_container *fc) {
	free(fc);
}

void hgt_file_container_flush(hgt_file_container *fc) {
	fflush(fc->p2);
	fflush(fc->p3);
	fflush(fc->p4);
	fflush(fc->cov);
	fflush(fc->ks);
	fflush(fc->t2);
    fflush(fc->t3);
    fflush(fc->t4);
    fflush(fc->q2);
    fflush(fc->q3);
    fflush(fc->q4);
}

void hgt_file_container_write_headers(hgt_file_container *fc) {
	fprintf(fc->p2, "#l\tp00\tp01\tp10\tp11\tp00 var\tp01 var\tp10 var\tp11 var\tsample n\tgenerations\n");
	fprintf(fc->p3, "#l\tp00\tp01\tp10\tp11\tp00 var\tp01 var\tp10 var\tp11 var\tsample n\tgenerations\n");
	fprintf(fc->p4, "#l\tp00\tp01\tp10\tp11\tp00 var\tp01 var\tp10 var\tp11 var\tsample n\tgenerations\n");
	fprintf(fc->cov, "#l\tscov\trcov\tpxpy\ttcov\tscov var\trcov var\tpxpy var\ttcov var\tsample n\tgenerations\n");
	fprintf(fc->ks, "#ks\tvd\tvar ks\tvar vd\tgenerations\n");
	fprintf(fc->t2, "#l\tt\tgeneration\n");
    fprintf(fc->t3, "#l\tt\tgeneration\n");
    fprintf(fc->t4, "#l\tt\tgeneration\n");
    fprintf(fc->q2, "#l\tt\tgeneration\n");
    fprintf(fc->q3, "#l\tt\tgeneration\n");
    fprintf(fc->q4, "#l\tt\tgeneration\n");
}