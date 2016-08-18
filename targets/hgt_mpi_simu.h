#ifndef hgt_mpi_simu_h
#define hgt_mpi_simu_h
#include "hgt_stat.h"
int check_mpi_error_code(int error_code, char * ops_type, char * parent_func);
FILE* create_file(char* prefix, char* appdix, char *fmt);
typedef struct _filecontainer file_container;
struct _filecontainer {
	FILE *cov;
};

typedef struct _hgt_simu_results hgt_simu_results;
struct _hgt_simu_results {
    int maxl;
    hgt_stat_meanvar_list *cs;
    hgt_stat_meanvar_list *cm;
    hgt_stat_meanvar_list *cr;
    hgt_stat_meanvar_list *ks;
    hgt_stat_meanvar_list *vd;
    hgt_stat_meanvar_list *cm2;
};
hgt_simu_results* hgt_simu_results_new(int maxl);
int hgt_simu_results_free(hgt_simu_results *r, int maxl);
int hgt_simu_results_increment(hgt_simu_results *r, hgt_cov_result *cov_result);

int cov_calc(hgt_simu_results *results, hgt_pop **ps, hgt_params *params, int cluster_num, int rank, int numprocs, gsl_rng *rng);
int write_cov(FILE *fp, hgt_simu_results *results, unsigned long gen, int cluster_num);
#endif