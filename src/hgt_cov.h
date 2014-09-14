//
//  hgt_cov.h
//  hgt
//
//  Created by Mingzhi Lin on 1/20/14.
//
//

#ifndef hgt_hgt_cov_h
#define hgt_hgt_cov_h

#include <gsl/gsl_rng.h>

typedef struct{
    unsigned long maxl;
    double ks;
    double vd;
    double *scov;
    double *rcov;
    double *pxpy;
    double *tcov;
} hgt_cov_result;

typedef int(*hgt_cov_sample_func)(unsigned long *a, unsigned long *b, unsigned long *c, unsigned long *d, unsigned long p_size, const gsl_rng *r);

hgt_cov_result *hgt_cov_result_alloc(unsigned long maxl);
void hgt_cov_result_free(hgt_cov_result *result);

int hgt_cov_result_calc_matrix(hgt_cov_result *r, short **matrix, unsigned long size, unsigned long length);

int hgt_cov_sample_p2(unsigned long *a, unsigned long *b, unsigned long *c, unsigned long *d, unsigned long p_size, const gsl_rng *r);
int hgt_cov_sample_p3(unsigned long *a, unsigned long *b, unsigned long *c, unsigned long *d, unsigned long p_size, const gsl_rng *r);
int hgt_cov_sample_p4(unsigned long *a, unsigned long *b, unsigned long *c, unsigned long *d, unsigned long p_size, const gsl_rng *r);
#endif
