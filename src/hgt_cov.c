//
//  hgt_cov.c
//  hgt
//
//  Created by Mingzhi Lin on 1/20/14.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
#include "hgt_cov.h"

hgt_cov_result *hgt_cov_result_alloc(unsigned long maxl) {
    hgt_cov_result *result = malloc(sizeof(hgt_cov_result));
    result->ks = 0;
    result->vd = 0;
    result->maxl = maxl;
    result->scov = malloc(maxl*sizeof(double));
    result->rcov = malloc(maxl*sizeof(double));
    result->pxpy = malloc(maxl*sizeof(double));
    result->tcov = malloc(maxl*sizeof(double));
    return result;
}

void hgt_cov_result_free(hgt_cov_result *r) {
    free(r->pxpy);
    free(r->scov);
    free(r->rcov);
    free(r->tcov);
    free(r);
}

// calculate covariance from a binary matrix.
int hgt_cov_result_calc_matrix(hgt_cov_result *r, short **matrix, unsigned long size, unsigned long length) {
    double * xs;
    double * xy;
    double * smx;
    double * smy;
    double * xvy;
    double * smxy;
    double * ks_arr;
    int i, j, k;
    double ks;
    
    unsigned long maxl = r->maxl;
    
    xs = calloc(length, sizeof(double));
    xy = calloc(maxl, sizeof(double));
    smx = calloc(maxl, sizeof(double));
    smy = calloc(maxl, sizeof(double));
    xvy = calloc(maxl, sizeof(double));
    smxy = calloc(maxl, sizeof(double));
    ks_arr = calloc(size, sizeof(double));
    
    for (i = 0; i < size; i ++) {
        ks = 0;
        for (j = 0; j < length; j ++) {
            if (matrix[i][j] == 1) {
                xs[j]++;
                ks ++;
                for (k = 0; k < maxl; k++) {
                    if (j+k < length) {
                        if (matrix[i][j+k] == 1) {
                            xy[k] ++;
                        }
                    } else {
                        if (matrix[i][j+k-length] == 1) {
                            xy[k] ++;
                        }
                    }
                }
            }
        }
        ks /= (double) length;
        ks_arr[i] = ks;
    }
    
    for (i = 0; i < length; i ++) {
        xs[i] /= (double) size;
    }
    
    for (i = 0; i < maxl; i++) {
        xy[i] = xy[i]/ ((double) size * (double) length);
    }
    
    for (i = 0; i < maxl; i ++) {
        for (j = 0; j < length; j ++) {
            if (i+j < length) {
                xvy[i] += xs[j] * xs[j+i];
                smx[i] += xs[j];
                smy[i] += xs[j+i];
            } else {
                xvy[i] += xs[j] * xs[j+i - length];
                smx[i] += xs[j];
                smy[i] += xs[j+i - length];
            }
        }
        xvy[i] /= (double) length;
        smx[i] /= (double) length;
        smy[i] /= (double) length;
    }
    
    for (i = 0; i < maxl; i ++) {
        smxy[i] = smx[i] * smy[i];
    }
    
    r->ks = gsl_stats_mean(ks_arr, 1, size);
    r->vd = gsl_stats_variance(ks_arr, 1, size);
    
    for (i = 0; i < maxl; i++) {
        r->scov[i] = xy[i] - xvy[i];
        r->rcov[i] = xvy[i] - smxy[i];
        r->pxpy[i] = smxy[i];
        r->tcov[i] = r->scov[i] + r->rcov[i] - r->vd;
    }
    
    free(xs);
    free(xy);
    free(smx);
    free(smy);
    free(xvy);
    free(smxy);
    free(ks_arr);
    
    return EXIT_SUCCESS;
}

int hgt_cov_sample_p2(unsigned long *a, unsigned long *b, unsigned long *c, unsigned long *d, unsigned long p_size, const gsl_rng *r) {
    *a = gsl_rng_uniform_int(r, p_size);
    *b = gsl_rng_uniform_int(r, p_size);
    while (*a == *b) {
        *b = gsl_rng_uniform_int(r, p_size);
    }
    *c = *a;
    *d = *b;
    return EXIT_SUCCESS;
}

int hgt_cov_sample_p3(unsigned long *a, unsigned long *b, unsigned long *c, unsigned long *d, unsigned long p_size, const gsl_rng *r) {
    *a = gsl_rng_uniform_int(r, p_size);
    *b = gsl_rng_uniform_int(r, p_size);
    while (*a == *b) {
        *b = gsl_rng_uniform_int(r, p_size);
    }
    *c = gsl_rng_uniform_int(r, p_size);
    while (*a == *c || *b == *c) {
        *c = gsl_rng_uniform_int(r, p_size);
    }
    *d = *a;
    return EXIT_SUCCESS;
}

int hgt_cov_sample_p4(unsigned long *a, unsigned long *b, unsigned long *c, unsigned long *d, unsigned long p_size, const gsl_rng *r) {
    *a = gsl_rng_uniform_int(r, p_size);
    *b = gsl_rng_uniform_int(r, p_size);
    while (*a == *b) {
        *b = gsl_rng_uniform_int(r, p_size);
    }
    *c = gsl_rng_uniform_int(r, p_size);
    while (*a == *c || *b == *c) {
        *c = gsl_rng_uniform_int(r, p_size);
    }
    *d = gsl_rng_uniform_int(r, p_size);
    while (*a == *d || *b == *d || *c == *d) {
        *d = gsl_rng_uniform_int(r, p_size);
    }
    return EXIT_SUCCESS;
}
