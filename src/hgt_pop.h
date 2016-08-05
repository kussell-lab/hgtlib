//
//  hgt_pop.h
//  hgt
//
//  Created by Mingzhi Lin on 1/18/14.
//
//


#ifndef hgt_hgt_pop_h
#define hgt_hgt_pop_h

#include "bstrlib.h"
#include "hgt_cov.h"
#include "hgt_genome.h"
#include "hgt_linkage.h"
#include "hgt_params.h"
#include <gsl/gsl_rng.h>
#include "hgt_pop_bsc.h"
#include "hgt_pop_define.h"

#define MORAN 0;
#define WRIGHT_FISHER 1;
#define LINEAR_SELECTION 2;
#define BSC 3;
#define CONSTANT_FRAG 0;
#define EXP_FRAG 1;

int hgt_pop_prune_linkages(hgt_pop *p);
int hgt_pop_calc_fitness(hgt_pop *p, double * fitness);

int hgt_pop_calc_coal_time(
    hgt_linkage **pop_linkages,
    unsigned size, 
    unsigned int sample_size,
    double *res,
    unsigned linkage_size, 
    hgt_linkage_find_time_func find_func,
    const gsl_rng *r);
typedef double (*hgt_pop_calc_most_recent_coal_func) (hgt_linkage **pop_linkages, unsigned size, unsigned long sample_size, double * res, unsigned linkage_size, const gsl_rng *r);
double hgt_pop_calc_most_recent_coal_time(hgt_linkage **pop_linkages, unsigned size, unsigned long sample_size, double * res, unsigned linkage_size, const gsl_rng *r);
double hgt_pop_calc_most_recent_ancestor_time(hgt_linkage **pop_linkages, unsigned size, unsigned long sample_size, double * res, unsigned linkage_size, const gsl_rng *r);

hgt_pop * hgt_pop_alloc(hgt_params *params, const gsl_rng * r);
hgt_pop * hgt_pop_copy(hgt_pop * p);
int hgt_pop_free(hgt_pop * r);
char *hgt_pop_to_json(hgt_pop *p, hgt_params *params);
double hgt_pop_get_time(hgt_pop *p);

double hgt_pop_mean_fitness(hgt_pop *p);

typedef double (*hgt_pop_sample_func)(hgt_pop *p, const gsl_rng *r);
double hgt_pop_sample_moran(hgt_pop *p, const gsl_rng *r);
double hgt_pop_sample_wf(hgt_pop *p, const gsl_rng *r);
double hgt_pop_sample_linear_selection(hgt_pop *p, const gsl_rng *r);

typedef double (*hgt_pop_coal_time_func)(unsigned long p_size, const gsl_rng *r);
double hgt_pop_coal_time_moran(unsigned long p_size, const gsl_rng *r);
double hgt_pop_coal_time_wf(unsigned long p_size, const gsl_rng *r);
double hgt_pop_coal_time_linear_selection(unsigned long p_size, const gsl_rng *r);
double hgt_pop_coal_time_bsc(unsigned long p_size, const gsl_rng *r);

typedef unsigned (*hgt_pop_frag_func)(hgt_params *params, const gsl_rng *r);
unsigned hgt_pop_frag_constant(hgt_params *params, const gsl_rng *r);
unsigned hgt_pop_frag_exp(hgt_params *params, const gsl_rng *r);
unsigned hgt_pop_frag_geom(hgt_params *params, const gsl_rng *r);
int hgt_pop_evolve(hgt_pop *p, 
                   hgt_params *params, 
                   hgt_pop_sample_func sample_f,
                   hgt_pop_frag_func frag_f, 
                   const gsl_rng *r);

double hgt_pop_calc_ks(hgt_pop *p);

int hgt_pop_calc_cov(hgt_cov_result *result, hgt_pop *p, unsigned sample, const gsl_rng* rng);
int hgt_pop_calc_cov_all(hgt_cov_result *result, hgt_pop *p);
int hgt_pop_calc_t2(hgt_pop *p, unsigned long sample_size, unsigned long * res, const gsl_rng *rng);

int hgt_pop_calc_dist(hgt_pop *p, double *ds1, double *ds2, unsigned long sample_size, hgt_cov_sample_func sample_func, const gsl_rng *r);

typedef int(*hgt_pop_calc_pxy_func)(double *pxy, unsigned long maxl, double *d1, double *d2, unsigned long len);
int hgt_pop_calc_pxy(double *pxy, unsigned int maxl, double *d1, double *d2, unsigned int len, int circular);
int hgt_pop_calc_pxy_fft(double *pxy, unsigned int maxl, double *d1, double *d2, unsigned int len, int circular);
int hgt_pop_linkage_prune_p(hgt_pop *p);

int hgt_pop_mutate(hgt_pop *p, hgt_params* params, const gsl_rng* r);


#endif
