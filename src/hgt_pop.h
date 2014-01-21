//
//  hgt_pop.h
//  hgt
//
//  Created by Mingzhi Lin on 1/18/14.
//
//

#include "hgt_cov.h"

#ifndef hgt_hgt_pop_h
#define hgt_hgt_pop_h

typedef struct {
    unsigned long size;     // population size
    unsigned long seq_len;  // genome length
    char ** genomes;        // genome sequences
    
    // cache
    unsigned long * survived;
    unsigned long * new_born;
    int cache_allocated;
} hgt_pop;

typedef struct {
    // population parameters
    unsigned long size;     // population size
    unsigned long seq_len;  // genome length
    double mu_rate;         // mutation rate
    double tr_rate;         // transfer rate
    unsigned long frag_len; // fragment length
    
    // sample parameters
    unsigned long generations;          // generations
    unsigned long sample_size;  // sample size
    unsigned long sample_time;  // sample time
    unsigned long replicates;         // replicates
    unsigned long maxl;         // maxl
    
    // output parameters
    char * prefix;
    
    // fitting paramters
    int fit_range; // fitting range
    int fit_flat;  // fitting flat
} hgt_pop_params;

hgt_pop * hgt_pop_alloc(unsigned long size, unsigned long seq_len, const gsl_rng * r);
hgt_pop * hgt_pop_copy(hgt_pop * p);
int hgt_pop_free(hgt_pop * r);

int hgt_pop_mutate_at(hgt_pop *p, unsigned long g, unsigned long s, const gsl_rng *r);
int hgt_pop_mutate(hgt_pop *p, const gsl_rng *r);
int hgt_pop_transfer_at(hgt_pop *p, 
                        unsigned long donor, 
                        unsigned long reciever, 
                        unsigned long frag_len, 
                        unsigned long start);
int hgt_pop_transfer(hgt_pop *p, unsigned long frag_len, const gsl_rng *r);

typedef int (*hgt_pop_sample_func)(hgt_pop *p, const gsl_rng *r);
int hgt_pop_sample_moran(hgt_pop *p, const gsl_rng *r);
int hgt_pop_sample_wf(hgt_pop *p, const gsl_rng *r);

typedef double (*hgt_pop_coal_time_func)(unsigned long p_size, const gsl_rng *r);
double hgt_pop_coal_time_moran(unsigned long p_size, const gsl_rng *r);
double hgt_pop_coal_time_wf(unsigned long p_size, const gsl_rng *r);

int hgt_pop_evolve(hgt_pop *p, 
                   hgt_pop_params *params, 
                   hgt_pop_sample_func sample_f, 
                   hgt_pop_coal_time_func c_time_f, 
                   const gsl_rng *r);

int hgt_pop_params_parse(hgt_pop_params *params, int argc, char **argv, char * progname);
int hgt_pop_params_free(hgt_pop_params *params);

double hgt_pop_calc_ks(hgt_pop *p);
int hgt_pop_calc_dist(hgt_pop *p, double *ds1, double *ds2, unsigned long sample_size, hgt_cov_sample_func sample_func, const gsl_rng *r);

typedef int(*hgt_pop_calc_pxy_func)(double **pxy, unsigned long maxl, double *d1, double *d2, unsigned long len);
int hgt_pop_calc_pxy(double **pxy, unsigned long maxl, double *d1, double *d2, unsigned long len, int circular);
int hgt_pop_calc_pxy_fft(double **pxy, unsigned long maxl, double *d1, double *d2, unsigned long len, int circular);
#endif
