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
#include <gsl/gsl_rng.h>

typedef struct {
    unsigned long size;     // population size
    unsigned long seq_len;  // genome length
    unsigned long generation;
    char ** genomes;        // genome sequences
    
    int ** transfer_hotspots; // transfer hotspots
    
    // cache
    unsigned long * survived;
    unsigned long * new_born;
    int cache_allocated;
} hgt_pop;

typedef struct {
    // population parameters
    unsigned long size;     // population size
    unsigned long seq_len;  // genome length
    
    // transfer parameters
    double tr_rate;                             // transfer rate
    unsigned long frag_len;                     // fragment length
    unsigned int tr_hotspot_num;                // number of hotspot in the genome
    unsigned long tr_hotspot_length;            // length of hotspot
    double tr_hotspot_ratio;                    // ratio of transfer rate
    unsigned long **tr_hotspots;                // locations of transfer hotspots
    
    // mutation parameters
    double mu_rate;                     // mutation rate
    unsigned int mu_hotspot_num;        // number of mutation hotspots
    unsigned long mu_hotspot_length;    // average hotspot length
    double mu_hotspot_ratio;            // ratio of mutation rate
    unsigned long **mu_hotspots;        // locations of mutation hotspots
    
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

hgt_pop * hgt_pop_alloc(hgt_pop_params *params, const gsl_rng * r);
hgt_pop * hgt_pop_copy(hgt_pop * p);
int hgt_pop_free(hgt_pop * r);
char *hgt_pop_to_json(hgt_pop *p, hgt_pop_params *params);

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
/*
 population evolves under exponentially distributed fragment sizes
 */
int hgt_pop_evolve_expon_frag(hgt_pop *p,
                              hgt_pop_params *params,
                              hgt_pop_sample_func sample_f,
                              hgt_pop_coal_time_func c_time_f,
                              const gsl_rng *r);

int hgt_pop_params_parse(hgt_pop_params *params, int argc, char **argv, char * progname);
int hgt_pop_params_free(hgt_pop_params *params);
int hgt_pop_params_printf(hgt_pop_params *params, FILE *stream);
hgt_pop_params *hgt_pop_params_alloc();

double hgt_pop_calc_ks(hgt_pop *p);

int hgt_pop_calc_cov(hgt_cov_result *result, hgt_pop *p, int sample, const gsl_rng* rng);

int hgt_pop_calc_dist(hgt_pop *p, double *ds1, double *ds2, unsigned long sample_size, hgt_cov_sample_func sample_func, const gsl_rng *r);

typedef int(*hgt_pop_calc_pxy_func)(double *pxy, unsigned long maxl, double *d1, double *d2, unsigned long len);
int hgt_pop_calc_pxy(double *pxy, unsigned long maxl, double *d1, double *d2, unsigned long len, int circular);
int hgt_pop_calc_pxy_fft(double *pxy, unsigned long maxl, double *d1, double *d2, unsigned long len, int circular);
#endif
