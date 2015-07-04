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

typedef struct hgt_pop_linkage hgt_pop_linkage;
typedef struct hgt_pop hgt_pop;
typedef struct hgt_pop_params hgt_pop_params;

struct hgt_pop {
    unsigned long size;     // population size
    unsigned long seq_len;  // genome length
    unsigned long generation;
    char ** genomes;        // genome sequences
    double *fitness;        // genome fitness
    hgt_pop_linkage ** linkages; // linkages.
    
    int ** transfer_hotspots; // transfer hotspots
    
    // cache
    unsigned long * survived;
    unsigned long * new_born;
    int cache_allocated;
};

struct hgt_pop_params {
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

    // fitness
    double fitness_scale;
    double fitness_shape;
    double b_mu_rate;
};

struct hgt_pop_linkage {
    int numChildren;
    unsigned long birthTime;
    hgt_pop_linkage * parent;
};

int hgt_pop_linkage_free(hgt_pop_linkage *l);
hgt_pop_linkage * hgt_pop_linkage_alloc();
typedef unsigned long hgt_pop_linkage_find_time_func(hgt_pop_linkage ** linkages, int size);
unsigned long hgt_pop_linkage_find_most_rescent_ancestor(hgt_pop_linkage ** linkages, int size);
unsigned long hgt_pop_linkage_find_most_rescent_coalescence(hgt_pop_linkage ** linkages, int size);
int hgt_pop_calc_coal_time(hgt_pop *p, 
    unsigned long sample_size, 
    unsigned long *res, 
    int linkage_size, 
    hgt_pop_linkage_find_time_func find_func,
    const gsl_rng *r);
int hgt_pop_calc_most_recent_coal_time(hgt_pop *p, unsigned long sample_size, unsigned long * res, int linkage_size, const gsl_rng *r);
int hgt_pop_calc_most_recent_ancestor_time(hgt_pop *p, unsigned long sample_size, unsigned long * res, int linkage_size, const gsl_rng *r);

hgt_pop * hgt_pop_alloc(hgt_pop_params *params, const gsl_rng * r);
hgt_pop * hgt_pop_copy(hgt_pop * p);
int hgt_pop_free(hgt_pop * r);
char *hgt_pop_to_json(hgt_pop *p, hgt_pop_params *params);

double hgt_pop_mean_fitness(hgt_pop *p);

typedef int (*hgt_pop_sample_func)(hgt_pop *p, const gsl_rng *r);
int hgt_pop_sample_moran(hgt_pop *p, const gsl_rng *r);
int hgt_pop_sample_wf(hgt_pop *p, const gsl_rng *r);

typedef double (*hgt_pop_coal_time_func)(unsigned long p_size, const gsl_rng *r);
double hgt_pop_coal_time_moran(unsigned long p_size, const gsl_rng *r);
double hgt_pop_coal_time_wf(unsigned long p_size, const gsl_rng *r);

typedef double (*hgt_pop_frag_func)(hgt_pop_params *params, const gsl_rng *r);
double hgt_pop_frag_constant(hgt_pop_params *params, const gsl_rng *r);
double hgt_pop_frag_exp(hgt_pop_params *params, const gsl_rng *r);

int hgt_pop_evolve(hgt_pop *p, 
                   hgt_pop_params *params, 
                   hgt_pop_sample_func sample_f, 
                   hgt_pop_coal_time_func c_time_f, 
                   hgt_pop_frag_func frag_f, 
                   const gsl_rng *r);

int hgt_pop_params_parse(hgt_pop_params *params, int argc, char **argv, char * progname);
int hgt_pop_params_free(hgt_pop_params *params);
int hgt_pop_params_printf(hgt_pop_params *params, FILE *stream);
hgt_pop_params *hgt_pop_params_alloc();

double hgt_pop_calc_ks(hgt_pop *p);

int hgt_pop_calc_cov(hgt_cov_result *result, hgt_pop *p, int sample, const gsl_rng* rng);
int hgt_pop_calc_cov_all(hgt_cov_result *result, hgt_pop *p);
int hgt_pop_calc_t2(hgt_pop *p, unsigned long sample_size, unsigned long * res, const gsl_rng *rng);

int hgt_pop_calc_dist(hgt_pop *p, double *ds1, double *ds2, unsigned long sample_size, hgt_cov_sample_func sample_func, const gsl_rng *r);

typedef int(*hgt_pop_calc_pxy_func)(double *pxy, unsigned long maxl, double *d1, double *d2, unsigned long len);
int hgt_pop_calc_pxy(double *pxy, unsigned long maxl, double *d1, double *d2, unsigned long len, int circular);
int hgt_pop_calc_pxy_fft(double *pxy, unsigned long maxl, double *d1, double *d2, unsigned long len, int circular);
#endif
