//
//  hgt_pop.h
//  hgt
//
//  Created by Mingzhi Lin on 1/18/14.
//
//

#ifndef hgt_hgt_pop_h
#define hgt_hgt_pop_h

typedef struct {
    unsigned long size;     // population size
    unsigned long seq_len;  // genome length
    char ** genomes;        // genome sequences
    
    // cache
    unsigned long * survived;
    unsigned long * new_born;
} hgt_pop;

hgt_pop * hgt_pop_alloc(unsigned long size, unsigned long seq_len, const gsl_rng * r);
hgt_pop * hgt_pop_copy(hgt_pop * p);
int hgt_pop_free(hgt_pop * r);

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

int hgt_pop_params_parse(hgt_pop_params *params, int argc, char **argv, char * progname);
int hgt_pop_params_free(hgt_pop_params *params);
#endif
