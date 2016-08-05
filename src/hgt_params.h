#ifndef HGT_HGT_PARAMS_H
#define HGT_HGT_PARAMS_H
#include <stdio.h>
typedef struct _hgt_params hgt_params;
struct _hgt_params {
    // population parameters
    unsigned int size;     // population size
    unsigned int seq_len;  // genome length

	// genome sequence parameters
	unsigned alphabet_size;
    
    // transfer parameters
    double tr_rate;                             // transfer rate
    unsigned int frag_len;						// fragment length
    unsigned int tr_hotspot_num;                // number of hotspot in the genome
    unsigned int tr_hotspot_length;				// length of hotspot
    double tr_hotspot_ratio;                    // ratio of transfer rate
    unsigned int **tr_hotspots;					// locations of transfer hotspots
	double tr_eff;								// transfer efficiency
    int tr_eff_len;                             // length of sequence for calculating distance
    
    // mutation parameters
    double mu_rate;                     // mutation rate
    unsigned int mu_hotspot_num;        // number of mutation hotspots
    unsigned int mu_hotspot_length;    // average hotspot length
    double mu_hotspot_ratio;            // ratio of mutation rate
    unsigned int **mu_hotspots;        // locations of mutation hotspots
    
    // sample parameters
    unsigned int generations;          // generations
    unsigned int sample_size;  // sample size
    unsigned int sample_time;  // sample time
    unsigned int replicates;   // replicates
    unsigned int maxl;         // maxl
	unsigned int sample_generations; // sample generations
    
    // output parameters
    char * prefix;
	int save_pop;
    int save_pxy;
    
    // fitting paramters
    int fit_range; // fitting range
    int fit_flat;  // fitting flat
    
    // fitness
    unsigned int fitness_size;
    double fitness_scale;
    double fitness_shape;
    double fitness_mutation_rate;   // fitness mutation rate
    int fitness_coupled;         // coupled with neutral genome.
    
    // linkage tracking size.
    int linkage_size;
    
    // reproduction model
    unsigned int reprodution;
    // fragment type
    unsigned int frag_type;
};

int hgt_params_parse(hgt_params *params, int argc, char **argv, char * progname);
int hgt_params_free(hgt_params *params);
int hgt_params_printf(hgt_params *params, FILE *stream);
hgt_params *hgt_params_alloc();

#endif
