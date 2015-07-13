#ifndef HGT_HGT_PARAMS_H
#define HGT_HGT_PARAMS_H
#include <stdio.h>
typedef struct _hgt_params hgt_params;
struct _hgt_params {
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
    unsigned int fitness_size;
    double fitness_scale;
    double fitness_shape;
    double b_mu_rate;
    
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