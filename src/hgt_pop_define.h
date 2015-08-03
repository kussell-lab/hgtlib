#ifndef hgt_hgt_pop_define_h
#define hgt_hgt_pop_define_h
typedef struct _hgt_pop hgt_pop;

struct _hgt_pop {
    unsigned int size;     // population size
    unsigned int seq_len;  // genome length
    unsigned int generation;
    double total_time;
    hgt_genome ** genomes;        // genome sequences
    hgt_linkage ** linkages; // linkages.
    hgt_linkage *** locus_linkages; // locus linkages.
    
    int ** transfer_hotspots; // transfer hotspots
    
    // cache
    unsigned int * survived;
    unsigned int * new_born;
    int cache_allocated;
    unsigned int linkage_size;
    unsigned int target_size; // target population size.
};
#endif