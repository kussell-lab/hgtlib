//
// hgt_pop_bsc.c
// hgt
//
// Created by Mingzhi Lin on 8/3/15
//
//
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "hgt_pop.h"
#include "hgt_genome.h"
#include "hgt_pop_internal.h"

// Based on the population model proposed by
// Jason Schweinsberg in "Dynamics of the evolving Bolthausen-Sznitman coalescent" paper.
double hgt_pop_sample_bsc(hgt_pop *p, const gsl_rng *r) {
    double current_time;
    current_time = hgt_pop_get_time(p);
    
    // reproduction.
    // random choose an individual to reproduce.
    unsigned long birth;
    int num_offsprings;
    birth = gsl_rng_uniform_int(r, p->size);
    num_offsprings = random_num_offsprings_bsc(p->size, r);
    // random choose num_offsprings to death.
    int *indices;
    indices = (int *) malloc((p->size - 1) * sizeof(int));
    int i, k;
    k = 0;
    for (i = 0; i < p->size; i++) {
        if (i != birth) {
            indices[k] = i;
            k++;
        }
    }
    int *deaths;
    deaths = (int *) malloc(num_offsprings * sizeof(int));
    gsl_ran_choose(r, deaths, num_offsprings, indices, p->size - 1, sizeof(int));
    for (i = 0; i < num_offsprings; i++) {
        // copy the genome.
        int d = deaths[i];
        hgt_genome_copy(p->genomes[d], p->genomes[birth]);
        // update linkage tracking.
        linkage_birth_dead(p->linkages, birth, d, current_time);
        unsigned l;
        for (l = 0; l < p->linkage_size; l++) {
            linkage_birth_dead(p->locus_linkages[l], birth, d, current_time);
        }
    }
    
    free(deaths);
    free(indices);
    
    double time;
    time = hgt_pop_coal_time_moran(p->size, r);
    increase_population_time(p, time);
    
    return time;
}

int random_num_offsprings_bsc(unsigned int n, const gsl_rng *r) {
    double u = gsl_rng_uniform_pos(r);
    double total = 0;
    
    int k;
    for (k = 1; k < n-1; k++) {
        double prop = (double) n / (double) (n - 1);
        prop /= ((double) (k * (k+1)));
        total += prop;
        if (total > u) {
            break;
        }
    }
    return k;
}