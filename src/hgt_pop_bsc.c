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
    
    // copy genomes.
    for (i = 0; i < num_offsprings; i++) {
        int d;
        d = deaths[i];
        hgt_genome_copy(p->genomes[d], p->genomes[birth]);
    }
    
    // update linkages.
    linkage_update_bsc(p->linkages, birth, deaths, num_offsprings, current_time);
    unsigned l;
    for (l = 0; l < p->linkage_size; l++) {
        linkage_update_bsc(p->locus_linkages[l], birth, deaths, num_offsprings, current_time);
    }
    
    
    free(deaths);
    free(indices);
    
    double time;
    time = hgt_pop_coal_time_bsc(p->size, r);
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

void linkage_update_bsc(hgt_linkage **linkages, int birth, int *death, int num, double birth_time) {
    hgt_linkage *parent;
    parent = linkages[birth];
    linkages[birth] = hgt_linkage_new(parent, birth_time);
    int i;
    for (i = 0; i < num; i++) {
        int d;
        d = death[i];
        hgt_linkage_free(linkages[d]);
        linkages[d] = hgt_linkage_new(parent, birth_time);
    }
}

double hgt_pop_coal_time_bsc(unsigned long size, const gsl_rng *r) {
    double time;
    time = gsl_ran_exponential(r, 2.0);
    return time;
}