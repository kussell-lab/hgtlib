//
// hgt_forward_coal_back.c
// hgt targets
//
// Created by Mingzhi Lin on 7/23/15.
//
//

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include "hgt_params.h"
#include "hgt_pop.h"
#include "hgt_utils.h"
#include "hgt_genome.h"
#include "hgt_linkage.h"
#include "hgt_forward_coal_back.h"


int main(int argc, char *argv[]) {
    
    // read hgt params.
    hgt_params *params = hgt_params_alloc();
    int exit_code;
    exit_code = hgt_params_parse(params, argc, argv, "hgt_forward_coal_back");
    if (exit_code == EXIT_FAILURE) {
        goto exit;
    }
    
    hgt_params_printf(params, stdout);
    
    // setting up random generator.
    const gsl_rng_type *T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    long int seed = (long int) time(NULL);
    gsl_rng *r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);
    
    // allocate populations.
    hgt_pop **ps;
    ps = hgt_utils_alloc_populations(params, 0, r);
    
    hgt_pop_sample_func sample_f;
    hgt_pop_coal_time_func coal_time_f;
    hgt_pop_frag_func frag_f;
    
    sample_f = hgt_pop_sample_moran;
    coal_time_f = hgt_pop_coal_time_moran;
    frag_f = hgt_pop_frag_constant;
    
    file_container *fc = create_file_container(params->prefix);
    
    int i;
    for (i = 0; i < params->sample_time; i++) {
        unsigned long current_generation;
        current_generation = (i + 1) * params->generations;
        hgt_utils_batch_evolve(ps, params->replicates, params, sample_f, coal_time_f, frag_f, r);
        sample(ps, params, 2, current_generation, fc, r);
    }
    
    close_file_container(fc);
    destroy_file_container(fc);
exit:
    
    return EXIT_SUCCESS;
}

file_container * create_file_container(char * prefix) {
    file_container *fc = malloc(sizeof(file_container));
    char * fn;
    asprintf(&fn, "%s.p2.txt", prefix);
    fc->p2 = fopen(fn, "w");
    
    free(fn);
    return fc;
}

int destroy_file_container(file_container *fc) {
    free(fc);
    return EXIT_SUCCESS;
}

int close_file_container(file_container *fc) {
    fclose(fc->p2);
    return EXIT_SUCCESS;
}

// randomly sample two linkage,
// calculate their coalescent time,
// and randomly add mutations in two genomes,
// calculate p2.
int sample(hgt_pop **ps, hgt_params *params, 
           int linkage_size, unsigned long gen, file_container *files,
           const gsl_rng *r) {
    
    hgt_linkage ** linkages = malloc(linkage_size * sizeof(hgt_linkage*));
    // allocate
    linkages = malloc(linkage_size * sizeof(hgt_linkage*));

    int i;
    for (i = 0; i < params->replicates; i++) {
        hgt_pop *p = ps[i];
        int s;
        for (s = 0; s < params->sample_size; ++s)
        {
            unsigned long coal_time;
            gsl_ran_choose(r, linkages, linkage_size, p->linkages, p->size, sizeof(hgt_linkage*));
            coal_time = hgt_linkage_find_most_rescent_ancestor_time(linkages, linkage_size);
            coal_evolve(params, linkage_size, params->seq_len, coal_time, gen, files, r);
        }
    }
    
    return EXIT_SUCCESS;
}

int coal_evolve(hgt_params *params, int size, int length, unsigned long time, unsigned long gen, file_container *files, const gsl_rng *r) {
    double mutation_rate = params->mu_rate;
    hgt_genome **genomes = create_genomes(size, length, r);
    mutate_genomes(genomes, size, mutation_rate, time, r);
    calc_pxy(genomes, size, params->maxl, files, gen);
    destroy_genomes(genomes, size);
    return EXIT_SUCCESS;
}

int calc_pxy(hgt_genome **genomes, int size, int maxl, file_container *files, unsigned long gen) {
    int calc_p2(hgt_genome *g1, hgt_genome *g2, int maxl, FILE *f, unsigned long gen);
    int i, j;
    for (i = 0; i < size; i++) {
        for (j = i+1; j < size; j++) {
            hgt_genome *g1, *g2;
            g1 = genomes[i];
            g2 = genomes[j];
            calc_p2(g1, g2, maxl, files->p2, gen);
        }
    }
    
    return EXIT_SUCCESS;
}

int calc_p2(hgt_genome *g1, hgt_genome *g2, int maxl, FILE *f, unsigned long gen) {
    int write_pxy(FILE *f, double *pxy, int maxl, unsigned long gen);
    int length = hgt_genome_get_seq_size(g1);
    int count = maxl * 4;
    double *pxy = malloc(count * sizeof(double));
    double *distance = malloc(length * sizeof(double));
    compare_genomes(g1, g2, distance, length);
    hgt_pop_calc_pxy_fft(pxy, maxl, distance, distance, length, 1);
    write_pxy(f, pxy, maxl, gen);
    free(distance);
    free(pxy);
    return EXIT_SUCCESS;
}

int write_pxy(FILE *f, double *pxy, int maxl, unsigned long gen) {
    int i, j;
    for (j = 0; j < maxl; j++) {
        fprintf(f, "%d\t", j);
        for (i = 0; i < 4; i++) {
            fprintf(f, "%g\t", pxy[4*j+i]);
        }
        fprintf(f, "%lu\n", gen);
    }
    fflush(f);
    return EXIT_SUCCESS;
}


hgt_genome ** create_genomes(int size, int length, const gsl_rng * r) {
    int random_seq(char *seq, int length, const gsl_rng *r);
    char *ancestor = malloc(length * sizeof(char));
    random_seq(ancestor, length, r);
    
    hgt_genome ** genomes = malloc(size * sizeof(hgt_genome*));
    int i;
    for (i = 0; i < size; i++) {
        genomes[i] = hgt_genome_new(ancestor, length, 1);
    }
    
    free(ancestor);
    
    return genomes;
}

int destroy_genomes(hgt_genome **genomes, int size) {
    int i;
    for (i = 0; i < size; i++) {
        hgt_genome_free(genomes[i]);
    }
    free(genomes);
    return EXIT_SUCCESS;
}

int mutate_genomes(hgt_genome **genomes, int size, double mutation_rate, unsigned long time, const gsl_rng *r) {
    int i;
    for (i = 0; i < size; i++) {
        hgt_genome *g = genomes[i];
        int length = hgt_genome_get_seq_size(g);
        double mu = mutation_rate * (double) length;
        int count = gsl_ran_poisson(r, mu);
        int c;
	for (c = 0; c < count; c++) {
            int pos = gsl_rng_uniform_int(r, length);
            hgt_genome_mutate(g, pos, r);
        }
    }
    return EXIT_SUCCESS;
}

int compare_genomes(hgt_genome *g1, hgt_genome *g2, double *distance, int length) {
    int i;
    char *s1 = hgt_genome_get_seq(g1);
    char *s2 = hgt_genome_get_seq(g2);
    for (i = 0; i < length; i++) {
        if (s1[i] != s2[i]) {
            distance[i] = 1.0;
        } else {
            distance[i] = 0.0;
        }
    }
    return EXIT_SUCCESS;
}
