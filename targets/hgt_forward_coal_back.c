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
#include "hgt_stat.h"

typedef struct _meanvar meanvar;
struct _meanvar {
    hgt_stat_mean* mean;
    hgt_stat_variance* var;
};

int main(int argc, char* argv[])
{
    // read hgt params.
    hgt_params* params = hgt_params_alloc();
    int exit_code;
    exit_code = hgt_params_parse(params, argc, argv, "hgt_forward_coal_back");
    if (exit_code == EXIT_FAILURE) {
        goto exit;
    }

    hgt_params_printf(params, stdout);

    // setting up random generator.
    const gsl_rng_type* T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    long int seed = (long int)time(NULL);
    gsl_rng* r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);

    // allocate populations.
    hgt_pop** ps;
    ps = hgt_utils_alloc_populations(params, 0, r);

    // specify process functions.
    hgt_pop_sample_func sample_f = hgt_pop_sample_moran;
    hgt_pop_coal_time_func coal_time_f = hgt_pop_coal_time_moran;
    hgt_pop_frag_func frag_f = hgt_pop_frag_constant;
    
    // prepare output files.
    file_container* fc = create_file_container(params->prefix);

    int i;
    for (i = 0; i < params->sample_time; i++) {
        unsigned long current_generation;
        current_generation = (i + 1) * params->generations;
        hgt_utils_batch_evolve(ps, params->replicates, params, sample_f,
            coal_time_f, frag_f, r);
        sample(ps, params, 2, coal_time_f, current_generation, fc, r);
    }

    close_file_container(fc);
    destroy_file_container(fc);
exit:

    return EXIT_SUCCESS;
}

file_container* create_file_container(char* prefix)
{
    FILE* create_file(char* prefix, char* appdix);
    file_container* fc = malloc(sizeof(file_container));
    fc->p2 = create_file(prefix, "p2");
    fc->t2 = create_file(prefix, "t2");
    return fc;
}

FILE* create_file(char* prefix, char* appdix)
{
    char* fn;
    asprintf(&fn, "%s.%s.txt", prefix, appdix);
    FILE* f;
    f = fopen(fn, "w");
    free(fn);
    return f;
}

int destroy_file_container(file_container* fc)
{
    free(fc);
    return EXIT_SUCCESS;
}

int close_file_container(file_container* fc)
{
    fclose(fc->p2);
    fclose(fc->t2);
    return EXIT_SUCCESS;
}

int flush_file_container(file_container* fc)
{
    fflush(fc->t2);
    fflush(fc->p2);
    return EXIT_SUCCESS;
}

// randomly sample two linkage,
// calculate their coalescent time,
// and randomly add mutations in two genomes,
// calculate p2.
int sample(hgt_pop** ps, hgt_params* params, int linkage_size, hgt_pop_coal_time_func coal_time_f,
    unsigned long gen, file_container* files, const gsl_rng* r)
{
    hgt_linkage** linkages = malloc(linkage_size * sizeof(hgt_linkage*));
    // allocate
    linkages = malloc(linkage_size * sizeof(hgt_linkage*));

    int i;
    for (i = 0; i < params->replicates; i++) {
        hgt_pop* p = ps[i];
        
        // create mean and variance containers.
        int count = params->maxl * 4;
        hgt_stat_meanvar_list *list = hgt_stat_meanvar_list_new(count);
        hgt_stat_meanvar *t2mv = hgt_stat_meanvar_new();
        int s;
        for (s = 0; s < params->sample_size; ++s) {
            // randomly choose two linear and measure their coalescent time.
            gsl_ran_choose(r, linkages, linkage_size, p->linkages, p->size, sizeof(hgt_linkage*));
            unsigned long time = hgt_linkage_find_most_rescent_ancestor_time(linkages, linkage_size);
            // update sample results.
            if (time > 0) {
                double coal_time = (double)(gen - time + 1) * coal_time_f(p->size, r);
                coal_evolve(params, linkage_size, params->seq_len, coal_time, list, r);
                update_t2(t2mv, coal_time);
            }
        }
        // write to files.
        write_pxy(files->p2, list, params->maxl, gen);
        write_t2(files->t2, t2mv, gen);
        flush_file_container(files);
        // destroy mean and variance containers.
        hgt_stat_meanvar_list_destroy(list);
        hgt_stat_meanvar_destroy(t2mv);
    }

    return EXIT_SUCCESS;
}

int coal_evolve(hgt_params* params, int size, int length, double time,
    hgt_stat_meanvar_list *list, const gsl_rng* r)
{
    double mutation_rate = params->mu_rate;
    hgt_genome** genomes = create_genomes(size, length, r);
    mutate_genomes(genomes, size, mutation_rate, time, r);
    calc_pxy(genomes, size, params->maxl, list);
    destroy_genomes(genomes, size);
    return EXIT_SUCCESS;
}

int calc_pxy(hgt_genome** genomes, int size, int maxl, hgt_stat_meanvar_list* list)
{
    int i, j;
    for (i = 0; i < size; i++) {
        for (j = i + 1; j < size; j++) {
            hgt_genome *g1, *g2;
            g1 = genomes[i];
            g2 = genomes[j];
            double* pxy = calc_p2(g1, g2, maxl);
            update_pxy(list, pxy);
            free(pxy);
        }
    }

    return EXIT_SUCCESS;
}

double* calc_p2(hgt_genome* g1, hgt_genome* g2, int maxl)
{
    int length = hgt_genome_get_seq_size(g1);
    double* pxy = (double*)malloc(maxl * 4 * sizeof(double));
    double* distance = compare_genomes(g1, g2);
    hgt_pop_calc_pxy_fft(pxy, maxl, distance, distance, length, 1);
    free(distance);
    return pxy;
}

int update_pxy(hgt_stat_meanvar_list* list, double* pxy)
{
    int i;
    for (i = 0; i < list->n; i++) {
        hgt_stat_meanvar_list_increment(list, i, pxy[i]);
    }
    return EXIT_SUCCESS;
}

int write_pxy(FILE* f, hgt_stat_meanvar_list *list, int maxl, unsigned long gen)
{
    int i, j;
    for (i = 0; i < maxl; i++) {
        fprintf(f, "%d\t", i);
        unsigned long n;
        for (j = 0; j < 4; j++) {
            int index = 4 * i + j;
            hgt_stat_meanvar *mv = hgt_stat_meanvar_list_get(list, index);
            double v = hgt_stat_mean_get(mv->mean);
            fprintf(f, "%g\t", v);
            n = hgt_stat_mean_get_n(mv->mean);
        }
        for (j = 0; j < 4; j++) {
            int index = 4 * i + j;
            hgt_stat_meanvar *mv = hgt_stat_meanvar_list_get(list, index);
            double v = hgt_stat_variance_get(mv->var);
            fprintf(f, "%g\t", v);
        }
        fprintf(f, "%lu\t%lu\n", n, gen);
    }
    return EXIT_SUCCESS;
}

int update_t2(hgt_stat_meanvar *mv, double t2) {
    hgt_stat_meanvar_increment(mv, t2);
    return EXIT_SUCCESS;
}

int write_t2(FILE* f, hgt_stat_meanvar *mv, unsigned long gen)
{
    double m = hgt_stat_mean_get(mv->mean);
    double v = hgt_stat_variance_get(mv->var);
    unsigned long n = hgt_stat_mean_get_n(mv->mean);
    fprintf(f, "%g\t%g\t%lu\t%lu\n", m, v, n, gen);
    return EXIT_SUCCESS;
}

hgt_genome** create_genomes(int size, int length, const gsl_rng* r)
{
    int random_seq(char* seq, int length, const gsl_rng* r);
    char* ancestor = hgt_genome_random_sequence(length, r);

    hgt_genome** genomes = malloc(size * sizeof(hgt_genome*));
    int i;
    for (i = 0; i < size; i++) {
        genomes[i] = hgt_genome_new(ancestor, length, 1);
    }

    free(ancestor);

    return genomes;
}

int destroy_genomes(hgt_genome** genomes, int size)
{
    int i;
    for (i = 0; i < size; i++) {
        hgt_genome_free(genomes[i]);
    }
    free(genomes);
    return EXIT_SUCCESS;
}

int mutate_genomes(hgt_genome** genomes, int size, double mutation_rate,
    double time, const gsl_rng* r)
{
    int i;
    for (i = 0; i < size; i++) {
        hgt_genome* g = genomes[i];
        int length = hgt_genome_get_seq_size(g);
        double mu = mutation_rate * (double)length * time;
        int count = gsl_ran_poisson(r, mu);
        int c;
        for (c = 0; c < count; c++) {
            int pos = gsl_rng_uniform_int(r, length);
            hgt_genome_mutate(g, pos, r);
        }
    }
    return EXIT_SUCCESS;
}

double* compare_genomes(hgt_genome* g1, hgt_genome* g2)
{
    int length = hgt_genome_get_seq_size(g1);
    double* distance = (double*)malloc(length * sizeof(double));
    char* s1 = hgt_genome_get_seq(g1);
    char* s2 = hgt_genome_get_seq(g2);
    int i;
    for (i = 0; i < length; i++) {
        if (s1[i] == s2[i]) {
            distance[i] = 0;
        }
        else {
            distance[i] = 1;
        }
    }
    return distance;
}
