//
//  hgt_utils.c
//  hgt
//
//  Created by Mingzhi Lin on 1/24/14.
//
//

#include <stdio.h>
#include <stdlib.h>
#include "hgt_utils.h"

// allocate HGT populations
// num: number of populations to be created
// size: population size
// seq_len: genome sequence length
// rank: identity for MPI process
// rs: gsl_rng random generators
hgt_pop ** hgt_utils_alloc_populations(hgt_pop_params *params, int rank, gsl_rng * rng) {
    int i;
    unsigned long num = params->replicates;
    hgt_pop ** pops;
    pops = malloc(num * sizeof(hgt_pop*));
    for (i = 0; i < num; i++) {
        pops[i] = hgt_pop_alloc(params, rng);
    }
    return pops;
}

void hgt_utils_free_populations(hgt_pop ** pops, int num) {
    int i;
    for (i = 0; i < num; i++) {
        hgt_pop_free(pops[i]);
    }
    free(pops);
}

// allocate gsl_rngs
// num: number of gsl_rng to be created
// rank: identity for MPI process
gsl_rng ** hgt_utils_alloc_gsl_rngs(int num, int rank) {
    gsl_rng ** rs;
    rs = malloc(num * sizeof(gsl_rng*));
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    int i, seed;
    for (i = 0; i < num; i++) {
        seed = rank * num + i;
        rs[i] = gsl_rng_alloc(T);
        gsl_rng_set(rs[i], seed);
    }
    return rs;
}

void hgt_utils_free_gsl_rngs(gsl_rng ** rs, int num) {
    int i;
    for (i = 0; i < num; i++) {
        gsl_rng_free(rs[i]);
    }
    free(rs);
}

hgt_stat_mean *** hgt_utils_alloc_stat_means(int num, int dimension) {
    hgt_stat_mean *** means;
    means = malloc(num * sizeof(hgt_stat_mean**));
    int i, j;
    for (i = 0; i < num; i++) {
        means[i] = malloc(dimension * sizeof(hgt_stat_mean*));
        for (j = 0; j < dimension; j++) {
            means[i][j] = hgt_stat_mean_alloc();
        }
    }
    return means;
}

void hgt_utils_free_stat_means(hgt_stat_mean *** means, int num, int dimension) {
    int i, j;
    for (i = 0; i < num; i++) {
        for (j = 0; j < dimension; j++) {
            hgt_stat_mean_free(means[i][j]);
        }
        free(means[i]);
    }
    free(means);
}

void hgt_utils_clean_stat_means(hgt_stat_mean *** means, int num, int dimension) {
    int i, j;
    for (i = 0; i < num; i++) {
        for (j = 0; j < dimension; j++) {
            hgt_stat_mean_clean(means[i][j]);
        }
    }
}

void hgt_utils_increment_stat_means(hgt_stat_mean *** means, double ** vals, int num, int dimension) {
    int i, j;
    for (i = 0; i < num; i++) {
        for (j = 0; j < dimension; j++) {
            hgt_stat_mean_increment(means[i][j], vals[i][j]);
        }
    }
}

hgt_stat_variance *** hgt_utils_alloc_stat_variance(int num, int dimension) {
    hgt_stat_variance *** vars;
    vars = malloc(num * sizeof(hgt_stat_variance**));
    int i, j;
    for (i = 0; i < num; i++) {
        vars[i] = malloc(dimension * sizeof(hgt_stat_variance*));
        for (j = 0; j < dimension; j++) {
            vars[i][j] = hgt_stat_variance_alloc();
            vars[i][j]->correct = 1;
        }
    }
    return vars;
}

void hgt_utils_clean_stat_variances(hgt_stat_variance *** vars, int num, int dimension) {
    int i, j;
    for (i = 0; i < num; i++) {
        for (j = 0; j < dimension; j++) {
            hgt_stat_variance_clean(vars[i][j]);
        }
    }
}

void hgt_utils_free_stat_variances(hgt_stat_variance *** vars, int num, int dimension) {
    int i, j;
    for (i = 0; i < num; i++) {
        for (j = 0; j < dimension; j++) {
            hgt_stat_variance_free(vars[i][j]);
        }
        free(vars[i]);
    }
    free(vars);
}

void hgt_utils_increment_stat_variances(hgt_stat_variance *** vars, double ** vals, int num, int dimension) {
    int i, j;
    for (i = 0; i < num; i++) {
        for (j = 0; j < dimension; j++) {
            hgt_stat_variance_increment(vars[i][j], vals[i][j]);
        }
    }
}

// do batch evoluation
int hgt_utils_batch_evolve_moran(hgt_pop ** pops, int num, hgt_pop_params * params, gsl_rng *rng) {
    int i, j;
    for (i = 0; i < num; i++) {
        for (j = 0; j < params->generations; j++) {
            hgt_pop_evolve(pops[i], params, hgt_pop_sample_moran, hgt_pop_coal_time_moran, hgt_pop_frag_constant, rng);
        }
    }
    return 0;
}

int hgt_utils_batch_evolve_moran_expon_frag(hgt_pop ** pops, int num, hgt_pop_params * params, gsl_rng *rng) {
    int i, j;
    for (i = 0; i < num; i++) {
        for (j = 0; j < params->generations; j++) {
            hgt_pop_evolve(pops[i], params, hgt_pop_sample_moran, hgt_pop_coal_time_moran, hgt_pop_frag_exp, rng);
        }
    }
    return EXIT_SUCCESS;
}

int hgt_utils_batch_evolve(hgt_pop **ps, int num, hgt_pop_params *params, hgt_pop_sample_func sample_func, hgt_pop_coal_time_func coal_time_func, gsl_rng *r) {
    int i, j;
    for (i = 0; i < num; i++) {
        for (j = 0; j < params->generations; j++) {
            hgt_pop_evolve(ps[i], params, sample_func, coal_time_func, hgt_pop_frag_constant, r);
        }
    }
    return EXIT_SUCCESS;
}

int hgt_utils_Roulette_Wheel_select(double * weights, int size, const gsl_rng *r) {
    int i;
    double total, v, accu;
    total = 0;
    for (i = 0; i < size; i++) {
        total += weights[i];
    }

    v = gsl_rng_uniform_pos(r);

    accu = 0;
    for (i = 0; i < size; i++) {
        accu += weights[i];
        if (accu/total >= v) {
           return i;
        }
    }

    return -1;
}
