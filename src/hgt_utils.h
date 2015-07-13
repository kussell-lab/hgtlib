//
//  hgt_utils.h
//  hgt
//
//  Created by Mingzhi Lin on 1/24/14.
//
//

#include <gsl/gsl_rng.h>
#include "hgt_stat.h"
#include "hgt_pop.h"
#include "hgt_params.h"

#ifndef hgt_hgt_utils_h
#define hgt_hgt_utils_h

hgt_pop ** hgt_utils_alloc_populations(hgt_params *params, int rank, gsl_rng *rng);
void hgt_utils_free_populations(hgt_pop ** pops, int num);

gsl_rng ** hgt_utils_alloc_gsl_rngs(int num, int rank);
void hgt_utils_free_gsl_rngs(gsl_rng ** rngs, int num);

hgt_stat_mean *** hgt_utils_alloc_stat_means(int num, int dimension);
void hgt_utils_free_stat_means(hgt_stat_mean *** means, int num, int dimension);
void hgt_utils_increment_stat_means(hgt_stat_mean *** means, double ** vals, int num, int dimension);
void hgt_utils_clean_stat_means(hgt_stat_mean *** means, int num, int dimension);

hgt_stat_variance *** hgt_utils_alloc_stat_variance(int num, int dimension);
void hgt_utils_free_stat_variances(hgt_stat_variance *** vars, int num, int dimension);
void hgt_utils_increment_stat_variances(hgt_stat_variance *** vars, double ** vals, int num, int dimension);
void hgt_utils_clean_stat_variances(hgt_stat_variance *** vars, int num, int dimension);

int hgt_utils_batch_evolve_moran(hgt_pop ** pops, int num, hgt_params * params, gsl_rng *rng);
int hgt_utils_batch_evolve_moran_expon_frag(hgt_pop ** pops, int num, hgt_params * params, gsl_rng *rng);
int hgt_utils_batch_evolve(hgt_pop **ps, int num, hgt_params *params, hgt_pop_sample_func sample_func, hgt_pop_coal_time_func coal_time_func, hgt_pop_frag_func frag_f, gsl_rng *r);

int hgt_utils_Roulette_Wheel_select(double * weights, int size, const gsl_rng *r);

#endif
