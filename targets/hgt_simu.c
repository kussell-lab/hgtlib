//
// hgt_simu.c
// hgt targets
// 
// Created by Mingzhi Lin on 7/31/15.
//
//

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "hgt_params.h"
#include "hgt_pop.h"
#include "hgt_file.h"

// global variables.
hgt_pop_sample_func SAMPLE_FUNC;
hgt_pop_frag_func FRAG_FUNC;

void specify_evolution_processes(hgt_params *params);
void evolve(hgt_params *params, hgt_pop *p, unsigned generation, const gsl_rng *r);
void sample(hgt_params *parmas, hgt_pop *p, unsigned sample_size, const gsl_rng *r, hgt_file_container *fc);
double *calc_p2(hgt_genome *g1, hgt_genome *g2);
double *calc_p3(hgt_genome *g1, hgt_genome *g2, hgt_genome *g3);
double *calc_px(double *distance1, double *distance2, int length); 
double *compare_genomes(hgt_genome *g1, hgt_genome *g2);
int main(int argc, char **argv) {
	// parse hgt parameters.
	hgt_params *params;
	params = hgt_params_alloc();
	int exit_code;
	exit_code = hgt_params_parse(params, argc, argv, "hgt_simu");
	if (exit_code == EXIT_FAILURE)
	{
		printf("ERROR: could not parse params.\n");
		exit(EXIT_FAILURE);
	}
	else {
		hgt_params_printf(params, stdout);
	}

	// set up random number generator.
	const gsl_rng_type *T;
	const gsl_rng *r;
	long int seed;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	seed = (long int)time(NULL);
	gsl_rng_set(r, seed);

	// allocate population.
	hgt_pop *p;
	p = hgt_pop_alloc(params, r);
	return EXIT_SUCCESS;

	// specify evolution processes.
	specify_evolution_processes(params);

	// prepare output files.
	hgt_file_container *fc;
	fc = hgt_file_container_create(params->prefix);

	// evolve to equilibrium.
	evolve(params, p, params->generations, r);

	// sample
	unsigned i;
	for ( i = 0; i < params->sample_time; i++)
	{
		evolve(params, p, params->sample_generations, r);

	}

}

void specify_evolution_processes(hgt_params *params) {
	hgt_pop_sample_func sample_f;
	hgt_pop_coal_time_func coal_time_f;
	hgt_pop_frag_func frag_f;
	switch (params->frag_type) {
	case 1:
		frag_f = hgt_pop_frag_exp;
		break;
	default:
		frag_f = hgt_pop_frag_constant;
		break;
	}
	switch (params->reprodution) {
	case 1:
		sample_f = hgt_pop_sample_wf;
		break;
	case 2:
		sample_f = hgt_pop_sample_linear_selection;
		break;
	default:
		sample_f = hgt_pop_sample_moran;
		coal_time_f = hgt_pop_coal_time_moran;
		break;
	}

	SAMPLE_FUNC = sample_f;
	FRAG_FUNC = frag_f;
}

void evolve(hgt_params *params, hgt_pop *p, unsigned generation, const gsl_rng *r) {
	unsigned g;
	for ( g = 0; g < generation; g++)
	{
		hgt_pop_evolve(p, params, SAMPLE_FUNC, FRAG_FUNC, r);
	}
}

void sample(hgt_params *params, hgt_pop *p, unsigned sample_size, const gsl_rng *r, hgt_file_container *fc) {
	unsigned i;
	for (i = 0; i < sample_size; i++)
	{

	}
}

double *calc_p2(hgt_genome *g1, hgt_genome *g2)
{
	double *distance;
	distance = compare_genome(g1, g2);
	int length;
	double *pxy;
	length = hgt_genome_get_seq_size(g1);
	pxy = calc_pxy(distance, distance, length);
	return pxy;
}

double *calc_pxy(double *distance1, double *distance2, int length) {
	int circular = 1;
	int dim = 4;
	// allocate results.
	double *pxy;
	pxy = (double *)malloc(length * dim * sizeof(double));
	hgt_pop_calc_pxy_fft(pxy, length, distance1, distance2, length, circular);
	return pxy;
}

double *compare_genome(hgt_genome *g1, hgt_genome *g2)
{
	// allocate an array of results.
	double *results;
	int length;
	length = hgt_genome_get_seq_size(g1);
	results = (double *)malloc(length * sizeof(double));

	// compare two sequences.
	char *s1, *s2;
	s1 = hgt_genome_get_seq(g1);
	s2 = hgt_genome_get_seq(g2);
	int i;
	for (i = 0; i < length; i++)
	{
		if (s1[i] == s2[i])
		{
			results[i] = 0;
		}
		else
		{
			results[i] = 1;
		}
	}

	return results;
}