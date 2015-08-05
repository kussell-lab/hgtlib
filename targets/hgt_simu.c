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
#include <gsl/gsl_randist.h>
#include <time.h>
#include "hgt_params.h"
#include "hgt_pop.h"
#include "hgt_file.h"
#include "hgt_stat.h"
#include "hgt_cov.h"
// global variables.
hgt_pop_sample_func SAMPLE_FUNC;
hgt_pop_frag_func FRAG_FUNC;

void specify_evolution_processes(hgt_params *params);
void evolve(hgt_params *params, hgt_pop *p, unsigned generation, const gsl_rng *r);
void sample(hgt_params *parmas, hgt_pop *p, unsigned sample_size, const gsl_rng *r, hgt_file_container *fc);
double *sample_pxy(hgt_pop *p, hgt_cov_sample_func sample_func, const gsl_rng *r);
void update_pxy(hgt_stat_meanvar_list *list, double *pxy);
double *calc_pxy(double *distance1, double *distance2, int length);
void write_pxy(FILE* f, hgt_stat_meanvar_list *list, int maxl, unsigned long gen);
double *calc_distance(hgt_genome *g1, hgt_genome *g2);
double *compare_genomes(hgt_genome *g1, hgt_genome *g2);
double *sample_t2(hgt_linkage_find_time_func linkage_find_func, hgt_pop *p, int linkage_size, const gsl_rng *r);
double sample_t2_one(hgt_linkage_find_time_func linkage_find_func, hgt_linkage **linkage, int size, int linkage_size, double current_time, const gsl_rng *r);
void write_t2(FILE *f, double *t2, int size, unsigned long gen);
void loadBar(int x, int n, int r, int w);
hgt_cov_result * calc_cov(hgt_pop *p, unsigned maxl, unsigned sample_size, const gsl_rng *r);
void update_cov(hgt_stat_meanvar_list *list, hgt_cov_result *result);
void write_cov(FILE *f, hgt_stat_meanvar_list *list, int maxl, unsigned long gen);
void update_ks(hgt_stat_meanvar_list *list, hgt_cov_result *result);
void write_ks(FILE *f, hgt_stat_meanvar_list *list, unsigned long gen);
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
	gsl_rng *r;
	long int seed;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	seed = (long int)time(NULL);
	gsl_rng_set(r, seed);

	// allocate population.
	hgt_pop *p;
	p = hgt_pop_alloc(params, r);

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
        sample(params, p, params->sample_size, r, fc);
        loadBar(i, params->sample_time, 100, 10);
	}
    
    // close files.
    hgt_file_container_close(fc);
    hgt_file_container_destroy(fc);
    
    hgt_pop_free(p);
    hgt_params_free(params);
    gsl_rng_free(r);
    
    return EXIT_SUCCESS;
}

// Process has done i out of n rounds,
// and we want a bar of width w and resolution r.
static void loadBar(int x, int n, int r, int w)
{
    // Only update r times.
    if ( x % (n/r +1) != 0 ) return;
    
    // Calculuate the ratio of complete-to-incomplete.
    float ratio = x/(float)n;
    int c= (int) ratio * w;
    
    // Show the percentage complete.
    printf("%3d%% [", (int)(ratio*100) );
    
    // Show the load bar.
    for (int x=0; x<c; x++)
        printf("=");
    
    for (int x=c; x<w; x++)
        printf(" ");
    
    // ANSI Control codes to go back to the
    // previous line and clear it.
    printf("]\n\033[F\033[J");
}

void specify_evolution_processes(hgt_params *params) {
	hgt_pop_sample_func sample_f;
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
    case 3:
        sample_f = hgt_pop_sample_bsc;
        break;
	default:
		sample_f = hgt_pop_sample_moran;
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
    unsigned int gen = p->generation;
    int count = params->seq_len * 4;
    hgt_stat_meanvar_list *p2_list = hgt_stat_meanvar_list_new(count);
    hgt_stat_meanvar_list *p3_list = hgt_stat_meanvar_list_new(count);
    hgt_stat_meanvar_list *p4_list = hgt_stat_meanvar_list_new(count);
	hgt_stat_meanvar_list *cov_list = hgt_stat_meanvar_list_new(params->maxl * 4);
	hgt_stat_meanvar_list *ks_list = hgt_stat_meanvar_list_new(2);
    unsigned i;
	for (i = 0; i < sample_size; i++)
	{
        // sample p2, p3, p4;
        double *pxy;
        pxy = sample_pxy(p, hgt_cov_sample_p2, r);
        update_pxy(p2_list, pxy);
        free(pxy);
        pxy = sample_pxy(p, hgt_cov_sample_p3, r);
        update_pxy(p3_list, pxy);
        free(pxy);
        pxy = sample_pxy(p, hgt_cov_sample_p4, r);
        update_pxy(p4_list, pxy);
        free(pxy);
        
        // sample t2;
        double *t2;
        int t2_size;
        t2 = sample_t2(hgt_linkage_find_most_rescent_ancestor_time, p, 2, r);
        t2_size = p->linkage_size + 1;
        write_t2(fc->t2, t2, t2_size, gen);
        free(t2);
        
        // sample t3;
        t2 = sample_t2(hgt_linkage_find_most_rescent_ancestor_time, p, 3, r);
        t2_size = p->linkage_size + 1;
        write_t2(fc->t3, t2, t2_size, gen);
        free(t2);
        
        // sample t4;
        t2 = sample_t2(hgt_linkage_find_most_rescent_ancestor_time, p, 4, r);
        t2_size = p->linkage_size + 1;
        write_t2(fc->t4, t2, t2_size, gen);
        free(t2);
        
        // sample q2;
        t2 = sample_t2(hgt_linkage_find_most_rescent_coalescence_time, p, 2, r);
        t2_size = p->linkage_size + 1;
        write_t2(fc->q2, t2, t2_size, gen);
        free(t2);
        
        // sample q3;
        t2 = sample_t2(hgt_linkage_find_most_rescent_coalescence_time, p, 3, r);
        t2_size = p->linkage_size + 1;
        write_t2(fc->q3, t2, t2_size, gen);
        free(t2);
        
        // sample q4;
        t2 = sample_t2(hgt_linkage_find_most_rescent_coalescence_time, p, 4, r);
        t2_size = p->linkage_size + 1;
        write_t2(fc->q4, t2, t2_size, gen);
        free(t2);
	}

	// sample and calculate covs.
	hgt_cov_result *cov_result;
	cov_result = calc_cov(p, params->maxl, sample_size, r);
	update_cov(cov_list, cov_result);
	update_ks(ks_list, cov_result);
    
    // write results to files.
    write_pxy(fc->p2, p2_list, params->maxl, gen);
    write_pxy(fc->p3, p3_list, params->maxl, gen);
    write_pxy(fc->p4, p4_list, params->maxl, gen);
	write_cov(fc->cov, cov_list, params->maxl, gen);
	write_ks(fc->ks, ks_list, gen);
    
    // destroy meanvar lists.
    hgt_stat_meanvar_list_destroy(p2_list);
    hgt_stat_meanvar_list_destroy(p3_list);
    hgt_stat_meanvar_list_destroy(p4_list);
	hgt_stat_meanvar_list_destroy(cov_list);
	hgt_stat_meanvar_list_destroy(ks_list);
    
    hgt_file_container_flush(fc);
}

double *sample_pxy(hgt_pop *p, hgt_cov_sample_func sample_func, const gsl_rng *r) {
    // sample two pairs of genomes.
    unsigned long a, b, c, d;
    sample_func(&a, &b, &c, &d, p->size, r);
    // calculate pairwise distance;
    double *distance1, *distance2;
    distance1 = calc_distance(p->genomes[a], p->genomes[b]);
    distance2 = calc_distance(p->genomes[c], p->genomes[d]);
    // compute pxy;
    double *pxy;
    int length;
    length = hgt_genome_get_seq_size(p->genomes[a]);
    pxy = calc_pxy(distance1, distance2, length);
    free(distance1);
    free(distance2);
    return pxy;
}

void update_pxy(hgt_stat_meanvar_list *list, double *pxy) {
    int i;
    for (i = 0 ; i < list->n; i++) {
        hgt_stat_meanvar_list_increment(list, i, pxy[i]);
    }
}

double *calc_distance(hgt_genome *g1, hgt_genome *g2) {
    double *distance;
    distance = compare_genomes(g1, g2);
    return distance;
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

void write_pxy(FILE* f, hgt_stat_meanvar_list *list, int maxl, unsigned long gen)
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
}

double *compare_genomes(hgt_genome *g1, hgt_genome *g2)
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

double *sample_t2(hgt_linkage_find_time_func linkage_find_func, hgt_pop *p, int linkage_size, const gsl_rng *r) {
    double current_time;
    current_time = hgt_pop_get_time(p);
    double *t2;
    t2 = (double *) malloc((p->linkage_size + 1) * sizeof(double));
    unsigned l;
    for (l = 0; l < p->linkage_size; l++) {
        hgt_linkage **linkages;
        linkages = p->locus_linkages[l];
        t2[l] = sample_t2_one(linkage_find_func, linkages, p->size, linkage_size, current_time, r);
    }
    t2[p->linkage_size] = sample_t2_one(linkage_find_func, p->linkages, p->size, linkage_size, current_time, r);
    return t2;
}

double sample_t2_one(hgt_linkage_find_time_func linkage_find_func, hgt_linkage **linkages, int size, int linkage_size, double current_time, const gsl_rng *r) {
    hgt_linkage ** chose = (hgt_linkage **) malloc(linkage_size * sizeof(hgt_linkage*));
    gsl_ran_choose(r, chose, linkage_size, linkages, size, sizeof(hgt_linkage*));
    double time, coal_time;
    time = linkage_find_func(chose, linkage_size);
    coal_time = current_time - time;
    free(chose);
    return coal_time;
}

void write_t2(FILE *f, double *t2, int size, unsigned long gen) {
    int i;
    for (i = 0; i < size; i++) {
        fprintf(f, "%d\t%g\t%lu\n", i, t2[i], gen);
    }
}

hgt_cov_result * calc_cov(hgt_pop *p, unsigned maxl, unsigned sample_size, const gsl_rng *r) {
	hgt_cov_result *result;
	result = hgt_cov_result_alloc(maxl);
	hgt_pop_calc_cov(result, p, sample_size, r); 
	return result;
}

void update_cov(hgt_stat_meanvar_list *list, hgt_cov_result *result) {
	int i;
	for (i = 0; i < result->maxl; i++)
	{
		hgt_stat_meanvar_list_increment(list, 4*i, result->scov[i]);
		hgt_stat_meanvar_list_increment(list, 4 * i + 1, result->rcov[i]);
		hgt_stat_meanvar_list_increment(list, 4 * i + 2, result->pxpy[i]);
		hgt_stat_meanvar_list_increment(list, 4 * i + 3, result->tcov[i]);
	}
}

void write_cov(FILE *f, hgt_stat_meanvar_list *list, int maxl, unsigned long gen) {
	int i, j, k, n;
	for (i = 0; i < maxl; i++)
	{
		fprintf(f, "%d\t", i);
		for (j = 0; j < 4; j++)
		{
			k = 4 * i + j;
			hgt_stat_meanvar *mv;
			mv = hgt_stat_meanvar_list_get(list, k);
			double m;
			m = hgt_stat_mean_get(mv->mean);
			fprintf(f, "%g\t", m);
		}

		for ( j = 0; j < 4; j++)
		{
			k = 4 * i + j;
			hgt_stat_meanvar *mv;
			mv = hgt_stat_meanvar_list_get(list, k);
			double v;
			v = hgt_stat_variance_get(mv->var);
			fprintf(f, "%g\t", v);
			n = hgt_stat_variance_get_n(mv->var);
		}
		fprintf(f, "%d\t%lu\n", n, gen);
	}
}

void update_ks(hgt_stat_meanvar_list *list, hgt_cov_result *result) {
	hgt_stat_meanvar_list_increment(list, 0, result->ks);
	hgt_stat_meanvar_list_increment(list, 1, result->vd);
}

void write_ks(FILE *f, hgt_stat_meanvar_list *list, unsigned long gen) 
{
	double ks, vd, ksvar, vdvar;
	int n;
	hgt_stat_meanvar *ksmv, *vdmv;
	ksmv = hgt_stat_meanvar_list_get(list, 0);
	vdmv = hgt_stat_meanvar_list_get(list, 1);
	ks = hgt_stat_mean_get(ksmv->mean);
	vd = hgt_stat_mean_get(vdmv->mean);
	ksvar = hgt_stat_variance_get(ksmv->var);
	vdvar = hgt_stat_variance_get(vdmv->var);
	n = hgt_stat_mean_get_n(ksmv->mean);
	fprintf(f, "%g\t%g\t%g\t%g\t%d\t%lu\n", ks, vd, ksvar, vdvar, n, gen);
}