//
//  hgt_pop.c
//  hgt
//
//  Created by Mingzhi Lin on 1/18/14.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include "hgt_genome.h"
#include "hgt_pop.h"
#include "hgt_corr.h"
#include "hgt_utils.h"
#include "hgt_pop_internal.h"

hgt_pop * hgt_pop_alloc(hgt_params *params, const gsl_rng * r) {
    unsigned i, j;

    hgt_pop * p = (hgt_pop *) malloc (sizeof(hgt_pop));
    
	// initialize population properities.
    p->size = params->size;
    p->target_size = params->size;
    p->seq_len = params->seq_len;
    p->generation = 0;
	p->total_time = 0;
    p->linkage_size = params->linkage_size;
	p->total_tr_count = 0;
	p->succ_tr_count = 0;

	if (params->fitness_mutation_rate > 0) {
		p->track_fitness = 1;
	} else {
		p->track_fitness = 0;
	}
    
	// create linkages for tracking.
    p->linkages = (hgt_linkage **) malloc(p->size * sizeof(hgt_linkage*));
    for (j = 0; j < p->size; ++j) {
        p->linkages[j] = hgt_linkage_new(NULL, p->total_time);
    }
    p->locus_linkages = (hgt_linkage ***) malloc(p->linkage_size * sizeof(hgt_linkage**));
    for (i = 0; i < p->linkage_size; i++) {
        p->locus_linkages[i] = (hgt_linkage**) malloc(p->size * sizeof(hgt_linkage*));
        for (j = 0; j < p->size; j++) {
            p->locus_linkages[i][j] = hgt_linkage_new(NULL, p->total_time);
        }
    }
    
	// create genomes.
    p->genomes = (hgt_genome **) malloc(p->size * sizeof(hgt_genome *));
	// set genome alphabet size beforing create genome sequences.
	hgt_genome_set_alphabet_size(params->alphabet_size);
    // random initilize genomes
	char * ancestor = hgt_genome_random_sequence(p->seq_len, r);
    for (i = 0; i < p->size; i ++) {
        p->genomes[i] = hgt_genome_new(ancestor, p->seq_len, params->fitness_size);
    }
    
    // define transfer hotspots
    if (params->tr_hotspot_num > 0) {
        params->tr_hotspots = malloc(params->tr_hotspot_num * sizeof(unsigned long*));
        unsigned int init_start = (unsigned int) gsl_rng_uniform_int(r, params->seq_len);
        unsigned int expected_space = (params->seq_len - params->tr_hotspot_num*params->tr_hotspot_length) / (params->tr_hotspot_num);
        for (i = 0; i < params->tr_hotspot_num; i++) {
            params->tr_hotspots[i] = malloc(2*sizeof(unsigned long));
            unsigned int start = init_start;
            unsigned int end = start + params->tr_hotspot_length;
            unsigned int space = expected_space;
            init_start = end + space;
            params->tr_hotspots[i][0] = start;
            params->tr_hotspots[i][1] = end;
        }
    }
    
    // define mutation hotspots
    if (params->mu_hotspot_num > 0) {
        params->mu_hotspots = malloc(params->mu_hotspot_num * sizeof(unsigned long*));
        unsigned int init_start = (unsigned int) gsl_rng_uniform_int(r, params->seq_len);
        unsigned int expected_space = (params->seq_len - params->mu_hotspot_num*params->mu_hotspot_length) / (params->mu_hotspot_num);
        for (i = 0; i < params->mu_hotspot_num; i++) {
            params->mu_hotspots[i] = malloc(2*sizeof(unsigned long));
            unsigned int start = init_start;
            unsigned int end = start + params->mu_hotspot_length;
            unsigned int space = expected_space;
            init_start = end + space;
            params->mu_hotspots[i][0] = start;
            params->mu_hotspots[i][1] = end;
        }
    }
    
    p->survived = NULL;
    p->new_born = NULL;
    p->cache_allocated = 0;
    
    free(ancestor);
    return p;
}

hgt_pop * hgt_pop_copy(hgt_pop * p) {
    hgt_pop * new_p = (hgt_pop *) calloc(1, sizeof(hgt_pop));
    new_p->size = p->size;
    new_p->seq_len = p->seq_len;
    new_p->genomes = (hgt_genome **) malloc(p->size * sizeof(hgt_genome *));
	unsigned i;
    for (i = 0; i < p->size; i++) {
        unsigned int seq_len = p->genomes[i]->seq_len;
        unsigned int fitness_size = p->genomes[i]->fitness_size;
        new_p->genomes[i] = hgt_genome_alloc(seq_len, fitness_size);
    }
    
    // initilize caches
    new_p->survived = (unsigned int *) calloc(p->size, sizeof(unsigned int));
    new_p->new_born = (unsigned int *) calloc(p->size, sizeof(unsigned int));
    
    return new_p;
}

int hgt_pop_free(hgt_pop * p) {
    unsigned i, j;
    for (i = 0; i < p->size; i ++) {
        hgt_genome_free(p->genomes[i]);
        hgt_linkage_free(p->linkages[i]);
    }

    for (i = 0; i < p->linkage_size; i++) {
        for (j = 0; j < p->size; j++) {
            hgt_linkage_free(p->locus_linkages[i][j]);
        }
        free(p->locus_linkages[i]);
    }

    free(p->genomes);
    free(p->linkages);
    free(p->locus_linkages);
    if (p->cache_allocated != 0) {
        free(p->survived);
        free(p->new_born);
    }
    free(p);
    
    return EXIT_SUCCESS;
}

char *hgt_pop_to_json(hgt_pop *p, hgt_params *params){
    char *c;
    bstring b = bfromcstr("");
    bformata(b, "{\n");
    bformata(b, "\"Size\": %ld,\n", params->size);
    bformata(b, "\"Length\": %ld,\n", params->seq_len);
    bformata(b, "\"MutationRate\": %g,\n", params->mu_rate);
    bformata(b, "\"TransferRate\": %g,\n", params->tr_rate);
    bformata(b, "\"FragLen\": %ld,\n", params->frag_len);
    bformata(b, "\"Generation\": %ld,\n", params->generations);
    bformata(b, "\"Genomes\": [\n");
    unsigned j;
    for (j = 0; j < p->size; j ++) {
        if (j < p->size - 1) {
            bformata(b, "\"%s\",\n", p->genomes[j]->seq);
        } else {
            bformata(b, "\"%s\"\n", p->genomes[j]->seq);
        }
    }
    bformata(b, "]\n}");
    c = bstr2cstr(b, '\n');
    bdestroy(b);
    
    return c;
}

int hgt_pop_fitness_mutate_step(hgt_pop *p, hgt_params *params, const gsl_rng *r) {
    unsigned int g, s, fitness_size;
    double fitness_scale;
    // randomly choose a cell.
    g = (unsigned int)gsl_rng_uniform_int(r, p->size);
    fitness_size = hgt_genome_get_fitness_size(p->genomes[g]);
    fitness_scale = params->fitness_scale;
    s = (unsigned int)gsl_rng_uniform_int(r, fitness_size);
    hgt_genome_fitness_mutate(p->genomes[g], s, fitness_scale);
    return EXIT_SUCCESS;
}

int hgt_pop_mutate(hgt_pop *p, hgt_params* params, const gsl_rng* r) {

    unsigned int g, s, hs_len, pos;
    double ratio;
    // randomly choose a genome for mutation
    g = (unsigned int)gsl_rng_uniform_int(r, p->size);
    
    if (params -> mu_hotspot_num > 0) { // we need to deal hotspots
        // first calculate total hotspot length, and its ratio to sequence length
        hs_len = params->mu_hotspot_num * params->mu_hotspot_length;
        ratio = (double) hs_len / (double) p->seq_len;
        
        // we need to decise a mutation is happened inside hotspot or outside
        // we flip a coin
        double v = gsl_rng_uniform(r);
        if (v < (1-ratio)/(1+(params->mu_hotspot_ratio - 1)*ratio)) { // outside hotspot
            // we randomly choose a position from sites outside hotspots
            pos = (unsigned int) gsl_rng_uniform_int(r, p->seq_len - hs_len);
            // and then search the absolute position on the genome
            s = search_region(pos, params->mu_hotspots, params->mu_hotspot_num, 0);
        } else {
            // similarly, we randomly choose a position from hotspots
            pos = (unsigned int) gsl_rng_uniform_int(r, hs_len);
            // and then search the absolute position on the genome
            s = search_region(pos, params->mu_hotspots, params->mu_hotspot_num, 1) % p->seq_len;
        }
    } else { // otherwise, we just randomly choose a site from the genome.
        s = (unsigned int) gsl_rng_uniform_int(r, p->seq_len);
    }
    
    // do mutation given a genome and a site
    hgt_genome_mutate(p->genomes[g], s, r);
    
    return EXIT_SUCCESS;
}

int hgt_pop_transfer_linkages(hgt_pop *p,
                        unsigned int donor,
                        unsigned int receiver,
                        unsigned int frag_len,
                        unsigned int start)
{

    unsigned int i, track_linkage;
    unsigned int s = p->seq_len/2;
    for (i = 0; i < p->linkage_size; i++) {
        track_linkage = 0;
        if (start <= s && s+i < start + frag_len) {
            track_linkage = 1;
        }

        if (track_linkage == 1) {
            linkage_transfer(p->locus_linkages[i], donor, receiver, p->total_time);
        }
    }

    return EXIT_SUCCESS;
}

int linkage_transfer(hgt_linkage ** linkages, int donor, int receiver, double birth_time) {
    hgt_linkage * parent;
    if (donor != receiver) {
        parent = linkages[donor]->parent;
		birth_time = linkages[donor]->birthTime;
        hgt_linkage_free(linkages[receiver]);
        linkages[receiver] = hgt_linkage_new(parent, birth_time);
    }
    return EXIT_SUCCESS;
}

void increase_population_time(hgt_pop *p, double time) {
	p->generation++;
	p->total_time += time;
}

double hgt_pop_get_time(hgt_pop *p) {
	return p->total_time;
}

double hgt_pop_sample_moran(hgt_pop *p, const gsl_rng *r) {
	// obtain current generation.
	double current_time;
	current_time = hgt_pop_get_time(p);
	
	// reproduction.
    unsigned int b, d;
    double * fitness;
    // randomly choose a going-death cell.
    d = (unsigned int) gsl_rng_uniform_int(r, p->size);
	if (p->track_fitness > 0) {
    	// randomly choose a going-birth one according to the fitness.
    	// first calculate the fitness for each cell.
    	fitness = (double *) malloc(p->size * sizeof(double));
    	hgt_pop_calc_fitness(p, fitness);
    	// randomly choose one proportional to the fitness.
    	b = (unsigned int) hgt_utils_Roulette_Wheel_select(fitness, p->size, r);
	} else {
		b = (unsigned int) gsl_rng_uniform_int(r, p->size);
	}
    // copy the genome and its fitness from b to d.
    if (b != d) {
        hgt_genome_copy(p->genomes[d], p->genomes[b]);
    }

	if (p->linkage_size > 0)
	{
		// update linkages.
		linkage_birth_dead(p->linkages, b, d, current_time);
		unsigned i;
		for (i = 0; i < p->linkage_size; i++) {
			hgt_linkage ** linkages;
			linkages = p->locus_linkages[i];
			linkage_birth_dead(linkages, b, d, current_time);
		}
	}
	
	if (p->track_fitness > 0) {
    	free(fitness);
	}

	double time;
	time = hgt_pop_coal_time_moran(p->size, r);
	increase_population_time(p, time);
    return time;
}

double hgt_pop_sample_linear_selection(hgt_pop *p, const gsl_rng *r) {
    hgt_genome ** current_genomes, ** new_genomes;
    double * normalized_fitness;
    unsigned int current_size;

    current_genomes = p->genomes;
    current_size = p->size;
    
    normalized_fitness = (double *) malloc(p->size * sizeof(double));
    hgt_pop_calc_fitness(p, normalized_fitness);
    
    unsigned int i, num_offspring, new_size;
    unsigned int * offsprings;
    offsprings = (unsigned int *) malloc(p->size * sizeof(unsigned int));
    new_size = 0;
    
    for (i = 0; i < current_size; i++) {
        num_offspring = gsl_ran_poisson(r, normalized_fitness[i]);
        offsprings[i] = num_offspring;
        new_size += num_offspring;
    }
    
    new_genomes = (hgt_genome **) malloc(new_size * sizeof(hgt_genome*));
    unsigned j, k;
    hgt_genome * g;
    k = 0;
    for (i = 0; i < current_size; i++) {
        num_offspring = offsprings[i];
        if (num_offspring > 0) {
            for (j = 0; j < num_offspring; j++) {
                if (j == 0) {
                    g = current_genomes[i];
                } else {
                    g = hgt_genome_alloc(current_genomes[i]->seq_len, current_genomes[i]->fitness_size);
                    hgt_genome_copy(g, current_genomes[i]);
                }
                new_genomes[k] = g;
                k++;
            }
        } else {
            hgt_genome_free(current_genomes[i]);
        }
    }
    
    update_linkages_wf(p, offsprings);
    p->size = new_size;
    p->genomes = new_genomes;
    
    free(current_genomes);
    free(normalized_fitness);
    free(offsprings);

	// increase population generation and total time;
	double time; 
	time = hgt_pop_coal_time_linear_selection(p->size, r);
	increase_population_time(p, time);

    return time;
}

int linkage_birth_dead(hgt_linkage **linkages, int birth, int dead, double birth_time) {
    hgt_linkage * parent;
    parent = linkages[birth];
    linkages[birth] = hgt_linkage_new(parent, birth_time);
    if (birth != dead) {
        hgt_linkage_free(linkages[dead]);
        linkages[dead] = hgt_linkage_new(parent, birth_time);
    }
    return EXIT_SUCCESS;
}

int update_linkages_wf(hgt_pop *p, unsigned int *offsprings) {
    unsigned int current_size = 0, new_size = 0;
    current_size = p->size;
    unsigned i = 0;
    for (i = 0; i < current_size; i++) {
        new_size+=offsprings[i];
    }
    
    hgt_linkage ** current_linkages = p->linkages;
    hgt_linkage ** new_linkages = (hgt_linkage **) malloc(new_size * sizeof(hgt_linkage*));
    linkage_linear_selection(current_linkages, offsprings, current_size, new_linkages, p->total_time);
    hgt_linkage_free_more(current_linkages,current_size);
    p->linkages = new_linkages;
    free(current_linkages);

    for (i = 0; i < p->linkage_size; i++) {
        current_linkages = p->locus_linkages[i];
        new_linkages = (hgt_linkage**) malloc(new_size * sizeof(hgt_linkage*));
        linkage_linear_selection(current_linkages, offsprings, current_size, new_linkages, p->total_time);
        p->locus_linkages[i] = new_linkages;
        hgt_linkage_free_more(current_linkages, current_size);
        free(current_linkages);
    }

    return EXIT_SUCCESS;
}

int linkage_linear_selection(hgt_linkage ** current_linkages, unsigned int * offsprings, unsigned int current_size, hgt_linkage ** new_linkages, double birth_time){
    unsigned i, num_offspring;
    
    hgt_linkage *l, *parent;
    unsigned j, k;
    k = 0;
    for (i = 0; i < current_size; i++) {
        parent = current_linkages[i];
        num_offspring = offsprings[i];
        if (num_offspring > 0) {
            for (j = 0; j < num_offspring; j++) {
                l = hgt_linkage_new(parent, birth_time);
                new_linkages[k] = l;
                k++;
            }
            if (parent->numChildren != num_offspring) {
                printf("not equal offsprings: %d vs %d\n!", parent->numChildren, num_offspring);
            }
        } else {
            parent->numChildren = 0;
        }
    }
    
    return EXIT_SUCCESS;
}

double hgt_pop_sample_wf(hgt_pop *p, const gsl_rng *r) {
    unsigned i, j, k;
    
    unsigned int *offsprings = (unsigned int*) calloc(p->size, sizeof(unsigned int));
    for (i = 0; i < p->size; i++) {
        offsprings[i] = 0;
    }
    
    for (i = 0; i < p->size; i ++) {
        j = (unsigned int) gsl_rng_uniform_int(r, p->size);
        offsprings[j]++;
    }
    
    hgt_genome **new_genomes = (hgt_genome **) malloc(p->size * sizeof(char*));
    hgt_genome **current_genomes = p->genomes;
    hgt_genome *g;
    k = 0;
    unsigned num_offspring;
    for (i = 0; i < p->size; i++) {
        num_offspring = offsprings[i];
        if (num_offspring > 0) {
            for (j = 0; j < num_offspring; j++) {
                if (j == 0) {
                    g = current_genomes[i];
                } else {
                    g = hgt_genome_alloc(current_genomes[i]->seq_len, current_genomes[i]->fitness_size);
                    hgt_genome_copy(g, current_genomes[i]);
                }
                new_genomes[k] = g;
                k++;
            }
        } else {
            hgt_genome_free(current_genomes[i]);
        }
    }
    
    update_linkages_wf(p, offsprings);
    p->genomes = new_genomes;
    free(offsprings);
    free(current_genomes);

	// increase population generation and total time.
	double time = hgt_pop_coal_time_wf(p->size, r);
	increase_population_time(p, time);

    return time;

}

unsigned hgt_pop_frag_constant(hgt_params *params, const gsl_rng *r) {
    return params->frag_len;
}

unsigned hgt_pop_frag_exp(hgt_params *params, const gsl_rng *r) {
    return (unsigned) gsl_ran_exponential(r, (double) params->frag_len);
}

unsigned hgt_pop_frag_geom(hgt_params *params, const gsl_rng *r) {
    double p = 1.0 / (double) params->frag_len;
    return (unsigned) gsl_ran_geometric(r, p);
}

int hgt_pop_evolve(hgt_pop *p, hgt_params *params, hgt_pop_sample_func sample_f, hgt_pop_frag_func frag_f, const gsl_rng *r) {
    double time = sample_f(p, r);
    // mutate fitness.
    double fitness_mutation_rate = params->fitness_mutation_rate * (double) p->size;
    if (fitness_mutation_rate > 0) {
		double mu;
		int count;
        mu = time * fitness_mutation_rate;
        count = gsl_ran_poisson(r, mu);
		int k;
        for (k = 0; k < count; k++) {
            hgt_pop_fitness_mutate_step(p, params, r);
        }
    }
    
    int weight_size = 2;
    double weights[2];
    weights[0] = params->mu_rate * (double) (p->seq_len * p->size);
    weights[1] = params->tr_rate * (double) (p->seq_len * p->size);
    double total = weights[0] + weights[1];
    double mu = total * time;
    int count = gsl_ran_poisson(r, mu);
    int k;
    for (k = 0; k < count; k++) {
        int g = gsl_rng_uniform_int(r, p->size);
        int choice = hgt_utils_Roulette_Wheel_select(weights, weight_size, r);
        int pos = gsl_rng_uniform_int(r, p->seq_len);
        if (choice == 0) {
            hgt_genome_mutate(p->genomes[g], pos, r);
        }
        else {
            int d = gsl_rng_uniform_int(r, p->size);
            if (d != g) {
                hgt_genome *receiver = p->genomes[g];
                hgt_genome *donor = p->genomes[d];
                int frag_len = frag_f(params, r);

				// calculate distance in the transfer region,
				// and determine whether they are homologous enough for transfer.
				int to_transfer = 1;
				if (params->tr_eff > 0)
				{
					to_transfer = 0;
					double distance, r1, v;
					distance = hgt_genome_distance(receiver, donor, pos, params->tr_eff_len);
					r1 = gsl_rng_uniform(r);
					v = exp(-distance / params->tr_eff);
					if (r1 < v)
					{
						to_transfer = 1;
					}
				}
				p->total_tr_count++;
				if (to_transfer) {
					hgt_genome_transfer(receiver, donor, pos, frag_len);
                    if (params->fitness_coupled > 0) {
                        hgt_genome_fitness_transfer(receiver, donor, pos, frag_len);
                    }

					if (params->linkage_size > 0) {
						hgt_pop_transfer_linkages(p, d, g, frag_len, pos);
					}
					p->succ_tr_count++;
				}
            }
        }
    }

//	unsigned g;
//	for (g = 0; g < p->size; g++) {
//		int weight_size = 2;
//		double weights[2];
//		weights[0] = params->mu_rate * (double)p->seq_len;
//		weights[1] = params->tr_rate * (double)p->seq_len;
//		double total = weights[0] + weights[1];
//		double mu = total * time;
//		int count = gsl_ran_poisson(r, mu);
//		int k;
//		for (k = 0; k < count; k++) {
//			int choice = hgt_utils_Roulette_Wheel_select(weights, weight_size, r);
//			int pos = gsl_rng_uniform_int(r, p->seq_len);
//			if (choice == 0) {
//				hgt_genome_mutate(p->genomes[g], pos, r);
//			}
//			else {
//				int d = gsl_rng_uniform_int(r, p->size);
//                if (d != g) {
//                    hgt_genome *receiver = p->genomes[g];
//                    hgt_genome *donor = p->genomes[d];
//                    int frag_len = frag_f(params, r);
//                    hgt_genome_transfer(receiver, donor, pos, frag_len);
//                    hgt_pop_transfer_linkages(p, d, g, frag_len, pos);
//                }
//			}
//		}
//
//	}

	if (p->generation % 10 == 0 && params->linkage_size > 0)
	{
		hgt_pop_prune_linkages(p);
	}
    
    return EXIT_SUCCESS;
}


double hgt_pop_coal_time_moran(unsigned long p_size, const gsl_rng *r) {
    double mu = 1.0 / (double) p_size;
    double time = gsl_ran_exponential(r, mu);
    return time;
}

double hgt_pop_coal_time_wf(unsigned long p_size, const gsl_rng *r) {
    return 1.0;
}

double hgt_pop_coal_time_linear_selection(unsigned long p_size, const gsl_rng *r) {
    return hgt_pop_coal_time_wf(p_size, r);
}

double hgt_pop_calc_ks(hgt_pop *p) {
    double ks;
    unsigned i, j, k;
    ks = 0;
    
    for (i = 0; i < p->size; i++) {
        for (j = i + 1; j < p->size; j++) {
            for (k = 0; k < p->seq_len; k++) {
                if (p->genomes[i]->seq[k] != p->genomes[j]->seq[k]) {
                    ks++;
                }
            }
        }
    }
    
    ks /= (double) (p->size * (p->size - 1)/2 * p->seq_len);
    
    return ks;
}

int hgt_pop_calc_dist(hgt_pop *p, double *ds1, double *ds2, unsigned long sample_size, hgt_cov_sample_func sample_func, const gsl_rng *r) {
    unsigned i, j;
    unsigned long a, b, c, d;
    for (i = 0; i < p->seq_len; i++) {
        ds1[i] = 0;
        ds2[i] = 0;
    }
    
    for (i = 0; i < sample_size; i++) {
        sample_func(&a, &b, &c, &d, p->size, r);
        for (j = 0; j < p->seq_len; j++) {
            if (p->genomes[a]->seq[j] != p->genomes[b]->seq[j]) {
                ds1[j]++;
            }
            if (p->genomes[c]->seq[j] != p->genomes[d]->seq[j]) {
                ds2[j]++;
            }
        }
    }
    
    for (i = 0; i < p->seq_len; i++) {
        ds1[i] /= (double) sample_size;
        ds2[i] /= (double) sample_size;
    }
    
    return EXIT_SUCCESS;
}

int hgt_pop_calc_pxy(double *pxy, unsigned int maxl, double *ds1, double *ds2, unsigned int len, int circular) {
    unsigned i, l;
    for (l = 0; l < maxl; l++) {
        for (i = 0; i < 4; i++) {
            pxy[l*4+i] = 0;
        }
    }
    
    int a, b;
    for (l = 0; l < maxl; l++) {
        unsigned long num = 0;
        for (i = 0; i < len; i++) {
            if (circular != 0) {
                a = i;
                b = (i-l+len) % len;
                num++;
            } else {
                if (i < l) {
                    continue;
                }
                a = i;
                b = i - l;
                num++;
            }
            pxy[l*4+3] += ds1[a] * ds2[b];
            pxy[l*4+3] += ds2[a] * ds1[b];
            
            pxy[l*4+2] += ds1[a] * (1-ds2[b]);
            pxy[l*4+2] += (1-ds2[a]) * ds1[b];
            
            pxy[l*4+1] += (1-ds1[a]) * ds2[b];
            pxy[l*4+1] += ds2[a]* (1-ds1[b]);
            
            pxy[l*4+0] += (1-ds1[a]) * (1-ds2[b]);
            pxy[l*4+0] += (1-ds2[a]) * (1-ds1[b]);
        }
        pxy[l*4+3] /= (double) 2*num;
        pxy[l*4+2] /= (double) 2*num;
        pxy[l*4+1] /= (double) 2*num;
        pxy[l*4+0] /= (double) 2*num;
    }
    
    return EXIT_SUCCESS;
}

int hgt_pop_calc_pxy_fft(double *pxy, unsigned int maxl, double *d1, double *d2, unsigned int len, int circular) {

    unsigned long fft_len;
    unsigned i, j, l;
    fft_len = next_power2(len);
	double **ds1, **ds2;
	ds1 = (double **)malloc(2 *sizeof(double*));
	ds2 = (double **)malloc(2 *sizeof(double*));
	for (i = 0; i < 2; i++) {
		ds1[i] = (double *)malloc(2 * fft_len*sizeof(double));
		ds2[i] = (double *)malloc(2 * fft_len*sizeof(double));
	}
	double *mask;
	mask = (double *)malloc(2 * fft_len * sizeof(double));
    for (i = 0; i < 2*fft_len; i++) {
        if (i < len) {
            ds1[1][i] = d1[i%len];
            ds2[1][i] = d2[i%len];
            ds1[0][i] = 1.0 - ds1[1][i];
            ds2[0][i] = 1.0 - ds2[1][i];
            mask[i] = 1;
        } else if (i < 2*len) {
            if (circular != 0) {
                ds1[1][i] = d1[i%len];
                ds2[1][i] = d2[i%len];
                ds1[0][i] = 1.0 - ds1[1][i];
                ds2[0][i] = 1.0 - ds2[1][i];
                mask[i] = 1;
            } else {
                ds1[1][i] = 0;
                ds2[1][i] = 0;
                ds1[0][i] = 0;
                ds2[0][i] = 0;
                mask[i] = 0;
            }
        } else {
            ds1[1][i] = 0;
            ds2[1][i] = 0;
            ds1[0][i] = 0;
            ds2[0][i] = 0;
            mask[i] = 0;
        }
    }
    
    hgt_corr_auto_fft(mask, fft_len);
    
	double *buf1, *buf2;
	buf1 = (double *)malloc(2 * fft_len*sizeof(double));
	buf2 = (double *)malloc(2 * fft_len*sizeof(double));
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 2; j++) {
            for (l = 0; l < 2*fft_len; l++) {
                buf1[l] = ds1[i][l];
                buf2[l] = ds2[j][l];
            }
            hgt_corr_fft(buf1, buf2, fft_len);
            for (l = 1; l < maxl; l++) {
                if (circular != 0) {
                    pxy[l*4+2*i+j] = (buf1[l] + buf1[2*fft_len-len+l]) / (mask[l]+mask[2*fft_len-len+l]);
                } else {
                    pxy[l*4+2*i+j] = buf1[l] / mask[l];
                }
                
            }
            pxy[2*i+j] = buf1[0] / mask[0];
            
            for (l = 0; l < 2*fft_len; l++) {
                buf2[l] = ds1[i][l];
                buf1[l] = ds2[j][l];
            }
            hgt_corr_fft(buf1, buf2, fft_len);
            for (l = 1; l < maxl; l++) {
                if (circular != 0) {
                    pxy[l*4+2*i+j] += (buf1[l] + buf1[2*fft_len-len+l]) / (mask[l]+mask[2*fft_len-len+l]);
                } else {
                    pxy[l*4+2*i+j] += buf1[l] / mask[l];
                }
                pxy[l*4+2*i+j] /= 2.0;
            }
            pxy[2*i+j] += buf1[0] / mask[0];
            pxy[2*i+j] /= 2.0;
        }
    }

	for (i = 0; i < 2; i++) {
		free(ds1[i]);
		free(ds2[i]);
	}
	free(ds1);
	free(ds2);
	free(mask);
	free(buf1);
	free(buf2);
    
    return EXIT_SUCCESS;
}

// Calcualte covariances using all individuals.
int hgt_pop_calc_cov_all(hgt_cov_result *result, hgt_pop *p) {
    // allocate a binary matrix
    unsigned long matrix_size = p->size * (p->size - 1) / 2;
    short  **matrix = malloc( matrix_size * sizeof(short*));
    
    unsigned i, j, k, h;
    
    k = 0;
    for (i = 0; i < p->size; i++) {
        for (j = i + 1; j < p->size; j++) {
            matrix[k] = malloc(p->seq_len*sizeof(short));
            for (h = 0; h < p->seq_len; h++) {
                if (p->genomes[i]->seq[h] != p->genomes[j]->seq[h]) {
                    matrix[k][h] = 1;
                } else {
                    matrix[k][h] = 0;
                }
            }
            k++;
        }
    }
    
    hgt_cov_result_calc_matrix(result, matrix, matrix_size, p->seq_len);
    
    for (i = 0; i < matrix_size; i++) {
        free(matrix[i]);
    }

    free(matrix);
    
    return EXIT_SUCCESS;
}

int hgt_pop_calc_cov(hgt_cov_result *result, hgt_pop *p, unsigned sample, const gsl_rng* rng) {
	// randomly choose a set of genome of sample size.
	if (sample > p->size) {
        sample = p->size;
    }
    hgt_genome **selected_genomes = (hgt_genome **)malloc(sample * sizeof(hgt_genome*));
	gsl_ran_choose(rng, selected_genomes, sample, p->genomes, p->size, sizeof(hgt_genome*));

	unsigned seq_len = hgt_genome_get_seq_size(selected_genomes[0]);

	// create binary matrix.
	unsigned sizeOfMatrix = sample * (sample - 1) / 2;
	short **matrix = malloc(sizeOfMatrix * sizeof(short*));

	unsigned i, j, k, h;
	h = 0;
	for (i = 0; i < sample; i++) {
		for (j = i + 1; j < sample; j++) {
			matrix[h] = (short *)malloc(seq_len * sizeof(short));
			char *seq_a, *seq_b;
			seq_a = hgt_genome_get_seq(selected_genomes[i]);
			seq_b = hgt_genome_get_seq(selected_genomes[j]);
			for (k = 0; k < seq_len; k++) {
				if (seq_a[k] != seq_b[k]) {
					matrix[h][k] = 1;
				}
				else {
					matrix[h][k] = 0;
				}
			}
			h++;
		}
	}

	free(selected_genomes);

	hgt_cov_result_calc_matrix(result, matrix, sample, p->seq_len);

	for (i = 0; i < sizeOfMatrix; i++) {
		free(matrix[i]);
	}
	free(matrix);

	return EXIT_SUCCESS;
}

int hgt_pop_calc_cov2(hgt_cov_result *result, hgt_pop *p, unsigned sample, const gsl_rng* rng) {
    // allocate a binary matrix
    short **matrix = malloc(sample*sizeof(short*));
    unsigned i, j;
    
    for (i = 0; i < sample; i++) {
        // randomly choose two distinct genomes for comparison.
		hgt_genome **selected_genomes = (hgt_genome **)malloc(2 * sizeof(hgt_genome *));
		gsl_ran_choose(rng, selected_genomes, 2, p->genomes, p->size, sizeof(hgt_genome*));
		unsigned seq_len = hgt_genome_get_seq_size(selected_genomes[0]);
        matrix[i] = (short *) malloc(seq_len * sizeof(short));
		for ( j = 0; j < seq_len; j++)
		{
			char *seq_a, *seq_b;
			seq_a = hgt_genome_get_seq(selected_genomes[0]);
			seq_b = hgt_genome_get_seq(selected_genomes[1]);
			if (seq_a[j] != seq_b[j])
			{
				matrix[i][j] = 1;
			}
			else {
				matrix[i][j] = 0;
			}
		}
		free(selected_genomes);
    }
    hgt_cov_result_calc_matrix(result, matrix, sample, p->seq_len);
    
    for (i = 0; i < sample; i++) {
        free(matrix[i]);
    }
    free(matrix);
    
    return EXIT_SUCCESS;
}

double hgt_pop_mean_fitness(hgt_pop *p) {
    double f;
    unsigned i;
    f = 0;
    for (i = 0; i < p->size; ++i){
        f += hgt_genome_get_fitness(p->genomes[i]);
    }

    f /= (double) p->size;

    return f;
}

int hgt_pop_prune_linkages(hgt_pop *p) {
    hgt_linkage_prune_more(p->linkages, p->size);
    unsigned i;
    for (i = 0; i < p->linkage_size; i++) {
        hgt_linkage_prune_more(p->locus_linkages[i], p->size);
    }
    return EXIT_SUCCESS;
}

int hgt_pop_calc_coal_time(hgt_linkage ** pop_linkages, unsigned size, 
    unsigned int sample_size,
    double *res, 
    unsigned linkage_size, 
    hgt_linkage_find_time_func find_func,
    const gsl_rng *r) {

    unsigned i, j;
	int *a, *b;
	a = (int *)malloc(linkage_size * sizeof(int));
	b = (int *)malloc(size * sizeof(int));
    for (i = 0; i < size; i++) {
        b[i] = i;
    }

    hgt_linkage ** linkages = (hgt_linkage**) malloc(linkage_size * sizeof(hgt_linkage*));

    for (i = 0; i < sample_size; i++) {
        gsl_ran_choose(r, a, linkage_size, b, size, sizeof(int));
        for (j = 0; j < linkage_size; j++) {
            linkages[j] = pop_linkages[a[j]];
        }
        res[i] = find_func(linkages, linkage_size);
    }

	free(a);
	free(b);
	free(linkages);
    return EXIT_SUCCESS;
}

double hgt_pop_calc_most_recent_coal_time(hgt_linkage ** pop_linkages, unsigned size, unsigned long sample_size, double * res, unsigned linkage_size, const gsl_rng *r) {
    return hgt_pop_calc_coal_time(pop_linkages, size, sample_size, res, linkage_size, hgt_linkage_find_most_rescent_coalescence_time, r);
}

double hgt_pop_calc_most_recent_ancestor_time(hgt_linkage ** pop_linkages, unsigned size, unsigned long sample_size, double * res, unsigned linkage_size, const gsl_rng *r) {
    return hgt_pop_calc_coal_time(pop_linkages, size, sample_size, res, linkage_size, hgt_linkage_find_most_rescent_ancestor_time, r);
}

int hgt_pop_calc_fitness(hgt_pop *p, double * fitness) {
    double f, mean_fitness, cpot, size_ratio;
    size_ratio = (double) p->size / (double) p->target_size;
    mean_fitness = hgt_pop_mean_fitness(p);
    cpot = mean_fitness - (1.0 - size_ratio);
    unsigned i;
    for (i = 0; i < p->size; i++) {
        f = hgt_genome_get_fitness(p->genomes[i]);
        fitness[i] = exp(f - cpot);
    }

    return EXIT_SUCCESS;
}

// Calculate coalescent time for each pair of genome.
int hgt_pop_coal_time_matrix(double **matrix, hgt_pop *p) {
    unsigned int i, j;
    hgt_linkage **linkages;
    linkages = (hgt_linkage **) malloc(2 * sizeof(hgt_linkage*));
    double current_time = hgt_pop_get_time(p);
    for (i = 0; i < p->size; i++) {
        for (j = 0; j < p->size; j++) {
            double time = 0;
            if (i != j) {
                linkages[0] = p->linkages[i];
                linkages[1] = p->linkages[j];
                time = hgt_linkage_find_most_rescent_coalescence_time(linkages, 2);
                time = current_time - time;
            }
            matrix[i][j] = time;
        }
    }
    
    return EXIT_SUCCESS;
}

/******** PRIVATE FUNCTIONS ***********/

unsigned long next_power2(unsigned int len) {
    unsigned int v;
    v = (unsigned int) pow(2, (int)(log(len)/log(2))+1);
    if (v % len == 0) {
        v = len;
    }
    return v;
}

unsigned int search_region(unsigned int pos, unsigned int** region, unsigned num_regions, int inside) {
    unsigned i;
    unsigned long c = pos;
    for (i = 0; i < num_regions; i++) {
        if (inside == 0) {
            if (c < region[i][0]) {
                return c;
            } else {
                c += (region[i][1] - region[i][0]);
            }
        } else {
            if (c + region[i][0] < region[i][1] ) {
                return c + region[i][0];
            } else {
                c -= (region[i][1] - region[i][0]);
            }
        }
    }
    return c;
}
