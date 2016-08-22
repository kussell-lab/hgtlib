//
//  hgt_mpi_simu2.c
//  hgt
//
//  Created by Mingzhi Lin on 8/21/16.
//
//

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
#include <time.h>
#include <math.h>
#include "hgt_pop.h"
#include "hgt_utils.h"
#include "bstrlib.h"
#include "hgt_params.h"
#include "hgt_mpi_simu.h"
#ifndef MPI_SUCCESS
    #define MPI_SUCCESS 0
#endif

char *pop_to_json(hgt_pop *p, hgt_params *params);
int check_mpi_error_code(int error_code, char *ops_type, char *parent_func);
char *pop_to_json(hgt_pop *p, hgt_params *params);
short int ** get_coal_rank_matrix(hgt_pop *p);
short int * get_coal_ranks(hgt_pop *p, int index);
double ** get_coal_time_matrix(hgt_pop *p);
double * get_coal_times(hgt_pop *p, int index);
int write_pops(hgt_pop **ps, hgt_params *params, int rank, int numprocs);
bstring to_json(hgt_pop *pop, hgt_params *params);

int main(int argc, char *argv[]) {
    int numprocs, rank, exit_code;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    hgt_params *params = hgt_params_alloc();
    exit_code = hgt_params_parse(params, argc, argv, "hgt_mpi_moran_const");
    if (exit_code == EXIT_FAILURE) {
        goto exit;
    }
    
    if (rank == 0) {
        hgt_params_printf(params, stdout);
    }
    
    const gsl_rng_type *T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    long int seed = (long int) time(NULL) + rank;
    gsl_rng *rng = gsl_rng_alloc(T);
    gsl_rng_set(rng, seed);
    hgt_pop **ps = hgt_utils_alloc_populations(params, rank, rng);
    
    unsigned i;
    time_t start, end;
    hgt_pop_sample_func sample_f;
    hgt_pop_frag_func frag_f;
    switch (params->frag_type) {
        case 1:
            frag_f = hgt_pop_frag_exp;
            break;
        case 2:
            frag_f = hgt_pop_frag_geom;
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
            if (rank == 0) {
                printf("using BSC model\n");
            }
            break;
        default:
            sample_f = hgt_pop_sample_moran;
            break;
    }
    
	start = clock();

	hgt_utils_batch_evolve(ps, params->replicates, params, sample_f, frag_f, rng);

	end = clock();

	if (rank == 0) {
		printf("%d, evolve using time = %ld sec\n", i, (end - start) / CLOCKS_PER_SEC);
	}
    
	start = clock();
	write_pops(ps, params, rank, numprocs);
	end = clock();
	if (rank == 0) {
		printf("%d, writting using time = %ld sec\n", i, (end - start) / CLOCKS_PER_SEC);
	}
    
    // hgt_utils_free_populations(ps, params->replicates);
    gsl_rng_free(rng);
    hgt_params_free(params);
    
exit:
    MPI_Finalize();
    
    return EXIT_SUCCESS;
}


int check_mpi_error_code(int error_code, char *ops_type, char *parent_func) {
	if (error_code != MPI_SUCCESS) {
		printf("error when %s with error code %d in func %s\n", ops_type, error_code, parent_func);
	}
    return EXIT_SUCCESS;
}

FILE* create_file(char* prefix, char* appdix, char *file_format)
{
	char fn[100];
	int cx = snprintf(fn, 100, "%s.%s.%s", prefix, appdix, file_format);
	if (!(cx >= 0 && cx < 100)) {
		printf("could not create file name!\n");
		exit(EXIT_FAILURE);
	}
	FILE* f;
	f = fopen(fn, "w");
	return f;
}

char *pop_to_json(hgt_pop *p, hgt_params *params) {
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
	for (j = 0; j < p->size; j++) {
		if (j < p->size - 1) {
			bformata(b, "\"%s\",\n", p->genomes[j]->seq);
		}
		else {
			bformata(b, "\"%s\"\n", p->genomes[j]->seq);
		}
	}
	bformata(b, "],\n");

	double **matrix = get_coal_time_matrix(p);
	bformata(b, "\"Ranks\": [\n");
	int k;
	for (j = 0; j < p->size; j++) {
		bformata(b, "[");
		for (k = 0; k < p->size; k++) {
			bformata(b, "%g", matrix[j][k]);
			if (k < p->size - 1) {
				bformata(b, ",");
			}
		}
		bformata(b, "]");
		if (j < p->size - 1) {
			bformata(b, ",");
		}
		bformata(b, "\n");
	}
	bformata(b, "]");
	bformata(b, "\n}");

	c = bstr2cstr(b, '\n');

	bdestroy(b);
	for (j = 0; j < p->size; j++) {
		free(matrix[j]);
	}
	free(matrix);

	return c;
}

short int ** get_coal_rank_matrix(hgt_pop *p) {

	short int ** matrix;
	matrix = (short int **)malloc(p->size * sizeof(short *));
	int i;
	for (i = 0; i < p->size; i++) {
		matrix[i] = get_coal_ranks(p, i);
	}


	return matrix;
}

double **get_coal_time_matrix(hgt_pop *p) {
	double **matrix;
	matrix = (double **)malloc(p->size * sizeof(double*));
	int i;
	for (i = 0; i < p->size; i++) {
		matrix[i] = get_coal_times(p, i);
	}
	return matrix;
}

double *get_coal_times(hgt_pop *p, int index) {
	double *times;
	times = (double *)malloc(p->size * sizeof(double));
	double current_time;
	current_time = hgt_pop_get_time(p);
	int i;
	for (i = 0; i < p->size; i++) {
		hgt_linkage *linkages[2];
		linkages[0] = p->linkages[index];
		linkages[1] = p->linkages[i];
		double time;
		if (i == index) {
			time = 0;
		}
		else {
			time = current_time - hgt_linkage_find_most_rescent_ancestor_time(linkages, 2);
		}
		times[i] = time;
	}
	return times;
}

short int * get_coal_ranks(hgt_pop *p, int index) {
	short int * ranks;
	ranks = (short int *)malloc(p->size * sizeof(short int));
	double *times;
	times = (double *)malloc(p->size * sizeof(double));
	double current_time;
	current_time = hgt_pop_get_time(p);
	int i;
	for (i = 0; i < p->size; i++) {
		hgt_linkage *linkages[2];
		linkages[0] = p->linkages[index];
		linkages[1] = p->linkages[i];
		double time;
		if (i == index) {
			time = 0;
		}
		else {
			time = current_time - hgt_linkage_find_most_rescent_ancestor_time(linkages, 2);
		}
		times[i] = time;
	}

	hgt_utils_pair *pairs;
	pairs = (hgt_utils_pair *)malloc(p->size * sizeof(hgt_utils_pair));
	for (i = 0; i < p->size; i++) {
		pairs[i].index = i;
		pairs[i].value = times[i];
	}

	qsort(pairs, p->size, sizeof(hgt_utils_pair), hgt_utils_compare);

	for (i = 0; i < p->size; i++) {
		ranks[i] = (short int) pairs[i].index;
	}

	free(pairs);
	free(times);
	return ranks;
}


int write_pops(hgt_pop **ps, hgt_params *params, int rank, int numprocs) {
	bstring to_json(hgt_pop *p, hgt_params *params);

	int dest, tag;
	dest = 0;
	tag = 0;

	FILE *fp = create_file(params->prefix, "populations", "json");
	int k;
	for (k = 0; k < params->replicates; k++) {
		hgt_pop *p = ps[k];
		bstring b = to_json(p, params);
		char *rc = bstr2cstr(b, '\n');
		if (rank != 0) {
			MPI_Send(rc, blength(b), MPI_CHAR, dest, tag, MPI_COMM_WORLD);
			free(rc);
		}
		else {
			MPI_Status status;
			int i;
			for (i = 0; i < numprocs; i++) {
				if (i != 0) {
					MPI_Probe(i, tag, MPI_COMM_WORLD, &status);
					int count;
					MPI_Get_count(&status, MPI_CHAR, &count);
					rc = malloc(count * sizeof(char));
					MPI_Recv(rc, count, MPI_CHAR, i, tag, MPI_COMM_WORLD, &status);
				}
				fprintf(fp, "%s", rc);
				free(rc);
			}
		}
		bdestroy(b);
	}

	fclose(fp);

	return EXIT_SUCCESS;
}

bstring to_json(hgt_pop *pop, hgt_params *params) {
	char *c;
	bstring b;
	b = bfromcstr("");
	c = pop_to_json(pop, params);
	bformata(b, "%s", c);
	bformata(b, "\n");
	free(c);
	return b;
}
