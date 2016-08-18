//
//  hgt_mpi_simu.c
//  hgt
//
//  Created by Mingzhi Lin on 8/18/16.
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

file_container *create_file_container(char *prefix) {
	file_container *fc;
	fc = (file_container *)malloc(sizeof(file_container));
	fc->cov = create_file(prefix, "cov", "txt");
	return fc;
}

void file_container_close(file_container *fc) {
	fclose(fc->cov);
}

void file_container_destroy(file_container *fc) {
	free(fc);
}

void file_container_flush(file_container *fc) {
	fflush(fc->cov);
}

void file_container_write_headers(file_container *fc) {
	fprintf(fc->cov, "l,m,v,n,t,g,c\n");
}

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

	file_container *fc;
    if (rank == 0) {
		fc = create_file_container(params->prefix);
		file_container_write_headers(fc);
    }
    
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
    
    for (i = 0; i < params->sample_time; i++) {
        unsigned long current_generation;
        current_generation = (i + 1) * params->generations;
        
        start = clock();

        hgt_utils_batch_evolve(ps, params->replicates, params, sample_f, frag_f, rng);
        
        end = clock();
        
        if (rank == 0) {
            printf("%d, evolve using time = %ld sec\n", i, (end - start)/CLOCKS_PER_SEC);
        }
        
        start = clock();
		
        if (params->sample_bias != 1) {
            params->cluster_num = 1;
        }
        int k;
        for (k = 0; k < params->cluster_num; k++) {
            hgt_simu_results *results;
            if (rank == 0) {
                results = hgt_simu_results_new(params->maxl);
            }

            cov_calc(results, ps, params, k, rank, numprocs, rng);
            if (rank == 0) {
                write_cov(fc->cov, results, current_generation, k);
                file_container_flush(fc);
            }

            if (rank == 0) {
                hgt_simu_results_free(results, params->maxl);
            }
        }
        
        
        end = clock();
        
        if (rank == 0) {
            printf("%d, calculation using time = %ld sec\n", i, (end - start)/CLOCKS_PER_SEC);
        }
        
    }
    
    
    if (rank == 0) {
        file_container_close(fc);
		file_container_destroy(fc);
    }
    
    // hgt_utils_free_populations(ps, params->replicates);
    gsl_rng_free(rng);
    hgt_params_free(params);
    
exit:
    MPI_Finalize();
    
    return EXIT_SUCCESS;
}

int cov_calc(hgt_simu_results *results, hgt_pop **ps, hgt_params *params, int cluster_num, int rank, int numprocs, gsl_rng *rng) {
	char *func_name = "cov_calc";
	int count, dest, tag, error_code;
    count = hgt_cov_result_length(params->maxl); 
    dest = 0;
    tag = 0;
    double *buf;
    hgt_cov_result *result ;
    result = hgt_cov_result_alloc(params->maxl);
    buf = malloc(count*sizeof(double));
	unsigned i;
    for (i = 0; i < params->replicates; i++) {
        // sample genomes and calculate correlations.
        if (params->sample_bias == 1) {
            int sample_size1, cluster_num1, sample_size2, cluster_num2;
            sample_size1 = params->sample_size;
            cluster_num1 = 1;
            sample_size2 = params->cluster_size;
            cluster_num2 = cluster_num;
            hgt_pop_calc_cov_bias(result, ps[i], sample_size1, cluster_num1, sample_size2, cluster_num2, rng);
        } else {
            hgt_pop_calc_cov(result, ps[i], params->sample_size, rng);
        }
        
        hgt_cov_result_to_array(buf, result, params->maxl);

        if (rank != 0) {
            error_code = MPI_Send(buf, count, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
			check_mpi_error_code(error_code, "MPI_Send", func_name);
        } else {
			int j;
            for (j = 0; j < numprocs; j++) {
                if (j != 0) {
                    error_code = MPI_Recv(buf, count, MPI_DOUBLE, j, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					check_mpi_error_code(error_code, "MPI_Recv", func_name);
                }

                hgt_cov_array_to_result(buf, result, params->maxl);
                hgt_simu_results_increment(results, result);
            }
        }
    }

    hgt_cov_result_free(result);
    free(buf);

    return EXIT_SUCCESS;
}

int write_cov(FILE *fp, hgt_simu_results *results, unsigned long gen, int cluster_num) {
    // write l, m, v, n, g
    int write_cov_one(FILE *fp, hgt_stat_meanvar **meanvars, int i, char *t, unsigned long gen, int cluster_num);
    write_cov_one(fp, results->ks->meanvars, 0, "Ks", gen, cluster_num);
    write_cov_one(fp, results->vd->meanvars, 0, "Vd", gen, cluster_num);
    int i;
    for (i = 0; i < results->maxl; i++) {
        write_cov_one(fp, results->cs->meanvars, i, "Cs", gen, cluster_num);
        write_cov_one(fp, results->cr->meanvars, i, "Cr", gen, cluster_num);
        write_cov_one(fp, results->cm->meanvars, i, "Cm", gen, cluster_num);
        write_cov_one(fp, results->cm2->meanvars, i, "Cm2", gen, cluster_num);
    }
    return EXIT_SUCCESS;
}

int write_cov_one(FILE *fp, hgt_stat_meanvar **meanvars, int i, char *t, unsigned long gen, int cluster_num) {
    fprintf(fp, "%d,%g,%g,%lu,%s,%lu,%d\n", i, hgt_stat_mean_get(meanvars[i]->mean), hgt_stat_variance_get(meanvars[i]->var), hgt_stat_mean_get_n(meanvars[i]->mean), t, gen, cluster_num);
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

hgt_simu_results * hgt_simu_results_new(int maxl) {
    hgt_simu_results *r = (hgt_simu_results *) malloc(sizeof(hgt_simu_results));
    r->maxl = maxl;
    r->cs = hgt_stat_meanvar_list_new(maxl);
    r->cr = hgt_stat_meanvar_list_new(maxl);
    r->cm = hgt_stat_meanvar_list_new(maxl);
    r->cm2 = hgt_stat_meanvar_list_new(maxl);
    r->ks = hgt_stat_meanvar_list_new(1);
    r->vd = hgt_stat_meanvar_list_new(1);
    return r;
}

int hgt_simu_results_free(hgt_simu_results *r, int maxl) {
    hgt_stat_meanvar_list_destroy(r->cs);
    hgt_stat_meanvar_list_destroy(r->cr);
    hgt_stat_meanvar_list_destroy(r->cm);
    hgt_stat_meanvar_list_destroy(r->cm2);
    hgt_stat_meanvar_list_destroy(r->ks);
    hgt_stat_meanvar_list_destroy(r->vd);
    free(r);

    return EXIT_SUCCESS;
}

int hgt_simu_results_increment(hgt_simu_results *r, hgt_cov_result *cov_result) {
    int i;
    for (i = 0; i < r->maxl; i++) {
        hgt_stat_meanvar_list_increment(r->cs, i, cov_result->scov[i]);
        hgt_stat_meanvar_list_increment(r->cm, i, cov_result->tcov[i]);
        hgt_stat_meanvar_list_increment(r->cr, i, cov_result->rcov[i]);
        hgt_stat_meanvar_list_increment(r->cm2, i, cov_result->tcov[i]/cov_result->ks);
    }

    hgt_stat_meanvar_list_increment(r->ks, 0, cov_result->ks);
    hgt_stat_meanvar_list_increment(r->vd, 0, cov_result->vd);

    return EXIT_SUCCESS;
}