//
//  hgt_mpi_moran_step.c
//  hgt
//
//  Created by Mingzhi Lin on 3/14/14.
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
#ifndef MPI_SUCCESS
    #define MPI_SUCCESS 0
#endif
FILE* create_file(char* prefix, char* appdix, char *fmt);
typedef struct _filecontainer file_container;
struct _filecontainer {
	FILE *p2;
	FILE *p3;
	FILE *p4;
	FILE *cov;
	FILE *ks;
	FILE *t2;
	FILE *q2;
};

file_container *create_file_container(char *prefix) {
	file_container *fc;
	fc = (file_container *)malloc(sizeof(file_container));
	fc->p2 = create_file(prefix, "p2", "txt");
	fc->p3 = create_file(prefix, "p3", "txt");
	fc->p4 = create_file(prefix, "p4", "txt");
	fc->cov = create_file(prefix, "cov", "txt");
	fc->ks = create_file(prefix, "ks", "txt");
	fc->t2 = create_file(prefix, "t2", "txt");
	fc->q2 = create_file(prefix, "q2", "txt");
	return fc;
}

void file_container_close(file_container *fc) {
	fclose(fc->p2);
	fclose(fc->p3);
	fclose(fc->p4);
	fclose(fc->cov);
	fclose(fc->ks);
	fclose(fc->t2);
	fclose(fc->q2);
}

void file_container_destroy(file_container *fc) {
	free(fc);
}

void file_container_flush(file_container *fc) {
	fflush(fc->p2);
	fflush(fc->p3);
	fflush(fc->p4);
	fflush(fc->cov);
	fflush(fc->ks);
	fflush(fc->t2);
	fflush(fc->q2);
}

void file_container_write_headers(file_container *fc) {
	fprintf(fc->p2, "#l\tp00\tp01\tp10\tp11\tp00 var\tp01 var\tp10 var\tp11 var\tsample n\tgenerations\n");
	fprintf(fc->p3, "#l\tp00\tp01\tp10\tp11\tp00 var\tp01 var\tp10 var\tp11 var\tsample n\tgenerations\n");
	fprintf(fc->p4, "#l\tp00\tp01\tp10\tp11\tp00 var\tp01 var\tp10 var\tp11 var\tsample n\tgenerations\n");
	fprintf(fc->cov, "#l\tscov\trcov\tpxpy\ttcov\tscov var\trcov var\tpxpy var\ttcov var\tsample n\tgenerations\n");
	fprintf(fc->ks, "#ks\tvd\tvar ks\tvar vd\tgenerations\n");
	fprintf(fc->t2, "#l\tt2\tt3\tt4\tt2_var\tt3_var\tt4_var\tn\tgeneration\n");
	fprintf(fc->q2, "#l\tt2\tt3\tt4\tt2_var\tt3_var\tt4_var\tn\tgeneration\n");
}

int update_t2(hgt_stat_mean ***t2means, hgt_stat_variance ***t2vars, double *buf, int dim, int max_linkage, int sample_size, double current_time);
int check_mpi_error_code(int error_code, char * ops_type, char * parent_func);

int main(int argc, char *argv[]) {
    int write_pops(hgt_pop **ps, hgt_params *params, int rank, int numprocs);
    int pxy_calc(hgt_stat_mean ***means, hgt_stat_variance ***vars, double *pxy, double *d1, double *d2,hgt_pop **ps, hgt_params *params, int rank, int numprocs, hgt_cov_sample_func sample_func, gsl_rng *r);
    int cov_calc(hgt_stat_mean ***means, hgt_stat_variance ***vars, hgt_pop **ps, hgt_params *params, int rank, int numprocs, gsl_rng *rng);
    int write_pxy(FILE * fp, unsigned  maxl, hgt_stat_mean ***means, hgt_stat_variance ***vars, unsigned long gen);
    int write_cov(FILE * fp, unsigned  maxl, hgt_stat_mean ***means, hgt_stat_variance ***vars, unsigned long gen);
    int write_ks(FILE * fp, unsigned  maxl, hgt_stat_mean ***means, hgt_stat_variance ***vars, unsigned long gen);
    int t2_calc(hgt_stat_mean *** t2means,
                hgt_stat_variance *** t2vars,
                hgt_params * params,
                hgt_pop ** ps,
                hgt_pop_calc_most_recent_coal_func coal_time_func,
                int rank,
                int numprocs,
                const gsl_rng * rng);
    int write_t2(FILE *fp, hgt_stat_mean ***means, hgt_stat_variance ***vars, int dim, int size, unsigned long gen);
    
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
        
    hgt_stat_mean ***p2means;
    hgt_stat_mean ***p3means;
    hgt_stat_mean ***p4means;
    hgt_stat_variance ***p2vars;
    hgt_stat_variance ***p3vars;
    hgt_stat_variance ***p4vars;
    hgt_stat_mean ***covmeans;
    hgt_stat_variance ***covvars;
    hgt_stat_mean ***t2means;
    hgt_stat_variance ***t2vars;
    hgt_stat_mean ***q2means;
    hgt_stat_variance ***q2vars;
    
    int linkage_dim = params->linkage_size + 1;
    if (rank == 0) {
        p2means = hgt_utils_alloc_stat_means(params->maxl, 4);
        p2vars = hgt_utils_alloc_stat_variance(params->maxl, 4);
        p3means = hgt_utils_alloc_stat_means(params->maxl, 4);
        p3vars = hgt_utils_alloc_stat_variance(params->maxl, 4);
        p4means = hgt_utils_alloc_stat_means(params->maxl, 4);
        p4vars = hgt_utils_alloc_stat_variance(params->maxl, 4);
        covmeans = hgt_utils_alloc_stat_means(params->maxl+1, 4);
        covvars = hgt_utils_alloc_stat_variance(params->maxl+1, 4);
        t2means = hgt_utils_alloc_stat_means(linkage_dim, 3);
        t2vars = hgt_utils_alloc_stat_variance(linkage_dim, 3);
        q2means = hgt_utils_alloc_stat_means(linkage_dim, 3);
        q2vars = hgt_utils_alloc_stat_variance(linkage_dim, 3);
    }

	file_container *fc;
    if (rank == 0) {
		fc = create_file_container(params->prefix);
		file_container_write_headers(fc);
    }
    
    double *pxy = malloc(params->maxl*4*params->sample_size*sizeof(double));
    double *d1 = malloc(params->seq_len*sizeof(double));
    double *d2 = malloc(params->seq_len*sizeof(double));
    
    unsigned i;
    time_t start, end;
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
        pxy_calc(p2means, p2vars, pxy, d1, d2, ps, params, rank, numprocs, hgt_cov_sample_p2, rng);
        pxy_calc(p3means, p3vars, pxy, d1, d2, ps, params, rank, numprocs, hgt_cov_sample_p3, rng);
        pxy_calc(p4means, p4vars, pxy, d1, d2, ps, params, rank, numprocs, hgt_cov_sample_p4, rng);
        cov_calc(covmeans, covvars, ps, params, rank, numprocs, rng);
        t2_calc(t2means, t2vars, params, ps, hgt_pop_calc_most_recent_ancestor_time, rank, numprocs, rng);
        t2_calc(q2means, q2vars, params, ps, hgt_pop_calc_most_recent_coal_time, rank, numprocs, rng);
        if (rank == 0) {
            write_pxy(fc->p2, params->maxl, p2means, p2vars, current_generation);
            write_pxy(fc->p3, params->maxl, p3means, p3vars, current_generation);
            write_pxy(fc->p4, params->maxl, p4means, p4vars, current_generation);
            write_cov(fc->cov, params->maxl, covmeans, covvars, current_generation);
            write_ks(fc->ks, params->maxl, covmeans, covvars, current_generation);
            write_t2(fc->t2, t2means, t2vars, linkage_dim, 3, current_generation);
            write_t2(fc->q2, q2means, q2vars, linkage_dim, 3, current_generation);
			file_container_flush(fc);
            
            hgt_utils_clean_stat_means(p2means, params->maxl, 4);
            hgt_utils_clean_stat_means(p3means, params->maxl, 4);
            hgt_utils_clean_stat_means(p4means, params->maxl, 4);
            hgt_utils_clean_stat_variances(p2vars, params->maxl, 4);
            hgt_utils_clean_stat_variances(p3vars, params->maxl, 4);
            hgt_utils_clean_stat_variances(p4vars, params->maxl, 4);
            hgt_utils_clean_stat_means(covmeans, params->maxl+1, 3);
            hgt_utils_clean_stat_variances(covvars, params->maxl+1, 3);
            hgt_utils_clean_stat_means(t2means, linkage_dim, 3);
            hgt_utils_clean_stat_variances(t2vars, linkage_dim, 3);
            hgt_utils_clean_stat_means(q2means, linkage_dim, 3);
            hgt_utils_clean_stat_variances(q2vars, linkage_dim, 3);
        }
        
        end = clock();
        
        if (rank == 0) {
            printf("%d, calculation using time = %ld sec\n", i, (end - start)/CLOCKS_PER_SEC);
        }
        
    }
    
    // write_pops(ps, params, rank, numprocs);
    
    if (rank == 0) {
        hgt_utils_free_stat_means(p2means, params->maxl, 4);
        hgt_utils_free_stat_means(p3means, params->maxl, 4);
        hgt_utils_free_stat_means(p4means, params->maxl, 4);
        hgt_utils_free_stat_variances(p2vars, params->maxl, 4);
        hgt_utils_free_stat_variances(p3vars, params->maxl, 4);
        hgt_utils_free_stat_variances(p4vars, params->maxl, 4);
        hgt_utils_free_stat_means(covmeans, params->maxl+1, 4);
        hgt_utils_free_stat_variances(covvars, params->maxl+1, 4);
        hgt_utils_free_stat_means(t2means, linkage_dim, 3);
        hgt_utils_free_stat_variances(t2vars, linkage_dim, 3);
        hgt_utils_free_stat_means(q2means, linkage_dim, 3);
        hgt_utils_free_stat_variances(q2vars, linkage_dim, 3);
		printf("rank %d: Free mean and variances!\n", rank);
        file_container_close(fc);
		file_container_destroy(fc);
        printf("Succesffully close all files!\n");
    }
    
    free(pxy);
    free(d1);
    free(d2);
    printf("rank %d: Free caches!\n", rank);
    
    hgt_utils_free_populations(ps, params->replicates);
    printf("rank %d: Free populations!\n", rank);
    gsl_rng_free(rng);
    hgt_params_free(params);
    printf("rank %d: Free rng and params!\n", rank);
    
exit:
    MPI_Finalize();
    
    return EXIT_SUCCESS;
}

int pxy_calc(hgt_stat_mean ***means, hgt_stat_variance ***vars, double *pxy, double *d1, double *d2,hgt_pop **ps, hgt_params *params, int rank, int numprocs, hgt_cov_sample_func sample_func, gsl_rng *r) {
	char * func_name = "pxy_calc";
	unsigned count, dest, tag, error_code;
    count = params->maxl * 4;
    dest = 0;
    tag = 0;
	unsigned i;
    for (i = 0; i < params->replicates; i++) {
		unsigned n;
        for (n = 0; n < params->sample_size; n++) {
            hgt_pop_calc_dist(ps[i], d1, d2, 1, sample_func, r);
            hgt_pop_calc_pxy_fft(pxy+(n*params->maxl*4), params->maxl, d1, d2, params->seq_len, 0);
        }
		unsigned l;
        for (l = 0; l < params->maxl; l++) {
			unsigned k;
            for (k = 0; k < 4; k++) {
                pxy[l*4+k] = gsl_stats_mean(pxy+l*4+k, params->maxl*4, params->sample_size);
            }
        }
        
        if (rank != 0) {
            error_code = MPI_Send(pxy, count, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
			check_mpi_error_code(error_code, "MPI_Send", func_name);
        } else {
			int j;
            for (j = 0; j < numprocs; j++) {
                if (j != 0) {
                    error_code = MPI_Recv(pxy, count, MPI_DOUBLE, j, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					check_mpi_error_code(error_code, "MPI_Recv", func_name);
                }
				unsigned l;
                for (l = 0; l < params->maxl; l++) {
					unsigned k;
                    for (k = 0; k < 4; k++) {
                        hgt_stat_mean_increment(means[l][k], pxy[l*4+k]);
                        hgt_stat_variance_increment(vars[l][k], pxy[l*4+k]);
                    }
                }
            }
        }
    }

	return EXIT_SUCCESS;
}

int cov_calc(hgt_stat_mean ***means, hgt_stat_variance ***vars, hgt_pop **ps, hgt_params *params, int rank, int numprocs, gsl_rng *rng) {
	char *func_name = "cov_calc";
	int count, dest, tag, error_code;
    count = params->maxl * 4 + 4; // scov + rcov + pxpy + ks + vd
    dest = 0;
    tag = 0;
    double *buf;
    hgt_cov_result *result ;
    result = hgt_cov_result_alloc(params->maxl);
    buf = malloc(count*sizeof(double));
	unsigned i;
    for (i = 0; i < params->replicates; i++) {
        // calculate covariance result
         hgt_pop_calc_cov(result, ps[i], params->sample_size, rng);
        // hgt_pop_calc_cov_all(result, ps[i]);
        // load buf the result
		unsigned j;
        for (j = 0; j < params->maxl; j++) {
            buf[4*j] = result->scov[j];
            buf[4*j+1] = result->rcov[j];
            buf[4*j+2] = result->pxpy[j];
            buf[4*j+3] = result->tcov[j];
            
        }
        buf[4*params->maxl] = result->ks;
        buf[4*params->maxl + 1] = result->vd;
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
				unsigned k;
                for (k = 0; k < params->maxl; k++) {
                    hgt_stat_mean_increment(means[k][0], buf[4*k]);
                    hgt_stat_mean_increment(means[k][1], buf[4*k+1]);
                    hgt_stat_mean_increment(means[k][2], buf[4*k+2]);
                    hgt_stat_mean_increment(means[k][3], buf[4*k+3]);
                    hgt_stat_variance_increment(vars[k][0], buf[4*k]);
                    hgt_stat_variance_increment(vars[k][1], buf[4*k+1]);
                    hgt_stat_variance_increment(vars[k][2], buf[4*k+2]);
                    hgt_stat_variance_increment(vars[k][3], buf[4*k+3]);
                }
                
                hgt_stat_mean_increment(means[params->maxl][0], buf[4*params->maxl]);
                hgt_stat_mean_increment(means[params->maxl][1], buf[4*params->maxl+1]);
                hgt_stat_variance_increment(vars[params->maxl][0], buf[4*params->maxl]);
                hgt_stat_variance_increment(vars[params->maxl][1], buf[4*params->maxl+1]);
            }
        }
    }
    hgt_cov_result_free(result);
    free(buf);
    return EXIT_SUCCESS;
}

int t2_calc(hgt_stat_mean *** t2means,
             hgt_stat_variance *** t2vars,
             hgt_params * params,
             hgt_pop ** ps,
             hgt_pop_calc_most_recent_coal_func coal_time_func,
             int rank,
             int numprocs,
             const gsl_rng * rng) {
    unsigned i, j, count, dim, dest, tag, linkage_sizes[3], max_linkage, error_code;
	char * func_name = "t2_calc";
    max_linkage = 3;
    dim = 1 + params->linkage_size;
    count = params->sample_size * dim * max_linkage;
    for (i = 0; i < max_linkage; i++) {
        linkage_sizes[i] = i + 2;
    }
    
    double * buf;
    buf = (double *) malloc((count + 1) * sizeof(double));
    
    dest = 0;
    tag = 0;
    
    for (i = 0; i < params->replicates; i++) {
        double current_time = hgt_pop_get_time(ps[i]);
		int k;
        for (k = 0; k < params->linkage_size; k++) {
            for (j = 0; j < max_linkage; j++) {
                
                coal_time_func(ps[i]->locus_linkages[k], ps[i]->size, params->sample_size, buf+(max_linkage * k + j) * params->sample_size, linkage_sizes[j], rng);
            }
        }
        
        for (j = 0; j < max_linkage; j++) {
            coal_time_func(ps[i]->linkages, ps[i]->size, params->sample_size, buf+(max_linkage * params->linkage_size + j) * params->sample_size, linkage_sizes[j], rng);
        }
        if (rank != 0) {
            error_code = MPI_Send(buf, count, MPI_UNSIGNED_LONG, dest, tag, MPI_COMM_WORLD);
			check_mpi_error_code(error_code, "MPI_Send", func_name);
        } else {
			int j;
            for (j = 0; j < numprocs; j++) {
                if (j != 0) {
                    // recieve data from worker nodes.
                    error_code = MPI_Recv(buf, count, MPI_UNSIGNED_LONG, j, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					check_mpi_error_code(error_code, "MPI_Recv", func_name);
                }
                update_t2(t2means, t2vars, buf, max_linkage, dim, params->sample_size, current_time);
            }
        }
    }
    free(buf);
    return EXIT_SUCCESS;
}

int update_t2(hgt_stat_mean ***t2means, hgt_stat_variance ***t2vars, double *buf, int max_linkage, int dim, int sample_size, double current_time) {
    int i, j, k;
    double v, t;
    for (i = 0; i < dim; i++) {
        for (j = 0; j < max_linkage; j++) {
            for (k = 0; k < sample_size; k++) {
                 v =buf[(i*max_linkage + j) * sample_size + k];
                 if (v > 0) {
                    t = current_time - v;
                    hgt_stat_mean_increment(t2means[i][j], (double)t);
                    hgt_stat_variance_increment(t2vars[i][j], (double)t);
                 }
            }
        }
    }
    
    return EXIT_SUCCESS;
}

int write_pops(hgt_pop **ps, hgt_params *params, int rank, int numprocs) {
    bstring to_json(hgt_pop **ps, hgt_params *params);
    
    int dest, tag;
    dest = 0;
    tag = 0;
    
    char *rc;
    bstring b = to_json(ps, params);
    rc = bstr2cstr(b, '\n');
    if (rank != 0) {
        MPI_Send(rc, blength(b), MPI_CHAR, dest, tag, MPI_COMM_WORLD);
        free(rc);
    } else {
		FILE *fp = create_file(params->prefix, "populations", "json");
       
        fprintf(fp, "[\n");
        MPI_Status status;
        
        int i;
        for (i = 0; i < numprocs; i++) {
            if (i != 0) {
                MPI_Probe(i, tag, MPI_COMM_WORLD, &status);
                int count;
                MPI_Get_count(&status, MPI_CHAR, &count);
                rc = malloc(count*sizeof(char));
                MPI_Recv(rc, count, MPI_CHAR, i, tag, MPI_COMM_WORLD, &status);
            }
            
            fprintf(fp, "%s\n", rc);
            if (i < numprocs - 1){
                fprintf(fp, ",\n");
            }
            free(rc);
        }
        
        fprintf(fp, "]\n");
        fclose(fp);
    }
    bdestroy(b);
    
    return EXIT_SUCCESS;
}

bstring to_json(hgt_pop ** ps, hgt_params * params) {
    unsigned i;
    char *c;
    bstring b;
    b = bfromcstr("");
    
    for (i = 0; i < params->replicates; i ++) {
        hgt_pop * pop = ps[i];
        c = hgt_pop_to_json(pop, params);
        bformata(b, "%s", c);
        free(c);
        if (i < params->replicates-1)
        {
            bformata(b, ",\n");
        }
    }
    
    return b;
}

int write_t2(FILE *fp, hgt_stat_mean ***means, hgt_stat_variance ***vars, int dim, int size, unsigned long gen) {
    int i, j;
    double v;
    for (i = 0; i < dim; i++) {
        fprintf(fp, "%d\t", i);
        for (j = 0; j < size; j++) {
            v = hgt_stat_mean_get(means[i][j]);
            fprintf(fp, "%g\t", v);
        }

        for (j = 0; j < size; j++) {
            v = hgt_stat_variance_get(vars[i][j]);
            fprintf(fp, "%g\t", v);
        }
        fprintf(fp, "%ld\t%lu\n", hgt_stat_mean_get_n(means[0][0]), gen);
    }
    fflush(fp);
    return EXIT_SUCCESS;
}

int write_pxy(FILE * fp, unsigned maxl, hgt_stat_mean ***means, hgt_stat_variance ***vars, unsigned long gen) {
    unsigned i, j;
    for (j = 0; j < maxl; j ++) {
        fprintf(fp, "%d\t", j);
        for (i = 0; i < 4; i++) {
            fprintf(fp, "%g\t", hgt_stat_mean_get(means[j][i]));
        }
        for (i = 0; i < 4; i++) {
            fprintf(fp, "%g\t", hgt_stat_variance_get(vars[j][i]));
        }
        fprintf(fp, "%ld\t%lu\n", hgt_stat_mean_get_n(means[j][0]), gen);
    }
    
    fflush(fp);
    return EXIT_SUCCESS;
}

int write_cov(FILE * fp, unsigned maxl, hgt_stat_mean ***means, hgt_stat_variance ***vars, unsigned long gen){
    unsigned i, j;
    for (i = 0; i < maxl; i++) {
        fprintf(fp, "%d\t", i);
        for (j = 0; j < 4; j++) {
            fprintf(fp, "%g\t", hgt_stat_mean_get(means[i][j]));
        }
        for (j = 0; j < 4; j++) {
            fprintf(fp, "%g\t", hgt_stat_variance_get(vars[i][j]));
        }
        fprintf(fp, "%ld\t%lu\n", hgt_stat_mean_get_n(means[i][0]), gen);
    }
    fflush(fp);
    return EXIT_SUCCESS;
}

int write_ks(FILE *fp, unsigned maxl, hgt_stat_mean ***means, hgt_stat_variance ***vars, unsigned long gen) {
    fprintf(fp, "%g\t%g\t%g\t%g\t%lu\t%lu\n", hgt_stat_mean_get(means[maxl][0]), hgt_stat_mean_get(means[maxl][1]), hgt_stat_variance_get(vars[maxl][0]), hgt_stat_variance_get(vars[maxl][1]), hgt_stat_mean_get_n(means[maxl][0]), gen);
    fflush(fp);
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
