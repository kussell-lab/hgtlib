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
#include "hgt_predict.h"

int main(int argc, char *argv[]) {
    int write_pops(hgt_pop **ps, hgt_pop_params *params, int rank, int numprocs);
    int pxy_calc(hgt_stat_mean ***means, hgt_stat_variance ***vars, double *pxy, double *d1, double *d2,hgt_pop **ps, hgt_pop_params *params, int rank, int numprocs, hgt_cov_sample_func sample_func, gsl_rng *r);
    int cov_calc(hgt_stat_mean ***means, hgt_stat_variance ***vars, hgt_pop **ps, hgt_pop_params *params, int rank, int numprocs, gsl_rng *rng);
    int write_pxy(FILE * fp, unsigned long maxl, hgt_stat_mean ***means, hgt_stat_variance ***vars, unsigned long gen);
    int write_cov(FILE * fp, unsigned long maxl, hgt_stat_mean ***means, hgt_stat_variance ***vars, unsigned long gen);
    int write_ks(FILE * fp, unsigned long maxl, hgt_stat_mean ***means, hgt_stat_variance ***vars, unsigned long gen);
    
    int numprocs, rank, exit_code;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    hgt_pop_params *params = hgt_pop_params_alloc();
    exit_code = hgt_pop_params_parse(params, argc, argv, "hgt_mpi_moran_const");
    if (exit_code == EXIT_FAILURE) {
        goto exit;
    }
    
    if (rank == 0) {
        hgt_pop_params_printf(params, stdout);
    }
    
    const gsl_rng_type *T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    int seed = rank;
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
    
    if (rank == 0) {
        p2means = hgt_utils_alloc_stat_means(params->maxl, 4);
        p2vars = hgt_utils_alloc_stat_variance(params->maxl, 4);
        p3means = hgt_utils_alloc_stat_means(params->maxl, 4);
        p3vars = hgt_utils_alloc_stat_variance(params->maxl, 4);
        p4means = hgt_utils_alloc_stat_means(params->maxl, 4);
        p4vars = hgt_utils_alloc_stat_variance(params->maxl, 4);
        covmeans = hgt_utils_alloc_stat_means(params->maxl+1, 3);
        covvars = hgt_utils_alloc_stat_variance(params->maxl+1, 3);
    }
    
    FILE * fp2; // output p2
    FILE * fp3; // output p3
    FILE * fp4; // output p4
    FILE * fpcov; // output cov
    FILE * fpks;
    if (rank == 0) {
        char * fn;
        asprintf(&fn, "%s.p2.txt", params->prefix);
        fp2 = fopen(fn, "w");
        fprintf(fp2, "#l\tp00\tp01\tp10\tp11\tp00 var\tp01 var\tp10 var\tp11 var\tsample n\tgenerations\n");
        
        asprintf(&fn, "%s.p3.txt", params->prefix);
        fp3 = fopen(fn, "w");
        fprintf(fp3, "#l\tp00\tp01\tp10\tp11\tp00 var\tp01 var\tp10 var\tp11 var\tsample n\tgenerations\n");
        
        asprintf(&fn, "%s.p4.txt", params->prefix);
        fp4 = fopen(fn, "w");
        fprintf(fp4, "#l\tp00\tp01\tp10\tp11\tp00 var\tp01 var\tp10 var\tp11 var\tsample n\tgenerations\n");
        
        asprintf(&fn, "%s.cov.txt", params->prefix);
        fpcov = fopen(fn, "w");
        fprintf(fpcov, "#l\tscov\trcov\tpxpy\tscov var\trcov var\tpxpy var\tsample n\tgenerations\n");
        
        asprintf(&fn, "%s.ks.txt", params->prefix);
        fpks = fopen(fn, "w");
        fprintf(fpks, "ks\tvd\tvar ks\tvar vd\tgenerations\n");
        
        free(fn);
    }
    
    double *pxy = malloc(params->maxl*4*params->sample_size*sizeof(double));
    double *d1 = malloc(params->seq_len*sizeof(double));
    double *d2 = malloc(params->seq_len*sizeof(double));
    
    int i;
    time_t start, end;
    for (i = 0; i < params->sample_time; i++) {
        start = clock();
        hgt_utils_batch_evolve_moran(ps, params->replicates, params, rng);
        pxy_calc(p2means, p2vars, pxy, d1, d2, ps, params, rank, numprocs, hgt_cov_sample_p2, rng);
        pxy_calc(p3means, p3vars, pxy, d1, d2, ps, params, rank, numprocs, hgt_cov_sample_p3, rng);
        pxy_calc(p4means, p4vars, pxy, d1, d2, ps, params, rank, numprocs, hgt_cov_sample_p4, rng);
        
        cov_calc(covmeans, covvars, ps, params, rank, numprocs, rng);
        
        if (rank == 0) {
            printf("Ks = %g, Expected = %g, Std Err = %g\n", hgt_stat_mean_get(p2means[0][3]), hgt_predict_ks_moran(params->size, params->mu_rate, params->tr_rate, params->frag_len), sqrt(hgt_stat_variance_get(p2vars[0][3])/(double)hgt_stat_variance_get_n(p2vars[0][3])));
            write_pxy(fp2, params->maxl, p2means, p2vars, (i+1)*params->generations);
            write_pxy(fp3, params->maxl, p3means, p3vars, (i+1)*params->generations);
            write_pxy(fp4, params->maxl, p4means, p4vars, (i+1)*params->generations);
            write_cov(fpcov, params->maxl, covmeans, covvars, (i+1)*params->generations);
            write_ks(fpks, params->maxl, covmeans, covvars, (i+1)*params->generations);
            
            hgt_utils_clean_stat_means(p2means, params->maxl, 4);
            hgt_utils_clean_stat_means(p3means, params->maxl, 4);
            hgt_utils_clean_stat_means(p4means, params->maxl, 4);
            hgt_utils_clean_stat_variances(p2vars, params->maxl, 4);
            hgt_utils_clean_stat_variances(p3vars, params->maxl, 4);
            hgt_utils_clean_stat_variances(p4vars, params->maxl, 4);
            hgt_utils_clean_stat_means(covmeans, params->maxl+1, 3);
            hgt_utils_clean_stat_variances(covvars, params->maxl+1, 3);
        }
        
        end = clock();
        if (rank == 0) {
            printf("rank %d, %d using time = %ld sec\n", rank, i, (end - start)/CLOCKS_PER_SEC);
        }
    }
    
    write_pops(ps, params, rank, numprocs);
    
    if (rank == 0) {
        
        hgt_utils_free_stat_means(p2means, params->maxl, 4);
        hgt_utils_free_stat_means(p3means, params->maxl, 4);
        hgt_utils_free_stat_means(p4means, params->maxl, 4);
        hgt_utils_free_stat_variances(p2vars, params->maxl, 4);
        hgt_utils_free_stat_variances(p3vars, params->maxl, 4);
        hgt_utils_free_stat_variances(p4vars, params->maxl, 4);
        hgt_utils_free_stat_means(covmeans, params->maxl+1, 3);
        hgt_utils_free_stat_variances(covvars, params->maxl+1, 3);
        fclose(fp2);
        fclose(fp3);
        fclose(fp4);
        fclose(fpcov);
        fclose(fpks);
    }
    
    free(pxy);
    free(d1);
    free(d2);
    
    hgt_utils_free_populations(ps, params->replicates);
    gsl_rng_free(rng);
    hgt_pop_params_free(params);
    
exit:
    MPI_Finalize();
    
    return EXIT_SUCCESS;
}

int pxy_calc(hgt_stat_mean ***means, hgt_stat_variance ***vars, double *pxy, double *d1, double *d2,hgt_pop **ps, hgt_pop_params *params, int rank, int numprocs, hgt_cov_sample_func sample_func, gsl_rng *r) {
    int i, j, l, k, n, count, dest, tag;
    count = params->maxl * 4;
    dest = 0;
    tag = 0;
    for (i = 0; i < params->replicates; i++) {
        for (n = 0; n < params->sample_size; n++) {
            hgt_pop_calc_dist(ps[i], d1, d2, 1, sample_func, r);
            hgt_pop_calc_pxy_fft(pxy+(n*params->maxl*4), params->maxl, d1, d2, params->seq_len, 1);
        }
        
        for (l = 0; l < params->maxl; l++) {
            pxy[l] = gsl_stats_mean(pxy+l, params->maxl*4, params->sample_size);
        }
        
        if (rank != 0) {
            MPI_Send(pxy, count, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        } else {
            for (j = 0; j < numprocs; j++) {
                if (j != 0) {
                    MPI_Recv(pxy, count, MPI_DOUBLE, j, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                
                for (l = 0; l < params->maxl; l++) {
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

int cov_calc(hgt_stat_mean ***means, hgt_stat_variance ***vars, hgt_pop **ps, hgt_pop_params *params, int rank, int numprocs, gsl_rng *rng) {
    int count, dest, tag;
    count = params->maxl * 3 + 3; // scov + rcov + pxpy + ks + vd
    dest = 0;
    tag = 0;
    double *buf;
    hgt_cov_result *result ;
    result = hgt_cov_result_alloc(params->maxl);
    buf = malloc((params->maxl*3+3)*sizeof(double));
    int i, j, k;
    for (i = 0; i < params->replicates; i++) {
        // calculate covariance result
        hgt_pop_calc_cov(result, ps[i], params->sample_size, rng);
        // load buf the result
        for (j = 0; j < params->maxl; j++) {
            buf[3*j] = result->scov[j];
            buf[3*j+1] = result->rcov[j];
            buf[3*j+2] = result->pxpy[j];
            
        }
        buf[3*params->maxl] = result->ks;
        buf[3*params->maxl + 1] = result->vd;
        
        if (rank != 0) {
            MPI_Send(buf, count, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
            
        } else {
            for (j = 0; j < numprocs; j++) {
                if (j != 0) {
                    MPI_Recv(buf, count, MPI_DOUBLE, j, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                
                for (k = 0; k < params->maxl; k++) {
                    hgt_stat_mean_increment(means[k][0], buf[3*k]);
                    hgt_stat_mean_increment(means[k][1], buf[3*k+1]);
                    hgt_stat_mean_increment(means[k][2], buf[3*k+2]);
                    hgt_stat_variance_increment(vars[k][0], buf[3*k]);
                    hgt_stat_variance_increment(vars[k][1], buf[3*k+1]);
                    hgt_stat_variance_increment(vars[k][2], buf[3*k+2]);
                }
                
                hgt_stat_mean_increment(means[params->maxl][0], buf[3*params->maxl]);
                hgt_stat_mean_increment(means[params->maxl][1], buf[3*params->maxl+1]);
                hgt_stat_variance_increment(vars[params->maxl][0], buf[3*params->maxl]);
                hgt_stat_variance_increment(vars[params->maxl][1], buf[3*params->maxl+1]);
            }
        }
    }
    hgt_cov_result_free(result);
    free(buf);
    return EXIT_SUCCESS;
}

int write_pops(hgt_pop **ps, hgt_pop_params *params, int rank, int numprocs) {
    bstring to_json(hgt_pop **ps, hgt_pop_params *params);
    
    int dest, tag;
    dest = 0;
    tag = 0;
    
    char *rc;
    bstring b = to_json(ps, params);
    rc = bstr2cstr(b, '\n');
    if (rank != 0) {
        MPI_Send(rc, blength(b), MPI_CHAR, dest, tag, MPI_COMM_WORLD);
        bdestroy(b);
        free(rc);
    } else {
        char *fn;
        FILE *fp;
        asprintf(&fn, "%s_populations.json", params->prefix);
        fp = fopen(fn, "w");
        
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
        free(fn);
    }
    
    return EXIT_SUCCESS;
}

bstring to_json(hgt_pop ** ps, hgt_pop_params * params) {
    int i;
    bstring b;
    b = bfromcstr("");
    
    for (i = 0; i < params->replicates; i ++) {
        hgt_pop * pop = ps[i];
        bformata(b, "%s", hgt_pop_to_json(pop, params));
        if (i < params->replicates-1)
        {
            bformata(b, ",\n");
        }
    }
    
    return b;
}

int write_pxy(FILE * fp, unsigned long maxl, hgt_stat_mean ***means, hgt_stat_variance ***vars, unsigned long gen) {
    int i, j;
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

int write_cov(FILE * fp, unsigned long maxl, hgt_stat_mean ***means, hgt_stat_variance ***vars, unsigned long gen){
    int i, j;
    for (i = 0; i < maxl; i++) {
        fprintf(fp, "%d\t", i);
        for (j = 0; j < 3; j++) {
            fprintf(fp, "%g\t", hgt_stat_mean_get(means[i][j]));
        }
        for (j = 0; j < 3; j++) {
            fprintf(fp, "%g\t", hgt_stat_variance_get(vars[i][j]));
        }
        fprintf(fp, "%ld\t%lu\n", hgt_stat_mean_get_n(means[i][0]), gen);
    }
    fflush(fp);
    return EXIT_SUCCESS;
}

int write_ks(FILE *fp, unsigned long maxl, hgt_stat_mean ***means, hgt_stat_variance ***vars, unsigned long gen) {
    fprintf(fp, "%g\t%g\t%g\t%g\t%lu\t%lu\n", hgt_stat_mean_get(means[maxl][0]), hgt_stat_mean_get(means[maxl][1]), hgt_stat_variance_get(vars[maxl][0]), hgt_stat_variance_get(vars[maxl][1]), hgt_stat_mean_get_n(means[maxl][0]), gen);
    return EXIT_SUCCESS;
}