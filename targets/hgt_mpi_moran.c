//
//  hgt_mpi_moran.c
//  hgt
//
//  Created by Mingzhi Lin on 1/21/14.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "hgt_pop.h"
#include "bstrlib.h"
#include "hgt_utils.h"
#include <time.h>

hgt_pop ** init(hgt_pop_params * params, int rank, gsl_rng * r);
gsl_rng ** setup_rng(hgt_pop_params * params, int rank);
void evolve(hgt_pop ** pop, hgt_pop_params * params, hgt_pop_sample_func sample_f, hgt_pop_coal_time_func c_time_f, int rank, gsl_rng ** rr);
bstring to_json(hgt_pop ** ps, hgt_pop_params * params);

int main(int argc, char * argv[]) {
    int numprocs, rank, i;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    hgt_pop_params *params = hgt_pop_params_alloc();
    hgt_pop_params_parse(params, argc, argv, "hgt_mpi_moran");
    hgt_pop_params_printf(params, stdout);
    
    gsl_rng **rr = hgt_utils_alloc_gsl_rngs(params->replicates, rank);
    hgt_pop **ps = hgt_utils_alloc_populations(params->replicates, params->size, params->seq_len, rank, rr);
    
    time_t start, end;
    start = clock();
    hgt_utils_batch_evolve_moran(ps, params->replicates, params, rr);
    end = clock();
    
    if (rank == 0) {
        printf("Finish evolution: %ld seconds\n", (end - start)/CLOCKS_PER_SEC);
    }
    
    char * rc;
    if (rank == 0) {
        // write populations into json file
        char * fn;
        FILE * fp;
        asprintf(&fn, "%s_populations.json", params->prefix);
        fp = fopen(fn, "w");
        
        fprintf(fp, "[\n");
        MPI_Status status;
        
        for (i = 0; i < numprocs; i++) {
            if (i == 0) {
                bstring b = to_json(ps, params);
                rc = bstr2cstr(b, '\n');
                bdestroy(b);
            } else {
                MPI_Probe(i, 0, MPI_COMM_WORLD, &status);
                int count;
                MPI_Get_count(&status, MPI_CHAR, &count);
                rc = malloc(count * sizeof(char));
                MPI_Recv(rc, count, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
            }
            
            fprintf(fp, "%s\n", rc);
            if (i < numprocs - 1)
            {
                fprintf(fp, ",\n");
            }
            free(rc);
        }
        fprintf(fp, "]\n");
        fclose(fp);
        free(fn);
    } else {
        bstring b = to_json(ps, params);
        rc = bstr2cstr(b, '\n');
        MPI_Send(rc, blength(b), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        bdestroy(b);
        free(rc);
    }
    
    return EXIT_SUCCESS;
}

hgt_pop ** init(hgt_pop_params * params, int rank, gsl_rng * r) {
    hgt_pop ** pop;
    pop = malloc(params->replicates * sizeof(hgt_pop *));
    int i;
    
    for (i = 0; i < params->replicates; i++) {
        pop[i] = hgt_pop_alloc(params->size, params->seq_len, r);
    }
    return pop;
}

gsl_rng ** setup_rng(hgt_pop_params * params, int rank) {
    gsl_rng ** rr;
    rr = malloc(params->replicates * sizeof(gsl_rng*));
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    int s, i;
    for (i = 0; i < params->replicates; i++) {
        s = rank * (int) params->replicates + i;
        rr[i] = gsl_rng_alloc(T);
        gsl_rng_set(rr[i], s);
    }
    return rr;
}

void evolve(hgt_pop ** pop, hgt_pop_params * params, hgt_pop_sample_func sample_f, hgt_pop_coal_time_func c_time_f, int rank, gsl_rng ** rr) {
    int i, j;
    for (i = 0; i < params->replicates; i ++) {
        for (j = 0; j < params->generations; j++) {
            hgt_pop_evolve(pop[i], params, sample_f, c_time_f, rr[i]);
        }
    }
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