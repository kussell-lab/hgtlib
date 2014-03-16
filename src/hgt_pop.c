//
//  hgt_pop.c
//  hgt
//
//  Created by Mingzhi Lin on 1/18/14.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <argtable2.h>
#include <math.h>
#include "hgt_pop.h"
#include "hgt_cov.h"
#include "hgt_corr.h"
#include "bstrlib.h"

const char DNA[5] = "ATGC\0";
const int NUM_DNA_CHAR = 4;

hgt_pop * hgt_pop_alloc(unsigned long size, unsigned long seq_len, const gsl_rng * r) {
    char * ancestor;
    int i;
    int random_seq(char * seq, unsigned long seq_len, const gsl_rng * r);
    
    hgt_pop * p = (hgt_pop *) malloc (sizeof(hgt_pop));
    
    p->size = size;
    p->seq_len = seq_len;
    
    p->genomes = (char **) malloc(size * sizeof(char *));
    
    // random initilize genomes
    ancestor = malloc((seq_len+1) * sizeof(char));
    random_seq(ancestor, seq_len, r);
    for (i = 0; i < size; i ++) {
        p->genomes[i] = malloc((seq_len+1) * sizeof(char));
        strcpy(p->genomes[i], ancestor);
    }
    
    p->survived = NULL;
    p->new_born = NULL;
    p->cache_allocated = 0;
    
    free(ancestor);
    return p;
}

hgt_pop * hgt_pop_copy(hgt_pop * p) {
    int i;
    hgt_pop * new_p = (hgt_pop *) calloc(1, sizeof(hgt_pop));
    new_p->size = p->size;
    new_p->seq_len = p->seq_len;
    new_p->genomes = (char **) malloc(p->size * sizeof(char *));
    for (i = 0; i < p->size; i++) {
        new_p->genomes[i] = malloc((p->seq_len+1) * sizeof(char));
        strcpy(new_p->genomes[i], p->genomes[i]);
    }
    
    // initilize caches
    new_p->survived = (unsigned long *) calloc(p->size, sizeof(unsigned long));
    new_p->new_born = (unsigned long *) calloc(p->size, sizeof(unsigned long));
    
    return new_p;
}

int hgt_pop_free(hgt_pop * p) {
    int i;
    for (i = 0; i < p->size; i ++) {
        free(p->genomes[i]);
    }
    free(p->genomes);
    if (p->cache_allocated != 0) {
        free(p->survived);
        free(p->new_born);
    }
    free(p);
    
    return EXIT_SUCCESS;
    
}

char *hgt_pop_to_json(hgt_pop *p, hgt_pop_params *params){
    bstring b = bfromcstr("");
    bformata(b, "{\n");
    bformata(b, "\"Size\": %ld,\n", params->size);
    bformata(b, "\"Length\": %ld,\n", params->seq_len);
    bformata(b, "\"MutationRate\": %g,\n", params->mu_rate);
    bformata(b, "\"TransferRate\": %g,\n", params->tr_rate);
    bformata(b, "\"FragLen\": %ld,\n", params->frag_len);
    bformata(b, "\"Generation\": %ld,\n", params->generations);
    bformata(b, "\"Genomes\": [\n");
    int j;
    for (j = 0; j < p->size; j ++) {
        if (j < p->size - 1) {
            bformata(b, "\"%s\",\n", p->genomes[j]);
        } else {
            bformata(b, "\"%s\"\n", p->genomes[j]);
        }
    }
    bformata(b, "]\n}\n");
    
    return bstr2cstr(b, '\n');
}

int hgt_pop_params_parse(hgt_pop_params *params, int argc, char **argv, char * progname){
    struct arg_int *size = arg_int1("n", "size", "<unsigned long>", "population size");
    struct arg_int *seq_len = arg_int1("l", "len", "<unsigned long>", "genome length");
    struct arg_int *frag_len = arg_int1("f", "frag", "<unsigned long>", "fragment length");
    struct arg_dbl *mu_rate = arg_dbl1("u", "mu_rate", "<double>", "mutation rate");
    struct arg_dbl *tr_rate = arg_dbl1("t", "tr_rate", "<double>", "transfer rate");
    struct arg_int *gen = arg_int1("g", "gen", "<unsigned long>", "generations");
    
    struct arg_int *spl_size = arg_int0("s", "sample_size", "<unsigned long>", "sample size");
    struct arg_int *spl_time = arg_int0("i", "sample_time", "<unsigned long>", "sample time");
    struct arg_int *repl = arg_int0("r", "replication", "<unsigned long>", "replication");
    struct arg_int *maxl = arg_int0("m", "maxl", "<unsigned long>", "maxl");
    
    struct arg_file *prefix = arg_file1("o", "output", "<output>", "prefix");

    struct arg_lit  *help    = arg_lit0(NULL,"help", "print this help and exit");
    struct arg_end  *end     = arg_end(20);

    void* argtable[] = {size, seq_len, frag_len, mu_rate, tr_rate, gen, spl_time, spl_size, repl, maxl, prefix, help, end};
    int nerrors;
    int exit_code = EXIT_SUCCESS;
    /* verify the argtable[] entries were allocated sucessfully */
    if (arg_nullcheck(argtable) != 0){
        /* NULL entries were detected, some allocations must have failed */
        printf("%s: insufficient memory\n",progname);
        exit_code = EXIT_FAILURE;
        goto exit;
    }

    /* Parse the command line as defined by argtable[] */
    nerrors = arg_parse(argc,argv,argtable);

    /* special case: '--help' takes precedence over error reporting */
    if (help->count > 0)
    {
        printf("Usage: %s", progname);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
        exit_code =  EXIT_FAILURE;
        goto exit;
    }

    /* If the parser returned any errors then display them and exit */
    if (nerrors > 0)
    {
        /* Display the error details contained in the arg_end struct.*/
        arg_print_errors(stdout,end,progname);
        printf("Try '%s --help' for more information.\n",progname);
        exit_code =  EXIT_FAILURE;
        goto exit;
    }

    // begin parse
    params->size = size->ival[0];
    params->seq_len = seq_len->ival[0];
    params->frag_len = frag_len->ival[0];
    params->mu_rate = mu_rate->dval[0];
    params->tr_rate = tr_rate->dval[0];
    params->generations = gen->ival[0];
    params->prefix = (char *) prefix->filename[0];
    
    if (maxl->count > 0) {
        params->maxl = maxl->ival[0];
    } else {
        params->maxl = params->seq_len;
    }

    if (spl_time->count > 0) {
        params->sample_time = spl_time->ival[0];
    } else {
        params->sample_time = 1;
    }
    
    if (spl_size->count > 0) {
        params->sample_size = spl_size->ival[0];
    } else {
        params->sample_size = 100;
    }
    
    if (repl->count > 0) {
        params->replicates = repl->ival[0];
    } else {
        params->replicates = 1;
    }
    
    exit:
        arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0])); 
    return exit_code;
}

int hgt_pop_params_free(hgt_pop_params *params){
    free(params);

    return EXIT_SUCCESS;
}

int hgt_pop_params_printf(hgt_pop_params *params, FILE *stream) {
    fprintf(stream, "population size = %lu\n", params->size);
    fprintf(stream, "genome length = %lu\n", params->seq_len);
    fprintf(stream, "mutation rate = %g\n", params->mu_rate);
    fprintf(stream, "transfer rate = %g\n", params->tr_rate);
    fprintf(stream, "transfer frag = %lu\n", params->frag_len);
    fprintf(stream, "generations = %lu\n", params->generations);
    fprintf(stream, "sample size = %lu\n", params->sample_size);
    fprintf(stream, "sample time = %lu\n", params->sample_time);
    fprintf(stream, "replicates = %lu\n", params->replicates);
    fprintf(stream, "prefix = %s\n", params->prefix);
    fprintf(stream, "maxl = %lu\n", params->maxl);
    return EXIT_SUCCESS;
}

hgt_pop_params *hgt_pop_params_alloc() {
    hgt_pop_params *params = calloc(1, sizeof(hgt_pop_params));
    return params;
}

int hgt_pop_mutate_at(hgt_pop *p, unsigned long g, unsigned long s, const gsl_rng *r) {
    int exit_code;
    exit_code = EXIT_SUCCESS;

    char c;
    c = DNA[gsl_rng_uniform_int(r, NUM_DNA_CHAR)];
    while (c == p->genomes[g][s]) {
        c = DNA[gsl_rng_uniform_int(r, NUM_DNA_CHAR)];
    }
    p->genomes[g][s] = c;

    return exit_code;
}

int hgt_pop_mutate(hgt_pop *p, const gsl_rng *r) {
    int exit_code;
    exit_code = EXIT_SUCCESS;
    unsigned long g, s;
    g = gsl_rng_uniform_int(r, p->size);
    s = gsl_rng_uniform_int(r, p->seq_len);
    hgt_pop_mutate_at(p, g, s, r);
    return exit_code;
}

int hgt_pop_transfer_at(hgt_pop *p, 
                        unsigned long donor, 
                        unsigned long receiver, 
                        unsigned long frag_len, 
                        unsigned long start)
{
    int exit_code;
    exit_code = EXIT_SUCCESS;

    if (start + frag_len < p->seq_len) {
        strncpy(p->genomes[receiver]+start, p->genomes[donor]+start, frag_len);
    } else {
        strcpy(p->genomes[receiver]+start, p->genomes[donor]+start);
        strncpy(p->genomes[receiver], p->genomes[donor], start + frag_len - p->seq_len);
    }
    return exit_code;
}

int hgt_pop_transfer(hgt_pop *p, unsigned long frag_len, const gsl_rng *r) {
    unsigned long donor, receiver, start;
    donor = gsl_rng_uniform_int(r, p->size);
    receiver = gsl_rng_uniform_int(r, p->size);
    if (donor != receiver) {
        start = gsl_rng_uniform_int(r, p->seq_len);
        hgt_pop_transfer_at(p, donor, receiver, frag_len, start);
    }

    return EXIT_SUCCESS;
}

int hgt_pop_sample_moran(hgt_pop *p, const gsl_rng *r) {
    unsigned long a, b;
    a = gsl_rng_uniform_int(r, p->size);
    b = gsl_rng_uniform_int(r, p->size);
    if (a != b) {
        strcpy(p->genomes[b], p->genomes[a]);
    }
    return EXIT_SUCCESS;
}

int hgt_pop_sample_wf(hgt_pop *p, const gsl_rng *r) {
    int i, j, k;
    if (p->cache_allocated == 0) {
        p->survived = (unsigned long *) calloc(p->size, sizeof(unsigned long));
        p->new_born = (unsigned long *) calloc(p->size, sizeof(unsigned long));
        p->cache_allocated = 1;
    } else {
        for (i = 0; i < p->size; i++) {
            p->survived[i] = 0;
        }
    }

    k = 0;
    for (i = 0; i < p->size; i ++) {
        j = gsl_rng_uniform_int(r, p->size);
        if (p->survived[j] == 1) {
            p->new_born[k] = j;
            k++;
        } else {
            p->survived[j] = 1;
        }
    }
    
    k = 0;
    for (i = 0; i < p->size; i ++) {
        if (p->survived[i] != 1) {
            strcpy(p->genomes[i], p->genomes[p->new_born[k]]);
            k++;
        }
    }

    return EXIT_SUCCESS;

}

int hgt_pop_evolve(hgt_pop *p, hgt_pop_params *params, hgt_pop_sample_func sample_f, hgt_pop_coal_time_func c_time_f, const gsl_rng *r) {
    double mu;
    unsigned long count, k;
    sample_f(p, r);
    mu = c_time_f(p->size, r) * (params->mu_rate + params->tr_rate) * (double) p->seq_len * (double) p->size;
    count = gsl_ran_poisson(r, mu);
    for (k = 0; k < count; k ++) {
        if (gsl_rng_uniform(r) < params->mu_rate/(params->mu_rate + params->tr_rate)) {
            hgt_pop_mutate(p, r);
        } else {
            hgt_pop_transfer(p, params->frag_len, r);
        }
    }
    p->generation++;
    return EXIT_SUCCESS;
}

// evolve with exponential distributed fragment sizes
int hgt_pop_evolve_expon_frag(hgt_pop *p, hgt_pop_params *params, hgt_pop_sample_func sample_f, hgt_pop_coal_time_func c_time_f, const gsl_rng *r) {
    double mu;
    unsigned long count, k, frag_len;
    sample_f(p, r);
    mu = c_time_f(p->size, r) * (params->mu_rate + params->tr_rate) * (double) p->seq_len * (double) p->size;
    count = gsl_ran_poisson(r, mu);
    for (k = 0; k < count; k ++) {
        if (gsl_rng_uniform(r) < params->mu_rate/(params->mu_rate + params->tr_rate)) {
            hgt_pop_mutate(p, r);
        } else {
            frag_len = gsl_ran_exponential(r, (double) params->frag_len);
            hgt_pop_transfer(p, frag_len, r);
        }
    }
    p->generation++;
    return EXIT_SUCCESS;
}


double hgt_pop_coal_time_moran(unsigned long p_size, const gsl_rng *r) {
    return gsl_ran_exponential(r, 1.0 / (double) p_size);
}

double hgt_pop_coal_time_wf(unsigned long p_size, const gsl_rng *r) {
    return gsl_ran_exponential(r, 1.0);
}

double hgt_pop_calc_ks(hgt_pop *p) {
    double ks;
    int i, j, k;
    ks = 0;
    
    for (i = 0; i < p->size; i++) {
        for (j = i + 1; j < p->size; j++) {
            for (k = 0; k < p->seq_len; k++) {
                if (p->genomes[i][k] != p->genomes[j][k]) {
                    ks++;
                }
            }
        }
    }
    
    ks /= (double) (p->size * (p->size - 1)/2 * p->seq_len);
    
    return ks;
}

int hgt_pop_calc_dist(hgt_pop *p, double *ds1, double *ds2, unsigned long sample_size, hgt_cov_sample_func sample_func, const gsl_rng *r) {
    int i, j;
    unsigned long a, b, c, d;
    for (i = 0; i < p->seq_len; i++) {
        ds1[i] = 0;
        ds2[i] = 0;
    }
    
    for (i = 0; i < sample_size; i++) {
        sample_func(&a, &b, &c, &d, p->size, r);
        for (j = 0; j < p->seq_len; j++) {
            if (p->genomes[a][j] != p->genomes[b][j]) {
                ds1[j]++;
            }
            if (p->genomes[c][j] != p->genomes[d][j]) {
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

int hgt_pop_calc_pxy(double *pxy, unsigned long maxl, double *ds1, double *ds2, unsigned long len, int circular) {
    int i, l;
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
                if (i-l < 0) {
                    continue;
                }
                a = i;
                b = i - l;
                num++;
            }
            pxy[l*4+3] += ds1[a] * ds2[b];
//            pxy[l][3] += ds2[a] * ds1[b];
            pxy[l*4+2] += ds1[a] * (1-ds2[b]);
//            pxy[l][2] += ds2[a] * (1-ds1[b]);
            pxy[l*4+1] += (1-ds1[a]) * ds2[b];
//            pxy[l][1] += (1-ds2[a]) * ds1[b];
            pxy[l*4+0] += (1-ds1[a]) * (1-ds2[b]);
//            pxy[l][0] += (1-ds2[a]) * (1-ds1[b]);
        }
        pxy[l*4+3] /= (double) num;
        pxy[l*4+2] /= (double) num;
        pxy[l*4+1] /= (double) num;
        pxy[l*4+0] /= (double) num;
    }
    
    return EXIT_SUCCESS;
}

int hgt_pop_calc_pxy_fft(double *pxy, unsigned long maxl, double *d1, double *d2, unsigned long len, int circular) {
    unsigned long next_power2(unsigned long len);
    unsigned long fft_len;
    int i, j, l;
    fft_len = next_power2(len);
    
    double ds1[2][2*fft_len];
    double ds2[2][2*fft_len];
    double mask[2*fft_len];
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
    
    double buf1[2*fft_len];
    double buf2[2*fft_len];
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
        }
    }
    
    return EXIT_SUCCESS;
}

int hgt_pop_calc_cov(hgt_cov_result *result, hgt_pop *p, int sample, const gsl_rng* rng) {
    // allocate a binary matrix
    short **matrix = malloc(sample*sizeof(short*));
    int i, j;
    unsigned long a, b;
    
    for (i = 0; i < sample; i++) {
        // randomly choose two distinct genomes for comparison.
        a = gsl_rng_uniform_int(rng, p->size);
        b = gsl_rng_uniform_int(rng, p->size);
        while (a == b) {
            b = gsl_rng_uniform_int(rng, p->size);
        }
        
        // do comparison to binary sequence.
        matrix[i] = malloc(p->seq_len*sizeof(short));
        for (j = 0; j < p->seq_len; j++) {
            if (p->genomes[a][j] != p->genomes[b][j]) {
                matrix[i][j] = 1;
            } else {
                matrix[i][j] = 0;
            }
        }
    }
    hgt_cov_result_calc_matrix(result, matrix, sample, p->seq_len);
    
    for (i = 0; i < sample; i++) {
        free(matrix[i]);
    }
    free(matrix);
    
    return EXIT_SUCCESS;
}

/******** PRIVATE FUNCTIONS ***********/
int random_seq(char * seq, unsigned long seq_len, const gsl_rng * r) {
    int i;
    for (i = 0; i < seq_len; i ++) {
        seq[i] = DNA[gsl_rng_uniform_int(r, 4)];
    }
    seq[i] = '\0';
    return 0;
}

unsigned long next_power2(unsigned long len) {
    unsigned long v;
    v = (unsigned long) pow(2, (int)(log(len)/log(2))+1);
    if (v % len == 0) {
        v = len;
    }
    return v;
}
