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
#include <argtable2.h>
#include "hgt_pop.h"

const char * DNA = "ATGC\0";

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
    
    // initilize caches
    p->survived = (unsigned long *) calloc(size, sizeof(unsigned long));
    p->new_born = (unsigned long *) calloc(size, sizeof(unsigned long));
    
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
    free(p->survived);
    free(p->new_born);
    free(p);
    
    return EXIT_SUCCESS;
    
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

    struct arg_lit  *help    = arg_lit0(NULL,"help", "print this help and exit");
    struct arg_end  *end     = arg_end(20);

    void* argtable[] = {size, seq_len, frag_len, mu_rate, tr_rate, gen, spl_time, spl_size, repl, help, end};
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
        exit_code =  EXIT_SUCCESS;
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

    if (spl_time->count > 0)
    {
        params->sample_time = spl_time->ival[0];
    }
    if (spl_size->count > 0)
    {
        params->sample_size = spl_size->ival[0];
    }
    if (repl->count > 0)
    {
        params->replicates = repl->ival[0];
    }
    
    exit:
        arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0])); 
    return exit_code;
}

int hgt_pop_params_free(hgt_pop_params *params){
    free(params);

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
