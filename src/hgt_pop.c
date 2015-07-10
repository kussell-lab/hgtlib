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
#include "hgt_utils.h"
#include "bstrlib.h"
#include "ini.h"
#include <math.h>

const char DNA[5] = "ATGC\0";
const int NUM_DNA_CHAR = 4;

unsigned long search_region(unsigned long pos, unsigned long** region, int num_regions, int inside);

hgt_pop * hgt_pop_alloc(hgt_pop_params *params, const gsl_rng * r) {
    char * ancestor;
    int i, j;
    int random_seq(char * seq, unsigned long seq_len, const gsl_rng * r);
    
    hgt_pop * p = (hgt_pop *) malloc (sizeof(hgt_pop));
    
    p->size = params->size;
    p->target_size = params->size;
    p->seq_len = params->seq_len;
    p->generation = 0;
    p->linkage_size = params->linkage_size;
    
    p->genomes = (char **) malloc(p->size * sizeof(char *));
    p->fitness = (double *) malloc(p->size * sizeof(double));
    p->linkages = (hgt_pop_linkage **) malloc(p->size * sizeof(hgt_pop_linkage*));
    p->locus_linkages = (hgt_pop_linkage ***) malloc(p->linkage_size * sizeof(hgt_pop_linkage**));
    for (i = 0; i < p->linkage_size; i++) {
        p->locus_linkages[i] = (hgt_pop_linkage**) malloc(p->size * sizeof(hgt_pop_linkage*));
        for (j = 0; j < p->size; j++) {
            p->locus_linkages[i][j] = hgt_pop_linkage_new(NULL, p->generation);
        }
    }
    
    // random initilize genomes
    ancestor = malloc((p->seq_len+1) * sizeof(char));
    random_seq(ancestor, p->seq_len, r);
    for (i = 0; i < p->size; i ++) {
        p->genomes[i] = malloc((p->seq_len+1) * sizeof(char));
        strcpy(p->genomes[i], ancestor);
        p->fitness[i] = 0;
        p->linkages[i] = hgt_pop_linkage_new(NULL, p->generation);
    }
    
    // define transfer hotspots
    if (params->tr_hotspot_num > 0) {
        params->tr_hotspots = malloc(params->tr_hotspot_num * sizeof(unsigned long*));
        unsigned long init_start = gsl_rng_uniform_int(r, params->seq_len);
        unsigned long expected_space = (params->seq_len - params->tr_hotspot_num*params->tr_hotspot_length) / (params->tr_hotspot_num);
        for (i = 0; i < params->tr_hotspot_num; i++) {
            params->tr_hotspots[i] = malloc(2*sizeof(unsigned long));
            unsigned long start = init_start;
            unsigned long end = start + params->tr_hotspot_length;
            unsigned long space = expected_space;
            init_start = end + space;
            params->tr_hotspots[i][0] = start;
            params->tr_hotspots[i][1] = end;
        }
    }
    
    // define mutation hotspots
    if (params->mu_hotspot_num > 0) {
        params->mu_hotspots = malloc(params->mu_hotspot_num * sizeof(unsigned long*));
        unsigned long init_start = gsl_rng_uniform_int(r, params->seq_len);
        unsigned long expected_space = (params->seq_len - params->mu_hotspot_num*params->mu_hotspot_length) / (params->mu_hotspot_num);
        for (i = 0; i < params->mu_hotspot_num; i++) {
            params->mu_hotspots[i] = malloc(2*sizeof(unsigned long));
            unsigned long start = init_start;
            unsigned long end = start + params->mu_hotspot_length;
            unsigned long space = expected_space;
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
    int i;
    hgt_pop * new_p = (hgt_pop *) calloc(1, sizeof(hgt_pop));
    new_p->size = p->size;
    new_p->seq_len = p->seq_len;
    new_p->genomes = (char **) malloc(p->size * sizeof(char *));
    new_p->fitness = (double *) malloc(p->size * sizeof(double));
    for (i = 0; i < p->size; i++) {
        new_p->genomes[i] = malloc((p->seq_len+1) * sizeof(char));
        strcpy(new_p->genomes[i], p->genomes[i]);
        new_p->fitness[i] = p->fitness[i];
    }
    
    // initilize caches
    new_p->survived = (unsigned long *) calloc(p->size, sizeof(unsigned long));
    new_p->new_born = (unsigned long *) calloc(p->size, sizeof(unsigned long));
    
    return new_p;
}

int hgt_pop_free(hgt_pop * p) {
    int i, j;
    for (i = 0; i < p->size; i ++) {
        free(p->genomes[i]);
        hgt_pop_linkage_free(p->linkages[i]);
    }

    for (i = 0; i < p->linkage_size; i++) {
        for (j = 0; j < p->size; j++) {
            hgt_pop_linkage_free(p->locus_linkages[i][j]);
        }
        free(p->locus_linkages[i]);
    }

    free(p->genomes);
    free(p->fitness);
    free(p->linkages);
    free(p->locus_linkages);
    if (p->cache_allocated != 0) {
        free(p->survived);
        free(p->new_born);
    }
    free(p);
    
    return EXIT_SUCCESS;
    
}

char *hgt_pop_to_json(hgt_pop *p, hgt_pop_params *params){
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
    int j;
    for (j = 0; j < p->size; j ++) {
        if (j < p->size - 1) {
            bformata(b, "\"%s\",\n", p->genomes[j]);
        } else {
            bformata(b, "\"%s\"\n", p->genomes[j]);
        }
    }
    bformata(b, "]\n}\n");
    c = bstr2cstr(b, '\n');
    bdestroy(b);
    
    return c;
}

// define handler for parsing configure file
static int hgt_pop_params_handler(void *params, const char* section, const char* name, const char* value) {
    hgt_pop_params* pconfig = (hgt_pop_params*) params;
    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
    if (MATCH("population", "size")) {
        pconfig->size = atoi(value);
    } else if (MATCH("population", "length")) {
        pconfig->seq_len = atoi(value);
    } else if (MATCH("population", "generations")) {
        pconfig->generations = atoi(value);
    } else if (MATCH("mutation", "rate")) {
        pconfig->mu_rate = atof(value);
    } else if (MATCH("mutation", "hotspot_number")){
        pconfig->mu_hotspot_num = atoi(value);
    } else if (MATCH("mutation", "hotspot_length")) {
        pconfig->mu_hotspot_length = atoi(value);
    } else if (MATCH("mutation", "hotspot_ratio")) {
        pconfig->mu_hotspot_ratio = atof(value);
    } else if (MATCH("transfer", "rate")) {
        pconfig->tr_rate = atof(value);
    } else if (MATCH("transfer", "fragment")) {
        pconfig->frag_len = atoi(value);
    } else if (MATCH("sample", "size")) {
        pconfig->sample_size = atoi(value);
    } else if (MATCH("sample", "time")) {
        pconfig->sample_time = atoi(value);
    } else if (MATCH("sample", "replicates")) {
        pconfig->replicates = atoi(value);
    } else if (MATCH("cov", "maxl")) {
        pconfig->maxl = atoi(value);
    } else if (MATCH("output", "prefix")) {
        pconfig->prefix = strdup(value);
    } else if (MATCH("transfer", "hotspot_number")) {
        pconfig->tr_hotspot_num = atoi(value);
    } else if (MATCH("transfer", "hotspot_length")) {
        pconfig->tr_hotspot_length = atoi(value);
    } else if (MATCH("transfer", "hotspot_ratio")) {
        pconfig->tr_hotspot_ratio = atoi(value);
    } else if (MATCH("fitness", "rate")) {
        pconfig->b_mu_rate = atof(value);
    } else if (MATCH("fitness", "scale")) {
        pconfig->fitness_scale = atof(value);
    } else if (MATCH("fitness", "shape")) {
        pconfig->fitness_shape = atof(value);
    } else {
        return 0;
    }
    return 1;
}

// use inih library to parse configure file, returning hgt_pop_params
int hgt_pop_params_parse_from_ini(hgt_pop_params *params, char *filename) {
    if (ini_parse(filename, hgt_pop_params_handler, params) < 0) {
        printf("Can't load %s\n", filename);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int hgt_pop_params_parse(hgt_pop_params *params, int argc, char **argv, char * progname){
    struct arg_int *size = arg_int1("n", "size", "<unsigned long>", "population size");
    struct arg_int *seq_len = arg_int1("l", "len", "<unsigned long>", "genome length");
    struct arg_int *frag_len = arg_int0("f", "frag", "<unsigned long>", "fragment length");
    struct arg_dbl *mu_rate = arg_dbl1("u", "mu_rate", "<double>", "mutation rate");
    struct arg_dbl *tr_rate = arg_dbl0("t", "tr_rate", "<double>", "transfer rate");
    struct arg_int *gen = arg_int1("g", "gen", "<unsigned long>", "generations");
    struct arg_dbl *b_mu_rate = arg_dbl0("z", "b_mu_rate", "<double>", "beneficial mutation rate");
    struct arg_dbl *fitness_scale = arg_dbl0("x", "fitness_scale", "<double>", "selection efficient scale");
    struct arg_dbl *fitness_shape = arg_dbl0("y", "fitness_shape", "<double>", "selection efficient shape");
    
    struct arg_int *spl_size = arg_int0("s", "sample_size", "<unsigned long>", "sample size");
    struct arg_int *spl_time = arg_int0("i", "sample_time", "<unsigned long>", "sample time");
    struct arg_int *repl = arg_int0("r", "replication", "<unsigned long>", "replication");
    struct arg_int *maxl = arg_int0("m", "maxl", "<unsigned long>", "maxl");
    
    struct arg_file *prefix = arg_file1("o", "output", "<output>", "prefix");
    
    struct arg_file *config = arg_file0("c", "config", "<output>", "configure file");
    struct arg_int *reproduction = arg_int0("a", "reproduction_model", "<unsigned int>", "reproduction model");
    struct arg_int *frag_type = arg_int0("b", "frag_type", "<unsigned int>", "fragment type");

    struct arg_lit  *help    = arg_lit0(NULL,"help", "print this help and exit");
    struct arg_end  *end     = arg_end(22);

    void* argtable[] = {size, seq_len, frag_len, mu_rate, tr_rate, gen, b_mu_rate, fitness_scale, fitness_shape, spl_time,
        spl_size, repl, maxl, prefix, config, reproduction, frag_type, help, end};
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
    
    /* check if a configure file is supplied, use config ini parser */
    if (config->count > 0) {
        exit_code = hgt_pop_params_parse_from_ini(params, (char *) config->filename[0]);
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
    params->b_mu_rate = b_mu_rate->dval[0];
    params->fitness_shape = fitness_shape->dval[0];
    params->fitness_scale = fitness_scale->dval[0];
    params->reprodution = reproduction->ival[0];
    params->frag_type = frag_type->ival[0];
    
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
    fprintf(stream, "transfer hotspot number = %d\n", params->tr_hotspot_num);
    fprintf(stream, "transfer hotspot length = %lu\n", params->tr_hotspot_length);
    fprintf(stream, "transfer hotspot ratio = %g\n", params->tr_hotspot_ratio);
    fprintf(stream, "mutation hotspot number = %d\n", params->mu_hotspot_num);
    fprintf(stream, "mutation hotspot length = %lu\n", params->mu_hotspot_length);
    fprintf(stream, "mutation hotspot ratio = %g\n", params->mu_hotspot_ratio);
    fprintf(stream, "fitness mutation rate = %g\n", params->b_mu_rate);
    fprintf(stream, "fitness scale = %g\n", params->fitness_scale);
    fprintf(stream, "fitness shape = %g\n", params->fitness_shape);
    fprintf(stream, "reproduction model = %d\n", params->reprodution);
    fprintf(stream, "frag type = %d\n", params->frag_type);
    return EXIT_SUCCESS;
}

hgt_pop_params *hgt_pop_params_alloc() {
    hgt_pop_params *params = calloc(1, sizeof(hgt_pop_params));
    return params;
}

int hgt_pop_mutate_at(hgt_pop *p,
                      unsigned long g,
                      unsigned long s,
                      const gsl_rng *r) {
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

int hgt_pop_mutate_step(hgt_pop *p, hgt_pop_params *params, const gsl_rng *r) {
    unsigned g;
    // randomly choose a cell.
    g = gsl_rng_uniform_int(r, p->size);
    p->fitness[g] += params->fitness_scale;

    return EXIT_SUCCESS;
}

int hgt_pop_mutate(hgt_pop *p, hgt_pop_params* params, const gsl_rng* r) {
    int hgt_pop_mutate_at(hgt_pop *p,
                          unsigned long g,
                          unsigned long s,
                          const gsl_rng *r);
    unsigned long g, s, hs_len, pos;
    double ratio;
    // randomly choose a genome for mutation
    g = gsl_rng_uniform_int(r, p->size);
    
    if (params -> mu_hotspot_num > 0) { // we need to deal hotspots
        // first calculate total hotspot length, and its ratio to sequence length
        hs_len = params->mu_hotspot_num * params->mu_hotspot_length;
        ratio = (double) hs_len / (double) p->seq_len;
        
        // we need to decise a mutation is happened inside hotspot or outside
        // we flip a coin
        double v = gsl_rng_uniform(r);
        if (v < (1-ratio)/(1+(params->mu_hotspot_ratio - 1)*ratio)) { // outside hotspot
            // we randomly choose a position from sites outside hotspots
            pos = gsl_rng_uniform_int(r, p->seq_len - hs_len);
            // and then search the absolute position on the genome
            s = search_region(pos, params->mu_hotspots, params->mu_hotspot_num, 0);
        } else {
            // similarly, we randomly choose a position from hotspots
            pos = gsl_rng_uniform_int(r, hs_len);
            // and then search the absolute position on the genome
            s = search_region(pos, params->mu_hotspots, params->mu_hotspot_num, 1) % p->seq_len;
        }
    } else { // otherwise, we just randomly choose a site from the genome.
        s = gsl_rng_uniform_int(r, p->seq_len);
    }
    
    // do mutation given a genome and a site
    hgt_pop_mutate_at(p, g, s, r);
    
    return EXIT_SUCCESS;
}

int hgt_pop_transfer(hgt_pop *p, hgt_pop_params* params, unsigned long frag_len, const gsl_rng* r) {
    int hgt_pop_transfer_at(hgt_pop *p,
                            unsigned long donor,
                            unsigned long receiver,
                            unsigned long frag_len,
                            unsigned long start);
    unsigned long s, donor, receiver, hs_len, pos;
    double ratio;
    // randomly choose a reciver
    receiver = gsl_rng_uniform_int(r, p->size);
    // randomly choose a donor.
    donor = gsl_rng_uniform_int(r, p->size);

    if (donor != receiver) {
        if (params -> tr_hotspot_num > 0) { // we need to deal hotspots
            // first calculate total hotspot length, and its ratio to sequence length
            hs_len = params->tr_hotspot_num * params->tr_hotspot_length;
            ratio = (double) hs_len / (double) p->seq_len;
            
            // we need to decise a mutation is happened inside hotspot or outside
            // we flip a coin
            double v = gsl_rng_uniform(r);
            if (v < (1-ratio)/(1+(params->tr_hotspot_ratio - 1)*ratio)) { // outside hotspot
                // we randomly choose a position from sites outside hotspots
                pos = gsl_rng_uniform_int(r, p->seq_len - hs_len);
                // and then search the absolute position on the genome
                s = search_region(pos, params->tr_hotspots, params->tr_hotspot_num, 0);
            } else {
                // similarly, we randomly choose a position from hotspots
                pos = gsl_rng_uniform_int(r, hs_len);
                // and then search the absolute position on the genome
                s = search_region(pos, params->tr_hotspots, params->tr_hotspot_num, 1) % p->seq_len;
            }
        } else { // otherwise, we just randomly choose a site from the genome.
            s = gsl_rng_uniform_int(r, p->seq_len);
        }
        
        // do transfer given a genome and a site
        hgt_pop_transfer_at(p, donor, receiver, frag_len, s);
    }
    
    return EXIT_SUCCESS;
}

int hgt_pop_transfer_fitness(hgt_pop *p, hgt_pop_params* params, unsigned long frag_len, const gsl_rng* r) {
    int hgt_pop_transfer_at(hgt_pop *p,
                            unsigned long donor,
                            unsigned long receiver,
                            unsigned long frag_len,
                            unsigned long start);
    unsigned long s, donor, receiver, hs_len, pos;
    double ratio, * fitness;
    // randomly choose a reciver
    receiver = gsl_rng_uniform_int(r, p->size);
    // randomly choose a donor according to fitnesses.
    fitness = (double *) malloc(p->size * sizeof(double));
    hgt_pop_calc_fitness(p, fitness);
    donor = hgt_utils_Roulette_Wheel_select(fitness, p->size, r);
    free(fitness);

    if (donor != receiver) {
        if (params -> tr_hotspot_num > 0) { // we need to deal hotspots
            // first calculate total hotspot length, and its ratio to sequence length
            hs_len = params->tr_hotspot_num * params->tr_hotspot_length;
            ratio = (double) hs_len / (double) params->seq_len;
            
            // we need to decise a mutation is happened inside hotspot or outside
            // we flip a coin
            double v = gsl_rng_uniform(r);
            if (v < (1-ratio)/(1+(params->tr_hotspot_ratio - 1)*ratio)) { // outside hotspot
                // we randomly choose a position from sites outside hotspots
                pos = gsl_rng_uniform_int(r, params->seq_len - hs_len);
                // and then search the absolute position on the genome
                s = search_region(pos, params->tr_hotspots, params->tr_hotspot_num, 0);
            } else {
                // similarly, we randomly choose a position from hotspots
                pos = gsl_rng_uniform_int(r, hs_len);
                // and then search the absolute position on the genome
                s = search_region(pos, params->tr_hotspots, params->tr_hotspot_num, 1) % params->seq_len;
            }
        } else { // otherwise, we just randomly choose a site from the genome.
            s = gsl_rng_uniform_int(r, params->seq_len);
        }
        
        // do transfer given a genome and a site
        hgt_pop_transfer_at(p, donor, receiver, frag_len, s);
    }

    free(fitness);
    
    return EXIT_SUCCESS;
}

int hgt_pop_transfer_at(hgt_pop *p,
                        unsigned long donor, 
                        unsigned long receiver, 
                        unsigned long frag_len, 
                        unsigned long start)
{
    int linkage_transfer(hgt_pop_linkage ** linkages, int donor, int receiver, unsigned long generation);
    int exit_code;
    exit_code = EXIT_SUCCESS;
    if (start + frag_len < p->seq_len) {
        strncpy(p->genomes[receiver]+start, p->genomes[donor]+start, frag_len);
    } else {
        strcpy(p->genomes[receiver]+start, p->genomes[donor]+start);
        strncpy(p->genomes[receiver], p->genomes[donor], start + frag_len - p->seq_len);
    }

    int i, track_linkage;
    for (i = 0; i < p->linkage_size; i++) {
        track_linkage = 0;
        if (start + frag_len < p->seq_len) {
            if (start == 0 && i < frag_len) {
                track_linkage = 1;
            }
        } else {
            if (i < start + frag_len - p->seq_len) {
                track_linkage = 1;
            }
        }

        if (track_linkage == 1) {
            linkage_transfer(p->locus_linkages[i], donor, receiver, p->generation);
        }
    }

    return exit_code;
}

int linkage_transfer(hgt_pop_linkage ** linkages, int donor, int receiver, unsigned long generation) {
    hgt_pop_linkage * parent;
    if (donor != receiver) {
        parent = linkages[donor]->parent;
        generation = linkages[donor]->birthTime;
        hgt_pop_linkage_free(linkages[receiver]);
        linkages[receiver] = hgt_pop_linkage_new(parent, generation);
    }
    return EXIT_SUCCESS;
}

int hgt_pop_sample_moran(hgt_pop *p, const gsl_rng *r) {
    int linkage_birth_dead(hgt_pop_linkage **linkages, int birth, int dead, unsigned long generation);
    unsigned long b, d;
    double * fitness;
    // increase population generation by 1.
    p->generation++;
    // randomly choose a going-death cell.
    d = gsl_rng_uniform_int(r, p->size);
    // randomly choose a going-birth one according to the fitness.
    // first calculate the fitness for each cell.
    fitness = (double *) malloc(p->size * sizeof(double));
    hgt_pop_calc_fitness(p, fitness);
    // randomly choose one proportional to the fitness.
    b = hgt_utils_Roulette_Wheel_select(fitness, p->size, r);
    // copy the genome and its fitness from b to d.
    if (b != d) {
        strcpy(p->genomes[d], p->genomes[b]);
        p->fitness[d] = p->fitness[b];
    }

    linkage_birth_dead(p->linkages, b, d, p->generation);
    int i;
    for (i = 0; i < p->linkage_size; i++) {
        hgt_pop_linkage ** linkages;
        linkages = p->locus_linkages[i];
        linkage_birth_dead(linkages, b, d, p->generation);
    }

    free(fitness);
    return EXIT_SUCCESS;
}

int hgt_pop_sample_linear_selection(hgt_pop *p, const gsl_rng *r) {
    int update_linkages_wf(hgt_pop *p, int *offsprings);
    char ** current_genomes, ** new_genomes;
    double * current_fitness, * normalized_fitness, * new_fitness;
    int current_size;
    
    p->generation++;
    current_genomes = p->genomes;
    current_fitness = p->fitness;
    current_size = p->size;
    
    normalized_fitness = (double *) malloc(p->size * sizeof(double));
    hgt_pop_calc_fitness(p, normalized_fitness);
    
    int i, num_offspring, new_size;
    int * offsprings;
    offsprings = (int *) malloc(p->size * sizeof(int));
    new_size = 0;
    
    for (i = 0; i < current_size; i++) {
        num_offspring = gsl_ran_poisson(r, normalized_fitness[i]);
        offsprings[i] = num_offspring;
        new_size += num_offspring;
    }
    
    new_genomes = (char **) malloc(new_size * sizeof(char*));
    new_fitness = (double *) malloc(new_size * sizeof(double));
    int j, k;
    char * g;
    k = 0;
    for (i = 0; i < current_size; i++) {
        num_offspring = offsprings[i];
        if (num_offspring > 0) {
            for (j = 0; j < num_offspring; j++) {
                if (j == 0) {
                    g = current_genomes[i];
                } else {
                    g = (char *) malloc(p->seq_len * sizeof(char));
                    strcpy(g, current_genomes[i]);
                }
                new_genomes[k] = g;
                new_fitness[k] = current_fitness[i];
                k++;
            }
        } else {
            free(current_genomes[i]);
        }
    }
    
    update_linkages_wf(p, offsprings);
    p->size = new_size;
    p->genomes = new_genomes;
    p->fitness = new_fitness;
    
    free(current_genomes);
    free(current_fitness);
    free(normalized_fitness);
    free(offsprings);
    return EXIT_SUCCESS;
}

int linkage_birth_dead(hgt_pop_linkage **linkages, int birth, int dead, unsigned long generation) {
    hgt_pop_linkage * parent;
    parent = linkages[birth];
    linkages[birth] = hgt_pop_linkage_new(parent, generation);
    if (birth != dead) {
        hgt_pop_linkage_free(linkages[dead]);
        linkages[dead] = hgt_pop_linkage_new(parent, generation);
    }
    return EXIT_SUCCESS;
}

int update_linkages_wf(hgt_pop *p, int *offsprings) {
    int linkage_linear_selection(hgt_pop_linkage ** current_linkages, int * offsprings, int current_size, hgt_pop_linkage ** new_linkages, unsigned long generation);
    int current_size = 0, new_size = 0;
    current_size = p->size;
    int i = 0;
    for (i = 0; i < current_size; i++) {
        new_size+=offsprings[i];
    }
    
    hgt_pop_linkage ** current_linkages = p->linkages;
    hgt_pop_linkage ** new_linkages = (hgt_pop_linkage **) malloc(new_size * sizeof(hgt_pop_linkage*));
    linkage_linear_selection(current_linkages, offsprings, current_size, new_linkages, p->generation);
    hgt_pop_linkages_free(current_linkages,current_size);
    p->linkages = new_linkages;
    
    for (i = 0; i < p->linkage_size; i++) {
        current_linkages = p->locus_linkages[i];
        new_linkages = (hgt_pop_linkage**) malloc(new_size * sizeof(hgt_pop_linkage*));
        linkage_linear_selection(current_linkages, offsprings, current_size, new_linkages, p->generation);
        p->locus_linkages[i] = new_linkages;
        hgt_pop_linkages_free(current_linkages, current_size);
    }
    
    return EXIT_SUCCESS;
}

int linkage_linear_selection(hgt_pop_linkage ** current_linkages, int * offsprings, int current_size, hgt_pop_linkage ** new_linkages, unsigned long generation){
    int i, num_offspring;
    
    hgt_pop_linkage *l, *parent;
    int j, k;
    k = 0;
    for (i = 0; i < current_size; i++) {
        parent = current_linkages[i];
        num_offspring = offsprings[i];
        if (num_offspring > 0) {
            for (j = 0; j < num_offspring; j++) {
                l = hgt_pop_linkage_new(parent, generation);
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

hgt_pop_linkage * hgt_pop_linkage_new(hgt_pop_linkage * parent, unsigned long birthTime) {
    hgt_pop_linkage * l;
    l = hgt_pop_linkage_alloc();
    l->birthTime = birthTime;
    l->parent = parent;
    if (parent) {
        parent->numChildren++;
    }
    return l;
}

int hgt_pop_sample_wf(hgt_pop *p, const gsl_rng *r) {
    p->generation++;
    int update_linkages_wf(hgt_pop *p, int *offsprings);
    int i, j, k;
    
    int *offsprings = (int*) calloc(p->size, sizeof(int));
    for (i = 0; i < p->size; i++) {
        offsprings[i] = 0;
    }
    
    for (i = 0; i < p->size; i ++) {
        j = gsl_rng_uniform_int(r, p->size);
        offsprings[j]++;
    }
    
    char **new_genomes = (char **) malloc(p->size * sizeof(char*));
    char **current_genomes = p->genomes;
    double *new_fitness = (double *) malloc(p->size * sizeof(double));
    double *current_fitness = p->fitness;
    char *g;
    k = 0;
    int num_offspring;
    for (i = 0; i < p->size; i++) {
        num_offspring = offsprings[i];
        if (num_offspring > 0) {
            for (j = 0; j < num_offspring; j++) {
                if (j == 0) {
                    g = current_genomes[i];
                } else {
                    g = (char *) malloc(p->seq_len * sizeof(char));
                    strcpy(g, current_genomes[i]);
                }
                new_genomes[k] = g;
                new_fitness[k] = current_fitness[i];
                k++;
            }
        } else {
            free(current_genomes[i]);
        }
    }
    
    update_linkages_wf(p, offsprings);
    p->genomes = new_genomes;
    p->fitness = new_fitness;
    free(offsprings);
    free(current_fitness);
    free(current_genomes);

    return EXIT_SUCCESS;

}

double hgt_pop_frag_constant(hgt_pop_params *params, const gsl_rng *r) {
    return (double) params->frag_len;
}

double hgt_pop_frag_exp(hgt_pop_params *params, const gsl_rng *r) {
    return gsl_ran_exponential(r, (double) params->frag_len);
}

int hgt_pop_evolve(hgt_pop *p, hgt_pop_params *params, hgt_pop_sample_func sample_f, hgt_pop_coal_time_func c_time_f, hgt_pop_frag_func frag_f, const gsl_rng *r) {
    double mu, total;
    unsigned long count, k, frag_len;
    int i;
    double weights[3];
    
    sample_f(p, r);

    weights[0] = params->mu_rate * (double) (p->seq_len * p->size);
    weights[1] = params->tr_rate * (double) (p->seq_len * p->size);
    weights[2] = params->b_mu_rate * (double) p->size;
    total = 0;
    for (i = 0; i < 3; i ++) {
        total += weights[i];
    }

    mu = c_time_f(p->size, r) * total;
    count = gsl_ran_poisson(r, mu);
    for (k = 0; k < count; k ++) {
        i = hgt_utils_Roulette_Wheel_select(weights, 3, r);
        if (i == 0) {
            hgt_pop_mutate(p, params, r);
        } else if (i == 1) {
            frag_len = frag_f(params, r);
            hgt_pop_transfer(p, params, params->frag_len, r);
        } else {
            hgt_pop_mutate_step(p, params, r);
        }
    }
    
    if (p->generation % 10 == 0) {
        hgt_pop_linkage_prune_p(p);
    }

    return EXIT_SUCCESS;
}


double hgt_pop_coal_time_moran(unsigned long p_size, const gsl_rng *r) {
    return gsl_ran_exponential(r, 1.0 / (double) p_size);
}

double hgt_pop_coal_time_wf(unsigned long p_size, const gsl_rng *r) {
    return gsl_ran_exponential(r, 1.0);
}

double hgt_pop_coal_time_linear_selection(unsigned long p_size, const gsl_rng *r) {
    return hgt_pop_coal_time_wf(p_size, r);
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
    
    return EXIT_SUCCESS;
}

// Calcualte covariances using all individuals.
int hgt_pop_calc_cov_all(hgt_cov_result *result, hgt_pop *p) {
    // allocate a binary matrix
    unsigned long matrix_size = p->size * (p->size - 1) / 2;
    short  **matrix = malloc( matrix_size * sizeof(short*));
    
    int i, j, k, h;
    
    k = 0;
    for (i = 0; i < p->size; i++) {
        for (j = i + 1; j < p->size; j++) {
            matrix[k] = malloc(p->seq_len*sizeof(short));
            for (h = 0; h < p->seq_len; h++) {
                if (p->genomes[i][h] != p->genomes[j][h]) {
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

double hgt_pop_mean_fitness(hgt_pop *p) {
    double f;
    int i;
    f = 0;
    for (i = 0; i < p->size; ++i){
        f += p->fitness[i];
    }

    f /= (double) p->size;

    return f;
}

int linkage_prune(hgt_pop_linkage *l) {
    hgt_pop_linkage *parent;
    parent = l->parent;
    if (parent) {
        if (parent->numChildren <= 1) {
            l->parent = parent->parent;
            free(parent);
            linkage_prune(l);
        }
    }
    return EXIT_SUCCESS;
}

int hgt_pop_linkages_prune(hgt_pop_linkage ** linkages, int size) {
    int linkage_prune(hgt_pop_linkage *l);
    int i;
    for (i = 0; i < size; i++) {
        linkage_prune(linkages[i]);
    }
    return EXIT_SUCCESS;
}

int hgt_pop_linkage_prune_p(hgt_pop *p) {
    hgt_pop_linkages_prune(p->linkages, p->size);
    int i;
    for (i = 0; i < p->linkage_size; i++) {
        hgt_pop_linkages_prune(p->locus_linkages[i], p->size);
    }
    return EXIT_SUCCESS;
}

int hgt_pop_linkages_free(hgt_pop_linkage ** linkages, int size) {
    int i;
    for (i = 0; i < size; i++) {
        hgt_pop_linkage_free(linkages[i]);
    }
    free(linkages);
    return EXIT_SUCCESS;
}

int hgt_pop_linkage_free(hgt_pop_linkage * l) {
    hgt_pop_linkage * parent;
    parent = l->parent;
    if (l->numChildren < 1) {
        if (parent) {
            // reduce the number of children of the parent by 1,
            // and try to free the parent.
            parent->numChildren--;
            hgt_pop_linkage_free(parent);
        }
        free(l);
    }
    
    
    return EXIT_SUCCESS;
}

hgt_pop_linkage * hgt_pop_linkage_alloc() {
    hgt_pop_linkage * l;
    l = (hgt_pop_linkage *) malloc(sizeof(hgt_pop_linkage));
    l->numChildren = 0;
    l->birthTime = 0;
    l->parent = NULL;
    return l;
}

unsigned long hgt_pop_linkage_find_most_rescent_ancestor(hgt_pop_linkage ** linkages, int size) {
    int i;
    int bad, found, maxIndex;
    unsigned long maxBirthTime, birthTime;
    hgt_pop_linkage * parent;
    // init variables.
    bad = 0;
    found = 1;
    parent = linkages[0]->parent;
    for (i = 0; i < size; ++i) {
        if (!linkages[i]->parent) {
            printf("born at %lu\n", linkages[i]->birthTime);
            bad = 1;
            break;
        } else {
            if (parent != linkages[i]->parent ) {
                found = 0;
            }
        }
    }

    if (bad == 1) {
        return 0;
    } else {
        if (found == 1) {
            return linkages[0]->birthTime;
        } else {
            // find the one that has max birth time.
            maxBirthTime = linkages[0]->birthTime;
            maxIndex = 0;
            for (i = 0; i < size; i++) {
                birthTime = linkages[i]->birthTime;
                if (maxBirthTime < birthTime) {
                    maxBirthTime = birthTime;
                    maxIndex = i;
                }
            }
            linkages[maxIndex] = linkages[maxIndex]->parent;
            return hgt_pop_linkage_find_most_rescent_ancestor(linkages, size);
        }
    }

}

unsigned long hgt_pop_linkage_find_most_rescent_coalescence(hgt_pop_linkage ** linkages, int size) {
    int i, j;
    int bad, maxIndex;
    unsigned long birthTime, maxBirthTime;
    
    bad = 0;
    for (i = 0; i < size; ++i) {
        if (!linkages[i]->parent) {
            bad = 1;
            break;
        }
    }

    if (bad == 1) {
        return 0;
    } else {
        maxBirthTime = 0;
        for (i = 0; i < size; i ++) {
            for (j = i+1; j < size; j++) {
                if (linkages[i]->parent == linkages[j]->parent) {
			birthTime = linkages[i]->birthTime;
                    if (maxBirthTime < birthTime) {
                        maxBirthTime = birthTime;
                    }
                }
            }
        }
        if (maxBirthTime > 0) {
            return maxBirthTime;
        } else {
            // find the one that has max birth time.
            maxBirthTime = linkages[0]->birthTime;
            maxIndex = 0;
            for (i = 0; i < size; i++) {
                if (maxBirthTime < linkages[i]->birthTime) {
                    maxBirthTime = linkages[i]->birthTime;
                    maxIndex = i;
                }
            }
            linkages[maxIndex] = linkages[maxIndex]->parent;
            return hgt_pop_linkage_find_most_rescent_coalescence(linkages, size);
        }
    }
}

int hgt_pop_calc_coal_time(hgt_pop_linkage ** pop_linkages, int size, 
    unsigned long sample_size, 
    unsigned long *res, 
    int linkage_size, 
    hgt_pop_linkage_find_time_func find_func,
    const gsl_rng *r) {

    int i, j, a[linkage_size], b[size];
    for (i = 0; i < size; i++) {
        b[i] = i;
    }

    hgt_pop_linkage * linkages[linkage_size];

    for (i = 0; i < sample_size; i++) {
        gsl_ran_choose(r, a, linkage_size, b, size, sizeof(int));
        for (j = 0; j < linkage_size; j++) {
            linkages[j] = pop_linkages[a[j]];
        }
        res[i] = find_func(linkages, linkage_size);
    }

    return EXIT_SUCCESS;
}

int hgt_pop_calc_most_recent_coal_time(hgt_pop_linkage ** pop_linkages, int size, unsigned long sample_size, unsigned long * res, int linkage_size, const gsl_rng *r) {
    return hgt_pop_calc_coal_time(pop_linkages, size, sample_size, res, linkage_size, hgt_pop_linkage_find_most_rescent_coalescence, r);
}

int hgt_pop_calc_most_recent_ancestor_time(hgt_pop_linkage ** pop_linkages, int size, unsigned long sample_size, unsigned long * res, int linkage_size, const gsl_rng *r) {
    return hgt_pop_calc_coal_time(pop_linkages, size, sample_size, res, linkage_size, hgt_pop_linkage_find_most_rescent_ancestor, r);
}

int hgt_pop_calc_fitness(hgt_pop *p, double * fitness) {
    double f, mean_fitness, cpot, size_ratio;
    size_ratio = (double) p->size / (double) p->target_size;
    mean_fitness = hgt_pop_mean_fitness(p);
    cpot = mean_fitness - (1.0 - size_ratio);
    int i;
    for (i = 0; i < p->size; i++) {
        f = p->fitness[i];
        fitness[i] = exp(f - cpot);
    }

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

unsigned long search_region(unsigned long pos, unsigned long** region, int num_regions, int inside) {
    int i;
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
