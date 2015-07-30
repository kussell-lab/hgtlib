#include "hgt_params.h"
#include "ini.h"
#include "argtable3.h"
#include <string.h>
#include <stdlib.h>

int hgt_params_check_default(hgt_params *params);

// define handler for parsing configure file
static int hgt_params_handler(void *params, const char* section, const char* name, const char* value) {
    hgt_params *params1 = (hgt_params*) params;
    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
    if (MATCH("population", "size")) {
        params1->size = (unsigned int) atoi(value);
    } else if (MATCH("population", "length")) {
        params1->seq_len = (unsigned int) atoi(value);
	} else if (MATCH("genome", "alphabet_size")) {
		params1->alphabet_size = (unsigned int)atoi(value);
	} else if (MATCH("population", "generations")) {
        params1->generations = (unsigned int) atoi(value);
    } else if (MATCH("population", "model")){
        params1->reprodution = (unsigned int ) atoi(value);
    } else if (MATCH("mutation", "rate")) {
        params1->mu_rate = atof(value);
    } else if (MATCH("mutation", "hotspot_number")){
        params1->mu_hotspot_num = (unsigned int) atoi(value);
    } else if (MATCH("mutation", "hotspot_length")) {
        params1->mu_hotspot_length = (unsigned int) atoi(value);
    } else if (MATCH("mutation", "hotspot_ratio")) {
        params1->mu_hotspot_ratio = atof(value);
    } else if (MATCH("transfer", "rate")) {
        params1->tr_rate = atof(value);
    } else if (MATCH("transfer", "fragment")) {
        params1->frag_len = (unsigned int) atoi(value);
    } else if (MATCH("transfer", "distribution")) {
        params1->frag_type = (unsigned int) atoi(value);
    } else if (MATCH("sample", "size")) {
        params1->sample_size = (unsigned int) atoi(value);
    } else if (MATCH("sample", "time")) {
        params1->sample_time = (unsigned int) atoi(value);
	} else if (MATCH("sample", "generations")) {
		params1->sample_generations = (unsigned int)atoi(value);
	} else if (MATCH("sample", "replicates")) {
        params1->replicates = (unsigned int) atoi(value);
    } else if (MATCH("linkage", "size")) {
        params1->linkage_size = (unsigned int) atoi(value);
    } else if (MATCH("cov", "maxl")) {
        params1->maxl = (unsigned int) atoi(value);
    } else if (MATCH("output", "prefix")) {
        params1->prefix = strdup(value);
    } else if (MATCH("transfer", "hotspot_number")) {
        params1->tr_hotspot_num = (unsigned int) atoi(value);
    } else if (MATCH("transfer", "hotspot_length")) {
        params1->tr_hotspot_length = (unsigned int) atoi(value);
    } else if (MATCH("transfer", "hotspot_ratio")) {
        params1->tr_hotspot_ratio = atoi(value);
    } else if (MATCH("fitness", "rate")) {
        params1->fitness_mutation_rate = atof(value);
    } else if (MATCH("fitness", "scale")) {
        params1->fitness_scale = atof(value);
    } else if (MATCH("fitness", "shape")) {
        params1->fitness_shape = atof(value);
    } else if (MATCH("fitness", "coupled")) {
        params1->fitness_coupled = (unsigned int) atoi(value);
    } else {
        return 0;
    }
    return 1;
}

// use inih library to parse configure file, returning hgt_params
int hgt_params_parse_from_ini(hgt_params *params, const char *filename) {
    if (ini_parse(filename, hgt_params_handler, params) < 0) {
        printf("Can't load %s\n", filename);
        return EXIT_FAILURE;
    }
    hgt_params_check_default(params);
    return EXIT_SUCCESS;
}

int hgt_params_parse(hgt_params *params, int argc, char **argv, char * progname){
    struct arg_int *size = arg_int0("n", "size", "<unsigned long>", "population size");
    struct arg_int *seq_len = arg_int0("l", "len", "<unsigned long>", "genome length");
	struct arg_int *alphabet_size = arg_int0("A", "alphabet_size", "<unsigned int>", "alphabet sizse");
    struct arg_int *frag_len = arg_int0("f", "frag", "<unsigned long>", "fragment length");
    struct arg_dbl *mu_rate = arg_dbl0("u", "mu_rate", "<double>", "mutation rate");
    struct arg_dbl *tr_rate = arg_dbl0("t", "tr_rate", "<double>", "transfer rate");
    struct arg_int *gen = arg_int0("g", "gen", "<unsigned long>", "generations");
    struct arg_dbl *fitness_mutation_rate = arg_dbl0("z", "b_mu_rate", "<double>", "beneficial mutation rate");
    struct arg_dbl *fitness_scale = arg_dbl0("x", "fitness_scale", "<double>", "selection efficient scale");
    struct arg_dbl *fitness_shape = arg_dbl0("y", "fitness_shape", "<double>", "selection efficient shape");
    struct arg_int *fitness_coupled = arg_int0("c", "fitness_coupled", "<int>", "fitness coupled");
    
    struct arg_int *spl_size = arg_int0("s", "sample_size", "<unsigned int>", "sample size");
    struct arg_int *spl_time = arg_int0("i", "sample_time", "<unsigned int>", "sample time");
	struct arg_int *spl_gen = arg_int0("G", "sample_generations", "<unsigned int>", "sample generations");
    struct arg_int *repl = arg_int0("r", "replication", "<unsigned int>", "replication");
    struct arg_int *maxl = arg_int0("m", "maxl", "<unsigned long>", "maxl");
    struct arg_int *linkage_size = arg_int0("k", "linkage_size", "<unsigned int>", "locus linkage size for tracking");
    
    struct arg_file *prefix = arg_file0("o", "output", "<output>", "prefix");
    
    struct arg_file *config = arg_file0("C", "config", "<output>", "configure file");
    struct arg_int *reproduction = arg_int0("a", "reproduction_model", "<unsigned int>", "reproduction model");
    struct arg_int *frag_type = arg_int0("b", "frag_type", "<unsigned int>", "fragment type");
    
    struct arg_lit  *help    = arg_lit0(NULL,"help", "print this help and exit");
    struct arg_end  *end     = arg_end(26);
    
    void* argtable[] = {size, seq_len, alphabet_size, frag_len, mu_rate, tr_rate, gen, fitness_mutation_rate, fitness_scale, fitness_shape, fitness_coupled, spl_time,
        spl_size, spl_gen, repl, maxl, linkage_size, prefix, config, reproduction, frag_type, help, end};
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
        const char * filename = config->filename[0];
        printf("Using config file: %s\n", filename);
        exit_code = hgt_params_parse_from_ini(params, filename);
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
    params->fitness_mutation_rate = fitness_mutation_rate->dval[0];
    params->fitness_shape = fitness_shape->dval[0];
    params->fitness_scale = fitness_scale->dval[0];
    params->fitness_coupled = fitness_coupled->ival[0];
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

	if (spl_gen->count > 0)
	{
		params->sample_generations = spl_gen->ival[0];
	}
	else {
		params->sample_generations = params->generations;
	}
    
    if (repl->count > 0) {
        params->replicates = repl->ival[0];
    } else {
        params->replicates = 1;
    }

	if (alphabet_size->count > 0) {
		params->alphabet_size = alphabet_size->ival[0];
	} else {
		params->alphabet_size = 4;
	}

    hgt_params_check_default(params);
    
exit:
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
    return exit_code;
}

int hgt_params_free(hgt_params *params){
    free(params);
    
    return EXIT_SUCCESS;
}

int hgt_params_check_default(hgt_params *params) {
    // check fitness size.
    if (params->fitness_coupled == 1) {
        params->fitness_size = params->seq_len;
    } else {
        params->fitness_size = 1;
    }

	// check sample generations.
	if (params->sample_generations <= 0) {
		params->sample_generations = params->generations;
	}

    return EXIT_SUCCESS;
}

int hgt_params_printf(hgt_params *params, FILE *stream) {
    fprintf(stream, "population size = %u\n", params->size);
    fprintf(stream, "genome length = %u\n", params->seq_len);
	fprintf(stream, "genome alphabet size = %u\n", params->alphabet_size);
    fprintf(stream, "mutation rate = %g\n", params->mu_rate);
    fprintf(stream, "transfer rate = %g\n", params->tr_rate);
    fprintf(stream, "transfer frag = %u\n", params->frag_len);
    fprintf(stream, "generations = %u\n", params->generations);
    fprintf(stream, "sample size = %u\n", params->sample_size);
    fprintf(stream, "sample time = %u\n", params->sample_time);
	fprintf(stream, "sample generations = %u\n", params->sample_generations);
    fprintf(stream, "replicates = %u\n", params->replicates);
    fprintf(stream, "prefix = %s\n", params->prefix);
    fprintf(stream, "maxl = %u\n", params->maxl);
    fprintf(stream, "transfer hotspot number = %d\n", params->tr_hotspot_num);
    fprintf(stream, "transfer hotspot length = %u\n", params->tr_hotspot_length);
    fprintf(stream, "transfer hotspot ratio = %g\n", params->tr_hotspot_ratio);
    fprintf(stream, "mutation hotspot number = %d\n", params->mu_hotspot_num);
    fprintf(stream, "mutation hotspot length = %u\n", params->mu_hotspot_length);
    fprintf(stream, "mutation hotspot ratio = %g\n", params->mu_hotspot_ratio);
    fprintf(stream, "fitness mutation rate = %g\n", params->fitness_mutation_rate);
    fprintf(stream, "fitness scale = %g\n", params->fitness_scale);
    fprintf(stream, "fitness shape = %g\n", params->fitness_shape);
    fprintf(stream, "fitness coupled = %d\n", params->fitness_coupled);
    fprintf(stream, "fitness size = %u\n", params->fitness_size);
    fprintf(stream, "reproduction model = %d\n", params->reprodution);
    fprintf(stream, "frag type = %d\n", params->frag_type);
    return EXIT_SUCCESS;
}

hgt_params *hgt_params_alloc() {
    hgt_params *params = calloc(1, sizeof(hgt_params));
    return params;
}
