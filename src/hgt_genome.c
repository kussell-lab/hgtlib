//
//  hgt_genome.c
//  hgt
//
//  Created by LinMingzhi on 7/12/15.
//
//
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "hgt_genome.h"

static int NUM_DNA_CHAR = 4;

int transfer_mutate(hgt_genome *receiver, hgt_genome *donor, int start, int end);
int hgt_genome_mutate_(hgt_genome *g, unsigned pos, unsigned char_max, const gsl_rng *r);

void hgt_genome_set_alphabet_size(unsigned size) {
	NUM_DNA_CHAR = size;
}

int hgt_genome_get_alphabet_size() {
	return NUM_DNA_CHAR;
}

hgt_genome * hgt_genome_alloc(unsigned int seq_len, unsigned int fitness_size) {
    hgt_genome *g = (hgt_genome *) malloc(sizeof(hgt_genome)) ;
    g->seq_len = seq_len;
    g->fitness_size = fitness_size;
    g->seq = (char *) malloc((seq_len+1) * sizeof(char));
    g->fitness = (double *) malloc(fitness_size * sizeof(double));
    g->fitness_score = 0;
    return g;
}

int hgt_genome_free(hgt_genome *g) {
    free(g->seq);
    free(g->fitness);
    free(g);
    return EXIT_SUCCESS;
}

int hgt_genome_copy(hgt_genome *des, hgt_genome *src) {
    memcpy(des->fitness, src->fitness, sizeof(double) * src->fitness_size);
    strcpy(des->seq, src->seq);
    des->fitness_score = src->fitness_score;
    return EXIT_SUCCESS;
}

hgt_genome * hgt_genome_new(char * seq, unsigned int seq_size, unsigned int fitness_size) {
    hgt_genome * g = hgt_genome_alloc(seq_size, fitness_size);
    strcpy(g->seq, seq);
    int i;
    for (i = 0; i < fitness_size; i++) {
        g->fitness[i] = 0;
    }
    
    return g;
}

int hgt_genome_mutate(hgt_genome *g, unsigned int pos, const gsl_rng *r) {
    return hgt_genome_mutate_(g, pos, NUM_DNA_CHAR, r);
}

int hgt_genome_mutate_(hgt_genome *g, unsigned int pos, unsigned char_max, const gsl_rng *r) {
	char random_c = (char) gsl_rng_uniform_int(r, char_max);
	while (random_c == g->seq[pos] )
	{
		random_c = (char)gsl_rng_uniform_int(r, char_max);
	}
	g->seq[pos] = random_c;
	return EXIT_SUCCESS;
}

int hgt_genome_fitness_mutate(hgt_genome *g, unsigned int pos, double delta) {
    if (pos < g->fitness_size) {
        g->fitness[pos] += delta;
        g->fitness_score += delta;
    }
    return EXIT_SUCCESS;
}

int hgt_genome_transfer(hgt_genome *receiver, hgt_genome *donor, unsigned int start, unsigned int frag_len) {
    int end;
    end = start + frag_len;
    if (end < receiver->seq_len) {
        strncpy(receiver->seq+start, donor->seq+start, frag_len);
    } else {
        strcpy(receiver->seq+start, donor->seq+start);
        strncpy(receiver->seq, donor->seq, end - receiver->seq_len);
    }
    return EXIT_SUCCESS;
}

int hgt_genome_fitness_transfer(hgt_genome *receiver, hgt_genome *donor, unsigned int start, unsigned int frag_len) {
    unsigned int end;

    end = start + frag_len;
    if (end < receiver->fitness_size) {
        transfer_mutate(receiver, donor, start, end);
    } else {
        transfer_mutate(receiver, donor, start, receiver->fitness_size);
        transfer_mutate(receiver, donor, 0, end - receiver->fitness_size);
    }
    return EXIT_SUCCESS;
}

int transfer_mutate(hgt_genome *receiver, hgt_genome *donor, int start, int end) {
    int i;
    double delta;
    for (i = start;i < end; i++) {
        delta = donor->fitness[i] - receiver->fitness[i];
        hgt_genome_fitness_mutate(receiver, i, delta);
    }

    return EXIT_SUCCESS;
}

double hgt_genome_get_fitness(hgt_genome *g) {
    return g->fitness_score;
}

char * hgt_genome_get_seq(hgt_genome *g) {
    char * seq;
    seq = g->seq;
    return seq;
}

unsigned int hgt_genome_get_seq_size(hgt_genome *g) {
    unsigned int size;
    size = g->seq_len;
    return size;
}

unsigned int hgt_genome_get_fitness_size(hgt_genome *g) {
    unsigned int size;
    size = g->fitness_size;
    return size;
}

char * hgt_genome_random_sequence(int length, const gsl_rng *r) {
    char * seq = (char *) malloc((length+1) * sizeof(char));
    int i;
    for (i = 0; i < length; i++) {
		seq[i] = (char)gsl_rng_uniform_int(r, NUM_DNA_CHAR);
    }
    seq[length] = '\0';
    return seq;
}