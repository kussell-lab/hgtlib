//
//  hgt_genome.c
//  hgt
//
//  Created by LinMingzhi on 7/12/15.
//
//
#include <stdlib.h>
#include <string.h>
#include "hgt_genome.h"
#include <gsl/gsl_rng.h>

const char DNA[5] = "ATGC\0";
const int NUM_DNA_CHAR = 4;

int transfer_mutate(hgt_genome *receiver, hgt_genome *donor, int start, int end);

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
    char c;
    c = DNA[gsl_rng_uniform_int(r, NUM_DNA_CHAR)];
    while ( c == g->seq[pos] ) {
        c = DNA[gsl_rng_uniform_int(r, NUM_DNA_CHAR)];
    }
    g->seq[pos] = c;
    
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
    unsigned int end, i;

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