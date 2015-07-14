//
//  hgt_pop_genome.h
//  hgt
//
//  Created by LinMingzhi on 7/12/15.
//
//

#ifndef __hgt__hgt_genome__
#define __hgt__hgt_genome__
#include <gsl/gsl_rng.h>
extern const char DNA[5];
extern const int NUM_DNA_CHAR;

typedef struct hgt_genome hgt_genome;
struct hgt_genome {
    char * seq;
    double * fitness;
    unsigned int seq_len;
    unsigned int fitness_size;

    double fitness_score;
};

hgt_genome * hgt_genome_alloc(unsigned int seq_len, unsigned int fitness_size);
hgt_genome * hgt_genome_new(char * seq, unsigned int seq_len, unsigned int fitness_size);
int hgt_genome_copy(hgt_genome *g1, hgt_genome *g2);
int hgt_genome_free(hgt_genome *g);
int hgt_genome_mutate(hgt_genome *g, unsigned int pos, const gsl_rng *r);
int hgt_genome_fitness_mutate(hgt_genome *g, unsigned int pos, double delta);
int hgt_genome_transfer(hgt_genome *receiver, hgt_genome *donor, unsigned int start, unsigned int frag_len);
int hgt_genome_fitness_transfer(hgt_genome *receiver, hgt_genome *donor, unsigned int start, unsigned int frag_len);
double hgt_genome_get_fitness(hgt_genome *g);
char * hgt_genome_get_seq(hgt_genome *g);
unsigned int hgt_genome_get_seq_size(hgt_genome *g);
unsigned int hgt_genome_get_fitness_size(hgt_genome *g);
#endif /* defined(__hgt__hgt_genome__) */
