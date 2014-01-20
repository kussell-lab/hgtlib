//
//  hgt_cov.c
//  hgt
//
//  Created by Mingzhi Lin on 1/20/14.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include "hgt_cov.h"

int hgt_cov_sample_p2(unsigned long *a, unsigned long *b, unsigned long *c, unsigned long *d, unsigned long p_size, const gsl_rng *r) {
    *a = gsl_rng_uniform_int(r, p_size);
    *b = gsl_rng_uniform_int(r, p_size);
    while (*a == *b) {
        *b = gsl_rng_uniform_int(r, p_size);
    }
    *c = *a;
    *d = *b;
    return EXIT_SUCCESS;
}

int hgt_cov_sample_p3(unsigned long *a, unsigned long *b, unsigned long *c, unsigned long *d, unsigned long p_size, const gsl_rng *r) {
    *a = gsl_rng_uniform_int(r, p_size);
    *b = gsl_rng_uniform_int(r, p_size);
    while (*a == *b) {
        *b = gsl_rng_uniform_int(r, p_size);
    }
    *c = gsl_rng_uniform_int(r, p_size);
    while (*a == *c || *b == *c) {
        *c = gsl_rng_uniform_int(r, p_size);
    }
    *d = *a;
    return EXIT_SUCCESS;
}

int hgt_cov_sample_p4(unsigned long *a, unsigned long *b, unsigned long *c, unsigned long *d, unsigned long p_size, const gsl_rng *r) {
    *a = gsl_rng_uniform_int(r, p_size);
    *b = gsl_rng_uniform_int(r, p_size);
    while (*a == *b) {
        *b = gsl_rng_uniform_int(r, p_size);
    }
    *c = gsl_rng_uniform_int(r, p_size);
    while (*a == *c || *b == *c) {
        *c = gsl_rng_uniform_int(r, p_size);
    }
    *d = gsl_rng_uniform_int(r, p_size);
    while (*a == *d || *b == *d || *c == *d) {
        *d = gsl_rng_uniform_int(r, p_size);
    }
    return EXIT_SUCCESS;
}
