//
//  hgt_cov.h
//  hgt
//
//  Created by Mingzhi Lin on 1/20/14.
//
//

#ifndef hgt_hgt_cov_h
#define hgt_hgt_cov_h

typedef int(*hgt_cov_sample_func)(unsigned long *a, unsigned long *b, unsigned long *c, unsigned long *d, unsigned long p_size, const gsl_rng *r);
int hgt_cov_sample_p2(unsigned long *a, unsigned long *b, unsigned long *c, unsigned long *d, unsigned long p_size, const gsl_rng *r);
int hgt_cov_sample_p3(unsigned long *a, unsigned long *b, unsigned long *c, unsigned long *d, unsigned long p_size, const gsl_rng *r);
int hgt_cov_sample_p4(unsigned long *a, unsigned long *b, unsigned long *c, unsigned long *d, unsigned long p_size, const gsl_rng *r);
#endif
