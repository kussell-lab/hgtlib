//
//  hgt_stat.h
//  hgt
//
//  Created by Mingzhi Lin on 1/24/14.
//
//

#ifndef hgt_hgt_stat_h
#define hgt_hgt_stat_h

typedef struct {
    double m1;
    double dev;
    double n_dev;
    unsigned long n;
} hgt_stat_first_moment;

typedef struct {
    hgt_stat_first_moment * m1;
    double m2;
} hgt_stat_second_moment;

typedef struct {
    hgt_stat_first_moment * m;
} hgt_stat_mean;

typedef struct {
    int correct;
    hgt_stat_second_moment * m2;
    
} hgt_stat_variance;

typedef struct {
    int correct;
    hgt_stat_variance * var;
} hgt_stat_standard_deviation;

// First moment
void hgt_stat_first_moment_increment(hgt_stat_first_moment * m, double v);
hgt_stat_first_moment * hgt_stat_first_moment_alloc();
double hgt_stat_first_moment_get(hgt_stat_first_moment * m);
unsigned long hgt_stat_first_moment_get_n(hgt_stat_first_moment * m);
void hgt_stat_first_moment_free(hgt_stat_first_moment * m);
void hgt_stat_first_moment_clean(hgt_stat_first_moment * m);
// Mean
hgt_stat_mean * hgt_stat_mean_alloc();
void hgt_stat_mean_increment(hgt_stat_mean * m, double v);
double hgt_stat_mean_get(hgt_stat_mean * m);
unsigned long hgt_stat_mean_get_n(hgt_stat_mean * m);
void hgt_stat_mean_free(hgt_stat_mean * m);
void hgt_stat_mean_clean(hgt_stat_mean * m);
// Second moment
hgt_stat_second_moment * hgt_stat_second_moment_alloc();
void hgt_stat_second_moment_increment(hgt_stat_second_moment * m, double v);
double hgt_stat_second_moment_get(hgt_stat_second_moment * m);
unsigned long hgt_stat_second_moment_get_n(hgt_stat_second_moment * m);
void hgt_stat_second_moment_free(hgt_stat_second_moment * m);
void hgt_stat_second_moment_clean(hgt_stat_second_moment * m);
// Variance
hgt_stat_variance * hgt_stat_variance_alloc();
void hgt_stat_variance_increment(hgt_stat_variance * var, double v);
double hgt_stat_variance_get(hgt_stat_variance * var);
unsigned long hgt_stat_variance_get_n(hgt_stat_variance * var);
void hgt_stat_variance_free(hgt_stat_variance * var);
void hgt_stat_variance_clean(hgt_stat_variance * var);

// Standard deviation
hgt_stat_standard_deviation *hgt_stat_standard_deviation_alloc();
void hgt_stat_standard_deviation_increment(hgt_stat_standard_deviation *std, double v);
double hgt_stat_standard_deviation_get(hgt_stat_standard_deviation *std);
unsigned long hgt_stat_standard_deviation_get_n(hgt_stat_standard_deviation *std);
void hgt_stat_standard_deviation_free(hgt_stat_standard_deviation *std);
void hgt_stat_standard_deviation_clean(hgt_stat_standard_deviation *std);

// Mean and Variance.
typedef struct {
    hgt_stat_mean *mean;
    hgt_stat_variance *var;
} hgt_stat_meanvar;
hgt_stat_meanvar *hgt_stat_meanvar_new();
void hgt_stat_meanvar_destroy(hgt_stat_meanvar *mv);
void hgt_stat_meanvar_increment(hgt_stat_meanvar *mv, double v);

// MeanVar list.
typedef struct{
    int n;
    hgt_stat_meanvar **meanvars;
} hgt_stat_meanvar_list;
hgt_stat_meanvar_list *hgt_stat_meanvar_list_new(int n);
void hgt_stat_meanvar_list_destroy(hgt_stat_meanvar_list *l);
void hgt_stat_meanvar_list_increment(hgt_stat_meanvar_list *l, int i, double v);
hgt_stat_meanvar *hgt_stat_meanvar_list_get(hgt_stat_meanvar_list *l, int i);

#endif
