//
//  hgt_stat.c
//  hgt
//
//  Created by Mingzhi Lin on 1/24/14.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hgt_stat.h"

/*
 * First moment
 */
void hgt_stat_first_moment_increment(hgt_stat_first_moment * m, double v) {
    if (m->n == 0) {
        m->m1 = 0;
    }
    m->n++;
    unsigned long n0 = m->n;
    m->dev = v - m->m1;
    m->n_dev = m->dev / (double) n0;
    m->m1 += m->n_dev;
}

hgt_stat_first_moment * hgt_stat_first_moment_alloc() {
    hgt_stat_first_moment * moment;
    moment = calloc(1, sizeof(hgt_stat_first_moment));
    return moment;
    
}

double hgt_stat_first_moment_get(hgt_stat_first_moment * m) {
    return m->m1;
}

unsigned long hgt_stat_first_moment_get_n(hgt_stat_first_moment * m) {
    return m->n;
}

void hgt_stat_first_moment_free(hgt_stat_first_moment * m) {
    free(m);
}

void hgt_stat_first_moment_clean(hgt_stat_first_moment * m) {
    m->dev = 0;
    m->m1 = 0;
    m->n = 0;
    m->n_dev = 0;
}

/*
 * Mean
 */

hgt_stat_mean * hgt_stat_mean_alloc() {
    hgt_stat_first_moment * moment;
    hgt_stat_mean * mean;
    
    moment = hgt_stat_first_moment_alloc();
    mean = calloc(1, sizeof(hgt_stat_mean));
    mean->m = moment;
    
    return mean;
}

void hgt_stat_mean_increment(hgt_stat_mean * m, double v) {
    hgt_stat_first_moment_increment(m->m, v);
}

double hgt_stat_mean_get(hgt_stat_mean * m) {
    return hgt_stat_first_moment_get(m->m);
}

unsigned long hgt_stat_mean_get_n(hgt_stat_mean * m) {
    return hgt_stat_first_moment_get_n(m->m);
}

void hgt_stat_mean_free(hgt_stat_mean * m) {
    hgt_stat_first_moment_free(m->m);
    free(m);
}

void hgt_stat_mean_clean(hgt_stat_mean * m) {
    hgt_stat_first_moment_clean(m->m);
}

/*
 * Second moment
 */

hgt_stat_second_moment * hgt_stat_second_moment_alloc() {
    hgt_stat_first_moment * m1;
    hgt_stat_second_moment * m2;
    m1 = hgt_stat_first_moment_alloc();
    m2 = calloc(1, sizeof(hgt_stat_second_moment));
    m2->m1 = m1;
    
    return m2;
}

void hgt_stat_second_moment_increment(hgt_stat_second_moment * m, double v) {
    if (m->m1->n < 1) {
        m->m1->m1 = 0;
        m->m2 = 0;
    }
    hgt_stat_first_moment_increment(m->m1, v);
    m->m2 += ((double) m->m1->n - 1.0) * m->m1->dev * m->m1->n_dev;
}

double hgt_stat_second_moment_get(hgt_stat_second_moment * m) {
    return m->m2;
}

unsigned long hgt_stat_second_moment_get_n(hgt_stat_second_moment * m) {
    return m->m1->n;
}

void hgt_stat_second_moment_free(hgt_stat_second_moment * m) {
    hgt_stat_first_moment_free(m->m1);
    free(m);
}

void hgt_stat_second_moment_clean(hgt_stat_second_moment * m) {
    hgt_stat_first_moment_clean(m->m1);
    m->m2 = 0;
}

/*
 * Variance
 */
hgt_stat_variance * hgt_stat_variance_alloc() {
    hgt_stat_second_moment * m2;
    hgt_stat_variance * var;
    
    m2 = hgt_stat_second_moment_alloc();
    var = calloc(1, sizeof(hgt_stat_variance));
    var->m2 = m2;
    
    return var;
}

void hgt_stat_variance_increment(hgt_stat_variance * var, double v) {
    hgt_stat_second_moment_increment(var->m2, v);
}

double hgt_stat_variance_get(hgt_stat_variance * var) {
    unsigned long n = hgt_stat_second_moment_get_n(var->m2);
    if (n == 0) {
        //return NAN;
        return -1;
    } else if (n == 1) {
        return 0;
    } else {
        if (var->correct == 0) {
            return hgt_stat_second_moment_get(var->m2) / (double) n;
        } else {
            return hgt_stat_second_moment_get(var->m2) / ((double) n - 1.0);
        }
    }
}

unsigned long hgt_stat_variance_get_n(hgt_stat_variance * var) {
    return hgt_stat_second_moment_get_n(var->m2);
}

void hgt_stat_variance_free(hgt_stat_variance * var) {
    hgt_stat_second_moment_free(var->m2);
    free(var);
}

void hgt_stat_variance_clean(hgt_stat_variance * var) {
    hgt_stat_second_moment_clean(var->m2);
}

hgt_stat_standard_deviation *hgt_stat_standard_deviation_alloc() {
    hgt_stat_standard_deviation *std = malloc(sizeof(hgt_stat_standard_deviation));
    std->var = hgt_stat_variance_alloc();
    
    return std;
}
void hgt_stat_standard_deviation_increment(hgt_stat_standard_deviation *std, double v){
    hgt_stat_variance_increment(std->var, v);
}
double hgt_stat_standard_deviation_get(hgt_stat_standard_deviation *std){
    std->var->correct = std->correct;
    return sqrt(hgt_stat_variance_get(std->var));
}
unsigned long hgt_stat_standard_deviation_get_n(hgt_stat_standard_deviation *std){
    return hgt_stat_variance_get_n(std->var);
}
void hgt_stat_standard_deviation_free(hgt_stat_standard_deviation *std){
    hgt_stat_variance_free(std->var);
    free(std);
}
void hgt_stat_standard_deviation_clean(hgt_stat_standard_deviation *std){
    hgt_stat_variance_clean(std->var);
}

