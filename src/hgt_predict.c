//
//  hgt_predict.c
//  hgt
//
//  Created by Mingzhi Lin on 1/21/14.
//
//

#include <stdio.h>
#include "hgt_predict.h"

double hgt_predict_ks_moran(unsigned long size, double mu_rate, double tr_rate, unsigned long frag_len) {
    return (size * mu_rate) / (1.0 + tr_rate * frag_len + 4.0/3.0 * size * mu_rate);
}

double hgt_predict_ks_wf(unsigned long size, double mu_rate, double tr_rate, unsigned long frag_len) {
    return (2 * size * mu_rate) / (1.0 + 2 * tr_rate * frag_len + 4.0/3.0 * 2 * size * mu_rate);
}