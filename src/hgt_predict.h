//
//  hgt_predict.h
//  hgt
//
//  Created by Mingzhi Lin on 1/21/14.
//
//

#ifndef hgt_hgt_predict_h
#define hgt_hgt_predict_h

typedef double (*hgt_predict_ks_func) (unsigned long size, double mu_rate, double tr_rate, unsigned long frag_len);
double hgt_predict_ks_moran(unsigned long size, double mu_rate, double tr_rate, unsigned long frag_len);
double hgt_predict_ks_wf(unsigned long size, double mu_rate, double tr_rate, unsigned long frag_len);

#endif
