#ifndef hgt_hgt_corr_h
#define hgt_hgt_corr_h

int hgt_corr_fft(double *f, double *g, unsigned long len);
int hgt_corr_brute_force(double *f, double *g, unsigned long len);
int hgt_corr_auto_fft(double *f, unsigned long len);

#endif
