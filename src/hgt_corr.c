#include <stdlib.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

/* 
 * hgt_corr_fft
 * calculate cross correlation using fft.
 */
int hgt_corr_fft(double *f, double *g, unsigned long len)
{
	unsigned long k, ft_len;
	double a, b, c, d;

	ft_len = 2 * len;
	gsl_fft_real_radix2_transform(f, 1, ft_len);
    gsl_fft_real_radix2_transform(g, 1, ft_len);
    
    for (k = 1; k < len; k++) {
        a = f[k];
        b = f[ft_len - k];
        c = g[k];
        d = g[ft_len - k];
        f[k] = (a * c + b * d);
        f[ft_len - k] = (b * c - a * d);
    }
    
    f[0] = f[0]*g[0];
	f[len] = f[len]*g[len];
	gsl_fft_halfcomplex_radix2_inverse(f, 1, ft_len);

	return EXIT_SUCCESS;
}

int hgt_corr_brute_force(double *f, double *g, unsigned long len)
{
	unsigned long i, k, ft_len;
	double *s;

	ft_len = 2 * len;
	s = calloc(ft_len, sizeof(double));

	for (i = 0; i < len; i++) {
		for (k = 0; k < len; k++) {
			s[i] += f[k] * g[(ft_len+k-i)%ft_len];
			s[ft_len-i] += g[k] * f[(ft_len+k-i)%ft_len];
		}
	}

	for (i = 0; i < ft_len; i++) {
		f[i] = s[i];
	}

	free(s);
	return EXIT_SUCCESS;
}
