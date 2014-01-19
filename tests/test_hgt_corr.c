#include <stdio.h>
#include <stdlib.h>
#include <check.h>
#include <math.h>
#include "hgt_corr.h"

START_TEST (test_hgt_corr_fft_auto)
{
	double data[10] = {
		0.1576,
	    0.9706,
	    0.9572,
	    0.4854,
	    0.8003,
	    0.1419,
	    0.4218,
	    0.9157,
	    0.7922,
	    0.9595
	};
	double expected[32] = {
		5.34401, 3.98031, 3.13718, 2.4438, 1.88223, 2.46069, 2.17929, \
1.83166, 1.05614, 0.151217,
		0,0,0,0,0,0,0,0,0,0,0,0,0,
	    0.151217, 1.05614, 1.83166, 2.17929, 2.46069, 1.88223, 2.4438, \
3.13718, 3.98031
	};

	int k;
	double *f, *g, *r, *q;
	f = (double *) calloc(32, sizeof(double));
	g = (double *) calloc(32, sizeof(double));
	r = (double *) calloc(32, sizeof(double));
	q = (double *) calloc(32, sizeof(double));
	for (k = 0; k < 10; k++)
	{
		f[k] = data[k];
		g[k] = data[k];
		r[k] = data[k];
		q[k] = data[k];
	}
	hgt_corr_fft(f, g, 16);
	hgt_corr_brute_force(r, q, 16);

	for (k = 0; k < 32; k++)
	{
//		printf("%d, %g, %g, %g\n", k, expected[k], f[k], r[k]);
		ck_assert_msg(fabs(f[k] - expected[k]) < 1E-5, "Was expecting %g, but got %g, at %d", expected[k], f[k], k);
		ck_assert_msg(fabs(r[k] - expected[k]) < 1E-5, "Was expecting %g, but got %g, at %d", expected[k], r[k], k);
		ck_assert_msg(fabs(r[k] - f[k]) < 1E-10, "Was expecting %g, but got %g, at %d", f[k], r[k], k);
	}
}
END_TEST

START_TEST (test_hgt_corr_fft_cros)
{
	double data1[10] = {
		0.6557,
	    0.0357,
	    0.8491,
	    0.9340,
	    0.6787,
	    0.7577,
	    0.7431,
	    0.3922,
	    0.6555,
	    0.1712
	};
	double data2[10] = {
		0.1576,
	    0.9706,
	    0.9572,
	    0.4854,
	    0.8003,
	    0.1419,
	    0.4218,
	    0.9157,
	    0.7922,
	    0.9595
	};
	double expected[32] = {
		3.41092, 3.86624, 3.40214, 2.79604, 3.00792, 2.27675, 1.87809, 
1.44342, 0.5537, 0.629144,
		0,0,0,0,0,0,0,0,0,0,0,0,0,
		0.0269811, 0.269474, 0.861912, 1.20833, 1.67127, 2.29295, 2.37102,
3.14141, 3.66636
	};
	int k;
	double *f, *g, *r, *q;
	f = (double *) calloc(32, sizeof(double));
	g = (double *) calloc(32, sizeof(double));
	r = (double *) calloc(32, sizeof(double));
	q = (double *) calloc(32, sizeof(double));
	for (k = 0; k < 10; k++)
	{
		f[k] = data2[k];
		g[k] = data1[k];
		r[k] = data2[k];
		q[k] = data1[k];
	}
	hgt_corr_fft(f, g, 16);
	hgt_corr_brute_force(r, q, 16);

	for (k = 0; k < 32; k++)
	{
//		printf("%d, %.10g, %.10g, %.10g\n", k, expected[k], f[k], r[k]);
		ck_assert_msg(fabs(f[k] - expected[k]) < 1E-5, "Was expecting %g, but got %g, at %d", expected[k], f[k], k);
		ck_assert_msg(fabs(r[k] - expected[k]) < 1E-5, "Was expecting %g, but got %g, at %d", expected[k], r[k], k);
		ck_assert_msg(fabs(r[k] - f[k]) < 1E-10, "Was expecting %g, but got %g, at %d", f[k], r[k], k);
	}
}
END_TEST

Suite *
hgt_corr_suite (void)
{
	Suite *s = suite_create("hgt_corr_fft");
	TCase *tc_core = tcase_create("Core");
	tcase_add_test(tc_core, test_hgt_corr_fft_auto);
	tcase_add_test(tc_core, test_hgt_corr_fft_cros);
	suite_add_tcase(s, tc_core);

	return s;
}

int main(void)
{
    int number_failed;
	Suite *s = hgt_corr_suite ();
	SRunner *sr = srunner_create (s);
	srunner_run_all (sr, CK_NORMAL);
	number_failed = srunner_ntests_failed (sr);
	srunner_free (sr);
	return (number_failed == 0) ? EXIT_SUCCESS : CK_FAILURE;
}
