//
//  test_hgt_cov.c
//  hgt
//
//  Created by Mingzhi Lin on 1/20/14.
//
//

#include <stdio.h>
#include <check.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include "hgt_cov.h"

const gsl_rng *RNG;

void init() {
    gsl_rng * alloc_rng();
    
    RNG = alloc_rng();
}

gsl_rng * alloc_rng() {
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    return r;
}

START_TEST(test_hgt_cov_sample_p2)
{
    int i, n;
    unsigned long p_size;
    p_size = 4;
    
    double count[p_size];
    for (i = 0; i < p_size; i++) {
        count[i] = 0;
    }
    
    unsigned long a, b, c, d;
    a = b = c = d = 0;
    
    n = 100000;
    p_size = 4;
    for (i = 0; i < n; i++) {
        hgt_cov_sample_p2(&a, &b, &c, &d, p_size, RNG);
        ck_assert_int_ne(a, b);
        ck_assert_int_eq(a, c);
        ck_assert_int_eq(b, d);
        count[a]++;
    }
    
    for (i = 0; i < p_size; i++) {
        count[i] /= (double)n;
        ck_assert_msg(fabs(count[i] - 1.0/(double)p_size) < 0.001, "was expecting %g, but got %g, at %d\n", 1.0/(double)p_size, count[i], i);
    }
}
END_TEST

START_TEST(test_hgt_cov_sample_p3)
{
    int i, n;
    unsigned long p_size;
    p_size = 4;
    
    double count[p_size];
    for (i = 0; i < p_size; i++) {
        count[i] = 0;
    }
    
    unsigned long a, b, c, d;
    a = b = c = d = 0;
    
    n = 100000;
    p_size = 4;
    for (i = 0; i < n; i++) {
        hgt_cov_sample_p3(&a, &b, &c, &d, p_size, RNG);
        ck_assert_int_ne(a, b);
        ck_assert_int_ne(a, c);
        ck_assert_int_eq(a, d);
        count[a]++;
    }
    
    for (i = 0; i < p_size; i++) {
        count[i] /= (double)n;
        ck_assert_msg(fabs(count[i] - 1.0/(double)p_size) < 0.001, "was expecting %g, but got %g, at %d\n", 1.0/(double)p_size, count[i], i);
    }
}
END_TEST

START_TEST(test_hgt_cov_sample_p4)
{
    int i, n;
    unsigned long p_size;
    p_size = 4;
    
    double count[p_size];
    for (i = 0; i < p_size; i++) {
        count[i] = 0;
    }
    
    unsigned long a, b, c, d;
    a = b = c = d = 0;
    
    n = 100000;
    p_size = 4;
    for (i = 0; i < n; i++) {
        hgt_cov_sample_p4(&a, &b, &c, &d, p_size, RNG);
        ck_assert_int_ne(a, b);
        ck_assert_int_ne(a, c);
        ck_assert_int_ne(a, d);
        ck_assert_int_ne(b, c);
        ck_assert_int_ne(b, d);
        ck_assert_int_ne(c, d);
        count[a]++;
    }
    
    for (i = 0; i < p_size; i++) {
        count[i] /= (double)n;
        ck_assert_msg(fabs(count[i] - 1.0/(double)p_size) < 0.01, "was expecting %g, but got %g, at %d\n", 1.0/(double)p_size, count[i], i);
    }
}
END_TEST

Suite *
hgt_cov_suite (void)
{
    Suite *s = suite_create("hgt_cov");
    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_hgt_cov_sample_p2);
    tcase_add_test(tc_core, test_hgt_cov_sample_p3);
    tcase_add_test(tc_core, test_hgt_cov_sample_p4);
    suite_add_tcase(s, tc_core);
    
    return s;
}

int main(void)
{
    init();
    int number_failed;
	Suite *s = hgt_cov_suite ();
	SRunner *sr = srunner_create (s);
	srunner_run_all (sr, CK_NORMAL);
	number_failed = srunner_ntests_failed (sr);
	srunner_free (sr);
	return (number_failed == 0) ? EXIT_SUCCESS : CK_FAILURE;
}


