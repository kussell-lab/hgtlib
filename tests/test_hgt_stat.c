//
//  test_hgt_stat.c
//  hgt
//
//  Created by Mingzhi Lin on 1/24/14.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <check.h>
#include "hgt_stat.h"

double MEAN;
double VARIANCE;
double STD;
double SECOND_MOMENT;
double *TEST_ARRAY;
int TEST_ARRAY_LENGTH;
void init() {
    MEAN = 12.404545454545455;
    VARIANCE = 10.00235930735931;
    SECOND_MOMENT = 210.04954545454547;
    STD = sqrt(VARIANCE);
    TEST_ARRAY_LENGTH = 22;
    TEST_ARRAY = malloc(TEST_ARRAY_LENGTH*sizeof(double));
    double testArray[22] = {12.5, 12.0, 11.8, 14.2, 14.9, 14.5, 21.0, 8.2, 10.3, 11.3,
        14.1, 9.9, 12.2, 12.0, 12.1, 11.0, 19.8, 11.0, 10.0, 8.8,
        9.0, 12.3};
    int i;
    for (i = 0; i < TEST_ARRAY_LENGTH; i++) {
        TEST_ARRAY[i] = testArray[i];
    }
    
}

START_TEST(test_hgt_stat_desc)
{
    hgt_stat_first_moment *m = hgt_stat_first_moment_alloc();
    hgt_stat_mean *mean = hgt_stat_mean_alloc();
    hgt_stat_second_moment *sm = hgt_stat_second_moment_alloc();
    hgt_stat_variance *var = hgt_stat_variance_alloc();
    hgt_stat_standard_deviation *std = hgt_stat_standard_deviation_alloc();
    
    var->correct = 1;
    std->correct = 1;
    
    ck_assert_int_eq(hgt_stat_first_moment_get_n(m), 0);
    ck_assert(fabs(hgt_stat_first_moment_get(m) - 0) < 1e-10);
    
    int i;
    for (i = 0; i < TEST_ARRAY_LENGTH; i++) {
        double v = TEST_ARRAY[i];
        hgt_stat_first_moment_increment(m, v);
        hgt_stat_mean_increment(mean, v);
        hgt_stat_second_moment_increment(sm, v);
        hgt_stat_variance_increment(var, v);
        hgt_stat_standard_deviation_increment(std, v);
    }
    
    ck_assert_int_eq(hgt_stat_first_moment_get_n(m), TEST_ARRAY_LENGTH);
    ck_assert_msg(fabs(hgt_stat_first_moment_get(m) - MEAN) < 1e-10, "was expecting mean to be %g, but got %g\n", MEAN, hgt_stat_first_moment_get(m));
    ck_assert_msg(fabs(hgt_stat_mean_get(mean) - MEAN) < 1e-10, "was expecting mean to be %g, but got %g\n", MEAN, hgt_stat_mean_get(mean));
    ck_assert_msg(fabs(hgt_stat_second_moment_get(sm) - SECOND_MOMENT) < 1e-10, "was expecting second moment to be %g, but got %g\n", SECOND_MOMENT, hgt_stat_second_moment_get(sm));
    ck_assert_msg(fabs(hgt_stat_variance_get(var) - VARIANCE) < 1e-10, "was expecting variance to be %g, but got %g\n", VARIANCE, hgt_stat_variance_get(var));
    ck_assert_msg(fabs(hgt_stat_standard_deviation_get(std) - STD) < 1e-10, "was expecting standard deviation to be %g, but got %g\n", STD, hgt_stat_standard_deviation_get(std));
}
END_TEST

Suite *
hgt_stat_suite (void) {
    Suite *s = suite_create("hgt_stat");
    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_hgt_stat_desc);
    suite_add_tcase(s, tc_core);
    
    return s;
}

int main(void){
    void init();
    init();
    int number_failed;
    Suite *s = hgt_stat_suite();
    SRunner *sr = srunner_create (s);
    srunner_run_all (sr, CK_NORMAL);
    number_failed = srunner_ntests_failed (sr);
    srunner_free (sr);
    return (number_failed == 0) ? EXIT_SUCCESS : CK_FAILURE;
}