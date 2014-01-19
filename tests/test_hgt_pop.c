//
//  test_hgt_pop.c
//  hgt
//
//  Created by Mingzhi Lin on 1/18/14.
//
//

#include <stdio.h>
#include <check.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "hgt_pop.h"

START_TEST (test_hgt_pop_alloc)
{
    unsigned long size, seq_len;
    hgt_pop * p;
    int i, k;
    
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    
    size = 10;
    seq_len = 10;
    
    p = hgt_pop_alloc(size, seq_len, r);
    
    ck_assert_msg(p->size == size, "Was expecting population size to be %lu, but got %lu", size, p->size);
    ck_assert_msg(p->seq_len == seq_len, "Was expecting genome length to be %lu, but got %lu", seq_len, p->seq_len);

    for (i = 0; i < size; i ++) {
        for (k = 0; k < seq_len; k ++) {
            ck_assert_msg(p->genomes[0][k] == p->genomes[i][k], "The nucl was not identical at genome %d at position %d", i, k);
        }
    }
}
END_TEST

START_TEST (test_hgt_pop_copy)
{
    unsigned long size, seq_len;
    hgt_pop * p, * old_p;
    int i, k;
    
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    
    size = 10;
    seq_len = 10;
    
    old_p = hgt_pop_alloc(size, seq_len, r);
    p = hgt_pop_copy(old_p);
    
    ck_assert_msg(p->size == size, "Was expecting population size to be %lu, but got %lu", size, p->size);
    ck_assert_msg(p->seq_len == seq_len, "Was expecting genome length to be %lu, but got %lu", seq_len, p->seq_len);

    for (i = 0; i < size; i ++) {
        for (k = 0; k < seq_len; k ++) {
            ck_assert_msg(p->genomes[0][k] == p->genomes[i][k], "The nucl was not identical at genome %d at position %d", i, k);
        }
    }
}
END_TEST

START_TEST (test_hgt_pop_params_parse)
{
    int argc, exit_code;
    char *argv1[] = {"test_hgt_pop_params_parse", 
                    "-n", "1000", 
                    "-l", "100", 
                    "-u", "1e-3", 
                    "-t", "1e-4", 
                    "-f", "50",
                    "-g", "10000"
                };
    argc = 13;
    hgt_pop_params * params = malloc(sizeof(hgt_pop_params));
    exit_code = hgt_pop_params_parse(params, argc, argv1, "test_hgt_pop_params_parse");
    
    ck_assert(exit_code == EXIT_SUCCESS);
    ck_assert_msg(params->size == 1000, "Was expecting population size to be %lu, but got %lu", 1000, params->size);
    ck_assert_msg(params->seq_len == 100, "Was expecting genome length to be %lu, but got %lu", 100, params->seq_len);
    ck_assert_msg(params->frag_len == 50, "Was expecting frag_len to be %lu, but got %lu", 50, params->frag_len);
    ck_assert_msg(fabs(params->mu_rate - 1e-3) < 1e-10, "Was expecting mu_rate to be %g, but got %g", 1e-3, params->mu_rate);
    ck_assert_msg(fabs(params->tr_rate - 1e-4) < 1e-10, "Was expecting tr_rate to be %g, but got %g", 1e-4, params->tr_rate);

    char *argv2[] = {"test_hgt_pop_params_parse", "-n", "1000", 
                    "-l", "100", 
                    "-u", "1e-3", 
                    "-t", "1e-4", 
                    "-f", "50",
                    "-g", "10000",
                    "-s", "1000",
                    "-i", "222",
                    "-r", "55"
                };
    argc = 19;
    hgt_pop_params_parse(params, argc, argv2, "test_hgt_pop_params_parse");
    ck_assert(exit_code == EXIT_SUCCESS);
    ck_assert(params->generations == 10000);
    ck_assert(params->sample_size == 1000);
    ck_assert(params->sample_time == 222);
    ck_assert(params->replicates == 55);

    char *argv3[] = {"test_hgt_pop_params_parse", "-n", "1000", 
                    "-l", "100", 
                    "-u", "1e-3", 
                    "-t", "1e-4", 
                    "-f", "50"
                };
    argc = 11;
    exit_code = hgt_pop_params_parse(params, argc, argv3, "test_hgt_pop_params_parse");
    ck_assert(exit_code == EXIT_FAILURE);
}   
END_TEST

Suite *
hgt_pop_suite (void)
{
    Suite *s = suite_create("hgt_pop");
    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_hgt_pop_alloc);
    tcase_add_test(tc_core, test_hgt_pop_copy);
    tcase_add_test(tc_core, test_hgt_pop_params_parse);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s = hgt_pop_suite ();
    SRunner *sr = srunner_create (s);
    srunner_run_all (sr, CK_NORMAL);
    number_failed = srunner_ntests_failed (sr);
    srunner_free (sr);
    return (number_failed == 0) ? EXIT_SUCCESS : CK_FAILURE;
}