//
//  test_hgt_pop.c
//  hgt
//
//  Created by Mingzhi Lin on 1/18/14.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <check.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
#include "hgt_pop.h"

hgt_pop *SMALL_P;
const gsl_rng *RNG;
const int SMALL_SIZE = 10;
const int SMALL_SEQ_LEN = 10;

void init() {
    gsl_rng * alloc_rng();

    RNG = alloc_rng();
    SMALL_P = hgt_pop_alloc(SMALL_SIZE, SMALL_SEQ_LEN, RNG);

}

gsl_rng * alloc_rng() {
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    return r;
}

START_TEST (test_hgt_pop_alloc)
{
    unsigned long size, seq_len;
    hgt_pop * p;
    int i, k;
    
    const gsl_rng * r = RNG;
    
    size = SMALL_SIZE;
    seq_len = SMALL_SEQ_LEN;
    
    p = hgt_pop_alloc(size, seq_len, r);
    
    ck_assert_msg(p->size == size, "Was expecting population size to be %lu, but got %lu", size, p->size);
    ck_assert_msg(p->seq_len == seq_len, "Was expecting genome length to be %lu, but got %lu", seq_len, p->seq_len);

    for (i = 0; i < size; i ++) {
        for (k = 0; k < seq_len; k ++) {
            ck_assert_msg(p->genomes[0][k] == p->genomes[i][k], "The nucl was not identical at genome %d at position %d", i, k);
        }
    }

    hgt_pop_free(p);
}
END_TEST

START_TEST (test_hgt_pop_copy)
{
    unsigned long size, seq_len;
    hgt_pop * p, * old_p;
    int i, k;
    
    const gsl_rng * r = RNG;
    
    size = SMALL_SIZE;
    seq_len = SMALL_SEQ_LEN;
    
    old_p = hgt_pop_alloc(size, seq_len, r);
    p = hgt_pop_copy(old_p);
    
    ck_assert_msg(p->size == size, "Was expecting population size to be %lu, but got %lu", size, p->size);
    ck_assert_msg(p->seq_len == seq_len, "Was expecting genome length to be %lu, but got %lu", seq_len, p->seq_len);

    for (i = 0; i < size; i ++) {
        for (k = 0; k < seq_len; k ++) {
            ck_assert_msg(p->genomes[0][k] == p->genomes[i][k], "The nucl was not identical at genome %d at position %d", i, k);
        }
    }

    hgt_pop_free(old_p);
    hgt_pop_free(p);
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

    hgt_pop_params_free(params);
}   
END_TEST

START_TEST (test_hgt_pop_mutate_at)
{
    unsigned long size, seq_len;
    hgt_pop * p;
    int n, i;
    unsigned long g, s;
    char c;
    int count[4];
    
    const gsl_rng *r = RNG;
    
    p = SMALL_P;
    size = SMALL_SIZE;
    seq_len = SMALL_SEQ_LEN;

    g = 1;
    s = 7;
    c = p->genomes[g][s];
    hgt_pop_mutate_at(p, g, s, r);
    ck_assert(p->genomes[g][s] != c);

    for (i = 0; i < 4; i++){
        count[i] = 0;
    }

    n = 100000;
    for (i = 0; i < n; i++) {
        hgt_pop_mutate_at(p, g, s, r);
        switch (p->genomes[g][s]) {
            case 'A':
                count[0]++;
                break;
            case 'C':
                count[1]++;
                break;
            case 'G':
                count[2]++;
                break;
            case 'T':
                count[3]++;
            default:
                break;
        }
    }

    for (i = 0; i < 4; i++) {
        ck_assert_msg(fabs((double)count[i]/(double)n  - 0.25 ) < 1e-3,
            "was expecting %g, but got %g, at %d", 0.25, (double)count[i]/(double)n, i
            );
    }
}
END_TEST

START_TEST (test_hgt_pop_mutate)
{
    unsigned long size, seq_len;
    hgt_pop * p;
    int n, i;
    unsigned long g, s;
    int count[4];
    
    const gsl_rng *r = RNG;
    
    size = SMALL_SIZE;
    seq_len = SMALL_SEQ_LEN;
    p = SMALL_P;

    n = 1000000;
    g = gsl_rng_uniform_int(r, p->size);
    s = gsl_rng_uniform_int(r, p->seq_len);
    for (i = 0; i < 4; i++){
        count[i] = 0;
    }

    for (i = 0; i < n; i++) {
        hgt_pop_mutate(p, r);
        switch (p->genomes[g][s]) {
            case 'A':
                count[0]++;
                break;
            case 'C':
                count[1]++;
                break;
            case 'G':
                count[2]++;
                break;
            case 'T':
                count[3]++;
            default:
                break;
        }
    }

    for (i = 0; i < 4; i++) {
        ck_assert_msg(fabs((double)count[i]/(double)n  - 0.25 ) < 1e-2,
            "was expecting %g, but got %g, at %d", 0.25, (double)count[i]/(double)n, i
            );
    }
}
END_TEST

START_TEST (test_hgt_pop_transfer_at)
{
    const gsl_rng *r = RNG;
    hgt_pop * p = SMALL_P;
    unsigned long donor, reciever, start, frag_len;
    int n, i, k, s;

    n = 1000;
    for (i = 0; i < n; i++) {
        hgt_pop_mutate(p, r);
    }

    for (i = 0; i < n; i++) {
        donor = gsl_rng_uniform_int(r, p->size);
        reciever = gsl_rng_uniform_int(r, p->size);
        while (donor == reciever) {
            reciever = gsl_rng_uniform_int(r, p->size);
        }
        start = gsl_rng_uniform_int(r, p->seq_len);
        frag_len = gsl_rng_uniform_int(r, p->seq_len);
        hgt_pop_transfer_at(p, donor, reciever, frag_len, start);
        for (k = 0; k < frag_len; k++) {
            s = (start + k) % p->size;
            ck_assert_msg(p->genomes[donor][s] == p->genomes[reciever][s], 
                "not identical at %d, between %lu (%s) and %lu (%s)", s, donor, p->genomes[donor][s], reciever, p->genomes[reciever][s]);
        }
    }


}
END_TEST

int evolve_and_calc_ks(double * ks, hgt_pop_params *params, hgt_pop_sample_func sample_f, hgt_pop_coal_time_func c_time_f, const gsl_rng *r) {
    int i, m;
    for (m = 0; m < params->replicates; m++) {
        hgt_pop *p = hgt_pop_alloc(params->size, params->seq_len, r);
        for (i = 0; i < params->generations; i++) {
            hgt_pop_evolve(p, params, sample_f, c_time_f, r);
        }
        
        ks[m] = hgt_pop_calc_ks(p);
        hgt_pop_free(p);
    }

    return EXIT_SUCCESS;
}

START_TEST (test_hgt_pop_evolve) {
    int evolve_and_calc_ks(double * ks, hgt_pop_params *params, hgt_pop_sample_func sample_f, hgt_pop_coal_time_func c_time_f, const gsl_rng *r);
    const gsl_rng *r = RNG;
    unsigned long size, seq_len;

    size = 10;
    seq_len = 100;

    hgt_pop_params *params = malloc(sizeof(hgt_pop_params));
    params->size = size;
    params->seq_len = seq_len;
    params->frag_len = 10;
    params->mu_rate = 0.1;
    params->replicates = 100;

    double expect, mean, sd;
    double ks[params->replicates];

    // test evolution without horizontal gene transfer
    params->tr_rate = 0;

    params->generations = 1000;
    evolve_and_calc_ks(ks, params, hgt_pop_sample_moran, hgt_pop_coal_time_moran, r);
    expect = (size * params->mu_rate) / (1.0 + 4.0/3.0 * size * params->mu_rate);
    mean = gsl_stats_mean(ks, 1, params->replicates);
    sd = gsl_stats_sd(ks, 1, params->replicates)/sqrt(params->size);
    ck_assert_msg(fabs(mean - expect) < sd, "was expecting %g, but got %g, and the standard deviation is %g\n", expect, mean, sd);
    
    params->generations = 100;
    evolve_and_calc_ks(ks, params, hgt_pop_sample_wf, hgt_pop_coal_time_wf, r);
    expect = (2 * size * params->mu_rate) / (1.0 + 4.0/3.0 * 2 * size * params->mu_rate);
    mean = gsl_stats_mean(ks, 1, params->replicates);
    sd = gsl_stats_sd(ks, 1, params->replicates)/sqrt(params->size);
    ck_assert_msg(fabs(mean - expect) < sd, "was expecting %g, but got %g, and the standard deviation is %g\n", expect, mean, sd);

    // test evolution with horizontal gene transfer
    params->tr_rate = 0.1;

    params->generations = 1000;
    evolve_and_calc_ks(ks, params, hgt_pop_sample_moran, hgt_pop_coal_time_moran, r);
    expect = (size * params->mu_rate) / (1.0 + params->tr_rate * params->frag_len + 4.0/3.0 * size * params->mu_rate);
    mean = gsl_stats_mean(ks, 1, params->replicates);
    sd = gsl_stats_sd(ks, 1, params->replicates)/sqrt(params->size);
    ck_assert_msg(fabs(mean - expect) < sd, "was expecting %g, but got %g, and the standard deviation is %g\n", expect, mean, sd);
    
    params->generations = 100;
    evolve_and_calc_ks(ks, params, hgt_pop_sample_wf, hgt_pop_coal_time_wf, r);
    expect = (2 * size * params->mu_rate) / (1.0 + 2 * params->tr_rate * params->frag_len + 4.0/3.0 * 2 * size * params->mu_rate);
    mean = gsl_stats_mean(ks, 1, params->replicates);
    sd = gsl_stats_sd(ks, 1, params->replicates)/sqrt(params->size);
    ck_assert_msg(fabs(mean - expect) < sd, "was expecting %g, but got %g, and the standard deviation is %g\n", expect, mean, sd);

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

    TCase *tc_mut = tcase_create("Mutate");
    tcase_add_test(tc_mut, test_hgt_pop_mutate_at);
    tcase_add_test(tc_mut, test_hgt_pop_mutate);
    suite_add_tcase(s, tc_mut);

    TCase *tc_tra = tcase_create("Transfer");
    tcase_add_test(tc_tra, test_hgt_pop_transfer_at);
    suite_add_tcase(s, tc_tra);

    TCase *tc_evl = tcase_create("Evolve");
    tcase_set_timeout(tc_evl, 100);
    tcase_add_test(tc_evl, test_hgt_pop_evolve);
    suite_add_tcase(s, tc_evl);

    return s;
}

int main(void)
{
    void init();

    init();
    int number_failed;
    Suite *s = hgt_pop_suite ();
    SRunner *sr = srunner_create (s);
    srunner_run_all (sr, CK_NORMAL);
    number_failed = srunner_ntests_failed (sr);
    srunner_free (sr);
    return (number_failed == 0) ? EXIT_SUCCESS : CK_FAILURE;
}
