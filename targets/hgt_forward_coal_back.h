#ifndef hgt_forward_coal_back_h
#define hgt_forward_coal_back_h
typedef struct _file_container file_container;
struct _file_container {
    FILE *coal_p2;
	FILE *forw_p2;
    FILE *t2;
    FILE *t2_spl;
};
file_container * create_file_container(char * prefix);
int destroy_file_container(file_container *fc);
int close_file_container(file_container *fc);
int coal_evolve(hgt_params *params, int size, int length, double time, hgt_stat_meanvar_list *list, const gsl_rng *r);
int calc_pxy(hgt_genome **genomes, int size, int maxl, hgt_stat_meanvar_list *list);
double * calc_p2(hgt_genome *g1, hgt_genome *g2, int maxl);
int update_pxy(hgt_stat_meanvar_list *list, double *pxy);
int write_pxy(FILE *f, hgt_stat_meanvar_list *list, int maxl, unsigned long gen);
hgt_genome ** create_genomes(int size, int length, const gsl_rng * r);
int destroy_genomes(hgt_genome **genomes, int size);
int mutate_genomes(hgt_genome **genomes, int size, double mutation_rate, double time, const gsl_rng *r);
double *compare_genomes(hgt_genome *g1, hgt_genome *g2);
int random_seq(char * seq, int seq_len, const gsl_rng * r);
int sample(hgt_pop **ps, hgt_params *params,
           int linkage_size, file_container *files,
           const gsl_rng *r);
int update_t2(hgt_stat_meanvar *mv, double t2);
int write_t2(FILE *f, hgt_stat_meanvar *mv, unsigned long gen);
int write_t2_all(FILE *f, double t2, unsigned int gen);
int flush_file_container(file_container *fc);
#endif