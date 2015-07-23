#ifndef hgt_forward_coal_back_h
#define hgt_forward_coal_back_h
typedef struct _file_container file_container;
struct _file_container {
    FILE *p2;
};
file_container * create_file_container(char * prefix);
int destroy_file_container(file_container *fc);
int close_file_container(file_container *fc);
int coal_evolve(hgt_params *params, int size, int length, unsigned long time, unsigned long gen, file_container *files, const gsl_rng *r);
int calc_pxy(hgt_genome **genomes, int size, int maxl, file_container *files, unsigned long gen);
int calc_p2(hgt_genome *g1, hgt_genome *g2, int maxl, FILE *f, unsigned long gen);
int write_pxy(FILE *f, double *pxy, int maxl, unsigned long gen);
hgt_genome ** create_genomes(int size, int length, const gsl_rng * r);
int destroy_genomes(hgt_genome **genomes, int size);
int mutate_genomes(hgt_genome **genomes, int size, double mutation_rate, unsigned long time, const gsl_rng *r);
int compare_genomes(hgt_genome *g1, hgt_genome *g2, double *distance, int length);
int random_seq(char * seq, int seq_len, const gsl_rng * r);
int sample(hgt_pop **ps, hgt_params *params,
           int linkage_size, unsigned long gen, file_container *files,
           const gsl_rng *r);
#endif