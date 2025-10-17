#ifndef BASIC_FUNCTIONS_H
#define BASIC_FUNCTIONS_H

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>

// Déclaration des fonctions
double randn();
double distance(double* p1, double* p2);
double distance2(double* p1, double* p2);
double norm(double* vec);
double norm2(double* vec);
void afficher_structure(double** R_matrix, int N);
double** allocate_matrix(size_t rows, size_t cols);
void free_if_allocated(void **ptr);
void free_matrix_if_allocated(double ***matrix, int rows);
void free_matrix_cube_if_allocated(double ***matrix, int rows, int cols);
void free_matrix_4D(double ****matrix, int Nm, int T_msd);
int* particules_outside(double** R, double* origine, double rayon, double epaisseur, int N);
typedef struct {double moyenne;double std;} Mesures; Mesures calcul_mesures (double ** R, int N);
typedef struct {double ete; double etes;} Mesures_temp; Mesures_temp calcul_mesures_temp (double ** R, int N, int t);
void enregistrement(FILE* fichier, double** R_matrix, int N, int t);
double dot_product(double *vec1, double *vec2) ;
void calculate_autocorrelation(double** Rbb, double* stock_correl, int Tp, int N);
bool is_in_list(int *list, int size, int num);
void generate_unique_random_numbers(int *list, int N, int count);
void regression_lineaire(double* x, double* y, int taille, double* alpha, double* logA);
double somme(double* data, int taille);
void print_time_remaining(int current_iteration, int total_iterations, time_t start_time);
void calculate_msd(double ****R_monomere_arrays,double ***stock_msd, int Tp, int Nm, int *list_monomere, int nbr_simu, int nbr_total_simu);
void plot_msd(double **stock_msd, int Tp, int Nm, int *list_monomer, const char *msd_filename, const char *msd_log_filename);
void plot_cols(double **array, int rows, int cols, int column_index, const char *filename, const char *title, const char *xlabel, const char *ylabel);
void plot_autocorrelation(double *time, double *autocorrelation, int length, const char *filename, const char *title, const char *xlabel, const char *ylabel);
void plot_columns_from_file(const char *filename, const char *title, const char *xlabel, const char *ylabel);
double calculate_angle(double *vec1, double *vec2);
void update_link_vectors(double **r_new, double **t_link, int N);
bool is_allocated(void *ptr);
void compteur_grands_deplacements(int N, int T, double **R, double **r_new, int compteur);
double calculate_mean_for_correlation(double **R, int i, int j);
void create_histogram(double *data, int data_size, int *bins, int bin_count, double min_value, double max_value);
int find_nearest_particles(double **R, int N, int particle);
bool is_present(int *liste, int taille, int valeur);
void echanger_deux_particules(double **R, int i, int j);
void reogarnisation(double **R, int N);
void four1(double* data, int n, int Isign);
void realft(double* vec_double, double* data, int Isign, int data_size);
void save_gyration(double* gyration_radius, int nbr_simu, int T_centre_de_masse, int periode_centre_de_masse);
void save_centre_de_masse(double** R_centre_de_masse, int N, int T_centre_de_masse, int nbr_simu, int periode_centre_de_masse);
void calculate_gyration_radius(double** R_centre_de_masse, double* gyration_radius, double** R, int N, int T_centre_de_masse, double** stock_cdm);
void creation_polymere_aleatoire(int N, double a, double **R);
#endif // MES_FONCTIONS_H
