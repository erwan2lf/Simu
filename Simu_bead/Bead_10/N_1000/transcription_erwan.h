#ifndef TRANSCRIPTION_H
#define TRANSCRIPTION_H

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "basic_functions.h"
#include "structures_depart.h"
#include "Plots.h"
#include "potentiels.h"
#include "movement.h"
#include <time.h>
#include <omp.h>
#include "config.h"
#include "file.h"
#include "simulation.h"
#include "neighborlist.h"
struct SimVars;
struct Config;


double*** matrix_rnap(int rows, int cols);
double*** initialisation_matrix_rnap_NAN(int rows, int cols);
double*** resize_matrix(double*** matrix, int rows, int cols, int new_rows, int new_cols);
void free_matrix_rnap_2(double*** matrix, int nb_rnap, int rnap_subunits);
void matrix_rnap_0(double*** matrix, int nb_rnap, int nb_subunit);
void copy_matrix(double** R_rnap, double** R_rnap_new, int nb_subunit) ;
void free_matrix_rnap(double*** matrix, int rows, int cols);
int selection_promoteur(double** R, int N, int position_promoteur);
int verification_position_accessible(int promoteur, double** R);
void initialisation_new_RNAP(double** R, double*** R_rnap, int nb_subunits, int promoteur, int rnap);

void enregistrement_RNAP(FILE* fichier, double** R, int N, double*** R_rnap, int nb_rnap, int t, int mono_transcrpt, int* positions_bille_rnap, int nb_subunits, int nb_rnap_initial, int* num_rnap, int sortie);
void enregistrement_RNAP_dyn(FILE* fichier, double** R, int N, double*** R_rnap, int nb_rnap, int t, int* positions_bille_rnap, int position_promoteur, int premiere_rnap, int derniere_rnap, int promoteur_size);
void enregistrement_RNAP_position(FILE* fichier, int nb_rnap, int t, int* positions_bille_rnap, double avancement_transcription);
void enregistrement_promoteur(FILE* fichier, double** R, int N, int nb_rnap, int t, int mono_transcrpt);
int recuperer_RNAP(const char* nom_fichier, int N, int nb_rnap, int nb_subunit, double** R, double*** R_rnap);

double recuperer_positions_rnap(const char* nom_fichier, int nb_rnap, int** positions, double avancement_transcription_recup);

double** mvt_brownian_harmonic_bending(int p, int pmin, int pmax, double** r_new, double** R, double a, double K, double dt);
void mvt_brownian_harmonic_force_RNAP(int p, int pmin, int pmax, int rnap, double** r_new, double** R, double alpha, double K, double Delta, FILE *test, int t, int periode_enregistrement_force, FILE* fichier_force_rnap, FILE* fichier_force_thermique);
double** polymere_brownian_motion_ring_v2(double** R, double a, double K, double dt, int N, double** r_new);
void mvt_brownian_harmonic_force(int p, int pmin, int pmax, double** r_new, double** R, double a, double K, double Delta, FILE* test, int t);
double** polymere_brownian_motion_ring(double** R, double a, double K, double dt, int N, double** r_new);
void polymere_brownian_motion_ring_force(double** R_rnap, double alpha, double K_rnap, double Delta, int N, int rnap, double** r_new, FILE* test, int t,  int periode_enregistrement_force, FILE* fichier_force_rnap, FILE* fichier_force_thermique);
double** bond_rnap_bead(double** R, double** R_rnap, double a_transpt, double K_transpt, double dt, double** r_new, int mono_transcrpt);
double** bond_rnap_bead_progressive_mvt(double** R, double** R_rnap, int rnap_subunits, double a_transpt, double K_transpt, double dt, double** r_new, int positions_bille_rnap, double dx_avancement_rnap, double a, double alpha, FILE* fichier_force_lea, int periode_enregistrement_force, int t);
void liaison_sup(double alpha, double K_transpt, double Delta, double** r_new, int p1, int p2, FILE* fichier_force_rnap_2, int t, int periode_enregistrement_force);
void lennard_jones_forces_rnap(double ***R_rnap, int nb_rnap, int rnap_subunits, double **R, int N, NeighborList_rnap **neighbor_lists, double epsilon, double sigma6, double sigma12, double sigma6rnap2, double sigma12rnap2, double cut_rnap2, double Delta, int t, FILE* test, int T, FILE* fichier_force_rnap, int periode_enregistrement_force);
int actualisation_position_bille_rnap(int positions_bille_rnap, double** R_rnap, double** R, int N, int N_pol, int* rattacher);
void LJ_ouverture_accessibilite(double **R, NeighborList *neighbor_lists, double epsilon, double sigma6, double sigma12, int* positions_bille_rnap, int nb_rnap) ;
double transcrire_nucleotides(double a, double avancement_transcription, int* positions_bille_rnap, int rnap);
void afficher_structure_rnap(double*** R_matrix, int nb_rnap, int nb_subunit) ;
void simu_RNAP_LJ_suite(int N, double** R, int T, int periode_enregistrement, double* parametres, const char* nom_fichier, const char* nom_fichier_rnap, int nb_rnap);
void  simu_RNAP_LJ(int N, double** R, int T, int periode_enregistrement, double* parametres, const char* nom_fichier, const char* nom_fichier_rnap, int nb_rnap, double K_bend, int attache, int plan, int bending, double alpha);
void promoteur_ouvert(double** R, int N, int promoteur, int promoteur_size, double promoteur_rayon_ouverture);
void  simu_RNAP_LJ_structure_gene(int N, double** R, int T, int periode_enregistrement, double* parametres, const char* nom_fichier, const char* nom_fichier_rnap, int nb_rnap, double K_bend, int attache, int plan, int bending, double alpha);
void simu_LJ_RNAP_erwan(const Config *cfg, SimVars *sv, int nbr_simu, const Files *f, int nb_rnap_2, int nb_rnap_initial_2, int rnap_subunits_2);
void augmenter_nombre_rnap(double**** R_rnap, int *nb_rnap, int nb_subunits);
void gonflement(double** R, int N, int *critere, double ecart_monomere);
void ecarter_deux_particules(double* vec1, double* vec2, double ecart_particule);
void creation_1_rnap_erwan(double** R, int nb_rnap, int* positions_bille_rnap, double*** R_rnap, int N, int rnap_subunits, int debut_segment, int fin_segment, int nb_rnap_initial);
void creation_grosse_bille(double** R, int nb_rnap, int* positions_bille_rnap, double*** R_rnap, int N, int debut_segment, int fin_segment, int nb_rnap_initial);
// void creation_RNAP_erwan(double** R, int nb_rnap, int* positions_bille_rnap, double*** R_rnap, int N, int promoteur, int rnap_subunits, int debut_segment, int fin_segment);


#endif // MES_FONCTIONS_H
