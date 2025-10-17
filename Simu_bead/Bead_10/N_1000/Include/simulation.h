#ifndef SIMULATION_H
#define SIMULATION_H

#include <time.h>        
#include "config.h"  
#include "file.h"    
#include "neighborlist.h" 


// Structure qui regroupe tout l'état/l’output de ta simu
typedef struct {
    // --- compteurs et chronométrage ---
    int     bidule, truc, prout, compteur_grand_deplacement;
    double  duree_boucle, duree_tot, temps_restant;
    int premier_passage;
    int compteur_train; 
    int variable_test_2;
    int variable_test_3;


    // --- paramètres physiques (copiés depuis cfg->parametres[]) ---
    double K, a, Delta;
    double rayon, epaisseur;
    double epsilon, sigma;
    double K_rnap, K_transpt, K_bend;
    double r_sphere, vitesse_rnap;

    // --- buffer rnap 
    int nb_rnap; // Nombre de rnap à un instant t
    int* positions_bille_rnap;

    double*** R_rnap;
    double*** R_rnap_new;

    double* avancement_transcription;

    int* compteur_mono_rnap;

    int* num_rnap;
    int sortie;



    // --- buffers pour corrélations, MSD, end-to-end, giration, etc. ---
    double *stock_correlation;
    double *stock_correlation_segment;
    double *time;
    double *log10_time;

    double **stock;
    double **Rbb_segment;
    double **Rbb_avant;
    double **Rbb_apres;
    int    *nombre_voisins;

    double **Rbb;
    double **R_segment;

    double **R_centre_de_masse;
    double **stock_cdm;
    double *gyration_radius;

    double ***stock_msd;
    double ****R_monomere_arrays;

    double **t_link;
    double **R;
    double **R_new;
    int *list_monomere;

    double **bending_forces;
    double **R_matrix;

    // --- Periode

    int T_enregistrement;
    int T_msd;
    int T_correlation;
    int T_endtoend;
    int T_centre_de_masse;
    int T_force;
    int T_voisin;


    int val_rnap; 


    
} SimVars;

// Prototypes des fonctions de simulation
void init_sim_vars      (SimVars *sv, const Config *cfg);
void cleanup_sim_vars(SimVars *sv, Config *cfg);
void ajouter_rnap(SimVars *sv, const Config *cfg, NeighborList *neighbor_lists, NeighborList_rnap*** p_neighbor_lists_rnap, int t);
void retirer_rnap(SimVars *sv, const Config *cfg, NeighborList *neighbor_list, NeighborList_rnap** neighbor_lists_rnap, int rnap, int prout, int t);
void retirer_rnap_2(SimVars *sv, const Config *cfg, NeighborList *neighbor_list, NeighborList_rnap*** p_neighbor_lists_rnap, int rnap, int prout, int t);
void enregistrement_data(SimVars *sv, const Config *cfg, const Files *f, int t, double endtoend, int nbr_simu);
void finalize_simulationaaaaa(SimVars *sv, const Config *cfg, const Files *f, double endtoend, int nbr_simu);
void f_equilibriate(SimVars *sv, const Config *cfg, const Files *f, NeighborList *neighbor_lists, NeighborList_rnap **neighbor_lists_rnap, int nbr_simu);
void calcul(SimVars *sv, const Config *cfg, const Files *f, NeighborList *neighbor_lists, NeighborList_rnap **neighbor_lists_rnap, int nbr_simu, double endtoend, int nb_rnap_2, int nb_rnap_initial_2, int rnap_subunits_2);
void thermic_noise_rnap(SimVars *sv, const Config *cfg, int rnap);

#endif // SIMULATION_H