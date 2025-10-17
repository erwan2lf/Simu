#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "structures_depart.h"
#include "Plots.h"
#include "basic_functions.h"
#include "potentiels.h"
#include "movement.h"
#include "transcription_erwan.h"
#include <time.h>
#include "config.h"
#include "simulation.h"
#include "file.h"


double timespec_to_sec(struct timespec t) {
    return t.tv_sec + t.tv_nsec / 1e9;
}


// Fonction principale
int main(int argc, char*argv[]){

    // Desactiver le buffering sur stdout 
    setvbuf(stdout, NULL, _IONBF, 0); // La fonction fprintf écrit dans les fichiers en temps réel.

    Config cfg = parse_config(argc, argv);
    SimVars sv = {0}; 
    Files f = {0};
    init_sim_vars(&sv, &cfg);
    open_simulation_files(&cfg, &f);

    int nb_rnap_initial_2 = cfg.nb_rnap_initial; 
    int nb_rnap_2 = sv.nb_rnap; 
    int rnap_subunits_2 = cfg.rnap_subunits;


    //Affichage des paramètres variables 
    printf("Nombre de RNAP     : %d\n", cfg.nb_rnap_initial);
    printf("Vitesse des RNAP   : %f\n", cfg.vitesse_rnap);
    printf("Seed               : %lu\n", cfg.seed); 
    printf("Nombre de monmères : %d \n", cfg.N);

    

    // int bin_count = 1000;
    // double min_value = - 5;
    // double max_value = 5; 
    // int *bins = (int*)malloc(bin_count * sizeof(int));
    // int number_data = 10000000;
    // double *data = (double*)malloc(number_data * sizeof(double)); 
    // for (int i = 0; i < number_data; i++) {
    //     data[i] = randn();
    // }

    // create_histogram(data, number_data, bins, bin_count, min_value, max_value);

    // free_if_allocated((void **)&bins);
    

    for ( int nbr_simu = 0; nbr_simu < cfg.nbr_total_simu; nbr_simu++){
        
        // creation_polymere_droit(cfg.N, cfg.a, cfg.ecart_train, sv.R);
        // creation_polymere_aleatoire(cfg.N, cfg.a, sv.R);
        //Splot_polymere(R, N);
        // printf("oui \n");
        char file_path[512];
        // snprintf(
        //     file_path,
        //     sizeof(file_path),
        //     "/Users/erwan/Documents/These/Cluster/Start/simulation_seed_%lu/brownian_LJ.lammpstrj",
        //     cfg.seed
        // );
        snprintf(
            file_path,
            sizeof(file_path),
            "/home/elefloch/Simulation/Simu/Start/simulation_seed_%lu/brownian_LJ.lammpstrj",
	    cfg.seed
        );
        printf("📂 Ouverture du fichier : %s\n", file_path);
        
        double** R_matrix = recuperer_derniere_structure(file_path, cfg.N);
        if (R_matrix == NULL) {
            fprintf(stderr, "Error: Could not read the structure from the file.\n");
            return 1;
        }
        for (int i = 0; i < cfg.N; i++) {
            for (int j = 0; j < 3; j++) {
                sv.R[i][j] = R_matrix[i][j];
            }
        }
        free_matrix_if_allocated(&R_matrix, cfg.N);
        // printf("non \n");
        //creation_fractal_globule(N, a, ecart, R);
        
        //creation_structure_knot(N, a, R);
        //Splot_polymere(R, N);
        /*if(critere==0){
            Delta = 0.001;
            R = simu(N, R, T, periode_enregistrement, periode_enregistrement_correlation, periode_enregistrement_endtoend, periode_enregistrement_msd, parametres, nom_fichier, stock_msd, Nm, list_monomere, stock_correlation, nbr_simu, nbr_total_simu, stock, K_bend, r_sphere, t_link, attache, confinement, plan, bending);
        }*/

        simu_LJ_RNAP_erwan(&cfg, &sv, nbr_simu, &f, nb_rnap_2, nb_rnap_initial_2, rnap_subunits_2);
}




for ( int i = 0; i < cfg.T_correlation; i++){
    sv.stock_correlation[i] = sv.stock_correlation[i]/cfg.nbr_total_simu;
    sv.stock_correlation_segment[i] = sv.stock_correlation_segment[i]/cfg.nbr_total_simu;
}

//plot_autocorrelation(time, cfg.stock_correlation, cfg.T_correlation, "autocorrelation_data_LJ.txt", "Autocorrelation", "Temps", "Autocorrelation");
//plot_autocorrelation(time, sv->stock_correlation_segment, cfg.T_correlation, "autocorrelation_data_segment_LJ.txt", "Autocorrelation segment", "Temps", "Autocorrelation segment");

cleanup_sim_vars(&sv, &cfg);
close_simulation_files(&cfg, &f);

return 0;

}
