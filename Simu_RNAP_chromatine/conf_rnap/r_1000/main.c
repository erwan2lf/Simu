#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include "include/config.h"
#include "include/simulation.h"
#include "transcription.h"
#include "file.h"


#define MAX_RNAP 50
#define MAX_RNAP_SUBUNITS 8 
#define DIM 3

#ifdef CLUSTER
    #include </home/elefloch/Simulation/MT/mt19937ar.h>
#else
    #include </Users/erwan/Documents/These/MTwister/mt19937ar.h>
#endif

double timespec_to_sec(struct timespec t) {
    return t.tv_sec + t.tv_nsec / 1e9;
}

int main(int argc, char*argv[])
{

    // print_banner();
    setvbuf(stdout, NULL, _IONBF, 0);

    static SimVars sv = {0};
    static Files f = {0};

    fflush(stdout);

    // 1) Lire config provenant des arguments
    Config cfg = {0};
    cfg = parse_config(argc, argv);

    int t_start = 0;
    // int checkpoint_found = load_checkpoint_metadata(&cfg, &t_start);


    cfg.resume_from_checkpoint = checkpoint_found;

    // 2) Allouer les buffers AVEC les bonnes dimensions cfg
    init_sim_vars(&sv, &cfg);

    
    // // 3) Charger uniquement les données dans les buffers déjà allocés
    // if (checkpoint_found)
    // {
    //     load_checkpoint_data(&sv, &cfg, &t_start);
    //     printf("✅ Reprise depuis t = %d\n", t_start);
    // }
    // else
    // {
    //     init_genrand(cfg.seed);
    //     printf("🚀 Démarrage neuf\n");
    // }

    init_genrand(cfg.seed);

    // for(int i = 0; i < cfg.N; i++){
    //     printf("R[%d][0] = %lf R[%d][1] = %lf R[%d][2] = %lf \n", i, sv.R[i][0], i, sv.R[i][1], i, sv.R[i][2]);
    // }

        

    printf("Premier nombre aléatoire : %.10f\n", genrand_real2());

    open_simulation_files(&cfg, &f);


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////              Création chromatine           ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // creation_polymere_droit(cfg.N, cfg.a, cfg.ecart_train, sv.R);
    //creation_polymere_aleatoire(cfg.N, cfg.a, sv.R);
    // creation_fractal_globule(N, a, ecart, R);

    if (cfg.resume_from_checkpoint == 0)
    {
        char file_path[512];
    
        #ifdef CLUSTER
            snprintf(
                file_path,
                sizeof(file_path),
                "/home/elefloch/Simulation/Start_simu/N_1000/Simulations/nb-rnap_0/simulation_seed_%lu/brownian_LJ.lammpstrj",
                cfg.seed
            );
        #else 
            snprintf(
                file_path,
                sizeof(file_path),
                "/Users/erwan/Documents/These/Cluster/Start/simulation_seed_%lu/brownian_LJ.lammpstrj",
                cfg.seed
            );
        #endif

        printf("📂 Ouverture du fichier : %s\n", file_path);

        double** R_matrix = recuperer_derniere_structure(file_path, cfg.N);
        if (R_matrix == NULL)
        {
            fprintf(stderr, "Error: Could not read the structure from the file.\n");
            return 1;
        }
        for (int i = 0; i < cfg.N; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                sv.R[i][j] = R_matrix[i][j];
            }
        }
        free_matrix_if_allocated(&R_matrix, cfg.N);

    }
    

    simu_LJ_RNAP_erwan(&cfg, &sv, &f, t_start);

    
    cleanup_sim_vars(&sv, &cfg);
    close_simulation_files(&cfg, &f);
    FILE *fin = fopen(".FINISHED", "w");
    if (fin) fprintf(fin, "DONE\n");
    fclose(fin);
    printf("=== SIMULATION TERMINATED ===\n");

    return 0;
}
