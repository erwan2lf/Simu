#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "config.h"


Config parse_config(int argc, char *argv[])
{
    Config cfg; 

    // Lecture des arguments (définis dans .bash)
    if (argc < 5)
    {
        fprintf(stderr, "Usage : %s nb_rnap vitesse_rnap seed\n", argv[0]); // Vérification du bon nombre d'arguments
        exit(1);
    }

    cfg.nb_rnap_initial = atoi(argv[1]);
    cfg.vitesse_rnap    = atof(argv[2]); 
    cfg.K_transpt = atof(argv[3]); 
    cfg.seed            = strtoul(argv[4], NULL, 10); 
    

    
    // Paramètres par défaut

    cfg.nbr_total_simu  = 1; 
    cfg.N               = 1000;
    cfg.a               = 1.0;
    cfg.alpha           = 1.0;
    cfg.K               = 1000.0;
    cfg.K_rnap          = 1000.0;
    // cfg.K_transpt       = 1000.0;
    cfg.K_bend          = 0.0;
    cfg.Delta           = 1e-4;
    cfg.epsilon         = 0.0024;
    cfg.epsilon_rnap    = 0.0024;
    cfg.sigma           = cfg.a;
    cfg.sigma6 = pow(cfg.sigma,6);
    cfg.sigma12 = pow(cfg.sigma,12);
    cfg.debut_segment   = 300;
    cfg.fin_segment     = 400;
    cfg.rnap_subunits   = 8;
    cfg.rayon_ecrantage_LJ_chrom = 2.0;
    cfg.rayon_ecrantage_LJ_rnap  = 2.0;
    cfg.ecart_train     = 2;
    cfg.attente_train   = 100000;
    cfg.Nm = cfg.N;
    cfg.r_sphere = 0;

    // Options fixes 
    cfg.attache = 0; // attache
    cfg.confinement = 0; // confinement
    cfg.plan = 0; // plan
    cfg.bending = 0; // bending
    cfg.critere = 1; // critere
    cfg.temperature = 1; // temperature
    cfg.equilibriate = 1; // Mise a l'équilibre du système avant calcul

    // 3.4. Durées et périodicités
    // cfg.T      = 1e7;


    int N_rec = 10000; 
    cfg.T = (int)round(((cfg.fin_segment - cfg.debut_segment) + cfg.ecart_train * (cfg.nb_rnap_initial)) / (cfg.vitesse_rnap * cfg.Delta) );
    printf("T0 = %d\n", cfg.T);
    cfg.T = cfg.T + cfg.T/10;
    printf("T1 = %d\n", cfg.T);
    int k = (cfg.T + N_rec - 1) / N_rec; 
    cfg.T =  k * N_rec;
    printf("Tf = %d\n", cfg.T);

    cfg.T_eq   = cfg.T/10;
    cfg.periode_enregistrement = k;  // periode_enregistrement
    printf("periode enregistrement = %d \n", cfg.periode_enregistrement);
    cfg.periode_msd = k;   // msd
    cfg.periode_correlation = (cfg.T + N_rec - 1)/N_rec;   // correlation
    cfg.periode_endtoend = (cfg.T + N_rec/10 - 1)/(N_rec/10);   // endtoend
    
   
    cfg.periode_centre_de_masse = cfg.T;       // cdm
    cfg.periode_voisin = cfg.T;       // voisins
    cfg.periode_force = cfg.T;       // force

    cfg.T_enregistrement = cfg.T / cfg.periode_enregistrement;
    printf("T_enregistrement = %d \n", cfg.T_enregistrement);
    cfg.T_msd = cfg.T / cfg.periode_msd;
    cfg.T_correlation = cfg.T / cfg.periode_correlation;
    cfg.T_endtoend = cfg.T / cfg.periode_endtoend;
    cfg.T_centre_de_masse = cfg.T / cfg.periode_centre_de_masse;
    cfg.T_voisin = cfg.T / cfg.periode_voisin;
    cfg.T_force = cfg.T / cfg.periode_force; 
    


    // --- RNAP
    cfg.a_rnap = cfg.alpha; 
    cfg.a_transpt = (cfg.alpha + 1)/2;
    cfg.mono_transcrpt = 10; 
    cfg.sigma_rnap = (cfg.alpha + 1)/2; 
    cfg.sigma6_rnap = pow(cfg.sigma_rnap, 6);
    cfg.sigma12_rnap = pow(cfg.sigma_rnap, 12);

    cfg.sigma_rnap2 = cfg.alpha; 
    cfg.sigma6_rnap2 = pow(cfg.sigma_rnap2, 6);
    cfg.sigma12_rnap2 = pow(cfg.sigma_rnap2, 12);

    cfg.dx_avancement_rnap = cfg.vitesse_rnap * cfg.Delta; 

    // --- Fichiers

    cfg.nom_fichier = "brownian_LJ.lammpstrj";
    cfg.nom_fichier_rnap = "brownian_rnap.lammpstrj";

    cfg.nom_fichier_equilibre = "brownian_LJ_equilibre.lammpstrj"; 
    cfg.nom_fichier_rnap_equilibre = "brownian_LJ_rnap_equilibre.lammpstrj";

    cfg.nom_test = "test.txt";
    cfg.nom_test2 = "test2.txt";
    cfg.nom_R_centre_de_masse = "nom_R_centre_de_masse.txt";

    cfg.nom_fichier_force = "fichier_force.txt";
    cfg.nom_fichier_force_rnap = "fichier_force_rnap.txt";
    cfg.nom_fichier_force_rnap_2 = "fichier_force_rnap2.txt";
    cfg.nom_fichier_force_thermique = "fichier_force_rnap2.txt";
    cfg.nom_fichier_force_rnap_LJ = "fichier_force_rnap_LJ.txt";
    cfg.nom_fichier_force_LJ = "fichier_force_LJ.txt";
    cfg.nom_fichier_force_lea = "fichier_force_lea.txt";

    cfg.nom_fichier_endtoend_segment = "endtoend_segment.txt";
    cfg.nom_fichier_endtoend_avant = "endtoend_avant.txt";
    cfg.nom_fichier_endtoend_apres = "endtoend_apres.txt";
    cfg.nom_fichier_voisin = "voisin.txt";
    cfg.nom_fichier_endtoend = "endtoend.txt";
    cfg.nom_fichier_correl_segment = "correl_segment.txt";

    cfg.nom_fichier_lamps = "brownian_LJ.lammpstrj";
    cfg.nom_fichier_rnap_lamps = "brownian_rnap.lammpstrj";

    cfg.nom_fichier_equilibre_lamps = "brownian_LJ_equilibre.lammpstrj"; 
    cfg.nom_fichier_rnap_equilibre_lamps = "brownian_LJ_rnap_equilibre.lammpstrj";


    return cfg;

}
