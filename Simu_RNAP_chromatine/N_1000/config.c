#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "config.h"

#define MAX_RNAP 50


Config parse_config(int argc, char *argv[])
{
    Config cfg = {0}; 

    

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////// Lecture des arguments (définis dans .bash) ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    if (argc < 5)
    {
        fprintf(stderr, "Usage : %s nb_rnap vitesse_rnap seed\n", argv[0]); // Vérification du bon nombre d'arguments
        exit(1);
    }

    cfg.nb_rnap_initial = atoi(argv[1]);
    cfg.vitesse_rnap    = atof(argv[2]); 
    cfg.K_transpt       = atof(argv[3]); 
    cfg.seed            = strtoul(argv[4], NULL, 10);

    // print_header("Lecture des arguments (définis dans .bash)");

    // printf("Nombre de RNAP : %d \n", cfg.nb_rnap_initial);
    // printf("Vitesse des RNAP en a/Tau_0 : %lf \n", cfg.vitesse_rnap);
    // printf("Constant de raideur RNAP - chromatine : %lf \n", cfg.K_transpt);
    // printf("Seed  : %lu \n", cfg.seed);

    


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////            Paramètres par défaut           ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    cfg.nbr_total_simu  = 1; 
    cfg.N               = 1000;
    cfg.a               = 1;
    cfg.alpha           = 1.0;
    cfg.K               = 1000.0;
    cfg.K_rnap          = 500.0;
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
    cfg.r_sphere = 0;
    cfg.r_conf = 10.77;
    cfg.epsilon_conf = 0.0024; 
    cfg.sigma_conf = cfg.a;

    // print_header("Paramètres par défaut");


    // printf("Nombre de monomères : %d\n", cfg.N);
    // printf("Taille des monomères : %lf\n", cfg.a);
    // printf("Taille des RNAP (en unité de a) : %lf\n", cfg.alpha); 
    // printf("Constante de raideur de la chromatine : %lf\n", cfg.K);
    // printf("Constante de raideur RNAP-RNAP : %lf\n", cfg.K_rnap);
    // printf("Module de courbure: %lf\n", cfg.K_bend);
    // printf("Pas de temps : %lf ", cfg.Delta);
    // printf("Intensité de LJ chromatine-chromatine : %lf\n",cfg.epsilon);
    // printf("Intensité de LJ RNAP-RNAP : %lf\n",cfg.epsilon_rnap);
    // printf("Sigma LJ (en unitté de a) : %lf\n", cfg.sigma); 
    // printf("Début du segment transcrit : %d \n", cfg.debut_segment); 
    // printf("Fin du segment transcrit : %d \n", cfg.debut_segment); 
    // printf("Nombre de sous unités par RNAP : %d\n", cfg.rnap_subunits); 
    // printf("Rayon écrantage LJ chromatine-chromatine : %lf \n", cfg.rayon_ecrantage_LJ_chrom);
    // printf("Rayon écrantage LJ RNAP-RNAO : %lf \n", cfg.rayon_ecrantage_LJ_rnap);
    // printf("Nombre de billes entre deux RNAP : %d\n", cfg.ecart_train);
    // printf("Taille de la sphère de confinement (si 0 pas de confinement) : %d\n",cfg.r_sphere);


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////               Options fixes                ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    cfg.attache = 0; // attache
    cfg.confinement = 1; // confinement
    cfg.plan = 0; // plan
    cfg.bending = 0; // bending
    cfg.critere = 1; // critere
    cfg.temperature = 1; // temperature
    cfg.equilibriate = 1; // Mise a l'équilibre du système avant calcul
    cfg.quench = 0;
        


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////           Durées et périodicitéss           ///////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // print_header("Durées et périodicitéss");

    int N_rec = 1000; 
    int k;

    if(cfg.nb_rnap_initial == 0)
    {

        cfg.T = (int)round(((cfg.fin_segment - cfg.debut_segment) + cfg.ecart_train * (MAX_RNAP)) / (cfg.vitesse_rnap * cfg.Delta) );
        cfg.T = cfg.T + cfg.T/10;
        // cfg.T = 13200000;
        printf("T0 = %d\n", cfg.T);
        cfg.T_eq = 0;

        k = (cfg.T + N_rec - 1) / N_rec; 
        cfg.periode_enregistrement = k;  // periode_enregistrement
        printf("periode enregistrement = %d \n", cfg.periode_enregistrement);
        cfg.periode_msd = k;   // msd
        cfg.periode_correlation = (cfg.T + N_rec - 1)/N_rec;   // correlation
        cfg.periode_endtoend = (cfg.T + N_rec/10 - 1)/(N_rec/10);   // endtoend

        cfg.periode_centre_de_masse = cfg.T;       // cdm
        cfg.periode_voisin = cfg.T;       // voisins
        cfg.periode_force = cfg.T;       // force
    }
    else
    {
        cfg.T = (int)round(((cfg.fin_segment - cfg.debut_segment) + cfg.ecart_train * (cfg.nb_rnap_initial)) / (cfg.vitesse_rnap * cfg.Delta) );
        printf("T0 = %d\n", cfg.T);
        cfg.T = cfg.T + cfg.T/10;
        printf("T1 = %d\n", cfg.T);
        k = (cfg.T + N_rec - 1) / N_rec; 
        cfg.T =  k * N_rec;
        cfg.T + 10000
        printf("Tf = %d\n", cfg.T);

        cfg.periode_enregistrement = k;  // periode_enregistrement
        printf("periode enregistrement = %d \n", cfg.periode_enregistrement);
        cfg.periode_msd = k;   // msd
        cfg.periode_correlation = (cfg.T + N_rec - 1)/N_rec;   // correlation
        cfg.periode_endtoend = (cfg.T + N_rec/10 - 1)/(N_rec/10);   // endtoend
        
    
        cfg.periode_centre_de_masse = cfg.T;       // cdm
        cfg.periode_voisin = cfg.T;       // voisins
        cfg.periode_force = cfg.T;       // force
    }
    
    if(cfg.quench == 1)
    {
        cfg.vitesse_rnap = 0;
    }
    
    

    
    // cfg.periode_enregistrement = k;  // periode_enregistrement
    // printf("periode enregistrement = %d \n", cfg.periode_enregistrement);
    // cfg.periode_msd = k;   // msd
    // cfg.periode_correlation = (cfg.T + N_rec - 1)/N_rec;   // correlation
    // cfg.periode_endtoend = (cfg.T + N_rec/10 - 1)/(N_rec/10);   // endtoend
    
   
    // cfg.periode_centre_de_masse = cfg.T;       // cdm
    // cfg.periode_voisin = cfg.T;       // voisins
    // cfg.periode_force = cfg.T;       // force

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


    // === Noms des fichiers ===
    snprintf(cfg.nom_fichier, sizeof(cfg.nom_fichier), "brownian_LJ.lammpstrj");
    snprintf(cfg.nom_fichier_rnap, sizeof(cfg.nom_fichier_rnap), "brownian_rnap.lammpstrj");

    snprintf(cfg.nom_fichier_equilibre, sizeof(cfg.nom_fichier_equilibre), "brownian_LJ_equilibre.lammpstrj");
    snprintf(cfg.nom_fichier_rnap_equilibre, sizeof(cfg.nom_fichier_rnap_equilibre), "brownian_LJ_rnap_equilibre.lammpstrj");

    snprintf(cfg.nom_test, sizeof(cfg.nom_test), "test.txt");
    snprintf(cfg.nom_test2, sizeof(cfg.nom_test2), "test2.txt");
    snprintf(cfg.nom_R_centre_de_masse, sizeof(cfg.nom_R_centre_de_masse), "centre_de_masse.txt");

    snprintf(cfg.nom_fichier_force, sizeof(cfg.nom_fichier_force), "fichier_force.txt");
    snprintf(cfg.nom_fichier_force_rnap, sizeof(cfg.nom_fichier_force_rnap), "fichier_force_rnap.txt");
    snprintf(cfg.nom_fichier_force_rnap_2, sizeof(cfg.nom_fichier_force_rnap_2), "fichier_force_rnap2.txt");
    snprintf(cfg.nom_fichier_force_thermique, sizeof(cfg.nom_fichier_force_thermique), "fichier_force_thermique.txt");
    snprintf(cfg.nom_fichier_force_rnap_LJ, sizeof(cfg.nom_fichier_force_rnap_LJ), "fichier_force_rnap_LJ.txt");
    snprintf(cfg.nom_fichier_force_LJ, sizeof(cfg.nom_fichier_force_LJ), "fichier_force_LJ.txt");
    snprintf(cfg.nom_fichier_force_lea, sizeof(cfg.nom_fichier_force_lea), "fichier_force_lea.txt");

    snprintf(cfg.nom_fichier_endtoend_segment, sizeof(cfg.nom_fichier_endtoend_segment), "endtoend_segment.txt");
    snprintf(cfg.nom_fichier_endtoend_avant, sizeof(cfg.nom_fichier_endtoend_avant), "endtoend_avant.txt");
    snprintf(cfg.nom_fichier_endtoend_apres, sizeof(cfg.nom_fichier_endtoend_apres), "endtoend_apres.txt");
    snprintf(cfg.nom_fichier_endtoend, sizeof(cfg.nom_fichier_endtoend), "endtoend.txt");
    snprintf(cfg.nom_fichier_voisin, sizeof(cfg.nom_fichier_voisin), "voisin.txt");
    snprintf(cfg.nom_fichier_correl_segment, sizeof(cfg.nom_fichier_correl_segment), "correl_segment.txt");

    // === fichiers LAMMPS ===
    snprintf(cfg.nom_fichier_lamps, sizeof(cfg.nom_fichier_lamps), "brownian_LJ.lammpstrj");
    snprintf(cfg.nom_fichier_rnap_lamps, sizeof(cfg.nom_fichier_rnap_lamps), "brownian_rnap.lammpstrj");

    snprintf(cfg.nom_fichier_equilibre_lamps, sizeof(cfg.nom_fichier_equilibre_lamps), "brownian_LJ_equilibre.lammpstrj");
    snprintf(cfg.nom_fichier_rnap_equilibre_lamps, sizeof(cfg.nom_fichier_rnap_equilibre_lamps), "brownian_LJ_rnap_equilibre.lammpstrj");


        ///// Restart de simu 

    cfg.resume_from_checkpoint = 0; 


    return cfg;

}

// === Codes de couleur ANSI ===
#define C_RESET   "\033[0m"
#define C_CYAN    "\033[36m"
#define C_YELLOW  "\033[1;33m"
#define C_BOLD    "\033[1m"

// === Fonction d’affichage ===
void print_header(const char *title)
{
    const int total_width = 180;   // largeur totale de la ligne
    const int inner_padding = 4;   // espaces autour du titre
    const char border_char = '/';  // caractère des bordures

    int title_len = strlen(title);
    int total_inner = title_len + 2 * inner_padding;

    // Si le titre est trop long pour rentrer dans la largeur
    if (total_inner >= total_width - 4) {
        printf(C_CYAN "%c%c %s %c%c\n" C_RESET,
               border_char, border_char, title, border_char, border_char);
        return;
    }

    int side_width = (total_width - total_inner) / 2;

    // Ligne supérieure
    printf(C_CYAN);
    for (int i = 0; i < total_width; i++) putchar(border_char);
    printf(C_RESET "\n");

    // Ligne du milieu avec titre centré et coloré
    printf(C_CYAN);
    for (int i = 0; i < side_width; i++) putchar(border_char);
    printf(C_RESET);

    printf(C_YELLOW "%*s%s%*s" C_RESET, inner_padding, "", title, inner_padding, "");

    printf(C_CYAN);
    for (int i = 0; i < side_width; i++) putchar(border_char);
    printf(C_RESET "\n");

    // Ligne inférieure
    printf(C_CYAN);
    for (int i = 0; i < total_width; i++) putchar(border_char);
    printf(C_RESET "\n");
}

#define C_RESET   "\033[0m"
#define C_BOLD    "\033[1m"
#define C_BLUE    "\033[38;5;39m"
#define C_CYAN    "\033[36m"
#define C_YELLOW  "\033[1;33m"
#define C_WHITE   "\033[97m"
#define C_GRAY    "\033[90m"

void print_banner(void) {
    printf("\n");
    printf(C_CYAN "//////////////////////////////////////////////////////////////////////////////////////////////////////////\n");
    printf(C_CYAN "////" C_RESET C_BOLD C_YELLOW "       Modeling of RNA Polymerase II Driven Chromatin Dynamics Simulation Code      " C_RESET C_CYAN "////\n");
    printf(C_CYAN "//////////////////////////////////////////////////////////////////////////////////////////////////////////\n" C_RESET);

    printf(C_WHITE "\n");
    printf("  📘 " C_BOLD "Code développé au " C_BLUE "Laboratoire de Physique Théorique (LPT - CNRS, Toulouse)" C_RESET "\n");
    printf("  👨‍🔬 " C_BOLD "Auteur : " C_YELLOW "Erwan Le Floch" C_RESET "\n");
    printf("  🧬 " C_BOLD "Projet de thèse : " C_WHITE "\"Modeling of RNA Polymerase II Driven Chromatin Dynamics\"" C_RESET "\n");
    printf("  🗓️  " C_BOLD "Dernière mise à jour : " C_GRAY __DATE__ " - " __TIME__ C_RESET "\n");
    printf("  ⚙️  " C_BOLD "Compilateur : " C_GRAY "gcc / g++ avec OpenMP et optimisations HPC" C_RESET "\n");
    printf("  🧩 " C_BOLD "Simulation : " C_WHITE "Chromatin polymer + RNAP motor complex (coarse-grained BD)" C_RESET "\n");

    printf("\n" C_CYAN "//////////////////////////////////////////////////////////////////////////////////////////////////////////" C_RESET "\n\n");
}


