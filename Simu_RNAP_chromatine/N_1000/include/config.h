#ifndef CONFIG_H
#define CONFIG_H

typedef struct {
    int nb_rnap_initial; 
    double vitesse_rnap;
    double K_transpt; // Contante de raideur chrom-rnap
    unsigned long seed; 

    // --- Paramètres chromatine
    int nbr_total_simu; 
    int N; // Nombre de monomère dans la chromatine dans chaque chaine
    int M; // = 2*N, nombre de monomères total
    double a; // diamètre d'une bille
    double K; // Constante de raideur chromatine
    double K_bend; //Module de courbure
    double Delta; // Dimensionless timestep
    double epsilon; // Potentiel de LJ Chrom
    double sigma; // sigma LJ
    double sigma6; 
    double sigma12;
    int Nm;
    double r_sphere; 

    // Confinement sphère
    double r_conf;
    double epsilon_conf; 
    double sigma_conf;

    // Conhésines

    int K_c; // Constante de raideur pour les ressorts entre les deux chaînes

    
    // --- Paramètres RNAP
    int rnap_subunits; // Nombre de sous unités d'une rnap
    double alpha; // rapport de taille entre bille de la chrom et bille rnap
    double K_rnap; // Constante de raideur RNAP
    double epsilon_rnap; // Potentiel de LJ RNAP
    double rayon_ecrantage_LJ_chrom; 
    double rayon_ecrantage_LJ_rnap; 
    int ecart_train; // Nombre de bille entre deux rnap
    int attente_train; 
    int debut_segment; // Début du segment transcrit
    int fin_segment; // Fin du segment transcrit
    double a_rnap; 
    double a_transpt; 
    int mono_transcrpt; 
    double sigma_rnap; 
    double sigma6_rnap; 
    double sigma12_rnap; 

    double sigma_rnap2; 
    double sigma6_rnap2; 
    double sigma12_rnap2;

    double dx_avancement_rnap;

    // Temps
    int T; // Temps total de simu
    int T_eq; // Temps de simu sans RNAP

    // Fichier .txt
    char nom_fichier[256];
    char nom_fichier_rnap[256];

    char nom_fichier_equilibre[256]; 
    char nom_fichier_rnap_equilibre[256];

    char nom_test[256]; 
    char nom_test2[256];

    char nom_R_centre_de_masse[256];

    char nom_fichier_force[256];
    char nom_fichier_force_rnap[256];
    char nom_fichier_force_rnap_2[256];
    char nom_fichier_force_thermique[256];
    char nom_fichier_force_rnap_LJ[256];
    char nom_fichier_force_LJ[256];
    char nom_fichier_force_lea[256];

    char nom_fichier_endtoend_segment[256];
    char nom_fichier_endtoend_avant[256];
    char nom_fichier_endtoend_apres[256]; 
    char nom_fichier_voisin[256]; 
    char nom_fichier_endtoend[256];
    char nom_fichier_correl_segment[256];

    char nom_fichier_lamps[256];
    char nom_fichier_rnap_lamps[256];

    char nom_fichier_equilibre_lamps[256]; 
    char nom_fichier_rnap_equilibre_lamps[256];


    // --- options 
    int attache;
    int confinement;
    int plan;
    int bending;
    int temperature;
    int critere;
    int equilibriate; 
    int quench;
    int cohesine;

    int freq_cohesine;

    // --- périodes d'enregistrement
    int periode_enregistrement;
    int periode_msd;
    int periode_correlation;
    int periode_endtoend;
    int periode_centre_de_masse;
    int periode_force;
    int periode_voisin; 

    int T_enregistrement;
    int T_msd;
    int T_correlation;
    int T_endtoend;
    int T_centre_de_masse;
    int T_force;
    int T_voisin;


    //// Restart de simu 
    int resume_from_checkpoint; 

} Config; 

Config parse_config(int argc, char *argv[]);
void print_header(const char *title);
void print_banner(void);

#endif // CONFIG_H

