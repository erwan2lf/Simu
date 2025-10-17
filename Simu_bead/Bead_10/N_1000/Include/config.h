#ifndef CONFIG_H
#define CONFIG_H

typedef struct {
    int nb_rnap_initial; 
    double vitesse_rnap;
    double K_transpt; // Contante de raideur chrom-rnap
    unsigned long seed; 

    // --- Paramètres chromatine
    int nbr_total_simu; 
    int N; // Nombre de monomère dans la chromatine
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
    const char* nom_fichier;
    const char* nom_fichier_rnap;

    const char* nom_fichier_equilibre; 
    const char* nom_fichier_rnap_equilibre;

    const char* nom_test; 
    const char* nom_test2;

    const char* nom_R_centre_de_masse;

    const char* nom_fichier_force;
    const char* nom_fichier_force_rnap;
    const char* nom_fichier_force_rnap_2;
    const char* nom_fichier_force_thermique;
    const char* nom_fichier_force_rnap_LJ;
    const char* nom_fichier_force_LJ; 
    const char* nom_fichier_force_lea;

    const char* nom_fichier_endtoend_segment;
    const char* nom_fichier_endtoend_avant;
    const char* nom_fichier_endtoend_apres; 
    const char* nom_fichier_voisin; 
    const char* nom_fichier_endtoend;
    const char* nom_fichier_correl_segment;

    const char* nom_fichier_lamps;
    const char* nom_fichier_rnap_lamps;

    const char* nom_fichier_equilibre_lamps; 
    const char* nom_fichier_rnap_equilibre_lamps;


    // --- options 
    int attache;
    int confinement;
    int plan;
    int bending;
    int temperature;
    int critere;
    int equilibriate; 

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

} Config; 

Config parse_config(int argc, char *argv[]);

#endif // CONFIG_H