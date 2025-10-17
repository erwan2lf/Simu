#include "file.h"
#include"config.h"
#include <stdlib.h>

void open_simulation_files(const Config *cfg, Files *f){

    f->fichier = fopen(cfg->nom_fichier, "w");
    f->fichier_equilibre = fopen(cfg->nom_fichier_equilibre, "w");
    f->fichier_rnap_equilibre = fopen(cfg->nom_fichier_rnap_equilibre,"w");

    f->test2 = fopen(cfg->nom_test, "w");
    f->test = fopen(cfg->nom_test2, "w");
    f->centre_de_masse = fopen(cfg->nom_R_centre_de_masse, "w");

    f->fichier_force = fopen(cfg->nom_fichier_force, "w");
    f->fichier_force_rnap = fopen(cfg->nom_fichier_force_rnap, "w");
    f->fichier_force_rnap_2 = fopen(cfg->nom_fichier_force_rnap_2, "w");
    f->fichier_force_thermique = fopen(cfg->nom_fichier_force_thermique, "w");
    f->fichier_force_rnap_LJ = fopen(cfg->nom_fichier_force_rnap_LJ, "w");
    f->fichier_force_LJ = fopen(cfg->nom_fichier_force_LJ, "w"); 
    f->fichier_force_lea = fopen(cfg->nom_fichier_force_lea, "w");

    f->fichier_endtoend_segment = fopen(cfg->nom_fichier_endtoend_segment, "w");
    f->fichier_endtoend_avant = fopen(cfg->nom_fichier_endtoend_avant, "w");
    f->fichier_endtoend_apres = fopen(cfg->nom_fichier_endtoend_apres, "w");
    f->fichier_voisin = fopen(cfg->nom_fichier_voisin,"w");
    f->fichier_correl_segment = fopen(cfg->nom_fichier_correl_segment,"w");
    f->fichier_endtoend = fopen(cfg->nom_fichier_endtoend,"w");

    f->param = fopen("param.txt", "w");
    f->fichier_rnap = fopen("rnap.txt", "w");

    if (!f->fichier || !f->fichier_equilibre || !f->fichier_rnap_equilibre || !f->test || !f->test2 || !f->centre_de_masse || !f->fichier_force || !f->fichier_force_rnap || !f->fichier_force_rnap_2 || !f->fichier_force_thermique || !f->fichier_force_LJ || !f->fichier_force_LJ || !f->fichier_force_lea || !f->fichier_endtoend_segment || !f->fichier_endtoend_avant || !f->fichier_endtoend_apres || !f->fichier_endtoend || !f->fichier_voisin || !f->fichier_correl_segment || !f->param || !f->fichier_rnap) {
        perror("Ouverture d'un fichier de simulation");
        exit(EXIT_FAILURE);
    }

}

void close_simulation_files(const Config *cfg, Files *f){
        
    fclose(f->fichier);
    fclose(f->fichier_equilibre);
    fclose(f->fichier_rnap_equilibre);

    fclose(f->test);
    fclose(f->test2);
    fclose(f->centre_de_masse);

    fclose(f->centre_de_masse);
    fclose(f->fichier_force_rnap);
    fclose(f->fichier_force_rnap_2);
    fclose(f->fichier_force_thermique);
    fclose(f->fichier_force_rnap_LJ);
    fclose(f->fichier_force_LJ);
    fclose(f->fichier_force_lea);

    fclose(f->fichier_endtoend);
    fclose(f->fichier_endtoend_avant);
    fclose(f->fichier_endtoend_segment);
    fclose(f->fichier_endtoend_apres);
    fclose(f->fichier_voisin);
    fclose(f->fichier_correl_segment);
    
}