#include "file.h"
#include"config.h"
#include <stdlib.h>

#include <sys/stat.h>   // mkdir
#include <sys/types.h>  // mode_t
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX_RNAP 50
#define MAX_RNAP_SUBUNITS 8 
#define DIM 3

/**
 * 📂 open_simulation_files
 * --------------------------------------------------------------------------------------------------
 * Crée un dossier "Resultats/" (s’il n’existe pas encore) et ouvre tous les fichiers de sortie 
 * de la simulation à l’intérieur de ce dossier.
 *
 * Chaque fichier est ouvert en mode écriture ("w") dans le répertoire :
 *    ./Resultats/
 * Exemple : ./Resultats/brownian_LJ.lammpstrj
 *
 * Si un fichier échoue à s’ouvrir, la fonction affiche une erreur explicite et termine le programme.
 * --------------------------------------------------------------------------------------------------
 */
void open_simulation_files(const Config *cfg, Files *f)
{
    //////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////  📁 Création du répertoire "Resultats"  ///////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////

    printf("\n--- DEBUG CONFIG open_simulation_files() ---\n");
    printf("cfg->nom_fichier ptr = %p\n", (void*)cfg->nom_fichier);
    printf("cfg->nom_fichier = '%s'\n\n", cfg->nom_fichier);

    const char *result_dir = "Resultats";

    struct stat st = {0};

    if (stat(result_dir, &st) == -1) {
        if (mkdir(result_dir, 0777) != 0) {
            perror("❌ Erreur : impossible de créer le dossier 'Resultats'");
            exit(EXIT_FAILURE);
        } else {
            printf("📂 Dossier 'Resultats' créé avec succès.\n");
        }
    } else {
        printf("📂 Dossier 'Resultats' déjà existant.\n");
    }

    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////  ✏️ Mode d'ouverture: "w" ou "a"  ////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////

    // Si reprise depuis un checkpoint → on AJOUTE dans les fichiers existants
    // const char *mode = (cfg->resume_from_checkpoint ? "a" : "w");
    const char *mode = "w";

    char path[512];
    
    // === Fichiers principaux ===

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->nom_fichier);
    f->fichier = fopen(path, mode);

   

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->nom_fichier_equilibre);
    f->fichier_equilibre = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->nom_fichier_rnap_equilibre);
    f->fichier_rnap_equilibre = fopen(path, mode);

    

    // === Tests et suivi ===
    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->nom_test);
    f->test2 = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->nom_test2);
    f->test = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->nom_R_centre_de_masse);
    f->centre_de_masse = fopen(path, mode);

    // === Forces ===
    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->nom_fichier_force);
    f->fichier_force = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->nom_fichier_force_rnap);
    f->fichier_force_rnap = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->nom_fichier_force_rnap_2);
    f->fichier_force_rnap_2 = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->nom_fichier_force_thermique);
    f->fichier_force_thermique = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->nom_fichier_force_rnap_LJ);
    f->fichier_force_rnap_LJ = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->nom_fichier_force_LJ);
    f->fichier_force_LJ = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->nom_fichier_force_lea);
    f->fichier_force_lea = fopen(path, mode);

    // === End-to-end, voisinage, corrélations ===
    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->nom_fichier_endtoend_segment);
    f->fichier_endtoend_segment = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->nom_fichier_endtoend_avant);
    f->fichier_endtoend_avant = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->nom_fichier_endtoend_apres);
    f->fichier_endtoend_apres = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->nom_fichier_endtoend);
    f->fichier_endtoend = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->nom_fichier_voisin);
    f->fichier_voisin = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->nom_fichier_correl_segment);
    f->fichier_correl_segment = fopen(path, mode);

    // === Paramètres et RNAP ===
    snprintf(path, sizeof(path), "%s/param.txt", result_dir);
    f->param = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/rnap.txt", result_dir);
    f->fichier_rnap = fopen(path, mode);

    

    //////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  ✅ Vérification d’ouverture  ////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////

#define CHECK_FILE(ptr, name) \
    if (!(ptr)) { fprintf(stderr, "❌ Erreur : impossible d’ouvrir '%s/%s'\n", result_dir, name); perror("fopen"); exit(EXIT_FAILURE); }

    CHECK_FILE(f->fichier, cfg->nom_fichier);
    CHECK_FILE(f->fichier_equilibre, cfg->nom_fichier_equilibre);
    CHECK_FILE(f->fichier_rnap_equilibre, cfg->nom_fichier_rnap_equilibre);
    CHECK_FILE(f->test, cfg->nom_test2);
    CHECK_FILE(f->test2, cfg->nom_test);
    CHECK_FILE(f->centre_de_masse, cfg->nom_R_centre_de_masse);
    CHECK_FILE(f->fichier_force, cfg->nom_fichier_force);
    CHECK_FILE(f->fichier_force_rnap, cfg->nom_fichier_force_rnap);
    CHECK_FILE(f->fichier_force_rnap_2, cfg->nom_fichier_force_rnap_2);
    CHECK_FILE(f->fichier_force_thermique, cfg->nom_fichier_force_thermique);
    CHECK_FILE(f->fichier_force_LJ, cfg->nom_fichier_force_LJ);
    CHECK_FILE(f->fichier_force_lea, cfg->nom_fichier_force_lea);
    CHECK_FILE(f->fichier_endtoend_segment, cfg->nom_fichier_endtoend_segment);
    CHECK_FILE(f->fichier_endtoend_avant, cfg->nom_fichier_endtoend_avant);
    CHECK_FILE(f->fichier_endtoend_apres, cfg->nom_fichier_endtoend_apres);
    CHECK_FILE(f->fichier_endtoend, cfg->nom_fichier_endtoend);
    CHECK_FILE(f->fichier_voisin, cfg->nom_fichier_voisin);
    CHECK_FILE(f->fichier_correl_segment, cfg->nom_fichier_correl_segment);
    CHECK_FILE(f->param, "param.txt");
    CHECK_FILE(f->fichier_rnap, "rnap.txt");

#undef CHECK_FILE

    printf("✅ Tous les fichiers ont été ouverts dans ./Resultats/ (mode = %s)\n", mode);
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


