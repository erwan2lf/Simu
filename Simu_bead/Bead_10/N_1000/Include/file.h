#ifndef FILE_H
#define FILE_H

#include <stdio.h>
#include "config.h"

typedef struct {
    FILE* fichier;
    FILE* fichier_equilibre;
    FILE* fichier_rnap_equilibre;

    FILE* test2;
    FILE* test;
    FILE* centre_de_masse;

    FILE* fichier_force;
    FILE* fichier_force_rnap;
    FILE* fichier_force_rnap_2;
    FILE* fichier_force_thermique;
    FILE* fichier_force_rnap_LJ;
    FILE* fichier_force_LJ; 
    FILE* fichier_force_lea;

    FILE* fichier_endtoend_segment;
    FILE* fichier_endtoend_avant;
    FILE* fichier_endtoend_apres;
    FILE* fichier_voisin;
    FILE* fichier_correl_segment;
    FILE* fichier_endtoend;

    FILE* param;
    FILE* fichier_rnap;
} Files;

void open_simulation_files(const Config *cfg, Files *f);
void close_simulation_files(const Config *cfg, Files *f);

#endif 