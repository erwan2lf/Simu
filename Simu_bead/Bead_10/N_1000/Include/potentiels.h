#ifndef POTENTIELS_H
#define POTENTIELS_H

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "basic_functions.h"
#include "neighborlist.h"  

// Déclaration des fonctions
double ForceLJ(double sigma, double epsilon, double x, double d);
void potentiel_barriere(double** R, double* origine, double rayon, double force_dep, double epaisseur, int N);
void force_bille_bille(double** R1,double** R2, double K_cohesine, double distance, int particule1, int particule2, int distance_cohesine_eq, double dt);
void build_neighbor_list(double **R, NeighborList *neighbor_lists, int N, int RCUT, int SKIN);
void lennard_jones_forces(double **R, NeighborList *neighbor_lists, int N, double epsilon, double sigma6, double sigma12, double Delta, int attache, int periode_enregistrement_force, FILE* fichier_force_LJ, int t);
void f_bending_forces(double **R, double **t_link, double **bending_forces, double K_bend, int N, int t);


#endif // MES_FONCTIONS_H
