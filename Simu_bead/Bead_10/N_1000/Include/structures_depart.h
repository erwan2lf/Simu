#ifndef STRUCTURES_DEPART_H
#define STRUCTURES_DEPART_H

#include <stdio.h>
#include <time.h>

// Déclaration des fonctions
double** creation_polymere(int N, double a, double ecart);
int find_max_timestep(const char* nom_fichier) ;
double** recuperer_derniere_structure(const char* nom_fichier, int N);
void creation_polymere_solenoide(int N, double a, double ecart, double epaisseur, double** R) ;
void creation_arn_polymerase(int N_arn, double a, double rayon_arn, double** R_arn);
void creation_polymere_droit(int N, double a, double ecart, double** R);
void creation_fractal_globule(int N, double a, double ecart, double** R);
void creation_structure_knot(int N, double a, double **R);
void creation_polymere_aleatoire(int N, double a, double **R);

#endif // MES_FONCTIONS_H
