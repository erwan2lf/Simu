#ifndef MOVEMENT_H
#define MOVEMENT_H
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "basic_functions.h"

// Déclaration des fonctions

void polymere_brownian_motion(double **R, double K, double Delta, int N,
                              double K_bend, double **bending_forces, int attache, int plan,
                              int t, FILE *test, int bending, int truc, int T,
                              FILE *fichier_force, int periode_enregistrement_force,
                              FILE *fichier_force_thermique, int temperature);
double** gaz_motion(double **R, int N, double ** r_new, double Delta, int plan, int attache);
void interchain_harmonic_spring(double **R, int N,
                                int *idx, int n_idx,
                                double Kc, double Delta, 
                                int t);
#endif // MES_FONCTIONS_H
