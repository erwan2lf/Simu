#include "transcription_erwan.h"
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <sys/stat.h>
#include <errno.h>
#include "structures_depart.h"
#include "Plots.h"
#include "basic_functions.h"
#include "potentiels.h"
#include "movement.h"
#include "config.h"
#include "simulation.h"
#include "neighborlist.h"   
#include "file.h"


#define BOX_SIZE 1e3
#define SCALE_POS 1000.0
#define PI 3.14159265358979323846



double*** matrix_rnap(int rows, int cols) {
    double*** matrix = (double***)malloc(rows * sizeof(double**));
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double**)malloc(cols * sizeof(double*));
        for (int j=0; j<cols ; j++) {
            matrix[i][j] = (double*)malloc(3 * sizeof(double));
            }
    }
    return matrix;
}

void free_matrix_rnap(double ***matrix, int rows, int cols) {
    if (matrix != NULL) {
        for (int i = 0; i < rows; i++) {
            if (matrix[i] != NULL) {
                for (int j = 0; j < cols; j++) {
                    if (matrix[i][j] != NULL) {
                        free(matrix[i][j]);
                        matrix[i][j] = NULL; // Prevent double free
                    }
                }
                free(matrix[i]);
                matrix[i] = NULL; // Prevent double free
            }
        }
        free(matrix);
        matrix = NULL; // Prevent double free
    }
}


double*** initialisation_matrix_rnap_NAN(int rows, int cols) {
    double*** matrix = (double***)malloc(rows * sizeof(double**));
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double**)malloc(cols * sizeof(double*));
        for (int j=0; j<cols ; j++) {
            matrix[i][j] = (double*)malloc(3 * sizeof(double));
            matrix[i][j][0] = NAN;
            matrix[i][j][1] = NAN;
            matrix[i][j][2] = NAN;

            }
    }
    return matrix; 
}

double*** resize_matrix(double*** matrix, int rows, int cols, int new_rows, int new_cols) {
    // Réallouer la mémoire pour les lignes
    matrix = (double***)realloc(matrix, new_rows * sizeof(double**));

    if (matrix == NULL) {
        fprintf(stderr, "Erreur de réallocation mémoire pour les lignes.\n");
        return NULL;
    }

    // Réinitialiser les nouvelles lignes si besoin
    for (int i = rows; i < new_rows; i++) {
        matrix[i] = (double**)malloc(new_cols * sizeof(double*));
        for (int j = 0; j < new_cols; j++) {
            matrix[i][j] = (double*)malloc(3 * sizeof(double)); // Nouvelle allocation pour les nouvelles sous-unités
        }
    }

    // Réallouer la mémoire pour les colonnes existantes si besoin
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double**)realloc(matrix[i], new_cols * sizeof(double*));
        for (int j = cols; j < new_cols; j++) {
            matrix[i][j] = (double*)malloc(3 * sizeof(double)); // Nouvelle allocation pour les nouvelles sous-unités
        }
    }
    return matrix;
}

void free_matrix_rnap_2(double*** matrix, int nb_rnap, int rnap_subunits) {
    for (int i = 0; i < nb_rnap; i++) {
        for (int j = 0; j < rnap_subunits; j++) {
            free(matrix[i][j]); // Libération de chaque sous-unité
        }
        free(matrix[i]); // Libération de chaque tableau de sous-unité
    }
    free(matrix); // Libération du tableau principal
}

void matrix_rnap_0(double*** matrix, int nb_rnap, int nb_subunit) {
    for (int rnap = 0; rnap < nb_rnap ; rnap ++){
        for (int subunit = 0 ; subunit <nb_subunit ; subunit ++){
            for (int coord = 0; coord < 3; coord++) {
                matrix[rnap][subunit][coord] = 0.0 ;
            }
        }
    }
}

void copy_matrix(double** R_rnap, double** R_rnap_new, int nb_subunit) {
        for (int subunit = 0; subunit < nb_subunit; subunit++) {
            R_rnap[subunit][0] = R_rnap_new[subunit][0];  // Copier chaque valeur
            R_rnap[subunit][1] = R_rnap_new[subunit][1];
            R_rnap[subunit][2] = R_rnap_new[subunit][2];

    }
}


int selection_promoteur(double** R, int N, int position_promoteur) {
    int new_promoteur ;
    int promoteur = position_promoteur ;
    int rnap_subunits = 8;
    double taille_pol = 5;
    double threshold = 2.0;  // Seuil de distance pour vérifier les collisions
    promoteur += 3;  // Initialement, on avance de 50
    int position_valide = 0;  // Indicateur de validité de la position
    while (!position_valide) {
        // Calcul des sous-unités de la RNAP autour du promoteur
        double x0 = R[promoteur][0], y0 = R[promoteur][1], z0 = R[promoteur][2];
        double vecteur1[3] = { R[promoteur + 2][0], R[promoteur + 2][1], R[promoteur + 2][2] };
        double vecteur2[3] = { R[promoteur - 2][0], R[promoteur - 2][1], R[promoteur - 2][2] };
        double vecteur_axe[3] = { vecteur1[0] - vecteur2[0], vecteur1[1] - vecteur2[1], vecteur1[2] - vecteur2[2] };
        double norme_axe = sqrt(vecteur_axe[0] * vecteur_axe[0] + vecteur_axe[1] * vecteur_axe[1] + vecteur_axe[2] * vecteur_axe[2]);
        double vecteur_normal[3] = { vecteur_axe[0] * cos(PI / 3), -vecteur_axe[1] * sin(PI / 3), 0 };
        vecteur_normal[2] = -(vecteur_normal[0] * vecteur_axe[0] + vecteur_normal[1] * vecteur_axe[1]) / vecteur_axe[2];
        double vecteur_tangent[3] = {
            vecteur_axe[1] * vecteur_normal[2] - vecteur_axe[2] * vecteur_normal[1],
            vecteur_axe[2] * vecteur_normal[0] - vecteur_axe[0] * vecteur_normal[2],
            vecteur_axe[0] * vecteur_normal[1] - vecteur_axe[1] * vecteur_normal[0]};
        double norme_tangent = sqrt(vecteur_tangent[0] * vecteur_tangent[0] + vecteur_tangent[1] * vecteur_tangent[1] + vecteur_tangent[2] * vecteur_tangent[2]);
        double norme_normal = sqrt(vecteur_normal[0] * vecteur_normal[0] + vecteur_normal[1] * vecteur_normal[1] + vecteur_normal[2] * vecteur_normal[2]);

        // Normalisation des vecteurs
        for (int j = 0; j < 3; j++) {
        vecteur_tangent[j] = vecteur_tangent[j] / norme_tangent * 5;
        vecteur_normal[j] = vecteur_normal[j] / norme_normal * 5;
        }

        position_valide = 1;  // On suppose que la position est valide
        // Vérification des sous-unités de la RNAP
        for (int p = 0; p < rnap_subunits; p++) {
            // Calcul de la position de la sous-unité p
            double xi = vecteur_normal[0] * cos((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2 + vecteur_tangent[0] * sin((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2;
            double yi = vecteur_normal[1] * cos((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2 + vecteur_tangent[1] * sin((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2;
            double zi = vecteur_normal[2] * cos((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2 + vecteur_tangent[2] * sin((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2;

            double position_subunit[3] = {x0 + xi, y0 + yi, z0 + zi};

            // Vérification si cette subbunit est proche d'une bille existante
            for (int i = 0; i < N; i++) {  // Parcourir les billes avant le promoteur
                double distance_bille = distance(position_subunit, R[i]);
                if (distance_bille < threshold) {
                    position_valide = 0;  // Position non valide
                    promoteur += 2;  // Avancer le promoteur de 2
                    break;  // Sortir de la boucle
                }
            }

            // Si la position est valide, on met à jour R_rnap
            if (position_valide) {
                new_promoteur = promoteur ;
            }
        }
    }
    // Enregistrer la position du promoteur pour cette RNAP
    printf("Promoteur final : %d\n", new_promoteur);
    return new_promoteur ;
}


// Ajouter une RNAP sur le promoteur
void add_rnap(double** R, int N, double ***R_rnap, int* positions_bill_rnap, int nb_rnap, int rnap_subunits, int debut_segment, int fin_segment){
    double taille_pol = 2; 
    double threshold = 1.0;  // Seuil de distance pour vérifier les collisions
    int promoteur = debut_segment; 
    int rnap = 0;
            // Calcul des sous-unités de la RNAP autour du promoteur
            double x0 = R[promoteur][0], y0 = R[promoteur][1], z0 = R[promoteur][2];
            double vecteur1[3] = {R[promoteur + 2][0], R[promoteur + 2][1], R[promoteur + 2][2]};
            double vecteur2[3] = {R[promoteur - 2][0], R[promoteur - 2][1], R[promoteur - 2][2]};
            double vecteur_axe[3] = {vecteur1[0] - vecteur2[0], vecteur1[1] - vecteur2[1], vecteur1[2] - vecteur2[2]};
            double norme_axe = sqrt(vecteur_axe[0] * vecteur_axe[0] + vecteur_axe[1] * vecteur_axe[1] + vecteur_axe[2] * vecteur_axe[2]);
            double vecteur_normal[3] = {0,-vecteur_axe[1], vecteur_axe[1]*vecteur_axe[1]/vecteur_axe[2]}; 
            
            double vecteur_tangent[3] = {
                vecteur_axe[1] * vecteur_normal[2] - vecteur_axe[2] * vecteur_normal[1],
                vecteur_axe[2] * vecteur_normal[0] - vecteur_axe[0] * vecteur_normal[2],
                vecteur_axe[0] * vecteur_normal[1] - vecteur_axe[1] * vecteur_normal[0]};
            double norme_tangent = sqrt(vecteur_tangent[0] * vecteur_tangent[0] + vecteur_tangent[1] * vecteur_tangent[1] + vecteur_tangent[2] * vecteur_tangent[2]);
            double norme_normal = sqrt(vecteur_normal[0] * vecteur_normal[0] + vecteur_normal[1] * vecteur_normal[1] + vecteur_normal[2] * vecteur_normal[2]);

            // Normalisation des vecteurs
            for (int j = 0; j < 3; j++) {
            vecteur_tangent[j] = vecteur_tangent[j] / norme_tangent*0.9;
            vecteur_normal[j] = vecteur_normal[j] / norme_normal*0.9;
            }

            // Vérification des sous-unités de la RNAP
            for (int p = 0; p < rnap_subunits; p++) {
                // Calcul de la position de la sous-unité p
                double xi = vecteur_normal[0] * cos((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2 + vecteur_tangent[0] * sin((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2;
                double yi = vecteur_normal[1] * cos((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2 + vecteur_tangent[1] * sin((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2;
                double zi = vecteur_normal[2] * cos((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2 + vecteur_tangent[2] * sin((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2;

                double position_subunit[3] = {x0 + xi, y0 + yi, z0 + zi};

                // Vérification si cette position est proche d'une bille existante
                for (int i = 0; i < N; i++) {  // Parcourir les billes avant le promoteur
                    double distance_bille = distance(position_subunit, R[i]);
                    if (distance_bille < threshold) {
                        ecarter_deux_particules(R[i],position_subunit,threshold);
                    }
                }
                // Si la position est valide, on met à jour R_rnap
                
                R_rnap[rnap][p][0] = x0 + xi;
                R_rnap[rnap][p][1] = y0 + yi;
                R_rnap[rnap][p][2] = z0 + zi;
                
            }
}



void creation_RNAP_erwan(double** R, int nb_rnap, int* positions_bille_rnap, double*** R_rnap, int N, int rnap_subunits, int debut_segment, int fin_segment){
    double taille_pol = 2; 
    double threshold = 1.0;  // Seuil de distance pour vérifier les collisions
    int promoteur = debut_segment;
    for (int rnap = 0; rnap < nb_rnap; rnap++){
        if(rnap==0){promoteur = debut_segment;}
        else{promoteur = positions_bille_rnap[rnap-1] + 2;}
        int position_valide = 0;  // Indicateur de validité de la position
       while (!position_valide) {
            // Calcul des sous-unités de la RNAP autour du promoteur
            double x0 = R[promoteur][0], y0 = R[promoteur][1], z0 = R[promoteur][2];
            double vecteur1[3] = {R[promoteur + 2][0], R[promoteur + 2][1], R[promoteur + 2][2]};
            double vecteur2[3] = {R[promoteur - 2][0], R[promoteur - 2][1], R[promoteur - 2][2]};
            double vecteur_axe[3] = {vecteur1[0] - vecteur2[0], vecteur1[1] - vecteur2[1], vecteur1[2] - vecteur2[2]};
            double norme_axe = sqrt(vecteur_axe[0] * vecteur_axe[0] + vecteur_axe[1] * vecteur_axe[1] + vecteur_axe[2] * vecteur_axe[2]);
            double vecteur_normal[3] = {0,-vecteur_axe[1], vecteur_axe[1]*vecteur_axe[1]/vecteur_axe[2]}; 
            
            double vecteur_tangent[3] = {
                vecteur_axe[1] * vecteur_normal[2] - vecteur_axe[2] * vecteur_normal[1],
                vecteur_axe[2] * vecteur_normal[0] - vecteur_axe[0] * vecteur_normal[2],
                vecteur_axe[0] * vecteur_normal[1] - vecteur_axe[1] * vecteur_normal[0]};
            double norme_tangent = sqrt(vecteur_tangent[0] * vecteur_tangent[0] + vecteur_tangent[1] * vecteur_tangent[1] + vecteur_tangent[2] * vecteur_tangent[2]);
            double norme_normal = sqrt(vecteur_normal[0] * vecteur_normal[0] + vecteur_normal[1] * vecteur_normal[1] + vecteur_normal[2] * vecteur_normal[2]);

            // Normalisation des vecteurs
            for (int j = 0; j < 3; j++) {
            vecteur_tangent[j] = vecteur_tangent[j] / norme_tangent*0.9;
            vecteur_normal[j] = vecteur_normal[j] / norme_normal*0.9;
            }

            position_valide = 1;  // On suppose que la position est valide
            // Vérification des sous-unités de la RNAP
            for (int p = 0; p < rnap_subunits; p++) {
                // Calcul de la position de la sous-unité p
                double xi = vecteur_normal[0] * cos((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2 + vecteur_tangent[0] * sin((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2;
                double yi = vecteur_normal[1] * cos((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2 + vecteur_tangent[1] * sin((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2;
                double zi = vecteur_normal[2] * cos((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2 + vecteur_tangent[2] * sin((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2;

                double position_subunit[3] = {x0 + xi, y0 + yi, z0 + zi};

                // Vérification si cette position est proche d'une bille existante
                for (int i = 0; i < N; i++) {  // Parcourir les billes avant le promoteur
                    double distance_bille = distance(position_subunit, R[i]);
                    if (distance_bille < threshold) {
                        ecarter_deux_particules(R[i],position_subunit,threshold);
                    }
                }
                // Si la position est valide, on met à jour R_rnap
                if (position_valide) {
                    R_rnap[rnap][p][0] = x0 + xi;
                    R_rnap[rnap][p][1] = y0 + yi;
                    R_rnap[rnap][p][2] = z0 + zi;
                }
            }
        }

        // Enregistrer la position du promoteur pour cette RNAP
        positions_bille_rnap[rnap] = promoteur;
        printf("Promoteur final pour la RNAP %d : %d\n", rnap, promoteur);

        // Vérification des distances pour le debug
        double dist_subunit_chrom = distance(R[promoteur], R_rnap[0][0]);
        double dist_subunit = distance(R_rnap[0][1], R_rnap[0][0]);
        printf("Distance entre deux subunits : %.2f \n", dist_subunit);
        printf("Distance entre une subunit et le monomere : %.2f \n", dist_subunit_chrom);
    }
}

void initialisation_new_RNAP(double** R, double*** R_rnap, int rnap_subunits, int promoteur, int rnap) {
    double taille_pol = 8;
    double threshold = 2.0;  // Seuil de distance pour vérifier les collisions

    // Calcul des sous-unités de la RNAP autour du promoteur

    double x0 = R[promoteur][0], y0 = R[promoteur][1], z0 = R[promoteur][2];
    double vecteur1[3] = { R[promoteur + 2][0], R[promoteur + 2][1], R[promoteur + 2][2] };
    double vecteur2[3] = { R[promoteur - 2][0], R[promoteur - 2][1], R[promoteur - 2][2] };
    double vecteur_axe[3] = { vecteur1[0] - vecteur2[0], vecteur1[1] - vecteur2[1], vecteur1[2] - vecteur2[2] };
    double norme_axe = sqrt(vecteur_axe[0] * vecteur_axe[0] + vecteur_axe[1] * vecteur_axe[1] + vecteur_axe[2] * vecteur_axe[2]);
    double vecteur_normal[3] = {0,-vecteur_axe[1], vecteur_axe[1]*vecteur_axe[1]/vecteur_axe[2]}; 
    double vecteur_tangent[3] = {
        vecteur_axe[1] * vecteur_normal[2] - vecteur_axe[2] * vecteur_normal[1],
        vecteur_axe[2] * vecteur_normal[0] - vecteur_axe[0] * vecteur_normal[2],
        vecteur_axe[0] * vecteur_normal[1] - vecteur_axe[1] * vecteur_normal[0]};
    double norme_tangent = sqrt(vecteur_tangent[0] * vecteur_tangent[0] + vecteur_tangent[1] * vecteur_tangent[1] + vecteur_tangent[2] * vecteur_tangent[2]);
    double norme_normal = sqrt(vecteur_normal[0] * vecteur_normal[0] + vecteur_normal[1] * vecteur_normal[1] + vecteur_normal[2] * vecteur_normal[2]);
    // Normalisation des vecteurs
    for (int j = 0; j < 3; j++) {
    vecteur_tangent[j] = vecteur_tangent[j] / norme_tangent * 1;
    vecteur_normal[j] = vecteur_normal[j] / norme_normal * 1;
    }
    // Vérification des sous-unités de la RNAP
    for (int p = 0; p < rnap_subunits; p++) {
        // Calcul de la position de la sous-unité p
        double xi = vecteur_normal[0] * cos((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2 + vecteur_tangent[0] * sin((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2;
        double yi = vecteur_normal[1] * cos((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2 + vecteur_tangent[1] * sin((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2;
        double zi = vecteur_normal[2] * cos((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2 + vecteur_tangent[2] * sin((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2;
        R_rnap[rnap][p][0] = x0 + xi;
        R_rnap[rnap][p][1] = y0 + yi;
        R_rnap[rnap][p][2] = z0 + zi;
    }
}


void print_position(FILE* f, int id, int type, double x, double y, double z, int num) {
    fprintf(f, "%d %d %lf %lf %lf %d\n", id, type, x / SCALE_POS, y / SCALE_POS, z / SCALE_POS, num);
}

void enregistrement_RNAP(FILE* fichier, double** R, int N, double*** R_rnap, int nb_rnap, int t,
                         int mono_transcrpt, int* positions_bille_rnap, int nb_subunits, int nb_rnap_initial, int* num_rnap, int sortie) {
    
    int total_atoms = N + nb_rnap_initial * (8 + 2);  // 8 subunits + 2 monomères par RNAP
    int count = N;
    int type = 3;

    // En-têtes LAMMPS
    fprintf(fichier, "ITEM: TIMESTEP\n%d\n", t);
    fprintf(fichier, "ITEM: NUMBER OF ATOMS\n%d\n", total_atoms);
    fprintf(fichier, "ITEM: BOX BOUNDS ss ss ss\n");
    fprintf(fichier, "%lf %lf\n", -BOX_SIZE, BOX_SIZE);
    fprintf(fichier, "%lf %lf\n", -BOX_SIZE, BOX_SIZE);
    fprintf(fichier, "%lf %lf\n", -BOX_SIZE, BOX_SIZE);
    fprintf(fichier, "ITEM: ATOMS id type xs ys zs\n");

    // 🔹 Monomères de la chaîne
    for (int i = 0; i < N; i++) {
        print_position(fichier, i, 1, R[i][0], R[i][1], R[i][2], -2);
    }

    // 🔸 RNAPs actifs
    for (int rnap = 0; rnap < nb_rnap; rnap++) {
        mono_transcrpt = positions_bille_rnap[rnap];

        // Deux monomères en type 2
        print_position(fichier, count++, 2, R[mono_transcrpt][0], R[mono_transcrpt][1], R[mono_transcrpt][2], num_rnap[rnap]);
        mono_transcrpt++;
        print_position(fichier, count++, 2, R[mono_transcrpt][0], R[mono_transcrpt][1], R[mono_transcrpt][2], num_rnap[rnap]);

        // Les 8 sous-unités
        for (int i = 0; i < nb_subunits; i++) {
            print_position(fichier, count++, type, R_rnap[rnap][i][0], R_rnap[rnap][i][1], R_rnap[rnap][i][2], num_rnap[rnap]);

            // if(sortie == 0){
            //     print_position(fichier, count++, type, R_rnap[rnap][i][0], R_rnap[rnap][i][1], R_rnap[rnap][i][2], num_rnap[rnap]);
            // }
            // else{
            //     print_position(fichier, count++, type, R_rnap[rnap][i][0], R_rnap[rnap][i][1], R_rnap[rnap][i][2], num_rnap[nb_rnap_initial-(rnap+1)]);
            // }
        }

        type++;  // type unique par RNAP
    }

    // ⚫ RNAPs inactifs (vides)
    for (int rnap = nb_rnap; rnap < nb_rnap_initial; rnap++) {
        // Deux monomères vides
        fprintf(fichier, "%d %d 0 0 0 %d\n", count++, type, -2);
        fprintf(fichier, "%d %d 0 0 0 %d\n", count++, type, -2);

        // Huit sous-unités vides (type 100 par défaut)
        for (int i = 0; i < nb_subunits; i++) {
            fprintf(fichier, "%d %d 0 0 0 %d\n", count++, 100, -1);
        }

        type++;
    }
}

void enregistrement_RNAP_dyn(FILE* fichier, double** R, int N, double*** R_rnap, int nb_rnap, int t, int* positions_bille_rnap, int position_promoteur, int premiere_rnap, int derniere_rnap, int promoteur_size) {
    int mono_transcrpt;
    double TT = 1e+3;
    fprintf(fichier, "ITEM: TIMESTEP\n");
    fprintf(fichier, "%d\n", t);
    fprintf(fichier, "ITEM: NUMBER OF ATOMS\n");
    fprintf(fichier, "%d\n", N+nb_rnap*(8+2)+promoteur_size);
    fprintf(fichier, "ITEM: BOX BOUNDS ss ss ss\n");
    fprintf(fichier, "%lf %lf\n", -TT, TT);
    fprintf(fichier, "%lf %lf\n", -TT, TT);
    fprintf(fichier, "%lf %lf\n", -TT, TT);
    fprintf(fichier, "ITEM: ATOMS id type xs ys zs\n");
    //printf("pb 222 \n ");
    for (int particle = 0; particle < N; particle++) {
        //printf("pb 3333 \n ");printf("before pb %d \n ", particle);
        fprintf(fichier, "%d 1 %lf %lf %lf\n", particle, R[particle][0] / 1000, R[particle][1] / 1000, R[particle][2] / 1000);
        //printf("pb %d \n ", particle);
                        }
    for (int bille_promoteur = 0 ; bille_promoteur<promoteur_size ; bille_promoteur++){
        fprintf(fichier, "%d 2 %lf %lf %lf\n", N+bille_promoteur, R[position_promoteur+bille_promoteur][0] / 1000, R[position_promoteur+bille_promoteur][1] / 1000, R[position_promoteur+bille_promoteur][2] / 1000);
        }


    //printf("pb stop");
    int count = N + promoteur_size;
    int type = 4 ;
    for (int rnap = 0; rnap< nb_rnap; rnap++){
        //type += rnap ;
        //mono_transcrpt+=50;
        mono_transcrpt = positions_bille_rnap[rnap];
        //printf("oui %d\n ", mono_transcrpt);
            if (rnap <derniere_rnap && rnap>= premiere_rnap) {
            fprintf(fichier, "%d %d %lf %lf %lf\n", count, 3,  R[mono_transcrpt][0] / 1000, R[mono_transcrpt][1] / 1000, R[mono_transcrpt][2] / 1000);
            mono_transcrpt+=1;count+=1;
            fprintf(fichier, "%d %d %lf %lf %lf\n", count, 3,  R[mono_transcrpt][0] / 1000, R[mono_transcrpt][1] / 1000, R[mono_transcrpt][2] / 1000);
            count+=1;
            for (int particle = 0; particle < 8; particle++) {//printf("%d %d\n", count+particle, type ); printf("%f", R_rnap[rnap][particle][0]);
                fprintf(fichier, "%d %d %lf %lf %lf\n", count+particle, type, R_rnap[rnap][particle][0] / 1000, R_rnap[rnap][particle][1] / 1000, R_rnap[rnap][particle][2] / 1000);
                }//printf("lll\n");
            count+=8; type+=1; //printf("pb 2");
                } else {
            fprintf(fichier, "%d %d %lf %lf %lf\n", count, 3, 0.0, 0.0, 0.0);
            mono_transcrpt+=1;count+=1;
            fprintf(fichier, "%d %d %lf %lf %lf\n", count, 3,   0.0, 0.0, 0.0);
            count+=1;
            for (int particle = 0; particle < 8; particle++) {//printf("%d %d\n", count+particle, type ); printf("%f", R_rnap[rnap][particle][0]);
                fprintf(fichier, "%d %d %lf %lf %lf\n", count+particle, type,  0.0, 0.0, 0.0);
                }//printf("lll\n");
            count+=8; type+=1; //printf("pb 2");


                }
        }
}

void enregistrement_RNAP_position(FILE* fichier, int nb_rnap, int t, int* positions_bille_rnap, double avancement_transcription) {
    double TT = 1e+3; int mono_transcrpt;
    fprintf(fichier, "TIMESTEP %d ", t);
    fprintf(fichier, "%f ", avancement_transcription);
    for (int rnap = 0; rnap< nb_rnap; rnap++){
        mono_transcrpt = positions_bille_rnap[rnap];
        fprintf(fichier, " %d ", mono_transcrpt);
        }
    fprintf(fichier, " \n ");
}

void enregistrement_promoteur(FILE* fichier, double** R, int N,  int nb_rnap, int t, int mono_transcrpt) {
    double TT = 1e+3;
    fprintf(fichier, "ITEM: TIMESTEP\n");
    fprintf(fichier, "%d\n", t);
    fprintf(fichier, "ITEM: NUMBER OF ATOMS\n");
    fprintf(fichier, "%d\n", N+nb_rnap);
    fprintf(fichier, "ITEM: BOX BOUNDS ss ss ss\n");
    fprintf(fichier, "%lf %lf\n", -TT, TT);
    fprintf(fichier, "%lf %lf\n", -TT, TT);
    fprintf(fichier, "%lf %lf\n", -TT, TT);
    fprintf(fichier, "ITEM: ATOMS id type xs ys zs\n");
    for (int particle = 0; particle < N; particle++) {
        fprintf(fichier, "%d 1 %lf %lf %lf\n", particle, R[particle][0] / 1000, R[particle][1] / 1000, R[particle][2] / 1000);
    }
    int count = N;
    int type = 3 ;
    for (int rnap = 0; rnap< nb_rnap; rnap++){
        fprintf(fichier, "%d %d %lf %lf %lf\n", count, 2,  R[mono_transcrpt][0] / 1000, R[mono_transcrpt][1] / 1000, R[mono_transcrpt][2] / 1000);
        count+=1;
        mono_transcrpt+=50;
        }
}

int recuperer_RNAP(const char* nom_fichier, int N, int nb_rnap, int nb_subunit, double** R, double*** R_rnap) {
    FILE* fichier = fopen(nom_fichier, "r");
    /// Verifications
    if (fichier == NULL) {
        printf("Erreur : Impossible d'ouvrir le fichier %s\n", nom_fichier);
        return 0;}
    /// Initialisation
    int t, nb_particles;
    int t_enregistrement;
    double x, y, z;
    char buffer[256]; // Pour stocker temporairement les lignes lues
    // Lire le fichier jusqu'à la fin pour récupérer la dernière structure
    /// fgets et fscanf font avancer le curseur sur la ligne suivante

    while (fgets(buffer, sizeof(buffer), fichier) != NULL) {

        // Identifier le début des données de structure
        if (strstr(buffer, "ITEM: ATOMS id type xs ys zs") != NULL) {

            // Lire le nombre de particules
            //fscanf(fichier, "%d", &nb_particles);
            // Lire les données de structure pour chaque particule
            for (int particle = 0; particle < N; particle++) {
                fscanf(fichier, "%*d %*d %lf %lf %lf", &x, &y, &z);
                //printf("%lf %lf %lf \n", x, y, z);
                R[particle][0] = x * 1000; // Convertir xs back to original scale
                R[particle][1] = y * 1000;
                R[particle][2] = z * 1000; //printf("%f ", z);
            }


            for (int rnap = 0 ; rnap<nb_rnap ; rnap ++) {
                for (int nucl_transcrit = 0 ; nucl_transcrit<3 ; nucl_transcrit++){ // 3 billes type 2 qui sont en bleu
                    fscanf(fichier, "%*d %*d %lf %lf %lf", &x, &y, &z);}

                for (int subunit = 0 ; subunit <nb_subunit ; subunit ++ ){
                    fscanf(fichier, "%*d %*d %lf %lf %lf", &x, &y, &z);
                    R_rnap[rnap][subunit][0] = x * 1000; // Convertir xs back to original scale
                    R_rnap[rnap][subunit][1] = y * 1000;
                    R_rnap[rnap][subunit][2] = z * 1000;

                    }
                }
            }
        // Lire le numéro de pas de temps

        if (strstr(buffer, "ITEM: TIMESTEP") != NULL) {
            fscanf(fichier, "%d", &t);
            t_enregistrement = t ;
        }
    }
    fclose(fichier); printf("ok ! \n");
    return t_enregistrement ;
}

double recuperer_positions_rnap(const char* nom_fichier, int nb_rnap, int** positions, double avancement_transcription_recup) {
    FILE* fichier = fopen(nom_fichier, "r");
    /// Vérification de l'ouverture du fichier
    if (fichier == NULL) {
        printf("Erreur : Impossible d'ouvrir le fichier %s\n", nom_fichier);
        return 0;
    }
    /// Initialisation
    char buffer[256]; // Pour stocker temporairement les lignes lues
    char ligne_rnap_contruite[256];
    int t; // Pour stocker le numéro du timestep
    int position ;
    double avancement ;
    /// Lire le fichier ligne par ligne
    //printf("0 \n");
    while (fgets(buffer, sizeof(buffer), fichier) != NULL) {
        // Vérifier si la ligne contient "TIMESTEP"
        if (strstr(buffer, " TIMESTEP") != NULL) {
        //printf("1 \n");
            // Lire les données de la ligne
            sscanf(buffer, " TIMESTEP %d %lf", &t, &avancement); // Extraire le numéro du timestep (facultatif)
            avancement_transcription_recup = avancement ;
            //printf("yeah \n");
            //printf("%f \n", avancement_transcription_recup);
            //printf("%f \n", *avancement_transcription_recup);
            char t_str[20];
            sprintf(t_str, "%d", t); // Convertir l'entier en chaîne
            int t_nombre_de_chiffres = strlen(t_str);

            char at_str[20];
            sprintf(at_str, "%f", avancement_transcription_recup); // Convertir l'entier en chaîne
            int at_nombre_de_chiffres = strlen(at_str)-3 ; /// -3 empirique ...

            //printf("%d %d %d %lf\n", t_nombre_de_chiffres, t, at_nombre_de_chiffres, avancement_transcription_recup);
            int count_string = 10 + t_nombre_de_chiffres + 1 + at_nombre_de_chiffres + 2;

            //sprintf(ligne_rnap_construite, "TIMESTEP %d", t);


            for (int rnap = 0; rnap < nb_rnap; rnap++) {
                //sprintf(ligne_construite, "%s %d", ligne_rnap_construite, t);
                // Utiliser sscanf pour lire les positions directement
                if (sscanf(buffer +  count_string, "%d", &position) == 1) {
                *positions[rnap] = position;
                    // La position a été lue avec succès
                    //printf("Position RNAP %d: %d\n", rnap , *positions[rnap]);
                    char rnap_str[20];
                    sprintf(rnap_str, "%d", position);
                    int rnap_nbre_chiffre = strlen(rnap_str);
                    count_string += rnap_nbre_chiffre + 2 ;

                } else {
                    printf("Erreur lors de la lecture de la position pour RNAP %d\n", rnap );
                }
            }


        }
    }
    fclose(fichier);
    printf("Recup : avancement %lf \n", avancement_transcription_recup);
    for (int rnap = 0; rnap < nb_rnap; rnap++) {
        printf("Position RNAP %d: %d\n", rnap , *positions[rnap]);
            }
    return avancement_transcription_recup ;
}
double** mvt_brownian_harmonic_bending(int p, int pmin, int pmax, double** r_new, double** R, double a, double K, double dt){
    double dmin, dmax, fmin, fmax;
    double F_alea, pot_harm[3];
    double bending_angle_target = 3.14/8  ;
    double bending_K = K*100;
    dmin = distance(R[pmin], R[p]);
    dmax = distance(R[pmax], R[p]);
    //printf("d %f ", dmin);
    double vecmin[3] = { R[pmin][0] - R[p][0], R[pmin][1] - R[p][1], R[pmin][2] - R[p][2] };
    double vecmax[3] = { R[pmax][0] - R[p][0], R[pmax][1] - R[p][1], R[pmax][2] - R[p][2] };
    double rsqmin =  dot_product(vecmin, vecmin);
    double norm_vecmin = sqrt(rsqmin);
    double rsqmax =  dot_product(vecmax, vecmax);
    double norm_vecmax = sqrt(rsqmax);
    double norm = dot_product(vecmin, vecmax) ;
    double cos_theta = norm/ (norm_vecmin * norm_vecmax);
    if (cos_theta>1.0){cos_theta = 1.0;}
    if (cos_theta<-1.0){cos_theta = -1.0;}
    //double theta = acos(cos_theta);
    double dtheta = 1/(sqrt(1-cos_theta*cos_theta));
    double bending_force =  - bending_K * (acos(cos_theta) - bending_angle_target) * dtheta ;
    //printf(" a %f ", acos(cos_theta) );
    for (int j = 0; j < 3; j++) {
        pot_harm[j] = K * ((R[pmax][j] - R[p][j]) / dmax * (dmax - a) + (R[pmin][j] - R[p][j]) / dmin * (dmin - a));
        F_alea = sqrt(2 * dt) * randn();
        r_new[p][j] += R[p][j] + F_alea + dt * pot_harm[j];
        fmin = 0;// bending_force * cos_theta / rsqmin * vecmin [j] - bending_force * 1 / (norm_vecmin*norm_vecmax) * vecmax [j] ;
        fmax = 0;//bending_force * cos_theta / rsqmax * vecmax [j] - bending_force * 1 / (norm_vecmin*norm_vecmax) * vecmin [j] ;
        r_new[p][j] += -(fmin + fmax)* dt;
        r_new[pmin][j] += fmin * dt ;
        r_new[pmax][j] += fmax * dt;
        }
        //printf(" d %.2f a %.2f a0 %.2f",dmin, acos(cos_theta), bending_angle_target);
    return r_new ;
}
double** polymere_brownian_motion_ring_v2(double** R, double a, double K, double dt, int N, double** r_new) {
    int p, pmin, pmax ;
    for (int i = 1; i < N - 1; i++) {
        p = i ; pmin= i - 1 ; pmax = i + 1 ;
        r_new = mvt_brownian_harmonic_bending(p, pmin, pmax, r_new, R, a, K, dt);
        }
    p = 0, pmax = 1, pmin = N-1 ;
    r_new = mvt_brownian_harmonic_bending(p, pmin, pmax, r_new, R, a, K, dt);
    p = N - 1; pmin = N - 2, pmax = 0;
    r_new = mvt_brownian_harmonic_bending(p, pmin, pmax, r_new, R, a, K, dt);
    return r_new ;
    }

void mvt_brownian_harmonic_force(int p, int pmin, int pmax, double** r_new, double** R, double a, double K, double Delta, FILE *test, int t){
    double dmin, dmax;
    double F_alea, pot_harm[3];
    dmin = distance(R[pmin], R[p]);
    dmax = distance(R[pmax], R[p]);
 
    for (int j = 0; j < 3; j++) {
        pot_harm[j] = K * ((R[pmax][j] - R[p][j]) * (1- 1/dmax) + (R[pmin][j] - R[p][j]) * (1 - 1/dmin));
        F_alea = sqrt(2 * Delta) * randn();
        r_new[p][j] += R[p][j] + F_alea + Delta * pot_harm[j];
        }
}

void mvt_brownian_harmonic_force_RNAP(int p, int pmin, int pmax, int rnap, double** r_new, double** R, double alpha, double K, double Delta, FILE *test, int t, int periode_enregistrement_force, FILE* fichier_force_rnap, FILE* fichier_force_thermique, int temperature){
    double dmin, dmax;
    double F_alea, pot_harm[3];
    dmin = distance(R[pmin], R[p]);
    dmax = distance(R[pmax], R[p]);
    double deplacement = 0.0;

    if(t % periode_enregistrement_force == 0)
    {
        fprintf(fichier_force_rnap, "%d %d ", rnap, p);
        fprintf(fichier_force_thermique, "%d %d ", rnap, p);
    }
    for(int j = 0; j < 3; j ++)
    {
        pot_harm[j] = K * ((R[pmax][j] - R[p][j]) * (1- alpha/dmax) + (R[pmin][j] - R[p][j]) * (1 - alpha/dmin));
        if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap, "%f ",pot_harm[j]);}
        deplacement = - r_new[p][j];
        r_new[p][j] += R[p][j] + Delta * pot_harm[j];
        deplacement += r_new[p][j]; 
        if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap, "%f ", deplacement);}
    }
    for (int j = 0; j < 3; j++)
    {
        if(temperature == 1)
        {
            F_alea = sqrt(2 * Delta) * randn()/1; 
        }
        else
        {
            F_alea = 0;
        }
        if(t % periode_enregistrement_force == 0){fprintf(fichier_force_thermique, "%f ", F_alea);}
        // F_alea = 0;
        deplacement = - r_new[p][j];
        r_new[p][j] += F_alea;
        deplacement += r_new[p][j]; 
        if(t % periode_enregistrement_force == 0){fprintf(fichier_force_thermique, "%f ", deplacement);}
    }
     if(t % periode_enregistrement_force == 0){fprintf(fichier_force_thermique, "\n");}
     if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap, "\n");}
}

void polymere_brownian_motion_ring_force(double** R_rnap, double alpha, double K_rnap, double Delta, int N, int rnap, double** r_new, FILE* test, int t,  int periode_enregistrement_force, FILE* fichier_force_rnap, FILE* fichier_force_thermique, int temperature) {

    int p, pmin, pmax;

    p = 0; pmax = 1; pmin = N-1;
    mvt_brownian_harmonic_force_RNAP(p, pmin, pmax, rnap, r_new, R_rnap, alpha, K_rnap, Delta, test, t, periode_enregistrement_force, fichier_force_rnap, fichier_force_thermique, temperature);

    for (int i = 1; i < N - 1; i++) {
        p = i ; pmin= i - 1 ; pmax = i + 1;
        mvt_brownian_harmonic_force_RNAP(p, pmin, pmax, rnap, r_new, R_rnap, alpha, K_rnap, Delta, test, t, periode_enregistrement_force, fichier_force_rnap, fichier_force_thermique, temperature);
        }

    p = N - 1; pmin = N - 2; pmax = 0;
    mvt_brownian_harmonic_force_RNAP(p, pmin, pmax, rnap, r_new, R_rnap, alpha, K_rnap, Delta, test, t, periode_enregistrement_force, fichier_force_rnap, fichier_force_thermique, temperature); 

}

double** polymere_brownian_motion_ring(double** R, double a, double K, double dt, int N, double** r_new) {
    double d1, d0, di, dk, da, db, f1, f2;
    double F_alea, pot_harm[3];
    double bending_angle_target = 2*PI/8  ;
    double bending_K = K*1000;
    for (int i = 1; i < N - 1; i++) {
        da = distance(R[i - 1], R[i]);
        db = distance(R[i + 1], R[i]);
        printf("d %f ", da);
        double veca[3] = { R[i - 1][0] - R[i][0], R[i - 1][1] - R[i][1], R[i - 1][2] - R[i][2] };
        double vecb[3] = { R[i + 1][0] - R[i][0], R[i + 1][1] - R[i][1], R[i + 1][2] - R[i][2] };
        double rsqa =  dot_product(veca, veca);
        double norm_veca = sqrt(rsqa);
        double rsqb =  dot_product(vecb, vecb);
        double norm_vecb = sqrt(rsqb);
        double norm = dot_product(veca, vecb) ;
        double cos_theta = norm/ (norm_veca * norm_vecb);
        if (cos_theta>1.0){cos_theta = 1.0;}
        if (cos_theta<-1.0){cos_theta = -1.0;}
        //double theta = acos(cos_theta);
        double dtheta = 1/(sqrt(1-cos_theta*cos_theta));
        double bending_force =  - bending_K * (acos(cos_theta) - bending_angle_target) * dtheta ;
        printf(" a %f ", acos(cos_theta) );
        for (int j = 0; j < 3; j++) {
            pot_harm[j] = K * ((R[i + 1][j] - R[i][j]) / db * (db - a) + (R[i - 1][j] - R[i][j]) / da * (da - a));
            F_alea = 0;// sqrt(2 * dt) * randn();
            r_new[i][j] += R[i][j] + F_alea + dt * pot_harm[j];
            f1 = bending_force * cos_theta / rsqa * veca [j] - bending_force * 1 / (norm_veca*norm_vecb) * vecb [j] ;
            f2 = bending_force * cos_theta / rsqb * vecb [j] - bending_force * 1 / (norm_veca*norm_vecb) * veca [j] ;
            r_new[i][j] += -(f1 + f2)* dt;
            r_new[i-1][j] += f1 * dt ;
            r_new[i+1][j] += f2 * dt;
                                    }
                                  }
                                  printf(" d  %f", da);
    //a=a/4;
    int p = 0, pmax = 1, pmin = N-1 ; // p3 = pmax
    di = distance(R[pmin], R[p]);
    dk = distance(R[pmax], R[p]);
    printf("d %f d %f \n", di, dk);
    double veca[3] = { R[pmax][0] - R[p][0], R[pmax][1] - R[p][1], R[pmax][2] - R[p][2] };
    double vecb[3] = { R[pmin][0] - R[p][0], R[pmin][1] - R[p][1], R[pmin][2] - R[p][2] };
    double rsqa =  dot_product(veca, veca);
    double norm_veca = sqrt(rsqa);
    double rsqb =  dot_product(vecb, vecb);
    double norm_vecb = sqrt(rsqb);
    double norm = dot_product(veca, vecb) ;
    double cos_theta = norm/ (norm_veca * norm_vecb);
    if (cos_theta>1.0){cos_theta = 1.0;}
    if (cos_theta<-1.0){cos_theta = -1.0;}
    //double theta = acos(cos_theta);
    double dtheta = 1/(sqrt(1-cos_theta*cos_theta));
    double bending_force =  - bending_K * (acos(cos_theta) - bending_angle_target) * dtheta ;
    for (int j = 0; j < 3; j++) {
        pot_harm[j] = K * (R[pmin][j] - R[p][j]) / di * (di - a) + K * (R[pmax][j] - R[p][j]) / dk * (dk - a);
        r_new[p][j] += R[p][j] + dt * pot_harm[j];// + sqrt(2 * dt) * randn();
        f1 = bending_force * cos_theta / rsqa * veca [j] - bending_force * 1 / (norm_veca*norm_vecb) * vecb [j] ;
        f2 = bending_force * cos_theta / rsqb * vecb [j] - bending_force * 1 / (norm_veca*norm_vecb) * veca [j] ;
        r_new[p][j] += -(f1 + f2)* dt;
        r_new[pmin][j] += f1 * dt ;
        r_new[pmax][j] += f2 * dt;   }
        printf(" a %f a0 %f \n ", acos(cos_theta), bending_angle_target);


    p = N - 1; pmin = N - 2, pmax = 0;
    di = distance(R[pmin], R[p]);
    dk = distance(R[pmax], R[p]);
    veca[0] = R[pmax][0] - R[p][0];
    veca[1] = R[pmax][1] - R[p][1];
    veca[2] = R[pmax][2] - R[p][2];
    vecb[0] = R[pmin][0] - R[p][0];
    vecb[1] = R[pmin][1] - R[p][1];
    vecb[2] = R[pmin][2] - R[p][2];
    rsqa =  dot_product(veca, veca);
    norm_veca = sqrt(rsqa);
    rsqb =  dot_product(vecb, vecb);
    norm_vecb = sqrt(rsqb);
    norm = dot_product(veca, vecb) ;
    cos_theta = norm/ (norm_veca * norm_vecb);
    if (cos_theta>1.0){cos_theta = 1.0;}
    if (cos_theta<-1.0){cos_theta = -1.0;}
    //double theta = acos(cos_theta);
     dtheta = 1/(sqrt(1-cos_theta*cos_theta));
     bending_force =  - bending_K * (acos(cos_theta) - bending_angle_target) * dtheta ;
    for (int j = 0; j < 3; j++) {
        //pot_harm[j] = K * ((R[i + 1][j] - R[i][j]) / db * (db - a) + (R[i - 1][j] - R[i][j]) / da * (da - a));
        pot_harm[j] = K * ((R[pmin ][j] - R[p][j]) / di * (di - a) + (R[pmax ][j] - R[p][j]) / dk * (dk - a));
        r_new[p][j] += R[p][j] + dt * pot_harm[j];// + sqrt(2 * dt) * randn();
        f1 = bending_force * cos_theta / rsqa * veca [j] - bending_force * 1 / (norm_veca*norm_vecb) * vecb [j] ;
        f2 = bending_force * cos_theta / rsqb * vecb [j] - bending_force * 1 / (norm_veca*norm_vecb) * veca [j] ;
        r_new[p][j] += -(f1 + f2)* dt;
        r_new[pmin][j] += f1 * dt ;
        r_new[pmax][j] += f2 * dt;
                                }
    return r_new ;
    }

double** bond_rnap_bead(double** R, double** R_rnap, double a_transpt, double K_transpt, double dt, double** r_new, int mono_transcrpt) {
    for (int subunit=0; subunit<8; subunit+=1){
        double dist = distance(R[mono_transcrpt], R_rnap[subunit]);
        for (int j = 0; j < 3; j++) {
             R_rnap[subunit][j] += dt * K_transpt * (R[mono_transcrpt][j] - R_rnap[subunit][j]) / dist * (dist - a_transpt);
             R[mono_transcrpt][j] += - dt * K_transpt * (R[mono_transcrpt][j] - R_rnap[subunit][j]) / dist * (dist - a_transpt);
                  }}
    return r_new ;
    }

// Version pondérée + lissée (smoothstep) pour continuité de vitesse
// Action/réaction symétrique entre RNAP (r_new) et les monomères (R)
double** bond_rnap_bead_progressive_mvt(
    double** R,                 // positions de la chromatine
    double** R_rnap,            // positions des sous-unités RNAP [8][3]
    double a_transpt,           // longueur à l'équilibre (si utile dans votre modèle)
    double K_transpt,           // raideur du couplage
    double Delta,               // pas de temps
    double** r_new,             // positions RNAP next (à mettre à jour)
    int mono_transcrpt,         // index du monomère courant p
    double dx_avancement_rnap,  // progression dans [0,1] entre p et p+1
    double a,                   // (non utilisé ici, laissé pour compat)
    double alpha,               // paramètre de votre facteur (1 - (alpha+1)/(2*dist))
    FILE* fichier_force_lea,    // log optionnel
    int periode_enregistrement_force,
    int t
) {
    const int NB_SUBUNITS = 8;
    const double eps = 1e-12;   // évite division par 0
    (void)a; // pas utilisé ici, mais conservé pour compat

    // Progression lissée pour éviter les discontinuités de vitesse
    // smoothstep: s -> s^2(3-2s)
    double s = dx_avancement_rnap;
    if (s < 0.0) s = 0.0;
    if (s > 1.0) s = 1.0;
    double s_eff = s * s * (3.0 - 2.0 * s);

    // Poids continus et lissés entre mono p et p+1
    double w0 = 1.0 - s_eff;   // poids appliqué au monomère p
    double w1 = s_eff;         // poids appliqué au monomère p+1

    // indices sécurité (assume mono_transcrpt+1 valide)
    int p   = mono_transcrpt;
    int p1  = mono_transcrpt + 1;

    if (t % periode_enregistrement_force == 0 && fichier_force_lea) {
        fprintf(fichier_force_lea, "t=%d RNAP-bond progressive (p=%d → p+1)\n", t, p);
        fprintf(fichier_force_lea, "weights: w0=%.6g w1=%.6g (s=%.6g s_eff=%.6g)\n", w0, w1, dx_avancement_rnap, s_eff);
    }

    for (int sub = 0; sub < NB_SUBUNITS; ++sub) {
        // Vecteurs vers p et p+1
        double d0x = R[p][0]  - R_rnap[sub][0];
        double d0y = R[p][1]  - R_rnap[sub][1];
        double d0z = R[p][2]  - R_rnap[sub][2];

        double d1x = R[p1][0] - R_rnap[sub][0];
        double d1y = R[p1][1] - R_rnap[sub][1];
        double d1z = R[p1][2] - R_rnap[sub][2];

        // Distances
        double dist0 = sqrt(d0x*d0x + d0y*d0y + d0z*d0z); if (dist0 < eps) dist0 = eps;
        double dist1 = sqrt(d1x*d1x + d1y*d1y + d1z*d1z); if (dist1 < eps) dist1 = eps;

        // Facteur "nucléosome" utilisé dans votre code: (1 - (alpha+1)/(2*dist))
        // Si vous préférez un ressort harmonique linéaire pur, remplacer f0/f1 par (dist - a_transpt)/dist
        double f0 = (1.0 - (alpha + 1.0) / (2.0 * dist0));
        double f1 = (1.0 - (alpha + 1.0) / (2.0 * dist1));

        // Force totale exercée par la chromatine sur la sous-unité (pondérée p/p+1)
        // F_subunit = K * [ w0 * d0 * f0 + w1 * d1 * f1 ]
        double Fx = K_transpt * ( w0 * d0x * f0 + w1 * d1x * f1 );
        double Fy = K_transpt * ( w0 * d0y * f0 + w1 * d1y * f1 );
        double Fz = K_transpt * ( w0 * d0z * f0 + w1 * d1z * f1 );

        // Mise à jour RNAP (r_new reçoit +F * Delta)
        r_new[sub][0] += Delta * Fx;
        r_new[sub][1] += Delta * Fy;
        r_new[sub][2] += Delta * Fz;

        // Réaction sur les monomères (opposée et pondérée)
        // p reçoit -w0 * F ; p+1 reçoit -w1 * F
        R[p][0]  -= Delta * w0 * Fx;
        R[p][1]  -= Delta * w0 * Fy;
        R[p][2]  -= Delta * w0 * Fz;

        R[p1][0] -= Delta * w1 * Fx;
        R[p1][1] -= Delta * w1 * Fy;
        R[p1][2] -= Delta * w1 * Fz;

        if (t % periode_enregistrement_force == 0 && fichier_force_lea) {
            fprintf(fichier_force_lea,
                " subunit=%d  dist0=%.4g dist1=%.4g  F=(%.4g,%.4g,%.4g)\n",
                sub, dist0, dist1, Fx, Fy, Fz
            );
        }
    }

    if (t % periode_enregistrement_force == 0 && fichier_force_lea) {
        fprintf(fichier_force_lea, "\n");
        fflush(fichier_force_lea);
    }

    return r_new;
}

void liaison_sup(double alpha, double K_transpt, double Delta, double** r_new, int p1, int p2, FILE* fichier_force_rnap_2, int t, int periode_enregistrement_force){ 
        double dist = distance(r_new[p2], r_new[p1]);
        double deplacement = 0.0; 
        if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap_2, "%d ",p1);}

        for(int j = 0; j < 3; j ++ ){
            if(t % periode_enregistrement_force == 0){
                fprintf(fichier_force_rnap_2, "%lf ", K_transpt * (r_new[p2][j] - r_new[p1][j]) * (1 - alpha/dist));
            }
            deplacement = - r_new[p1][j];
            r_new[p1][j] +=   Delta * K_transpt * (r_new[p2][j] - r_new[p1][j]) * (1 - alpha/dist);
            deplacement += r_new[p1][j];
            if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap_2,"%lf ", deplacement);}
            
        }
        if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap_2, "\n");}
        if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap_2, "%d ",p2);}


        for (int j = 0; j < 3; j++) {
            if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap_2, "%lf ", - K_transpt * (r_new[p2][j] - r_new[p1][j]) * (1 - alpha/dist));}
            deplacement = - r_new[p2][j];
            r_new[p2][j] += - Delta * K_transpt * (r_new[p2][j] - r_new[p1][j]) * (1 - alpha/dist);
            deplacement += r_new[p2][j];
            if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap_2, "%lf ", deplacement);}
        }
        if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap_2, "\n");}
}




void lennard_jones_forces_rnap(double ***R_rnap, int nb_rnap, double **R, int N, NeighborList_rnap **neighbor_lists, double epsilon, double sigma6, double sigma12, double sigma6rnap2, double sigma12rnap2, double cut_rnap2, double Delta, int t, FILE* test, int T, FILE* fichier_force_rnap, int periode_enregistrement_force){
    if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap,"TIMESTEP: %d\n",t);}
    double deplacement = 0.0;
    for (int rnap = 0; rnap < nb_rnap; rnap++) {

        if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap,"RNAP: %d\n",rnap);}

        for (int subunit = 0; subunit<8 ; subunit++){

            if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap,"SUBUNIT: %d\n",subunit);}

            for (int k = 0; k < neighbor_lists[rnap][subunit].count; k++) {

                int j = neighbor_lists[rnap][subunit].neighbors[k];
                if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap,"NEIGHBOR: %d ",j);}
                double *Ri = R_rnap[rnap][subunit], *Rj = R[j];
                
                double dx = Ri[0] - Rj[0];
                double dy = Ri[1] - Rj[1];
                double dz = Ri[2] - Rj[2];

                double r2 = dx * dx + dy * dy + dz * dz;


                double r8 = r2 * r2 * r2 * r2; 
                double r14 = r8 * r2 * r2 * r2; 
                // double f = 4 * epsilon * (12 * sigma12/r14 - 6 * sigma6 / r8);
                double f = 4  * (12 * epsilon * sigma12/r14);

                double e = 300;
                if (f > e) {f = e;}


                if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap,"%f ", f * dx);}
                deplacement = - Ri[0];
                Ri[0] += Delta * f * dx;
                deplacement += Ri[0];
                if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap,"%f ", deplacement);}
                

                if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap,"%f ", f * dy);}
                deplacement = - Ri[1];
                Ri[1] += Delta * f * dy;
                deplacement += Ri[1];
                if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap,"%f ", deplacement);}
               
                if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap,"%f ", f * dz);}
                deplacement = - Ri[2];
                Ri[2] += Delta * f * dz;
                deplacement += Ri[2];
                if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap,"%f\n", deplacement);}
                
                Rj[0] -= Delta * f * dx;
                Rj[1] -= Delta * f * dy;
                Rj[2] -= Delta * f * dz;
            }
            for (int subunit2 = 0; subunit2<8 ; subunit2++) {
                double epsilon_test_att = epsilon; 
                double epsilon_test_rep = 1.5 * epsilon;
                if (subunit == subunit2) {continue;}
                double *Ri = R_rnap[rnap][subunit], *Rj = R_rnap[rnap][subunit2];
                double dx = Ri[0] - Rj[0];
                double dy = Ri[1] - Rj[1];
                double dz = Ri[2] - Rj[2];
                double r2 = dx * dx + dy * dy + dz * dz;
                if (sqrt(r2)>cut_rnap2) {continue;}
                

                double r8 = r2 * r2 * r2 * r2; 
                double r14 = r8 * r2 * r2 * r2; 
                double f = 4  * (12 * epsilon_test_rep * sigma12/r14 - 6 * epsilon_test_att * sigma6 / r8);
                

                double e = 300;
                if (f > e) {f = e;}

                Ri[0] += Delta * f * dx;
                Ri[1] += Delta * f * dy;
                Ri[2] += Delta * f * dz;
                Rj[0] -= Delta * f * dx;
                Rj[1] -= Delta * f * dy;
                Rj[2] -= Delta * f * dz;
            }
        } 
    }
}

int actualisation_position_bille_rnap(int positions_bille_rnap, double** R_rnap, double** R, int N, int N_pol, int* rattacher) {
    int bille_testee, dmin;
    int bille_supposee = positions_bille_rnap;
    double barycentre_rnap[3] = {0, 0, 0};
    double d ;

    for (int p =0 ; p<N_pol ; p++){
        barycentre_rnap[0]+= R_rnap[p][0] ;
        barycentre_rnap[1]+= R_rnap[p][1] ;
        barycentre_rnap[2]+= R_rnap[p][2] ;
        }
    barycentre_rnap[0]= barycentre_rnap[0]/N_pol;
    barycentre_rnap[1]= barycentre_rnap[1]/N_pol;
    barycentre_rnap[2]= barycentre_rnap[2]/N_pol;

    dmin = fabs(distance(barycentre_rnap, R[positions_bille_rnap]));
    bille_testee = positions_bille_rnap-1 ;//printf("1 %d \n", bille_testee);
    d = fabs(distance(barycentre_rnap, R[bille_testee]));
    if (d<dmin){bille_supposee=bille_testee;dmin = d;}


    bille_testee = positions_bille_rnap-2 ;
    d = fabs(distance(barycentre_rnap, R[bille_testee]));
    if (d<dmin){bille_supposee=bille_testee;dmin = d;}


    for (int p = positions_bille_rnap+1 ; p<positions_bille_rnap+2 ; p++){
        d = fabs(distance(barycentre_rnap, R[p]));
        if (d<dmin){bille_supposee=p;dmin = d;}
        }

    if (dmin>30) { printf("sauté yes \n");
    bille_supposee = positions_bille_rnap ;
    * rattacher = 1 ;
     }
    return bille_supposee;
}



void LJ_ouverture_accessibilite(double **R, NeighborList *neighbor_lists, double epsilon, double sigma6, double sigma12, int* positions_bille_rnap, int nb_rnap) {
    int rnap ;
    for (int i = 0; i < nb_rnap; i++) {
        rnap = positions_bille_rnap[i];
        for (int k = 0; k < neighbor_lists[rnap].count; k++) {
            int j = neighbor_lists[rnap].neighbors[k];
            double dx = R[rnap][0] - R[j][0];
            double dy = R[rnap][1] - R[j][1];
            double dz = R[rnap][2] - R[j][2];
            double r2 = dx * dx + dy * dy + dz * dz;
            //double r2 = distance(R[i], R[j])
            double r8 = r2 * r2 * r2 * r2;
            double r14 = r8 * r2 * r2* r2;
            //double y = 4 * epsilon * (-12 * pow(sigma, 12) * x / pow(d, 14) + 6 * pow(sigma, 6) * x / pow(d, 8))
            double f = 4 * epsilon * (12 * sigma12 / r14 - 6*sigma6 / r8);//*x /////////////*24
            //printf("%f %f \n", f, sqrt(r2));
            double e = 10;
            if (f > e) {f = e;}
            if (f < -e) {f = -e;}
            R[rnap][0] += f * dx;
            R[rnap][1] += f * dy;
            R[rnap][2] += f * dz;
            R[j][0] -= f * dx;
            R[j][1] -= f * dy;
            R[j][2] -= f * dz;
        }
    }
}



void rattacher_RNAP(double** R,  int nb_rnap, double** r_rnap, int rnap_subunits, int rnap_position){
    double taille_pol = 15;
    double x0 = R[rnap_position][0], y0 = R[rnap_position][1], z0 = R[rnap_position][2];
    double vecteur1[3] = { R[rnap_position + 2][0], R[rnap_position + 2][1], R[rnap_position + 2][2] };
    double vecteur2[3] = { R[rnap_position - 2][0], R[rnap_position - 2][1], R[rnap_position - 2][2] };
    double vecteur_axe[3] = { vecteur1[0] - vecteur2[0], vecteur1[1] - vecteur2[1], vecteur1[2] - vecteur2[2] };
    double norme_axe = sqrt(vecteur_axe[0] * vecteur_axe[0] + vecteur_axe[1] * vecteur_axe[1] + vecteur_axe[2] * vecteur_axe[2]);
    double vecteur_normal[3] = { vecteur_axe[0] * cos(PI / 3), -vecteur_axe[1] * sin(PI / 3), 0 };
    vecteur_normal[2] = -(vecteur_normal[0] * vecteur_axe[0] + vecteur_normal[1] * vecteur_axe[1]) / vecteur_axe[2];
    double vecteur_tangent[3] = {
        vecteur_axe[1] * vecteur_normal[2] - vecteur_axe[2] * vecteur_normal[1],
        vecteur_axe[2] * vecteur_normal[0] - vecteur_axe[0] * vecteur_normal[2],
        vecteur_axe[0] * vecteur_normal[1] - vecteur_axe[1] * vecteur_normal[0]};
    double norme_tangent = sqrt(vecteur_tangent[0] * vecteur_tangent[0] + vecteur_tangent[1] * vecteur_tangent[1] + vecteur_tangent[2] * vecteur_tangent[2]);
    double norme_normal = sqrt(vecteur_normal[0] * vecteur_normal[0] + vecteur_normal[1] * vecteur_normal[1] + vecteur_normal[2] * vecteur_normal[2]);

    // Normalisation des vecteurs
    for (int j = 0; j < 3; j++) {
    vecteur_tangent[j] = vecteur_tangent[j] / norme_tangent * 5;
    vecteur_normal[j] = vecteur_normal[j] / norme_normal * 5;
    }
    for (int p = 0; p < rnap_subunits; p++) {
        // Calcul de la position de la sous-unité p
        double xi = vecteur_normal[0] * cos((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2 + vecteur_tangent[0] * sin((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2;
        double yi = vecteur_normal[1] * cos((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2 + vecteur_tangent[1] * sin((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2;
        double zi = vecteur_normal[2] * cos((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2 + vecteur_tangent[2] * sin((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2;

        r_rnap[p][0] = x0 + xi;
        r_rnap[p][1] = y0 + yi;
        r_rnap[p][2] = z0 + zi;
        }
}

double transcrire_nucleotides(double a, double avancement_transcription, int* positions_bille_rnap, int rnap) {
    // Vérifie si l'avancement de la transcription dépasse le seuil
    if (a - avancement_transcription < 0.00001) {
        // Affiche un message de transcription
        printf("+1 nucle transcrit rnap %d\n", rnap);

        // Affiche la position initiale
        printf("init position %ls\n", positions_bille_rnap);

        // Incrémente la position de la bille
        positions_bille_rnap[rnap] += 1;

        // Réinitialise l'avancement de la transcription
        avancement_transcription = 0;

        // Affiche la nouvelle position après incrémentation
        printf("then position %ls\n", positions_bille_rnap);
    }
    return avancement_transcription ;
}

void afficher_structure_rnap(double*** R_matrix, int nb_rnap, int nb_subunit) {
    for (int i = 0; i < nb_rnap; i++) {
    for (int subunit =0 ; subunit <nb_subunit ; subunit ++){
            printf("Particule %d %d : x = %lf, y = %lf, z = %lf\n", i, subunit,  R_matrix[i][subunit][0], R_matrix[i][subunit][1], R_matrix[i][subunit][2]);
    } }
}



void promoteur_ouvert(double** R, int N, int promoteur, int promoteur_size, double promoteur_rayon_ouverture) {
    double d;
    for (int bille_prom = promoteur+1 ; bille_prom < promoteur + promoteur_size -1 ;  bille_prom ++){

        for (int bille = 0 ; bille < promoteur; bille ++ ){
            d = distance (R[bille_prom], R[bille]);
            if (d < promoteur_rayon_ouverture ) {
            //printf("yeah %d\n", bille);
                R[bille][0] += (R[bille][0]- R[bille_prom][0])/d*(promoteur_rayon_ouverture-d) ;
                R[bille][1] += (R[bille][1]- R[bille_prom][1])/d*(promoteur_rayon_ouverture-d) ;
                R[bille][2] += (R[bille][2]- R[bille_prom][2])/d*(promoteur_rayon_ouverture-d) ;
            }

         }
        for (int bille = bille_prom + promoteur_size + 1 ; bille < N-1; bille ++ ){
            d = distance (R[bille_prom], R[bille]);
            if (d < promoteur_rayon_ouverture  ) {

                R[bille][0] += (R[bille][0]- R[bille_prom][0])/d*(promoteur_rayon_ouverture-d) ;
                R[bille][1] += (R[bille][1]- R[bille_prom][1])/d*(promoteur_rayon_ouverture-d) ;
                R[bille][2] += (R[bille][2]- R[bille_prom][2])/d*(promoteur_rayon_ouverture-d) ;
            }

         }
    }
}

void simu_LJ_RNAP_erwan(const Config *cfg, SimVars *sv, int nbr_simu, const Files *f){


    clock_t start_2, end_2; 
    clock_t start, end; double duree_boucle ; double duree_tot = 0 ; double temps_restant ;

    double mesures[2];

    int taille_domaine = 50 ; int close_part;

    int t = 0;
    double **R = sv->R; 
    double **R_new = sv->R_new;

    ///////////////// Neighborlist ////////////////////
    NeighborList_rnap** neighbor_lists_rnap = allocate_neighbor_list_rnap(sv->nb_rnap, cfg->rnap_subunits);
    NeighborList *neighbor_lists = malloc(cfg->N * sizeof(NeighborList));  // Allocation mémoire pour les listes de voisins

    for (int i = 0; i < cfg->N; i++) {
        neighbor_lists[i].neighbors = malloc(10 * sizeof(int));
        neighbor_lists[i].capacity = 10;
        }

    printf("nb_rnap_initial = %d ",cfg->nb_rnap_initial);

    // creation_RNAP_erwan(R, nb_rnap, positions_bille_rnap, R_rnap, N, rnap_subunits, debut_segment, fin_segment);

    ///////////////// Fichier et enregistrement  ////////////////
    if(cfg->nb_rnap_initial > 0){
        enregistrement_RNAP(f->fichier, R, cfg->N, sv->R_rnap, sv->nb_rnap, t, cfg->mono_transcrpt, sv->positions_bille_rnap, cfg->rnap_subunits, cfg->nb_rnap_initial,sv->num_rnap, sv->sortie);
        enregistrement_RNAP_position(f->fichier_rnap, sv->nb_rnap, t, sv->positions_bille_rnap, sv->avancement_transcription[0]);
    }
    
    

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////// Mise à l'équilibre////////// //////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    f_equilibriate(sv, cfg, f, neighbor_lists, neighbor_lists_rnap, nbr_simu);
    sv->num_rnap[0] = 0;
    creation_1_rnap_erwan(R, sv->nb_rnap, sv->positions_bille_rnap, sv->R_rnap, cfg->N, cfg->rnap_subunits, cfg->debut_segment, cfg->fin_segment, cfg->nb_rnap_initial);

    
    double endtoend = 0.0;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////// Vrai Boucle ////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
    sv->premier_passage = 1;

    calcul(sv, cfg, f, neighbor_lists, neighbor_lists_rnap, nbr_simu, endtoend);

    endtoend = endtoend / cfg->T_endtoend;
    finalize_simulationaaaaa(sv, cfg, f, endtoend, nbr_simu);
    for (int i = 0; i < cfg->N; i++) {free(neighbor_lists[i].neighbors);} free(neighbor_lists);

    return;
}


void gonflement(double** R, int N, int *critere, double ecart_monomere){
    FILE* fichier = fopen("glonflement.lammpstrj", "w");

    if( fichier == NULL){
        perror("Erreur lors de l'ouverture de fichier gonflement");
        exit(EXIT_FAILURE);
    }

    while(*critere != 1){
        int count = 0;
        double vecteur[3];
        double vecteur_normalise[3];
        double delta;
        for(int i = 0; i < N; i++){
            for(int j = i +2; j < N; j++){
                if(distance(R[i],R[j])<ecart_monomere){
                    count +=1;
                    vecteur[0] = R[j][0] - R[i][0];
                    vecteur[1] = R[j][1] - R[i][1];
                    vecteur[2] = R[j][2] - R[i][2];

                    vecteur_normalise[0] = vecteur[0] / distance(R[i],R[j]);
                    vecteur_normalise[1] = vecteur[1] / distance(R[i],R[j]);
                    vecteur_normalise[2] = vecteur[2] / distance(R[i],R[j]);

                    delta = ecart_monomere - distance(R[i],R[j]);

                    R[i][0] -= delta * vecteur_normalise[0];
                    R[i][1] -= delta * vecteur_normalise[1];
                    R[i][2] -= delta * vecteur_normalise[2];

                    R[j][0] += delta * vecteur_normalise[0];
                    R[j][1] += delta * vecteur_normalise[1];
                    R[j][2] += delta * vecteur_normalise[2];
                }
            }
        }
        if(count<10){*critere=1;}
    }
    enregistrement(fichier,R,N,0);
    printf("critere = %d \n",*critere);
   
}

void ecarter_deux_particules(double* vec1, double* vec2, double ecart_particule){

    double vecteur[3] = {vec2[0]-vec1[0], vec2[1]-vec1[1], vec2[2]-vec1[2]};
    double norme = distance(vec1,vec2);
    double vecteur_normalise[3] = {vecteur[0]/norme, vecteur[1]/norme, vecteur[2]/norme};
    double delta = ecart_particule - norme;

    vec1[0] -= delta * vecteur_normalise[0];
    vec1[1] -= delta * vecteur_normalise[1];
    vec1[2] -= delta * vecteur_normalise[2];

    vec2[0] += delta * vecteur_normalise[0];
    vec2[1] += delta * vecteur_normalise[1];
    vec2[2] += delta * vecteur_normalise[2];
}

void creation_1_rnap_erwan(double** R, int nb_rnap, int* positions_bille_rnap, double*** R_rnap, int N, int rnap_subunits, int debut_segment, int fin_segment, int nb_rnap_initial){
    if(nb_rnap_initial > 0){
        double taille_pol = 2; 
        double threshold = 1.0;  // Seuil de distance pour vérifier les collisions
        int promoteur = debut_segment;
        int rnap = 0;
            if(rnap==0){promoteur = debut_segment;}
            int position_valide = 0;  // Indicateur de validité de la position
        while (!position_valide) {
                // Calcul des sous-unités de la RNAP autour du promoteur
                double x0 = R[promoteur][0], y0 = R[promoteur][1], z0 = R[promoteur][2];
                double vecteur1[3] = {R[promoteur + 2][0], R[promoteur + 2][1], R[promoteur + 2][2]};
                double vecteur2[3] = {R[promoteur - 2][0], R[promoteur - 2][1], R[promoteur - 2][2]};
                double vecteur_axe[3] = {vecteur1[0] - vecteur2[0], vecteur1[1] - vecteur2[1], vecteur1[2] - vecteur2[2]};
                double norme_axe = sqrt(vecteur_axe[0] * vecteur_axe[0] + vecteur_axe[1] * vecteur_axe[1] + vecteur_axe[2] * vecteur_axe[2]);
                double vecteur_normal[3] = {0,-vecteur_axe[1], vecteur_axe[1]*vecteur_axe[1]/vecteur_axe[2]}; 
                
                double vecteur_tangent[3] = {
                    vecteur_axe[1] * vecteur_normal[2] - vecteur_axe[2] * vecteur_normal[1],
                    vecteur_axe[2] * vecteur_normal[0] - vecteur_axe[0] * vecteur_normal[2],
                    vecteur_axe[0] * vecteur_normal[1] - vecteur_axe[1] * vecteur_normal[0]};
                double norme_tangent = sqrt(vecteur_tangent[0] * vecteur_tangent[0] + vecteur_tangent[1] * vecteur_tangent[1] + vecteur_tangent[2] * vecteur_tangent[2]);
                double norme_normal = sqrt(vecteur_normal[0] * vecteur_normal[0] + vecteur_normal[1] * vecteur_normal[1] + vecteur_normal[2] * vecteur_normal[2]);

                // Normalisation des vecteurs
                for (int j = 0; j < 3; j++) {
                vecteur_tangent[j] = vecteur_tangent[j] / norme_tangent*0.9;
                vecteur_normal[j] = vecteur_normal[j] / norme_normal*0.9;
                }

                position_valide = 1;  // On suppose que la position est valide
                // Vérification des sous-unités de la RNAP
                for (int p = 0; p < rnap_subunits; p++) {
                    // Calcul de la position de la sous-unité p
                    double xi = vecteur_normal[0] * cos((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2 + vecteur_tangent[0] * sin((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2;
                    double yi = vecteur_normal[1] * cos((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2 + vecteur_tangent[1] * sin((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2;
                    double zi = vecteur_normal[2] * cos((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2 + vecteur_tangent[2] * sin((p + 1) * 2 * PI / rnap_subunits) * taille_pol / 2;

                    double position_subunit[3] = {x0 + xi, y0 + yi, z0 + zi};

                    // Vérification si cette position est proche d'une bille existante
                    for (int i = 0; i < N; i++) {  // Parcourir les billes avant le promoteur
                        double distance_bille = distance(position_subunit, R[i]);
                        if (distance_bille < threshold) {
                            ecarter_deux_particules(R[i],position_subunit,threshold);
                        }
                    }
                    // Si la position est valide, on met à jour R_rnap
                    if (position_valide) {
                        R_rnap[rnap][p][0] = x0 + xi;
                        R_rnap[rnap][p][1] = y0 + yi;
                        R_rnap[rnap][p][2] = z0 + zi;
                    }
                }
            }

            // Enregistrer la position du promoteur pour cette RNAP
            positions_bille_rnap[rnap] = promoteur;
            //printf("Promoteur final pour la RNAP %d : %d\n", rnap, promoteur);

            // Vérification des distances pour le debug
            double dist_subunit_chrom = distance(R[promoteur], R_rnap[0][0]);
            double dist_subunit = distance(R_rnap[0][1], R_rnap[0][0]);
            //printf("Distance entre deux subunits : %.2f \n", dist_subunit);
            //printf("Distance entre une subunit et le monomere : %.2f \n", dist_subunit_chrom);
    }
}


void lennard_jones_forces_rnap_rnap(double ***R_rnap,
                                    int nb_rnap,
                                    double epsilon,
                                    double sigma6,
                                    double sigma12,
                                    double cut_rnap2,
                                    double Delta)
{

    double epsilon_rep = 10 * epsilon; 
    double epsilon_att = epsilon; 

    for (int rnap_i = 0; rnap_i < nb_rnap; rnap_i++) {

        for (int sub_i = 0; sub_i < 8; sub_i++) {

            double *Ri = R_rnap[rnap_i][sub_i];

            // --- Interactions RNAP_i ↔ RNAP_j (j > i pour éviter les doublons) ---
            for (int rnap_j = rnap_i + 1; rnap_j < nb_rnap; rnap_j++) {
                for (int sub_j = 0; sub_j < 8; sub_j++) {

                    double *Rj = R_rnap[rnap_j][sub_j];

                    double dx = Ri[0] - Rj[0];
                    double dy = Ri[1] - Rj[1];
                    double dz = Ri[2] - Rj[2];

                    double r2 = dx * dx + dy * dy + dz * dz;
                    if (r2 > cut_rnap2 || r2 == 0.0)
                        continue;

                    // Pré-calculs pour le potentiel LJ
                    double r8  = r2 * r2 * r2 * r2;
                    double r14 = r8 * r2 * r2 * r2;

                    // Force de Lennard-Jones classique : F = -dU/dr
                    double f = 4.0 * (12.0 * epsilon_rep * sigma12 / r14
                                    - 6.0 * epsilon_att * sigma6 / r8);

                    // Saturation (pour éviter les explosions)
                    const double f_max = 600.0;
                    if (f > f_max)
                        f = f_max;

                    // Application symétrique de la force
                    Ri[0] += Delta * f * dx;
                    Ri[1] += Delta * f * dy;
                    Ri[2] += Delta * f * dz;

                    Rj[0] -= Delta * f * dx;
                    Rj[1] -= Delta * f * dy;
                    Rj[2] -= Delta * f * dz;
                }
            }

            // --- Auto-interactions entre sous-unités du même RNAP ---
            for (int sub_j = 0; sub_j < 8; sub_j++) {
                if (sub_i == sub_j) continue;

                double *Rj = R_rnap[rnap_i][sub_j];
                double dx = Ri[0] - Rj[0];
                double dy = Ri[1] - Rj[1];
                double dz = Ri[2] - Rj[2];
                double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 > cut_rnap2 || r2 == 0.0)
                    continue;

                double r8  = r2 * r2 * r2 * r2;
                double r14 = r8 * r2 * r2 * r2;

                // Force modifiée : plus répulsive entre sous-unités d’un même RNAP
                double f = 4.0 * (12.0 * epsilon_rep * sigma12 / r14
                                - 6.0 * epsilon_att * sigma6 / r8);

                const double f_max = 300.0;
                if (f > f_max)
                    f = f_max;

                Ri[0] += Delta * f * dx;
                Ri[1] += Delta * f * dy;
                Ri[2] += Delta * f * dz;

                Rj[0] -= Delta * f * dx;
                Rj[1] -= Delta * f * dy;
                Rj[2] -= Delta * f * dz;
            }
        }
    }
}
