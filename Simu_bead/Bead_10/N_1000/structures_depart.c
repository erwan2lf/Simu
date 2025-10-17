#include "structures_depart.h"
#include "basic_functions.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define M_PI 3.14159265

// Fonction pour créer un polymère
double** creation_polymere(int N, double a, double ecart) {
    double** R = (double**)malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) {
        R[i] = (double*)malloc(3 * sizeof(double));}
    for (int i = 0; i < N; i++) {
        R[i][1] = ecart;
        R[i][2] = 0;}
    for (int i = 0; i < N; i++) {R[i][0] = -N * a / 2 + i * a;}

    //for (int part=0; part<N; part++) {printf("%d %1f %1f %1f \n",part, R[part][0], R[part][1],R[part][2] );}
    return R; }

int find_max_timestep(const char* nom_fichier) {
    FILE* fichier = fopen(nom_fichier, "r");
    if (fichier == NULL) {
        printf("Erreur : Impossible d'ouvrir le fichier %s\n", nom_fichier);
        return -1;
    }

    int max_timestep = -1;
    int stock = -1;
    char buffer[256];
    while (fgets(buffer, sizeof(buffer), fichier) != NULL) {

        if (strstr(buffer, "ITEM: TIMESTEP") != NULL) {
            int timestep;
            fscanf(fichier, "%d", &timestep);
            if (timestep > max_timestep) {
                stock = max_timestep;
                max_timestep = timestep;
            }
        }
    }

    fclose(fichier);
    printf("%d %d \n", max_timestep, stock);
    return stock;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double** recuperer_derniere_structure(const char* nom_fichier, int N) {
    FILE* fichier = fopen(nom_fichier, "r");
    if (fichier == NULL) {
        printf("Erreur : Impossible d'ouvrir le fichier %s\n", nom_fichier);
        return NULL;
    }

    int t;
    double x, y, z;
    double** R_matrix = (double**)malloc(N * sizeof(double*));
    if (R_matrix == NULL) {
        printf("Erreur : Échec de l'allocation de mémoire pour R_matrix\n");
        fclose(fichier);
        return NULL;
    }

    for (int i = 0; i < N; i++) {
        R_matrix[i] = (double*)malloc(3 * sizeof(double));
        if (R_matrix[i] == NULL) {
            printf("Erreur : Échec de l'allocation de mémoire pour R_matrix[%d]\n", i);
            for (int j = 0; j < i; j++) free(R_matrix[j]);
            free(R_matrix);
            fclose(fichier);
            return NULL;
        }
    }

    char buffer[256];
    while (fgets(buffer, sizeof(buffer), fichier) != NULL) {
        if (strstr(buffer, "ITEM: ATOMS id type xs ys zs") != NULL) {
            for (int i = 0; i < N; i++) {
                if (fscanf(fichier, "%*d %*d %lf %lf %lf", &x, &y, &z) != 3) {
                    printf("Erreur : Lecture échouée à la particule %d\n", i);
                    for (int j = 0; j < N; j++) free(R_matrix[j]);
                    free(R_matrix);
                    fclose(fichier);
                    return NULL;
                }
                R_matrix[i][0] = x * 1000.0;
                R_matrix[i][1] = y * 1000.0;
                R_matrix[i][2] = z * 1000.0;
            }
        }
        if (strstr(buffer, "ITEM: TIMESTEP") != NULL) {
            fscanf(fichier, "%d", &t);
        }
    }

    fclose(fichier);
    return R_matrix;
}

void creation_polymere_solenoide(int N, double a, double ecart, double epaisseur, double** R) {
    double distance_bout_a_bout = 30*20;// 30billes de taille 20
    double i_pas = distance_bout_a_bout/N ;
    double th = 2 * M_PI / 300.0; // 10 points sur le cercle
    double theta = 3 * M_PI / 2.0;
    double rayon = a / 2.0 / sin(th / 2.0);
    for (int i = 0; i < N; ++i) {
        R[i][0] = i * i_pas + ecart;
        R[i][1] = rayon * cos(theta);
        R[i][2] = rayon * sin(theta);
        theta += th;
    }
    // Ajustement de Z pour l'épaisseur et l'écart
    double min_Z = R[0][2];
    for (int i = 1; i < N; ++i) {
        if (R[i][2] < min_Z) {
            min_Z = R[i][2];
        }
    }
    for (int i = 0; i < N; ++i) {
        R[i][2] = R[i][2] - min_Z - epaisseur + ecart;
    }

    // Centrage de X
    double max_X = R[0][0];
    for (int i = 1; i < N; ++i) {
        if (R[i][0] > max_X) {
            max_X = R[i][0];
        }
    }
    for (int i = 0; i < N; ++i) {
        R[i][0] = R[i][0] - max_X / 2.0;
    }
}

void creation_polymere_droit(int N, double a, double ecart, double** R){
    for ( int i = 0; i < N; i ++){
        R[i][0] = 0;
        R[i][1] = 0;
        R[i][2] = i * a + ecart;
    }
}

void creation_fractal_globule(int N, double a, double ecart, double** R){
   for (int i = 0; i < N; i++){
    /*if (i<(int)N/3){
        R[i][0]=0;
        R[i][1]=0;
        R[i][2]=i*a+ecart;
    }
    if(i>=(int)N/3 && i<=(int)2*N/3){
        R[i][0]=0;
        R[i][1]=(i-N/3)*a+ecart;
        R[i][2]=(N/3)*a+ecart;
    }
    if ( i > (int)2*N/3){
        R[i][0]=0;
        R[i][1]=(N/3)*a+ecart;
        R[i][2]=-(i-N+1)*a+ecart;
    }*/
   if (i<(int)N/2){
        R[i][0]=0;
        R[i][1]=0;
        R[i][2]=i*a+ecart;
    }
    if ( i >= (int)N/2){
        R[i][0]=0;
        R[i][1]=a;
        R[i][2]=-(i-N+1)*a+ecart;
    }
   }

   }

void creation_structure_knot(int N, double a, double **R){
    for(int i = 0; i < N; i++){
        if(i <= (N-24)/2){
            R[i][0] = 0; 
            R[i][1] = 0;
            R[i][2] = i*a;
        }
        
        if( i >= (N+24)/2){
            R[i][0] = 1; 
            R[i][1] = 0; 
            R[i][2] = -(i-N+1)*a;
        }
    }

    R[(N-24)/2+1][0] = 0;
    R[(N-24)/2+1][1] = 0;
    R[(N-24)/2+1][2] = R[(N-24)/2][2] +1;

    R[(N-24)/2+2][0] = 0;
    R[(N-24)/2+2][1] = 1;
    R[(N-24)/2+2][2] = R[(N-24)/2][2] +1;

    R[(N-24)/2+3][0] = 0;
    R[(N-24)/2+3][1] = 2;
    R[(N-24)/2+3][2] = R[(N-24)/2][2] +1;

    R[(N-24)/2+4][0] = 1;
    R[(N-24)/2+4][1] = 2;
    R[(N-24)/2+4][2] = R[(N-24)/2][2] +1;

    R[(N-24)/2+5][0] = 2;
    R[(N-24)/2+5][1] = 2;
    R[(N-24)/2+5][2] = R[(N-24)/2][2] +1;

    R[(N-24)/2+6][0] = 3;
    R[(N-24)/2+6][1] = 2;
    R[(N-24)/2+6][2] = R[(N-24)/2][2] +1;

    R[(N-24)/2+7][0] = 3;
    R[(N-24)/2+7][1] = 2;
    R[(N-24)/2+7][2] = R[(N-24)/2][2];

    R[(N-24)/2+8][0] = 3;
    R[(N-24)/2+8][1] = 2;
    R[(N-24)/2+8][2] = R[(N-24)/2][2] - 1;

    R[(N-24)/2+9][0] = 3;
    R[(N-24)/2+9][1] = 1;
    R[(N-24)/2+9][2] = R[(N-24)/2][2] - 1;

    R[(N-24)/2+10][0] = 2;
    R[(N-24)/2+10][1] = 1;
    R[(N-24)/2+10][2] = R[(N-24)/2][2] - 1;

    R[(N-24)/2+11][0] = 1;
    R[(N-24)/2+11][1] = 1;
    R[(N-24)/2+11][2] = R[(N-24)/2][2] - 1;

    R[(N-24)/2+12][0] = 1;
    R[(N-24)/2+12][1] = 1;
    R[(N-24)/2+12][2] = R[(N-24)/2][2];

    R[(N-24)/2+13][0] = 1;
    R[(N-24)/2+13][1] = 1;
    R[(N-24)/2+13][2] = R[(N-24)/2][2] + 1;

    R[(N-24)/2+14][0] = 1;
    R[(N-24)/2+14][1] = 1;
    R[(N-24)/2+14][2] = R[(N-24)/2][2] + 2;

    R[(N-24)/2+15][0] = 1;
    R[(N-24)/2+15][1] = 2;
    R[(N-24)/2+15][2] = R[(N-24)/2][2] + 2;

    R[(N-24)/2+16][0] = 1;
    R[(N-24)/2+16][1] = 3;
    R[(N-24)/2+16][2] = R[(N-24)/2][2] + 2;

    R[(N-24)/2+17][0] = 2;
    R[(N-24)/2+17][1] = 3;
    R[(N-24)/2+17][2] = R[(N-24)/2][2] + 2;

    R[(N-24)/2+18][0] = 2;
    R[(N-24)/2+18][1] = 3;
    R[(N-24)/2+18][2] = R[(N-24)/2][2] + 1;

    R[(N-24)/2+19][0] = 2;
    R[(N-24)/2+19][1] = 3;
    R[(N-24)/2+19][2] = R[(N-24)/2][2];

    R[(N-24)/2+20][0] = 2;
    R[(N-24)/2+20][1] = 2;
    R[(N-24)/2+20][2] = R[(N-24)/2][2];

    R[(N-24)/2+21][0] = 2;
    R[(N-24)/2+21][1] = 1;
    R[(N-24)/2+21][2] = R[(N-24)/2][2];

    R[(N-24)/2+22][0] = 2;
    R[(N-24)/2+22][1] = 0;
    R[(N-24)/2+22][2] = R[(N-24)/2][2];

    R[(N-24)/2+23][0] = 1;
    R[(N-24)/2+23][1] = 0;
    R[(N-24)/2+23][2] = R[(N-24)/2][2];

}
