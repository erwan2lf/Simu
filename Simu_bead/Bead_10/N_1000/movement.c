#include "movement.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "basic_functions.h"

double** polymere_brownian_motion(double** R, double K, double Delta, int N, double** r_new,
                                   double K_bend, double **bending_forces, int attache, int plan,
                                   int t, FILE *test, int bending, int truc, int T,
                                   FILE *fichier_force, int periode_enregistrement_force,
                                   FILE* fichier_force_thermique, int temperature) {
    double d1, d0, di;
    double F_alea, pot_harm[3];
    double deplacement = 0.0;
    double epsilon_dist = 1e-8;

    double total_force_thermique = 0.0;
    double total_force_harmonique = 0.0;

    if (t % periode_enregistrement_force == 0) {
        fprintf(fichier_force, "Timestep %d\n", t);
    }

    // === Cas extrémité gauche ===
    if (attache == 0) {
        int p = 0, p2 = 1;

        if (t % periode_enregistrement_force == 0) {
            fprintf(fichier_force, " %d", p);
        }

        di = distance(R[p2], R[p]);

        for (int j = 0; j < 3; j++) {
            deplacement = -R[p][j];
            pot_harm[j] = (di > epsilon_dist) ? K * (R[p2][j] - R[p][j]) * (1.0 - 1.0 / di) : 0.0;
            r_new[p][j] = R[p][j] + Delta * pot_harm[j];
            deplacement += r_new[p][j];

            if (t % periode_enregistrement_force == 0) {
                fprintf(fichier_force, " %lf %lf", pot_harm[j], deplacement);
            }
        }

        double Fh_norm2 = 0.0;
        for (int j = 0; j < 3; j++) Fh_norm2 += pot_harm[j] * pot_harm[j];
        total_force_harmonique += sqrt(Fh_norm2);
        if(temperature == 1){
            double Fvec2 = 0.0;
            for (int j = 0; j < 3; j++) {
                F_alea = sqrt(2 * Delta) * randn();
                Fvec2 += F_alea * F_alea / (2 * Delta);
                deplacement = -r_new[p][j];
                r_new[p][j] += F_alea;
                deplacement += r_new[p][j];
            }
            total_force_thermique += sqrt(Fvec2);
        }
        

        if (t % periode_enregistrement_force == 0) {
            fprintf(fichier_force, "\n");
        }
    }

    // === Cas du milieu ===
    for (int i = 1; i < N - 1; i++) {
        if (t % periode_enregistrement_force == 0) {
            fprintf(fichier_force, " %d", i);
        }

        d1 = distance(R[i - 1], R[i]);
        d0 = distance(R[i + 1], R[i]);

        for (int j = 0; j < 3; j++) {
            deplacement = -R[i][j];
            pot_harm[j] = 0.0;

            if (d0 > epsilon_dist)
                pot_harm[j] += K * (R[i + 1][j] - R[i][j]) * (1.0 - 1.0 / d0);
            if (d1 > epsilon_dist)
                pot_harm[j] += K * (R[i - 1][j] - R[i][j]) * (1.0 - 1.0 / d1);

            r_new[i][j] = R[i][j] + Delta * pot_harm[j];
            deplacement += r_new[i][j];

            if (t % periode_enregistrement_force == 0) {
                fprintf(fichier_force, " %lf %lf", pot_harm[j], deplacement);
            }
        }

        double Fh_norm2 = 0.0;
        for (int j = 0; j < 3; j++) Fh_norm2 += pot_harm[j] * pot_harm[j];
        total_force_harmonique += sqrt(Fh_norm2);

        if(temperature == 1){
            double Fvec2 = 0.0;
            for (int j = 0; j < 3; j++) {
                F_alea = sqrt(2 * Delta) * randn();
                Fvec2 += F_alea * F_alea / (2 * Delta);

                deplacement = -r_new[i][j];
                r_new[i][j] += F_alea;
                deplacement += r_new[i][j];

                if (bending == 1) {
                    r_new[i][j] += Delta * bending_forces[i][j];
                }
            }
            total_force_thermique += sqrt(Fvec2);
        }
        

        if (plan == 1 && r_new[i][2] < 0) {
            r_new[i][2] = -r_new[i][2];
        }

        if (t % periode_enregistrement_force == 0) {
            fprintf(fichier_force, "\n");
        }
    }

    // === Cas extrémité droite ===
    if (attache == 0) {
        int p = N - 1, p2 = N - 2;

        if (t % periode_enregistrement_force == 0) {
            fprintf(fichier_force, " %d", p);
        }

        di = distance(R[p2], R[p]);

        for (int j = 0; j < 3; j++) {
            deplacement = -R[p][j];
            pot_harm[j] = (di > epsilon_dist) ? K * (R[p2][j] - R[p][j]) * (1.0 - 1.0 / di) : 0.0;
            r_new[p][j] = R[p][j] + Delta * pot_harm[j];
            deplacement += r_new[p][j];

            if (t % periode_enregistrement_force == 0) {
                fprintf(fichier_force, " %lf %lf", pot_harm[j], deplacement);
            }
        }

        double Fh_norm2 = 0.0;
        for (int j = 0; j < 3; j++) Fh_norm2 += pot_harm[j] * pot_harm[j];
        total_force_harmonique += sqrt(Fh_norm2);

        if(temperature == 1){
             double Fvec2 = 0.0;
            for (int j = 0; j < 3; j++) {
                F_alea = sqrt(2 * Delta) * randn();
                Fvec2 += F_alea * F_alea / ( 2 * Delta);

                deplacement = -r_new[p][j];
                r_new[p][j] += F_alea;
                deplacement += r_new[p][j];
            }
            total_force_thermique += sqrt(Fvec2);
        }
       
        if (t % periode_enregistrement_force == 0) {
            fprintf(fichier_force, "\n");
        }
    }

    // === Enregistrement des forces moyennes ===
    if (t % periode_enregistrement_force == 0) {
        double moyenne_thermique = total_force_thermique / N;
        double moyenne_harmonique = total_force_harmonique / N;

        fprintf(fichier_force, "# Moyennes de forces (thermique | harmonique) : %lf %lf\n",
                moyenne_thermique, moyenne_harmonique);
        fprintf(fichier_force, "# Moyennes (csv) : %d,%lf,%lf\n", t, moyenne_thermique, moyenne_harmonique);
    }

    return r_new;
}

void confinement_sphere(double **R, int N, double r_sphere){
    double origine[3] = {0,0,0};
    for (int i = 0; i < N; i++){
        if(distance(R[i], origine) > r_sphere){
            for (int j = 0; j < 3; j++){
                if (fabs(R[i][j]) > r_sphere){
                    if (R[i][j] < 0){
                        R[i][j] = - r_sphere - (R[i][j] + r_sphere);
                    }
                    else{
                        R[i][j] = r_sphere - (R[i][j] - r_sphere);
                    }
                }
            }
        }
    }
}

double** gaz_motion(double **R, int N, double ** r_new, double Delta, int plan, int attache){

    if(attache == 1){
        for(int i = 1; i < N-1; i++){
            for(int j = 0; j < 3; j++){
            r_new[i][j] = R[i][j] + sqrt(2 * Delta) * randn();
            }
        }
    }
    if(attache == 0){
        for(int i = 0; i < N; i++){
            for(int j = 0; j < 3; j++){
                r_new[i][j] = R[i][j] + sqrt(2 * Delta) * randn();
            }
        }
    }
    

    if(plan == 1){
        for(int i = 0; i < N; i++){
            if (r_new[i][2] < 0) {
                r_new[i][2] = - r_new[i][2];
            }
        }
    }

    return r_new;
}