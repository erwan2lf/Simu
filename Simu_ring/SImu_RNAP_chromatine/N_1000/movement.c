#include "movement.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "basic_functions.h"

/**
 * @brief Simulates one Brownian dynamics step for a polymer chain.
 *
 * This function updates the positions of all monomers in the polymer
 * according to harmonic spring forces, bending forces (optional),
 * and thermal noise. The update rule follows an overdamped Langevin
 * dynamics discretized with Euler–Maruyama scheme.
 *
 * Forces and displacements can optionally be written to files for debugging.
 *
 * @param R             Current monomer positions [N][3]
 * @param K             Spring constant between adjacent monomers
 * @param Delta         Integration timestep
 * @param N             Number of monomers in the chain
 * @param r_new         Output array for updated positions [N][3]
 * @param K_bend        Bending stiffness coefficient
 * @param bending_forces Precomputed bending forces [N][3]
 * @param attache        If 0 → both ends free, if 1 → ends fixed
 * @param plan           If 1 → reflective boundary in z<0
 * @param t              Current simulation time step
 * @param test           Optional debug file (unused here)
 * @param bending        If 1 → apply bending forces
 * @param truc           Unused / unclear parameter (to be removed or clarified)
 * @param T              Total simulation time (unused)
 * @param fichier_force  File handle for recording harmonic/thermal forces
 * @param periode_enregistrement_force Period (in steps) between force recordings
 * @param fichier_force_thermique Optional file for thermal forces (currently unused)
 * @param temperature    If 1 → include thermal noise
 *
 * @return Pointer to the updated position array r_new
 */
void polymere_brownian_motion(double **R, double K, double Delta, int N,
                              double K_bend, double **bending_forces, int attache, int plan,
                              int t, FILE *test, int bending, int truc, int T,
                              FILE *fichier_force, int periode_enregistrement_force,
                              FILE *fichier_force_thermique, int temperature)
{
    double pot_harm[3];
    const double eps = 1e-8;

    double total_force_thermique = 0.0;
    double total_force_harmonique = 0.0;

    if (t % periode_enregistrement_force == 0)
        fprintf(fichier_force, "Timestep %d\n", t);

    // --- Cas extrémité gauche ---
    if (attache == 0) {
        int p = 0, p2 = 1;
        double di = distance(R[p2], R[p]);

        for (int j = 0; j < 3; j++) {
            pot_harm[j] = (di > eps) ? K * (R[p2][j] - R[p][j]) * (1.0 - 1.0 / di) : 0.0;
            R[p][j] = R[p][j] + Delta * pot_harm[j];
        }

        total_force_harmonique += sqrt(pot_harm[0]*pot_harm[0] + pot_harm[1]*pot_harm[1] + pot_harm[2]*pot_harm[2]);

        if (temperature == 1) {
            double Fvec2 = 0.0;
            for (int j = 0; j < 3; j++) {
                double F_alea = sqrt(2 * Delta) * randn();
                R[p][j] += F_alea;
                Fvec2 += F_alea * F_alea / (2 * Delta);
            }
            total_force_thermique += sqrt(Fvec2);
        }
    }

    // --- Cas du milieu ---
    for (int i = 1; i < N - 1; i++) {
        double d1 = distance(R[i - 1], R[i]);
        double d0 = distance(R[i + 1], R[i]);

        for (int j = 0; j < 3; j++) {
            pot_harm[j] = 0.0;
            if (d0 > eps)
                pot_harm[j] += K * (R[i + 1][j] - R[i][j]) * (1.0 - 1.0 / d0);
            if (d1 > eps)
                pot_harm[j] += K * (R[i - 1][j] - R[i][j]) * (1.0 - 1.0 / d1);

            R[i][j] = R[i][j] + Delta * pot_harm[j];
            if (bending == 1)
                R[i][j] += Delta * bending_forces[i][j];
        }

        total_force_harmonique += sqrt(pot_harm[0]*pot_harm[0] + pot_harm[1]*pot_harm[1] + pot_harm[2]*pot_harm[2]);

        if (temperature == 1) {
            double Fvec2 = 0.0;
            for (int j = 0; j < 3; j++) {
                double F_alea = sqrt(2 * Delta) * randn();
                R[i][j] += F_alea;
                Fvec2 += F_alea * F_alea / (2 * Delta);
            }
            total_force_thermique += sqrt(Fvec2);
        }

        // Réflexion sur un plan z=0
        if (plan == 1 && R[i][2] < 0)
            R[i][2] = -R[i][2];
    }

    // --- Cas extrémité droite ---
    if (attache == 0) {
        int p = N - 1, p2 = N - 2;
        double di = distance(R[p2], R[p]);

        for (int j = 0; j < 3; j++) {
            pot_harm[j] = (di > eps) ? K * (R[p2][j] - R[p][j]) * (1.0 - 1.0 / di) : 0.0;
            R[p][j] = R[p][j] + Delta * pot_harm[j];
        }

        total_force_harmonique += sqrt(pot_harm[0]*pot_harm[0] + pot_harm[1]*pot_harm[1] + pot_harm[2]*pot_harm[2]);

        if (temperature == 1) {
            double Fvec2 = 0.0;
            for (int j = 0; j < 3; j++) {
                double F_alea = sqrt(2 * Delta) * randn();
                R[p][j] += F_alea;
                Fvec2 += F_alea * F_alea / (2 * Delta);
            }
            total_force_thermique += sqrt(Fvec2);
        }
    }

    // --- Enregistrement des forces moyennes ---
    if (t % periode_enregistrement_force == 0) {
        double moyenne_thermique = total_force_thermique / N;
        double moyenne_harmonique = total_force_harmonique / N;
        fprintf(fichier_force, "# Moyennes (thermique | harmonique): %lf %lf\n", moyenne_thermique, moyenne_harmonique);
        fprintf(fichier_force, "# CSV: %d,%lf,%lf\n", t, moyenne_thermique, moyenne_harmonique);
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