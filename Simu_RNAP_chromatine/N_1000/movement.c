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
    (void)K_bend;
    (void)t;
    (void)test;
    (void)truc;
    (void)T;
    (void)fichier_force;
    (void)periode_enregistrement_force;
    (void)fichier_force_thermique;

    const double eps = 1e-8;
    const double eps2 = eps * eps;
    const int do_temp = (temperature == 1);
    const int do_bending = (bending == 1);
    const int do_plan = (plan == 1);
    const double sigma_th = do_temp ? sqrt(2.0 * Delta) : 0.0;

    // tableau temporaire des forces harmoniques
    double (*F)[3] = calloc((size_t)N, sizeof(*F));
    if (F == NULL) {
        perror("calloc F");
        exit(EXIT_FAILURE);
    }

    // ----------------------------------------
    // 1) calcul des forces par liaison
    // ----------------------------------------
    for (int i = 0; i < N - 1; i++) {
        double *Ri = R[i];
        double *Rj = R[i + 1];

        double dx = Rj[0] - Ri[0];
        double dy = Rj[1] - Ri[1];
        double dz = Rj[2] - Ri[2];
        double r2 = dx * dx + dy * dy + dz * dz;

        if (r2 <= eps2) continue;

        double r = sqrt(r2);
        double coef = K * (1.0 - 1.0 / r);

        double fx = coef * dx;
        double fy = coef * dy;
        double fz = coef * dz;

        F[i][0]     += fx;
        F[i][1]     += fy;
        F[i][2]     += fz;

        F[i + 1][0] -= fx;
        F[i + 1][1] -= fy;
        F[i + 1][2] -= fz;
    }

    // ----------------------------------------
    // 2) mise à jour des positions
    // ----------------------------------------
    for (int i = 0; i < N; i++) {
        double *Ri = R[i];

        // si attache==1, on ne bouge pas les extrémités
        if (attache == 1 && (i == 0 || i == N - 1)) {
            continue;
        }

        Ri[0] += Delta * F[i][0];
        Ri[1] += Delta * F[i][1];
        Ri[2] += Delta * F[i][2];

        if (do_bending && i > 0 && i < N - 1) {
            Ri[0] += Delta * bending_forces[i][0];
            Ri[1] += Delta * bending_forces[i][1];
            Ri[2] += Delta * bending_forces[i][2];
        }

        if (do_temp) {
            Ri[0] += sigma_th * randn();
            Ri[1] += sigma_th * randn();
            Ri[2] += sigma_th * randn();
        }

        if (do_plan && i > 0 && i < N - 1 && Ri[2] < 0.0) {
            Ri[2] = -Ri[2];
        }
    }

    free(F);
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