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
    const int do_temp = (temperature == 1);
    const int do_bending = (bending == 1);
    const int do_plan = (plan == 1);
    const double sigma_th = do_temp ? sqrt(2.0 * Delta) : 0.0;

    // ----------------------------
    // Extrémité gauche
    // ----------------------------
    if (attache == 0) {
        int p = 0, p2 = 1;
        double *Rp = R[p];
        double *Rq = R[p2];

        double dx = Rq[0] - Rp[0];
        double dy = Rq[1] - Rp[1];
        double dz = Rq[2] - Rp[2];
        double r2 = dx*dx + dy*dy + dz*dz;

        if (r2 > eps * eps) {
            double r = sqrt(r2);
            double coef = K * (1.0 - 1.0 / r);

            Rp[0] += Delta * coef * dx;
            Rp[1] += Delta * coef * dy;
            Rp[2] += Delta * coef * dz;
        }

        if (do_temp) {
            Rp[0] += sigma_th * randn();
            Rp[1] += sigma_th * randn();
            Rp[2] += sigma_th * randn();
        }
    }

    // ----------------------------
    // Milieu
    // ----------------------------
    for (int i = 1; i < N - 1; i++) {
        double *Ri   = R[i];
        double *Rim1 = R[i - 1];
        double *Rip1 = R[i + 1];

        double fx = 0.0, fy = 0.0, fz = 0.0;

        // contribution du voisin i+1
        {
            double dx = Rip1[0] - Ri[0];
            double dy = Rip1[1] - Ri[1];
            double dz = Rip1[2] - Ri[2];
            double r2 = dx*dx + dy*dy + dz*dz;

            if (r2 > eps * eps) {
                double r = sqrt(r2);
                double coef = K * (1.0 - 1.0 / r);
                fx += coef * dx;
                fy += coef * dy;
                fz += coef * dz;
            }
        }

        // contribution du voisin i-1
        {
            double dx = Rim1[0] - Ri[0];
            double dy = Rim1[1] - Ri[1];
            double dz = Rim1[2] - Ri[2];
            double r2 = dx*dx + dy*dy + dz*dz;

            if (r2 > eps * eps) {
                double r = sqrt(r2);
                double coef = K * (1.0 - 1.0 / r);
                fx += coef * dx;
                fy += coef * dy;
                fz += coef * dz;
            }
        }

        Ri[0] += Delta * fx;
        Ri[1] += Delta * fy;
        Ri[2] += Delta * fz;

        if (do_bending) {
            Ri[0] += Delta * bending_forces[i][0];
            Ri[1] += Delta * bending_forces[i][1];
            Ri[2] += Delta * bending_forces[i][2];
        }

        if (do_temp) {
            Ri[0] += sigma_th * randn();
            Ri[1] += sigma_th * randn();
            Ri[2] += sigma_th * randn();
        }

        if (do_plan && Ri[2] < 0.0) {
            Ri[2] = -Ri[2];
        }
    }

    // ----------------------------
    // Extrémité droite
    // ----------------------------
    if (attache == 0) {
        int p = N - 1, p2 = N - 2;
        double *Rp = R[p];
        double *Rq = R[p2];

        double dx = Rq[0] - Rp[0];
        double dy = Rq[1] - Rp[1];
        double dz = Rq[2] - Rp[2];
        double r2 = dx*dx + dy*dy + dz*dz;

        if (r2 > eps * eps) {
            double r = sqrt(r2);
            double coef = K * (1.0 - 1.0 / r);

            Rp[0] += Delta * coef * dx;
            Rp[1] += Delta * coef * dy;
            Rp[2] += Delta * coef * dz;
        }

        if (do_temp) {
            Rp[0] += sigma_th * randn();
            Rp[1] += sigma_th * randn();
            Rp[2] += sigma_th * randn();
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