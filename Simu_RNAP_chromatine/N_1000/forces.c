#include "forces.h"
#include "basic_functions.h"   /* randn() */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ================================================================
 * CHROMATINE — ressort harmonique
 * ================================================================ */
void accumulate_spring_forces(double **R, double (*F)[3],
                              double K, int N)
{
    const double eps2 = 1e-16;

    for (int i = 0; i < N - 1; i++) {
        double dx = R[i+1][0] - R[i][0];
        double dy = R[i+1][1] - R[i][1];
        double dz = R[i+1][2] - R[i][2];

        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 <= eps2) continue;

        double r    = sqrt(r2);
        // double coef = K * (1.0 - 1.0 / r);
        double coef = K * (r - 1.0) / r;

        F[i][0]   += coef * dx;
        F[i][1]   += coef * dy;
        F[i][2]   += coef * dz;

        F[i+1][0] -= coef * dx;
        F[i+1][1] -= coef * dy;
        F[i+1][2] -= coef * dz;
    }
}

/* ================================================================
 * CHROMATINE — Lennard-Jones chrom-chrom
 * ================================================================ */
void accumulate_lj_forces(double **R, double (*F)[3],
                          NeighborList *neighbor_lists, int N,
                          double epsilon, double sigma6, double sigma12,
                          int attache)
{
    const double fmax = 300.0;
    const double c12  = 48.0 * epsilon * sigma12;
    const double c6   = 24.0 * epsilon * sigma6;

    for (int i = 1; i < N - 1; i++) {
        double *Ri = R[i];

        for (int k = 0; k < neighbor_lists[i].count; k++) {
            int j = neighbor_lists[i].neighbors[k];
            if (j <= i) continue;
            double *Rj = R[j];

            double dx = Ri[0] - Rj[0];
            double dy = Ri[1] - Rj[1];
            double dz = Ri[2] - Rj[2];

            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 == 0.0) continue;

            double inv_r2  = 1.0 / r2;
            double inv_r6  = inv_r2 * inv_r2 * inv_r2;
            double inv_r8  = inv_r6 * inv_r2;
            double inv_r14 = inv_r8 * inv_r6;

            double f = c12 * inv_r14 - c6 * inv_r8;
            if (f >  fmax) f =  fmax;
            if (f < -fmax) f = -fmax;

            double fx = f * dx;
            double fy = f * dy;
            double fz = f * dz;

            if (attache == 1) {
                if (i != 0 && i != N-1) {
                    F[i][0] += fx; F[i][1] += fy; F[i][2] += fz;
                }
                if (j != 0 && j != N-1) {
                    F[j][0] -= fx; F[j][1] -= fy; F[j][2] -= fz;
                }
            } else {
                F[i][0] += fx; F[i][1] += fy; F[i][2] += fz;
                F[j][0] -= fx; F[j][1] -= fy; F[j][2] -= fz;
            }
        }
    }
}

/* ================================================================
 * CHROMATINE — confinement sphérique
 * ================================================================ */
void accumulate_conf_forces(double **R, double (*F)[3], int N,
                            const double cdm[3],
                            double r_conf, double sigma_conf,
                            double epsilon_conf)
{
    const double fmax = 300.0;

    for (int i = 0; i < N; i++) {
        double dx = R[i][0] - cdm[0];
        double dy = R[i][1] - cdm[1];
        double dz = R[i][2] - cdm[2];

        double d = sqrt(dx*dx + dy*dy + dz*dz);
        if (d < 1e-12) continue;

        double d_wall = r_conf - d;
        if (d_wall >= sigma_conf || d_wall <= 0.0) continue;

        double inv   = sigma_conf / d_wall;
        double inv6  = inv*inv*inv*inv*inv*inv;
        double inv12 = inv6 * inv6;

        double f_mag = epsilon_conf * (12.0 * inv12 - 6.0 * inv6) / d_wall;
        if (f_mag >  fmax) f_mag =  fmax;
        if (f_mag < -fmax) f_mag = -fmax;

        double nx = -dx / d;
        double ny = -dy / d;
        double nz = -dz / d;

        F[i][0] += f_mag * nx;
        F[i][1] += f_mag * ny;
        F[i][2] += f_mag * nz;
    }
}

/* ================================================================
 * RNAP — ressort en anneau (topologie circulaire)
 * ================================================================ */
void accumulate_ring_forces(double **R_rnap, double (*F_rnap)[3],
                            double K_rnap, double alpha, int nsub)
{
    const double eps2 = 1e-16;

    for (int i = 0; i < nsub; i++) {
        int i_next = (i + 1) % nsub;

        double dx = R_rnap[i_next][0] - R_rnap[i][0];
        double dy = R_rnap[i_next][1] - R_rnap[i][1];
        double dz = R_rnap[i_next][2] - R_rnap[i][2];

        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 <= eps2) continue;

        double r    = sqrt(r2);
        // double coef = K_rnap * (1.0 - alpha / r);
        double coef = K_rnap * (r - 1.0) / r;

        F_rnap[i][0]      += coef * dx;
        F_rnap[i][1]      += coef * dy;
        F_rnap[i][2]      += coef * dz;

        F_rnap[i_next][0] -= coef * dx;
        F_rnap[i_next][1] -= coef * dy;
        F_rnap[i_next][2] -= coef * dz;
    }
}

/* ================================================================
 * RNAP — liaisons diagonales (nsub=8)
 * Paires : (0,4), (1,5), (2,6), (3,7)
 * ================================================================ */
void accumulate_diagonal_forces(double **R_rnap, double (*F_rnap)[3],
                                double K_rnap, double a_transpt, int nsub)
{
    if (nsub != 8) return;

    const double l0   = 3.0 * a_transpt;
    const double eps2 = 1e-16;
    const int pairs[4][2] = {{0,4}, {1,5}, {2,6}, {3,7}};

    for (int p = 0; p < 4; p++) {
        int i = pairs[p][0];
        int j = pairs[p][1];

        double dx = R_rnap[j][0] - R_rnap[i][0];
        double dy = R_rnap[j][1] - R_rnap[i][1];
        double dz = R_rnap[j][2] - R_rnap[i][2];

        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 <= eps2) continue;

        double r    = sqrt(r2);
        // double coef = 2.0 * K_rnap * (1.0 - l0 / r);
        double coef = 2.0 * K_rnap * (r - 1.0) / r;

        F_rnap[i][0] += coef * dx;
        F_rnap[i][1] += coef * dy;
        F_rnap[i][2] += coef * dz;

        F_rnap[j][0] -= coef * dx;
        F_rnap[j][1] -= coef * dy;
        F_rnap[j][2] -= coef * dz;
    }
}

/* ================================================================
 * RNAP ↔ CHROMATINE — couplage progressif (smoothstep)
 * ================================================================ */
void accumulate_bond_rnap_chrom(const Config *cfg,
                                double **R,
                                double **R_rnap,
                                double (*F_chrom)[3],
                                double (*F_rnap)[3],
                                int mono,
                                double avancement)
{
    const int    NB_SUB    = cfg->rnap_subunits;
    const double K_transpt = cfg->K_transpt;
    const double alpha     = cfg->alpha;
    const double eps       = 1e-12;

    if (mono < 0 || mono >= cfg->N - 1) return;

    int p  = mono;
    int p1 = mono + 1;

    double s = avancement;
    if (s < 0.0) s = 0.0;
    if (s > 1.0) s = 1.0;
    double s_eff = s * s * (3.0 - 2.0 * s);

    double w0 = 1.0 - s_eff;
    double w1 = s_eff;

    for (int sub = 0; sub < NB_SUB; sub++) {
        double d0x = R[p][0]  - R_rnap[sub][0];
        double d0y = R[p][1]  - R_rnap[sub][1];
        double d0z = R[p][2]  - R_rnap[sub][2];

        double d1x = R[p1][0] - R_rnap[sub][0];
        double d1y = R[p1][1] - R_rnap[sub][1];
        double d1z = R[p1][2] - R_rnap[sub][2];

        double dist0 = sqrt(d0x*d0x + d0y*d0y + d0z*d0z);
        double dist1 = sqrt(d1x*d1x + d1y*d1y + d1z*d1z);
        if (dist0 < eps) dist0 = eps;
        if (dist1 < eps) dist1 = eps;

        double f0 = (1.0 - (alpha + 1.0) / (2.0 * dist0));
        double f1 = (1.0 - (alpha + 1.0) / (2.0 * dist1));

        double Fx = K_transpt * (w0 * d0x * f0 + w1 * d1x * f1);
        double Fy = K_transpt * (w0 * d0y * f0 + w1 * d1y * f1);
        double Fz = K_transpt * (w0 * d0z * f0 + w1 * d1z * f1);

        F_rnap[sub][0] += Fx;
        F_rnap[sub][1] += Fy;
        F_rnap[sub][2] += Fz;

        F_chrom[p][0]  -= w0 * Fx;
        F_chrom[p][1]  -= w0 * Fy;
        F_chrom[p][2]  -= w0 * Fz;

        F_chrom[p1][0] -= w1 * Fx;
        F_chrom[p1][1] -= w1 * Fy;
        F_chrom[p1][2] -= w1 * Fz;
    }
}

/* ================================================================
 * RNAP ↔ CHROMATINE — LJ répulsif (via neighbor list)
 * ================================================================ */
void accumulate_lj_rnap_chrom(const Config *cfg,
                               double **R_rnap,
                               double **R,
                               double (*F_rnap)[3],
                               double (*F_chrom)[3],
                               NeighborList_rnap **neighbor_lists,
                               int rnap_id,
                               double epsilon,
                               double sigma6, double sigma12,
                               double cut2)
{
    const double fmax = 300.0;
    const double c12  = 4.0 * 12.0 * epsilon * sigma12;

    for (int sub = 0; sub < cfg->rnap_subunits; sub++) {
        double *Ri = R_rnap[sub];

        for (int k = 0; k < neighbor_lists[rnap_id][sub].count; k++) {
            int j = neighbor_lists[rnap_id][sub].neighbors[k];
            double *Rj = R[j];

            double dx = Ri[0] - Rj[0];
            double dy = Ri[1] - Rj[1];
            double dz = Ri[2] - Rj[2];

            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 == 0.0 || r2 > cut2) continue;

            double r8  = r2 * r2 * r2 * r2;
            double r14 = r8 * r2 * r2 * r2;

            double f = c12 / r14;
            if (f > fmax) f = fmax;

            double fx = f * dx;
            double fy = f * dy;
            double fz = f * dz;

            F_rnap[sub][0] += fx;
            F_rnap[sub][1] += fy;
            F_rnap[sub][2] += fz;

            F_chrom[j][0] -= fx;
            F_chrom[j][1] -= fy;
            F_chrom[j][2] -= fz;
        }
    }
}

/* ================================================================
 * RNAP ↔ RNAP — LJ inter et intra
 * ================================================================ */
void accumulate_lj_rnap_rnap(const Config *cfg,
                              double ***R_rnap,
                              double (**F_rnap_all)[3],
                              int nb_rnap,
                              int *l_rnap,
                              double epsilon,
                              double sigma6, double sigma12,
                              double cut2)
{
    const int    nsub       = cfg->rnap_subunits;
    const double fmax_inter = 1000.0;
    const double fmax_intra = 300.0;

    const double c12_inter = 48.0 * epsilon * sigma12;
    const double c6_inter  = 24.0 * epsilon * sigma6;

    const double c12_intra = 4.0 * 12.0 * 1.5 * epsilon * sigma12;
    const double c6_intra  = 4.0 *  6.0 * 1.0 * epsilon * sigma6;

    for (int ri = 0; ri < nb_rnap; ri++) {
        if (l_rnap[ri] < 0) continue;

        for (int si = 0; si < nsub; si++) {
            double *Ri = R_rnap[ri][si];

            /* Inter-RNAP */
            for (int rj = ri + 1; rj < nb_rnap; rj++) {
                if (l_rnap[rj] < 0) continue;

                for (int sj = 0; sj < nsub; sj++) {
                    double *Rj = R_rnap[rj][sj];

                    double dx = Ri[0] - Rj[0];
                    double dy = Ri[1] - Rj[1];
                    double dz = Ri[2] - Rj[2];

                    double r2 = dx*dx + dy*dy + dz*dz;
                    if (r2 == 0.0 || r2 > cut2) continue;

                    double inv_r2  = 1.0 / r2;
                    double inv_r6  = inv_r2 * inv_r2 * inv_r2;
                    double inv_r8  = inv_r6 * inv_r2;
                    double inv_r14 = inv_r8 * inv_r6;

                    double f = c12_inter * inv_r14 - c6_inter * inv_r8;
                    if (f > fmax_inter) f = fmax_inter;

                    double fx = f * dx;
                    double fy = f * dy;
                    double fz = f * dz;

                    F_rnap_all[ri][si][0] += fx;
                    F_rnap_all[ri][si][1] += fy;
                    F_rnap_all[ri][si][2] += fz;

                    F_rnap_all[rj][sj][0] -= fx;
                    F_rnap_all[rj][sj][1] -= fy;
                    F_rnap_all[rj][sj][2] -= fz;
                }
            }

            /* Intra-RNAP (sj > si) */
            for (int sj = si + 1; sj < nsub; sj++) {
                double *Rj = R_rnap[ri][sj];

                double dx = Ri[0] - Rj[0];
                double dy = Ri[1] - Rj[1];
                double dz = Ri[2] - Rj[2];

                double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 == 0.0 || r2 > cut2) continue;

                double inv_r2  = 1.0 / r2;
                double inv_r6  = inv_r2 * inv_r2 * inv_r2;
                double inv_r8  = inv_r6 * inv_r2;
                double inv_r14 = inv_r8 * inv_r6;

                double f = c12_intra * inv_r14 - c6_intra * inv_r8;
                if (f > fmax_intra) f = fmax_intra;

                double fx = f * dx;
                double fy = f * dy;
                double fz = f * dz;

                F_rnap_all[ri][si][0] += fx;
                F_rnap_all[ri][si][1] += fy;
                F_rnap_all[ri][si][2] += fz;

                F_rnap_all[ri][sj][0] -= fx;
                F_rnap_all[ri][sj][1] -= fy;
                F_rnap_all[ri][sj][2] -= fz;
            }
        }
    }
}

/* ================================================================
 * MISE À JOUR EULER-MARUYAMA
 * ================================================================ */
void euler_maruyama_update(double **R, double (*F)[3], int N,
                           double Delta, int temperature,
                           int attache, int plan)
{
    const double sigma_th = temperature ? sqrt(2.0 * Delta) : 0.0;

    for (int i = 0; i < N; i++) {

        if (attache == 1 && (i == 0 || i == N-1))
            continue;

        R[i][0] += Delta * F[i][0];
        R[i][1] += Delta * F[i][1];
        R[i][2] += Delta * F[i][2];

        if (temperature) {
            R[i][0] += sigma_th * randn();
            R[i][1] += sigma_th * randn();
            R[i][2] += sigma_th * randn();
        }

        if (plan && i > 0 && i < N-1 && R[i][2] < 0.0)
            R[i][2] = -R[i][2];
    }
}


/* ================================================================
 * Mise à jour Euler-Maruyama — RNAP (tableau plat)
 * ================================================================ */
void euler_maruyama_update_flat(double **R_rnap, double (*F)[3], int nsub,
                                double Delta, int temperature)
{
    const double sigma_th = temperature ? sqrt(2.0 * Delta) : 0.0;

    for (int s = 0; s < nsub; s++) {
        R_rnap[s][0] += Delta * F[s][0];
        R_rnap[s][1] += Delta * F[s][1];
        R_rnap[s][2] += Delta * F[s][2];

        if (temperature) {
            R_rnap[s][0] += sigma_th * randn();
            R_rnap[s][1] += sigma_th * randn();
            R_rnap[s][2] += sigma_th * randn();
        }
    }
}