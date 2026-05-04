#ifndef FORCES_H
#define FORCES_H

/*
 * forces.h
 * --------
 * Schéma Euler-Maruyama correct :
 *
 *   1. Accumuler toutes les forces sur R(t) dans F
 *   2. R(t+dt) = R(t) + dt*F + sqrt(2dt)*eta   (une seule passe)
 *
 * Chromatine :
 *   accumulate_spring_forces()    → F_chrom
 *   accumulate_lj_forces()        → F_chrom
 *   accumulate_conf_forces()      → F_chrom
 *
 * RNAP (appelé pour chaque RNAP active dans la boucle) :
 *   accumulate_ring_forces()      → F_rnap
 *   accumulate_diagonal_forces()  → F_rnap
 *   accumulate_bond_rnap_chrom()  → F_rnap + F_chrom (réaction)
 *   accumulate_lj_rnap_chrom()    → F_rnap + F_chrom (réaction)
 *   accumulate_lj_rnap_rnap()     → F_rnap de toutes les RNAP actives
 *
 * Mise à jour :
 *   euler_maruyama_update()       → pour R et R_rnap
 */

#include <stdio.h>
#include "config.h"
#include "neighborlist.h"

/* ================================================================
 * CHROMATINE
 * ================================================================ */

/* Forces de ressort harmonique (chaîne linéaire)
 * Accumule dans F[N][3] — ne modifie pas R. */
void accumulate_spring_forces(double **R, double (*F)[3],
                              double K, int N);

/* Forces de Lennard-Jones chromatine-chromatine
 * Accumule dans F[N][3] — ne modifie pas R. */
void accumulate_lj_forces(double **R, double (*F)[3],
                          NeighborList *neighbor_lists, int N,
                          double epsilon, double sigma6, double sigma12,
                          int attache);

/* Force de confinement sphérique (potentiel LJ répulsif)
 * Accumule dans F[N][3] — ne modifie pas R. */
void accumulate_conf_forces(double **R, double (*F)[3], int N,
                            const double cdm[3],
                            double r_conf, double sigma_conf,
                            double epsilon_conf);

/* ================================================================
 * RNAP — forces internes (anneau + diagonales)
 * ================================================================ */

/* Forces de ressort en anneau entre sous-unités RNAP consécutives
 * (topologie circulaire : 0-1-2-...-N-1-0)
 * Accumule dans F_rnap[nsub][3] — ne modifie pas R_rnap. */
void accumulate_ring_forces(double **R_rnap, double (*F_rnap)[3],
                            double K_rnap, double alpha, int nsub);

/* Liaisons diagonales hardcodées pour nsub=8 : (0,4),(1,5),(2,6),(3,7)
 * Accumule dans F_rnap[nsub][3] — ne modifie pas R_rnap. */
void accumulate_diagonal_forces(double **R_rnap, double (*F_rnap)[3],
                                double K_rnap, double a_transpt, int nsub);

/* ================================================================
 * RNAP — couplage RNAP ↔ chromatine
 * ================================================================ */

/* Couplage progressif RNAP↔chromatine (bond_rnap_bead_progressive_mvt)
 * - Accumule sur F_rnap[nsub][3]  (force sur la RNAP)
 * - Accumule sur F_chrom[N][3]    (réaction sur la chromatine)
 * Ne modifie ni R ni R_rnap. */
void accumulate_bond_rnap_chrom(const Config *cfg,
                                double **R,
                                double **R_rnap,
                                double (*F_chrom)[3],
                                double (*F_rnap)[3],
                                int mono,
                                double avancement);

/* LJ répulsif RNAP↔chromatine (via neighbor list)
 * - Accumule sur F_rnap[nsub][3]
 * - Accumule sur F_chrom[N][3]    (réaction)
 * Ne modifie ni R ni R_rnap. */
void accumulate_lj_rnap_chrom(const Config *cfg,
                               double **R_rnap,
                               double **R,
                               double (*F_rnap)[3],
                               double (*F_chrom)[3],
                               NeighborList_rnap **neighbor_lists,
                               int rnap_id,
                               double epsilon,
                               double sigma6, double sigma12,
                               double cut2);

/* ================================================================
 * RNAP — interactions RNAP↔RNAP
 * ================================================================ */

/* LJ entre sous-unités de RNAPs différentes + auto-interactions intra-RNAP
 * Accumule dans F_rnap_all[MAX_RNAP][nsub][3]
 * (tableau 3D alloué dans SimVars, indexé par [rnap][sub][dim])
 * Ne modifie pas R_rnap. */
void accumulate_lj_rnap_rnap(const Config *cfg,
                              double ***R_rnap,
                              double (**F_rnap_all)[3],
                              int nb_rnap,
                              int *l_rnap,
                              double epsilon,
                              double sigma6, double sigma12,
                              double cut2);

/* ================================================================
 * MISE À JOUR EULER-MARUYAMA
 * ================================================================ */

/* R[i](t+dt) = R[i](t) + dt*F[i] + sqrt(2dt)*eta[i]
 * Fonctionne pour la chromatine (N monomères) ET pour une RNAP (nsub sous-unités).
 * - attache == 1 : extrémités i=0 et i=N-1 fixes (chromatine uniquement)
 * - plan   == 1 : réflexion si R[i][2] < 0 */
void euler_maruyama_update(double **R, double (*F)[3], int N,
                           double Delta, int temperature,
                           int attache, int plan, double gamma_fric);

void euler_maruyama_update_flat(double **R_rnap, double (*F)[3], int nsub,
                                double Delta, int temperature, double gamma_fric);

#endif /* FORCES_H */
