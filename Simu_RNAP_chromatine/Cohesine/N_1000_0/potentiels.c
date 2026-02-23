#include "potentiels.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "basic_functions.h"
#include "neighborlist.h" 
double ForceLJ(double sigma, double epsilon, double x, double d) {
    double epsilon_rep = epsilon;
    double epsilon_att = 2 * epsilon;
    double y = 4  * (-12 * epsilon_rep * pow(sigma, 12) * x / pow(d, 14) + 6 * epsilon_att * pow(sigma, 6) * x / pow(d, 8));
    double e = 500;
    if (y > e) {y = e;}
    if (y < -e) {y = -e;}
    return y;}




void potentiel_barriere(double** R, double* origine, double rayon, double force_dep, double epaisseur, int N) {
    for (int part = 0; part < N; ++part) {
        double distance_origine = sqrt(pow(R[part][0] - origine[0], 2) + pow(R[part][1] - origine[1], 2) );

        if (distance_origine > rayon && fabs(R[part][2] - origine[2]) < epaisseur) {
            double D = distance_origine;
            R[part][0] += force_dep / D * (origine[0] - R[part][0]);
            R[part][1] += force_dep / D * (origine[1] - R[part][1]);
            R[part][2] += force_dep / D * (origine[2] - R[part][2]);
        }
    }
}


void force_bille_bille(double** R1,double** R2, double K_cohesine, double dist, int particule1, int particule2, int distance_cohesine_eq, double dt){
        /*R[particule1][0] += norme_force / distance * (R[particule2][0] - R[particule1][0]);
        R[particule1][1] += norme_force / distance * (R[particule2][1] - R[particule1][1]);
        R[particule1][2] += norme_force / distance * (R[particule2][2] - R[particule1][2]);
        R[particule2][0] += - norme_force / distance * (R[particule2][0] - R[particule1][0]);
        R[particule2][1] += - norme_force / distance * (R[particule2][1] - R[particule1][1]);
        R[particule2][2] += - norme_force / distance * (R[particule2][2] - R[particule1][2]);*/
        for (int j = 0; j < 3; j++) {
             R1[particule1][j] += dt * K_cohesine * (R2[particule2][j] - R1[particule1][j]) / dist * (dist - distance_cohesine_eq);
             R2[particule2][j] += - dt * K_cohesine * (R2[particule2][j] - R1[particule1][j]) / dist * (dist - distance_cohesine_eq);
                                    }}


static inline int is_chain_endpoint(int idx, int Nchain) {
    return (idx == 0) || (idx == Nchain - 1) || (idx == Nchain) || (idx == 2*Nchain - 1);
}


void lennard_jones_forces(double **R, NeighborList *neighbor_lists, int N,
                          double epsilon, double sigma6, double sigma12,
                          double Delta, int attache,
                          int periode_enregistrement_force,
                          FILE* fichier_force_LJ, int t) {

    double epsilon_att = epsilon; 
    double epsilon_rep = epsilon; 
    double deplacement = 0.0; 
    double e = 300;

    double total_force = 0.0;
    int nb_paires = 0;

    if (t % periode_enregistrement_force == 0) {
        fprintf(fichier_force_LJ, "TIMESTEP: %d\n", t);
    }

    for (int i = 1; i < N - 1; i++) {
        if (t % periode_enregistrement_force == 0) {
            fprintf(fichier_force_LJ, "Particule: %d\n", i);
        }

        for (int k = 0; k < neighbor_lists[i].count; k++) {
            int j = neighbor_lists[i].neighbors[k];
            if (j <= i) continue;

            if (t % periode_enregistrement_force == 0) {
                fprintf(fichier_force_LJ, "NEIGHBOR: %d ", j);
            }

            double dx = R[i][0] - R[j][0];
            double dy = R[i][1] - R[j][1];
            double dz = R[i][2] - R[j][2];

            double r2 = dx * dx + dy * dy + dz * dz;
            double r8 = r2 * r2 * r2 * r2;
            double r14 = r8 * r2 * r2 * r2;

            double f = 4 * (12 * epsilon_rep * sigma12 / r14 - 6 * epsilon_att * sigma6 / r8);
            if (f > e) f = e;

            double fx = f * dx;
            double fy = f * dy;
            double fz = f * dz;
            double norm_f = sqrt(fx * fx + fy * fy + fz * fz);

            total_force += norm_f;
            nb_paires++;

            // Application des forces
            if (attache == 1) {
                if (i != 0 && i != N - 1) {
                    R[i][0] += Delta * fx;
                    R[i][1] += Delta * fy;
                    R[i][2] += Delta * fz;
                }
                if (j != 0 && j != N - 1) {
                    R[j][0] -= Delta * fx;
                    R[j][1] -= Delta * fy;
                    R[j][2] -= Delta * fz;
                }
            }

            if (attache == 0) {
                if (t % periode_enregistrement_force == 0) fprintf(fichier_force_LJ, "%f ", fx);
                deplacement = -R[i][0];
                R[i][0] += Delta * fx;
                deplacement += R[i][0];
                if (t % periode_enregistrement_force == 0) fprintf(fichier_force_LJ, "%f ", deplacement);

                if (t % periode_enregistrement_force == 0) fprintf(fichier_force_LJ, "%f ", fy);
                deplacement = -R[i][1];
                R[i][1] += Delta * fy;
                deplacement += R[i][1];
                if (t % periode_enregistrement_force == 0) fprintf(fichier_force_LJ, "%f ", deplacement);

                if (t % periode_enregistrement_force == 0) fprintf(fichier_force_LJ, "%f ", fz);
                deplacement = -R[i][2];
                R[i][2] += Delta * fz;
                deplacement += R[i][2];
                if (t % periode_enregistrement_force == 0) fprintf(fichier_force_LJ, "%f ", deplacement);

                R[j][0] -= Delta * fx;
                R[j][1] -= Delta * fy;
                R[j][2] -= Delta * fz;
            }

            if (t % periode_enregistrement_force == 0) {
                fprintf(fichier_force_LJ, "\n");
            }
        }
    }

    if (t % periode_enregistrement_force == 0 && nb_paires > 0) {
        double force_moyenne = total_force / nb_paires;
        fprintf(fichier_force_LJ, "# Force moyenne LJ : %lf\n", force_moyenne);
    }
}

void lennard_jones_forces_2chains(double **R, NeighborList *neighbor_lists,
                                  int Nchain,  // N par chaîne
                                  double epsilon, double sigma6, double sigma12,
                                  double Delta, int attache,
                                  int periode_enregistrement_force,
                                  FILE* fichier_force_LJ, int t)
{
    const int Ntot = 2 * Nchain;
    double epsilon_att = epsilon;
    double epsilon_rep = epsilon;
    double deplacement = 0.0;
    double emax = 300.0;

    double total_force = 0.0;
    int nb_paires = 0;

    if (t % periode_enregistrement_force == 0) {
        fprintf(fichier_force_LJ, "TIMESTEP: %d\n", t);
    }

    // ⚠️ boucle sur TOUTES les particules
    for (int i = 0; i < Ntot; i++) {

        if (t % periode_enregistrement_force == 0) {
            fprintf(fichier_force_LJ, "Particule: %d\n", i);
        }

        for (int k = 0; k < neighbor_lists[i].count; k++) {
            int j = neighbor_lists[i].neighbors[k];
            if (j <= i) continue; // éviter double comptage

            double dx = R[i][0] - R[j][0];
            double dy = R[i][1] - R[j][1];
            double dz = R[i][2] - R[j][2];

            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < 1e-12) continue; // sécurité (évite /0)

            double r4  = r2*r2;
            double r8  = r4*r4;
            double r14 = r8*r4*r2;   // r^(8+4+2)=r^14

            double f = 4.0 * (12.0 * epsilon_rep * sigma12 / r14
                            - 6.0  * epsilon_att * sigma6  / r8);
            if (f > emax) f = emax;
            if (f < -emax) f = -emax;

            double fx = f * dx;
            double fy = f * dy;
            double fz = f * dz;

            total_force += sqrt(fx*fx + fy*fy + fz*fz);
            nb_paires++;

            // Application forces (Euler overdamped)
            if (attache == 1) {
                if (!is_chain_endpoint(i, Nchain)) {
                    R[i][0] += Delta * fx;
                    R[i][1] += Delta * fy;
                    R[i][2] += Delta * fz;
                }
                if (!is_chain_endpoint(j, Nchain)) {
                    R[j][0] -= Delta * fx;
                    R[j][1] -= Delta * fy;
                    R[j][2] -= Delta * fz;
                }
            } else { // attache == 0 : tout bouge
                // tes prints/debug peuvent rester
                R[i][0] += Delta * fx;  R[i][1] += Delta * fy;  R[i][2] += Delta * fz;
                R[j][0] -= Delta * fx;  R[j][1] -= Delta * fy;  R[j][2] -= Delta * fz;
            }

            if (t % periode_enregistrement_force == 0) {
                fprintf(fichier_force_LJ, "NEIGHBOR: %d %f %f %f\n", j, fx, fy, fz);
            }
        }
    }

    if (t % periode_enregistrement_force == 0 && nb_paires > 0) {
        fprintf(fichier_force_LJ, "# Force moyenne LJ : %lf\n", total_force / nb_paires);
    }
}


void f_bending_forces(double **R, double **t_link, double **bending_forces, double K_bend, int N, int t){
    double ddum1, ddum2, ddum3, ddum4, ddum5, ddum6, angle;
    for (int i = 0; i < N; i++){
        if ( i > 2 && i < N-2){
            ddum1 = norm(t_link[i]) * norm(t_link[i+1]);
            ddum2 = norm(t_link[i-1]) * norm(t_link[i]);
            ddum3 = norm(t_link[i-2]) * norm(t_link[i-1]);
            ddum4 = dot_product(t_link[i], t_link[i+1]);
            ddum5 = dot_product(t_link[i-1], t_link[i]);
            ddum6 = dot_product(t_link[i-2], t_link[i-1]);
            for(int j = 0; j < 3; j++){
                bending_forces[i][j] = K_bend * ((-t_link[i-2][j] + (t_link[i-1][j]*ddum6)/norm2(t_link[i-1]))/ddum3+((t_link[i-1][j])*(ddum5/norm2(t_link[i-1])+1)-(t_link[i][j])*(ddum5/norm2(t_link[i])+1))/ddum2+(-t_link[i][j]*ddum4/norm2(t_link[i])+t_link[i+1][j])/ddum1);
            }
            

        }
            
            else if(i == 0){
                ddum3 = norm(t_link[i]) * norm(t_link[i+1]);
                ddum6 = dot_product(t_link[i], t_link[i+1]);

                angle = calculate_angle(t_link[i], t_link[i+1]);

                for(int j = 0; j < 3; j++){
                    bending_forces[i][j] = K_bend * (cos(angle) - sin(angle)/(sqrt(ddum3*ddum3/(ddum6*ddum6)-1)))*(t_link[i+1][j] + t_link[i][j]*ddum6/norm2(t_link[i]))/ddum3;
                    //printf("%d %d %f ",t,i,sqrt(ddum3*ddum3/(ddum6*ddum6)-1));
                }
                //printf("\n");
            }

            else if(i==1){
                ddum2= norm(t_link[i-1])*norm(t_link[i]);
                ddum3= norm(t_link[i])*norm(t_link[i+1]);
                ddum5= dot_product(t_link[i-1],t_link[i]);
                ddum6= dot_product(t_link[i],t_link[i+1]);
                angle = calculate_angle(t_link[i], t_link[i+1]);

                for(int j = 0; j < 3; j++){
                    bending_forces[i][j] = K_bend * ( cos(angle) - sin(angle)/sqrt(ddum2*ddum2/(ddum5*ddum5)-1) * (t_link[i-1][j]*(ddum5/norm2(t_link[i-1])+1)-t_link[i][j]*(ddum5/norm2(t_link[i])+1))/ddum2 + (-t_link[i][j]*ddum6/norm2(t_link[i])+t_link[i+1][j])/ddum3);
                    
                }

            }

            else if(i == 2){
                ddum1 = norm(t_link[2])*norm(t_link[3]);
                ddum2 = norm(t_link[1])*norm(t_link[2]);
                ddum3 = norm(t_link[0])*norm(t_link[1]); 
                ddum4 = dot_product(t_link[2],t_link[3]);
                ddum5 = dot_product(t_link[1],t_link[2]);
                ddum6 = dot_product(t_link[0],t_link[1]); 

                angle = calculate_angle(t_link[2], t_link[3]);
                
                
                for(int j = 0; j < 3; j++){
                    bending_forces[i][j] = K_bend * (cos(angle) - sin(angle)/sqrt(ddum3*ddum3/(ddum6*ddum6)-1) * (-t_link[0][j]+t_link[1][j]*ddum6/norm2(t_link[1]))/ddum3 + (t_link[1][j]*(ddum5/norm2(t_link[1])+1)-t_link[2][j]*(ddum5/norm2(t_link[2])+1))/ddum2 + (-t_link[2][j]*ddum4/norm2(t_link[2])+t_link[3][j])/ddum1);
            
                }


            }

            else if(i == N-2){
                ddum1 = norm(t_link[N-4])*norm(t_link[N-3]);
                ddum2 = norm(t_link[N-2])*norm(t_link[N-3]);
                ddum4 = dot_product(t_link[N-4],t_link[N-3]);
                ddum5 = dot_product(t_link[N-3],t_link[N-2]);

                for(int j = 0; j < 3; j++){
                    bending_forces[i][j] =  -K_bend*((-t_link[N-4][j]+t_link[N-3][j]*ddum4/norm2(t_link[N-3]))/ddum1 + (t_link[N-3][j]*(ddum5*(norm2(t_link[N-3])+1))-t_link[N-2][j]*(ddum5/norm2(t_link[N-2])+1))/ddum2);
                    
                }
            }

            else if( i == N-1){
                ddum1 = norm(t_link[N-2])*norm(t_link[N-3]);
                ddum4 = dot_product(t_link[N-2],t_link[N-3]);
                
                for(int j = 0; j < 3; j++){
                    bending_forces[i][j] = K_bend * (-t_link[N-3][j]+t_link[N-2][j]*ddum4/norm2(t_link[N-2]))/ddum1;
                }
            }


        }
    }

    #include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

static inline size_t pair_index(int a, int b, int N) {
    return (size_t)a * (size_t)N + (size_t)b;
}

void debug_neighbor_double_count_strict(NeighborList *nl, int N, int t) {
    // counts[a*N + b] = nb fois qu'on a vu l'arête orientée a->b
    uint8_t *counts = (uint8_t*)calloc((size_t)N * (size_t)N, sizeof(uint8_t));
    if (!counts) {
        fprintf(stderr, "debug_neighbor_double_count_strict: alloc failed (N=%d)\n", N);
        return;
    }

    long long directed_edges = 0;
    long long symmetric_pairs = 0;      // nb de paires (i,j) telles que i->j et j->i existent
    long long multi_edges = 0;          // nb d'arêtes vues >1 fois (même orientation)
    long long self_edges = 0;           // i==j

    // 1) Remplir counts et détecter doublons orientés
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < nl[i].count; ++k) {
            int j = nl[i].neighbors[k];
            directed_edges++;

            if (j == i) self_edges++;

            size_t idx = pair_index(i, j, N);
            if (counts[idx] == 255) {
                // saturé
            } else {
                counts[idx]++;
            }
            if (counts[idx] > 1) {
                multi_edges++; // même (i->j) ajouté plusieurs fois -> bug de construction de NL
            }
        }
    }

    // 2) Mesurer la symétrie : pour chaque paire non orientée (i<j)
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            uint8_t cij = counts[pair_index(i, j, N)];
            uint8_t cji = counts[pair_index(j, i, N)];
            if (cij > 0 && cji > 0) symmetric_pairs++;
        }
    }

    fprintf(stderr,
        "[t=%d] STRICT neighbor check:\n"
        "  directed_edges(total i->j) = %lld\n"
        "  symmetric_pairs(i<j with both directions present) = %lld\n"
        "  multi_edges(duplicate same direction) = %lld\n"
        "  self_edges(i==j) = %lld\n",
        t, directed_edges, symmetric_pairs, multi_edges, self_edges
    );

    // Interprétation :
    // - si symmetric_pairs est grand (proche du nb de paires uniques), ta NL est "complète" (i->j et j->i)
    // - si symmetric_pairs ~ 0, ta NL est "demi" (typiquement seulement j>i)
    // - si multi_edges > 0 : bug (même voisin répété)
    // - self_edges > 0 : bug (un point est son propre voisin)

    // Bonus : afficher quelques exemples de doublons orientés
    if (multi_edges > 0) {
        int printed = 0;
        fprintf(stderr, "  Examples of duplicate directed edges (i->j count>1):\n");
        for (int i = 0; i < N && printed < 10; ++i) {
            for (int j = 0; j < N && printed < 10; ++j) {
                uint8_t c = counts[pair_index(i, j, N)];
                if (c > 1) {
                    fprintf(stderr, "    %d -> %d : %u times\n", i, j, (unsigned)c);
                    printed++;
                }
            }
        }
    }

    free(counts);
}



