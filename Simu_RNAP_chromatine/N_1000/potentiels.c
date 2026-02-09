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


 // Rayon de coupure pour l'interaction de Lennard-Jones
 // Marge pour réduire la fréquence de mise à jour des listes de voisinage

// Structure pour représenter la liste de voisins pour chaque particule


// void lennard_jones_forces(double **R, NeighborList *neighbor_lists, int N,
//                           double epsilon, double sigma6, double sigma12,
//                           double Delta, int attache,
//                           int periode_enregistrement_force,
//                           FILE* fichier_force_LJ, int t) {

//     double epsilon_att = epsilon; 
//     double epsilon_rep = epsilon; 
//     double deplacement = 0.0; 
//     double e = 300;

//     double total_force = 0.0;
//     int nb_paires = 0;

//     if (t % periode_enregistrement_force == 0) {
//         fprintf(fichier_force_LJ, "TIMESTEP: %d\n", t);
//     }

//     for (int i = 1; i < N - 1; i++) {
//         if (t % periode_enregistrement_force == 0) {
//             fprintf(fichier_force_LJ, "Particule: %d\n", i);
//         }

//         for (int k = 0; k < neighbor_lists[i].count; k++) {
//             int j = neighbor_lists[i].neighbors[k];
//             if (j < i ){continue;}

//             if (t % periode_enregistrement_force == 0) {
//                 fprintf(fichier_force_LJ, "NEIGHBOR: %d ", j);
//             }

//             double dx = R[i][0] - R[j][0];
//             double dy = R[i][1] - R[j][1];
//             double dz = R[i][2] - R[j][2];

//             double r2 = dx * dx + dy * dy + dz * dz;
//             double r8 = r2 * r2 * r2 * r2;
//             double r14 = r8 * r2 * r2 * r2;

//             double f = 4 * (12 * epsilon_rep * sigma12 / r14 - 6 * epsilon_att * sigma6 / r8);
//             if (f > e) f = e;

//             double fx = f * dx;
//             double fy = f * dy;
//             double fz = f * dz;
//             double norm_f = sqrt(fx * fx + fy * fy + fz * fz);

//             total_force += norm_f;
//             nb_paires++;

//             // Application des forces
//             if (attache == 1) {
//                 if (i != 0 && i != N - 1) {
//                     R[i][0] += Delta * fx;
//                     R[i][1] += Delta * fy;
//                     R[i][2] += Delta * fz;
//                 }
//                 if (j != 0 && j != N - 1) {
//                     R[j][0] -= Delta * fx;
//                     R[j][1] -= Delta * fy;
//                     R[j][2] -= Delta * fz;
//                 }
//             }

//             if (attache == 0) {
//                 if (t % periode_enregistrement_force == 0) fprintf(fichier_force_LJ, "%f ", fx);
//                 deplacement = -R[i][0];
//                 R[i][0] += Delta * fx;
//                 deplacement += R[i][0];
//                 if (t % periode_enregistrement_force == 0) fprintf(fichier_force_LJ, "%f ", deplacement);

//                 if (t % periode_enregistrement_force == 0) fprintf(fichier_force_LJ, "%f ", fy);
//                 deplacement = -R[i][1];
//                 R[i][1] += Delta * fy;
//                 deplacement += R[i][1];
//                 if (t % periode_enregistrement_force == 0) fprintf(fichier_force_LJ, "%f ", deplacement);

//                 if (t % periode_enregistrement_force == 0) fprintf(fichier_force_LJ, "%f ", fz);
//                 deplacement = -R[i][2];
//                 R[i][2] += Delta * fz;
//                 deplacement += R[i][2];
//                 if (t % periode_enregistrement_force == 0) fprintf(fichier_force_LJ, "%f ", deplacement);

//                 R[j][0] -= Delta * fx;
//                 R[j][1] -= Delta * fy;
//                 R[j][2] -= Delta * fz;
//             }

//             if (t % periode_enregistrement_force == 0) {
//                 fprintf(fichier_force_LJ, "\n");
//             }
//         }
//     }

//     if (t % periode_enregistrement_force == 0 && nb_paires > 0) {
//         double force_moyenne = total_force / nb_paires;
//         fprintf(fichier_force_LJ, "# Force moyenne LJ : %lf\n", force_moyenne);
//     }
// }

void lennard_jones_forces(
    double **R,
    NeighborList *neighbor_lists,
    int N,
    double epsilon,   
    double sigma6,
    double sigma12,
    double Delta,
    int attache, // 1: ne bouge pas les extrémités, 0: bouge tout
    int periode_enregistrement_force,
    FILE *fichier_force_LJ,
    int t
) {
    double epsilon_rep = 2 * epsilon; 
    double epsilon_att = 0.2 * epsilon;

    const double r2_min = 1e-12;      // évite division par 0
    const double fcoef_max = 1000.0;   // clamp symétrique sur le coef scalaire

    // Accumulation des forces
    double (*F)[3] = calloc((size_t)N, sizeof(*F));
    if (!F) {
        fprintf(stderr, "Erreur: allocation forces LJ impossible\n");
        return;
    }

    double total_pair_force = 0.0;
    int nb_paires = 0;

    const int do_log = (periode_enregistrement_force > 0) &&
                       (t % periode_enregistrement_force == 0);

    if (do_log) {
        fprintf(fichier_force_LJ, "TIMESTEP: %d\n", t);
    }

    // Boucle sur toutes les particules (inclut 0 et N-1)
    for (int i = 0; i < N; i++) {
        if (do_log) {
            fprintf(fichier_force_LJ, "Particule: %d\n", i);
        }

        for (int k = 0; k < neighbor_lists[i].count; k++) {
            int j = neighbor_lists[i].neighbors[k];

            // éviter double comptage
            if (j <= i) continue;
            if (j < 0 || j >= N) continue; // sécurité

            if (do_log) {
                fprintf(fichier_force_LJ, "NEIGHBOR: %d ", j);
            }

            double dx = R[i][0] - R[j][0];
            double dy = R[i][1] - R[j][1];
            double dz = R[i][2] - R[j][2];

            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < r2_min) {
                if (do_log) fprintf(fichier_force_LJ, "SKIP r2~0\n");
                continue;
            }

            // r^8 et r^14 (comme dans ta version)
            double r4  = r2 * r2;
            double r8  = r4 * r4;          // r^8
            double r14 = r8 * r4 * r2;     // r^(8+4+2) = r^14

            // coef scalaire tel que F_vec = fcoef * r_vec
            // F = 4 * [ 12 eps_rep sigma^12 / r^14 - 6 eps_att sigma^6 / r^8 ] * r_vec
            double fcoef = 4.0 * (12.0 * epsilon_rep * sigma12 / r14
                                - 6.0  * epsilon_att * sigma6  / r8);

            // clamp symétrique
            if (fcoef >  fcoef_max) fcoef =  fcoef_max;
            if (fcoef < -fcoef_max) fcoef = -fcoef_max;

            double fx = fcoef * dx;
            double fy = fcoef * dy;
            double fz = fcoef * dz;

            double norm_f = sqrt(fx*fx + fy*fy + fz*fz);
            total_pair_force += norm_f;
            nb_paires++;

            // accumulation action-réaction
            F[i][0] += fx;  F[i][1] += fy;  F[i][2] += fz;
            F[j][0] -= fx;  F[j][1] -= fy;  F[j][2] -= fz;

            if (do_log) {
                // tu peux logger ce que tu veux ici (force paire)
                fprintf(fichier_force_LJ, "fx=%g fy=%g fz=%g |F|=%g\n", fx, fy, fz, norm_f);
            }
        }
    }

    // Application des déplacements (Euler explicite)
    for (int i = 0; i < N; i++) {
        if (attache == 1 && (i == 0 || i == N-1)) continue;

        double oldx = R[i][0], oldy = R[i][1], oldz = R[i][2];
        R[i][0] += Delta * F[i][0];
        R[i][1] += Delta * F[i][1];
        R[i][2] += Delta * F[i][2];

        if (do_log) {
            fprintf(fichier_force_LJ,
                    "MOVE i=%d dX=%g dY=%g dZ=%g\n",
                    i, R[i][0]-oldx, R[i][1]-oldy, R[i][2]-oldz);
        }
    }

    if (do_log && nb_paires > 0) {
        fprintf(fichier_force_LJ, "# Force moyenne LJ (par paire) : %g\n",
                total_pair_force / nb_paires);
    }

    free(F);
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



