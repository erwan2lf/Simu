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


void lennard_jones_forces(double **R, NeighborList *neighbor_lists, int N,
                          double epsilon, double sigma6, double sigma12,
                          double Delta, int attache,
                          int periode_enregistrement_force,
                          FILE *fichier_force_LJ, int t)
{
    const double epsilon_att = epsilon;
    const double epsilon_rep = epsilon;
    const double fmax = 300.0;

    const double c12 = 48.0 * epsilon_rep * sigma12;
    const double c6  = 24.0 * epsilon_att * sigma6;

    const int do_log = (periode_enregistrement_force > 0) &&
                       (t % periode_enregistrement_force == 0);

    double total_force = 0.0;
    int nb_paires = 0;

    if (do_log) {
        fprintf(fichier_force_LJ, "TIMESTEP: %d\n", t);
    }

    for (int i = 1; i < N - 1; i++) {
        double *Ri = R[i];

        if (do_log) {
            fprintf(fichier_force_LJ, "Particule: %d\n", i);
        }

        for (int k = 0; k < neighbor_lists[i].count; k++) {
            int j = neighbor_lists[i].neighbors[k];

            // À garder si la liste de voisins est symétrique :
            // if (j <= i) continue;

            double *Rj = R[j];

            if (do_log) {
                fprintf(fichier_force_LJ, "NEIGHBOR: %d ", j);
            }

            double dx = Ri[0] - Rj[0];
            double dy = Ri[1] - Rj[1];
            double dz = Ri[2] - Rj[2];

            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 == 0.0) continue;  // sécurité

            double inv_r2  = 1.0 / r2;
            double inv_r4  = inv_r2 * inv_r2;
            double inv_r6  = inv_r4 * inv_r2;
            double inv_r8  = inv_r6 * inv_r2;
            double inv_r14 = inv_r8 * inv_r4 * inv_r2;

            double f = c12 * inv_r14 - c6 * inv_r8;
            if (f > fmax) f = fmax;

            double fx = f * dx;
            double fy = f * dy;
            double fz = f * dz;

            if (do_log) {
                double norm_f = fabs(f) * sqrt(r2);
                total_force += norm_f;
                nb_paires++;
            }

            if (attache == 1) {
                if (i != 0 && i != N - 1) {
                    Ri[0] += Delta * fx;
                    Ri[1] += Delta * fy;
                    Ri[2] += Delta * fz;
                }
                if (j != 0 && j != N - 1) {
                    Rj[0] -= Delta * fx;
                    Rj[1] -= Delta * fy;
                    Rj[2] -= Delta * fz;
                }
            } else {
                if (do_log) fprintf(fichier_force_LJ, "%f ", fx);
                double deplacement = Delta * fx;
                Ri[0] += deplacement;
                if (do_log) fprintf(fichier_force_LJ, "%f ", deplacement);

                if (do_log) fprintf(fichier_force_LJ, "%f ", fy);
                deplacement = Delta * fy;
                Ri[1] += deplacement;
                if (do_log) fprintf(fichier_force_LJ, "%f ", deplacement);

                if (do_log) fprintf(fichier_force_LJ, "%f ", fz);
                deplacement = Delta * fz;
                Ri[2] += deplacement;
                if (do_log) fprintf(fichier_force_LJ, "%f ", deplacement);

                Rj[0] -= Delta * fx;
                Rj[1] -= Delta * fy;
                Rj[2] -= Delta * fz;
            }

            if (do_log) {
                fprintf(fichier_force_LJ, "\n");
            }
        }
    }

    if (do_log && nb_paires > 0) {
        double force_moyenne = total_force / nb_paires;
        fprintf(fichier_force_LJ, "# Force moyenne LJ : %lf\n", force_moyenne);
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



