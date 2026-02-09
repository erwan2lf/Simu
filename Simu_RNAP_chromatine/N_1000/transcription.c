#include "transcription.h"
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <sys/stat.h>
#include <errno.h>
#include "structures_depart.h"
#include "basic_functions.h"
#include "potentiels.h"
#include "movement.h"
#include "config.h"
#include "simulation.h"
#include "neighborlist.h"   
#include "file.h"


#define BOX_SIZE 1e3
#define SCALE_POS 1000.0
#define PI 3.14159265358979323846
#define MAX_RNAP 50
#define MAX_RNAP_SUBUNITS 8 
#define DIM 3



double*** allocate_matrix_3D(int rows, int cols)
{
/*******************************************************************************************************
 * 📦 allocate_3D_matrix
 * -----------------------------------------------------------------------------------------------------
 * • Description :
 *   Alloue dynamiquement une matrice tridimensionnelle de type `double***`.
 *   Cette fonction est générique et peut être utilisée pour représenter n’importe quelle structure 3D :
 *   coordonnées spatiales (x, y, z), tenseurs, champs de forces, vitesses, etc.
 *
 *   La structure résultante correspond à :
 *       matrix[i][j][k]
 *   où :
 *     - `i` ∈ [0, dim1 - 1]  → premier axe  (ex. : index du polymère ou RNAP)
 *     - `j` ∈ [0, dim2 - 1]  → deuxième axe (ex. : sous-unité, monomère, ou segment)
 *     - `k` ∈ [0, dim3 - 1]  → troisième axe (ex. : coordonnée x/y/z, ou composante vectorielle)
 *
 * • Arguments :
 *   - dim1 : taille du premier axe
 *   - dim2 : taille du deuxième axe
 *   - dim3 : taille du troisième axe
 *
 * • Valeur de retour :
 *   - Retourne un pointeur `double***` vers la matrice allouée.
 *
 * • Remarques :
 *   - Chaque niveau est alloué indépendamment via `malloc()`.
 *   - Les valeurs ne sont **pas initialisées** (utiliser `calloc` si nécessaire).
 *   - Doit être libérée avec la fonction complémentaire `free_3D_matrix()`.
 *
 * • Exemple d’utilisation :
 *     double*** positions = allocate_3D_matrix(N_particles, N_sites, 3);
 *     positions[i][j][0] = 1.0;   // x
 *     positions[i][j][1] = 2.0;   // y
 *     positions[i][j][2] = 3.0;   // z
 *
 *******************************************************************************************************/
    double*** matrix = (double***)malloc(rows * sizeof(double**));
    for (int i = 0; i < rows; i++)
    {
        matrix[i] = (double**)malloc(cols * sizeof(double*));
        for (int j=0; j<cols ; j++)
        {
            matrix[i][j] = (double*)malloc(3 * sizeof(double));
        }
    }
    return matrix;
}

/*******************************************************************************************************
 * 🧹 free_3D_matrix
 * -----------------------------------------------------------------------------------------------------
 * Libère une matrice tridimensionnelle allouée avec `allocate_3D_matrix()`.
 *******************************************************************************************************/
void free_matrix_3D(double ***matrix, int rows, int cols)
{
    if (matrix != NULL)
    {
        for (int i = 0; i < rows; i++)
        {
            if (matrix[i] != NULL)
            {
                for (int j = 0; j < cols; j++)
                {
                    if (matrix[i][j] != NULL)
                    {
                        free(matrix[i][j]);
                        matrix[i][j] = NULL; // Prevent double free
                    }
                }
                free(matrix[i]);
                matrix[i] = NULL; // Prevent double free
            }
        }
        free(matrix);
        matrix = NULL; // Prevent double free
    }
}

/**
 * 🧾 print_position
 * ------------------------------------------------------------------------------
 * Écrit une ligne au format LAMMPS :
 *   id type xs ys zs
 *
 * Paramètres :
 *   - f : fichier ouvert (ex. .lammpstrj)
 *   - id : identifiant de particule (index unique)
 *   - type : type de particule (1 = monomère, 2 = RNAP, etc.)
 *   - x, y, z : coordonnées réelles en unités du code
 *
 * Les coordonnées sont divisées par SCALE_POS pour conversion d’échelle.
 * ------------------------------------------------------------------------------
 */
void print_position(FILE *f, int id, int type, double x, double y, double z, int marqueur)
{
    if (!f)
        return;
    fprintf(f, "%d %d %lf %lf %lf %d\n", id, type, x / SCALE_POS, y / SCALE_POS, z / SCALE_POS, marqueur);
}


/**
 * 🧬 enregistrement_RNAP (version avec gestion de l_rnap == -2)
 * ------------------------------------------------------------------------------------------------
 * Enregistre les coordonnées de la chromatine + RNAPs dans un format compatible LAMMPS.
 *
 * Chaque frame contient :
 *   • Les N monomères de la chaîne (type = 1)
 *   • Pour CHAQUE RNAP potentielle (0..nb_rnap_initial-1) :
 *       - l_rnap[i] >= 1 → RNAP active (positions réelles)
 *       - l_rnap[i] == -1 → inactive (jamais entrée, position 0,0,0)
 *       - l_rnap[i] == -2 → sortie définitive, position 0,0,0, type = 99
 *
 * ⚙️ Structure d'une RNAP :
 *   2 monomères sur la chaîne (type = 2)
 *   + nb_subunits sous-unités (type = 3 + i)
 *
 * Le nombre total de particules est constant : N + nb_rnap_initial × (nb_subunits + 2)
 * ------------------------------------------------------------------------------------------------
 */
void enregistrement_RNAP(FILE *fichier,
                        double **R,
                        int N,
                        double ***R_rnap,
                        int nb_rnap,
                        int t,
                        int *positions_bille_rnap,
                        int nb_subunits,
                        int nb_rnap_initial,
                        int *l_rnap, 
                        SimVars *sv,
                        const Config *cfg)
{
    if (!fichier)
    {
        fprintf(stderr, "❌ Erreur : fichier RNAP non ouvert.\n");
        return;
    }

    // === Nombre total constant ===
    int total_atoms = N + nb_rnap_initial * (nb_subunits + 2);

    // === En-têtes LAMMPS ===
    fprintf(fichier, "ITEM: TIMESTEP\n%d\n", t);
    fprintf(fichier, "ITEM: NUMBER OF ATOMS\n%d\n", total_atoms);
    fprintf(fichier, "ITEM: BOX BOUNDS ss ss ss\n");
    fprintf(fichier, "%lf %lf\n", -BOX_SIZE, BOX_SIZE);
    fprintf(fichier, "%lf %lf\n", -BOX_SIZE, BOX_SIZE);
    fprintf(fichier, "%lf %lf\n", -BOX_SIZE, BOX_SIZE);
    fprintf(fichier, "ITEM: ATOMS id type xs ys zs\n");


    // === 1️⃣ Monomères de la chaîne ===
    for (int i = 0; i < N; i++)
    {
        print_position(fichier, i, 1, R[i][0], R[i][1], R[i][2], sv->is_rnap[i]);
    }
    
    // === 2️⃣ RNAPs (actives, inactives, sorties) ===
    int count = N;
    int type = 3;

    for (int rnap = 0; rnap < nb_rnap_initial; rnap++)
    {
        int state = l_rnap[rnap];
        int mono_pos = positions_bille_rnap[rnap];
        if (mono_pos < 0 || mono_pos >= N - 1)
            mono_pos = 0; // sécurité

        // --- RNAP active ---
        if (state >= 1)
        {
            // Monomères liés à la RNAP
            print_position(fichier, count++, 2, R[mono_pos][0], R[mono_pos][1], R[mono_pos][2], 0);
            print_position(fichier, count++, 2, R[mono_pos + 1][0], R[mono_pos + 1][1], R[mono_pos + 1][2], 0);

            // Sous-unités
            for (int s = 0; s < nb_subunits; s++)
                print_position(fichier, count++, type,
                               R_rnap[rnap][s][0],
                               R_rnap[rnap][s][1],
                               R_rnap[rnap][s][2], 0);
        }

        // --- RNAP inactive temporairement ---
        else if (state == -1)
        {
            for (int i = 0; i < nb_subunits + 2; i++)
                print_position(fichier, count++, 2, 0.0, 0.0, 0.0, 0);
        }

        // --- RNAP sortie définitive ---
        else if (state == -2)
        {
            for (int i = 0; i < nb_subunits + 2; i++)
                print_position(fichier, count++, 99, 0.0, 0.0, 0.0, 0);
        }

        type++;
    }

    fflush(fichier); // sécurité
}




void polymere_brownian_motion_ring_force(double** R_rnap, double alpha, double K_rnap, double Delta, int N, int rnap, FILE* test, int t,  int periode_enregistrement_force, FILE* fichier_force_rnap, FILE* fichier_force_thermique, int temperature) {

    int p, pmin, pmax;

    p = 0; pmax = 1; pmin = N-1;
    mvt_brownian_harmonic_force_RNAP(p, pmin, pmax, rnap, R_rnap, alpha, K_rnap, Delta, test, t, periode_enregistrement_force, fichier_force_rnap, fichier_force_thermique, temperature);

    for (int i = 1; i < N - 1; i++) {
        p = i ; pmin= i - 1 ; pmax = i + 1;
        mvt_brownian_harmonic_force_RNAP(p, pmin, pmax, rnap, R_rnap, alpha, K_rnap, Delta, test, t, periode_enregistrement_force, fichier_force_rnap, fichier_force_thermique, temperature);
        }

    p = N - 1; pmin = N - 2; pmax = 0;
    mvt_brownian_harmonic_force_RNAP(p, pmin, pmax, rnap, R_rnap, alpha, K_rnap, Delta, test, t, periode_enregistrement_force, fichier_force_rnap, fichier_force_thermique, temperature); 

}



void mvt_brownian_harmonic_force_RNAP(int p, int pmin, int pmax, int rnap, double** R, double alpha, double K, double Delta, FILE *test, int t, int periode_enregistrement_force, FILE* fichier_force_rnap, FILE* fichier_force_thermique, int temperature){
    double dmin, dmax;
    double F_alea, pot_harm[3];
    dmin = distance(R[pmin], R[p]);
    dmax = distance(R[pmax], R[p]);
    double deplacement = 0.0;

    if(t % periode_enregistrement_force == 0)
    {
        fprintf(fichier_force_rnap, "%d %d ", rnap, p);
        fprintf(fichier_force_thermique, "%d %d ", rnap, p);
    }
    for(int j = 0; j < 3; j ++)
    {
        pot_harm[j] = K * ((R[pmax][j] - R[p][j]) * (1- alpha/dmax) + (R[pmin][j] - R[p][j]) * (1 - alpha/dmin));
        if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap, "%f ",pot_harm[j]);}
        deplacement = - R[p][j];
        R[p][j] += Delta * pot_harm[j];
        deplacement += R[p][j]; 
        if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap, "%f ", deplacement);}
    }
    for (int j = 0; j < 3; j++)
    {
        if(temperature == 1)
        {
            F_alea = sqrt(2 * Delta) * randn(); 
        }
        else
        {
            F_alea = 0;
        }
        if(t % periode_enregistrement_force == 0){fprintf(fichier_force_thermique, "%f ", F_alea);}
        // F_alea = 0;
        deplacement = - R[p][j];
        R[p][j] += F_alea;
        deplacement += R[p][j]; 
        if(t % periode_enregistrement_force == 0){fprintf(fichier_force_thermique, "%f ", deplacement);}
    }
     if(t % periode_enregistrement_force == 0){fprintf(fichier_force_thermique, "\n");}
     if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap, "\n");}
}


void liaison_sup(double alpha, double K_transpt, double Delta, double** r, int p1, int p2, FILE* fichier_force_rnap_2, int t, int periode_enregistrement_force){ 
        double dist = distance(r[p2], r[p1]);
        double deplacement = 0.0; 
        if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap_2, "%d ",p1);}

        for(int j = 0; j < 3; j ++ ){
            if(t % periode_enregistrement_force == 0){
                fprintf(fichier_force_rnap_2, "%lf ", K_transpt * (r[p2][j] - r[p1][j]) * (1 - alpha/dist));
            }
            deplacement = - r[p1][j];
            r[p1][j] +=   Delta * K_transpt * (r[p2][j] - r[p1][j]) * (1 - alpha/dist);
            deplacement += r[p1][j];
            if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap_2,"%lf ", deplacement);}
            
        }
        if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap_2, "\n");}
        if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap_2, "%d ",p2);}


        for (int j = 0; j < 3; j++) {
            if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap_2, "%lf ", - K_transpt * (r[p2][j] - r[p1][j]) * (1 - alpha/dist));}
            deplacement = - r[p2][j];
            r[p2][j] += - Delta * K_transpt * (r[p2][j] - r[p1][j]) * (1 - alpha/dist);
            deplacement += r[p2][j];
            if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap_2, "%lf ", deplacement);}
        }
        if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap_2, "\n");}
}

void centre_masse_rnap(
    double **R_rnap,
    int n_sub,
    double com[3]
)
{
    com[0] = 0.0;
    com[1] = 0.0;
    com[2] = 0.0;

    for (int i = 0; i < n_sub; i++) {
        com[0] += R_rnap[i][0];
        com[1] += R_rnap[i][1];
        com[2] += R_rnap[i][2];
    }

    double inv = 1.0 / (double)n_sub;
    com[0] *= inv;
    com[1] *= inv;
    com[2] *= inv;
}

int monomere_plus_proche(
    double **R,
    int N,
    const double pos[3],
    int p_ref,
    int W
)
{
    int i_min = p_ref;
    double d2_min = 1e300;

    int i0 = p_ref - W;
    int i1 = p_ref + W;
    if (i0 < 0) i0 = 0;
    if (i1 > N - 1) i1 = N - 1;

    for (int i = i0; i <= i1; i++) {
        double dx = R[i][0] - pos[0];
        double dy = R[i][1] - pos[1];
        double dz = R[i][2] - pos[2];
        double d2 = dx*dx + dy*dy + dz*dz;

        if (d2 < d2_min) {
            d2_min = d2;
            i_min = i;
        }
    }

    return i_min;
}

void force_rnap_along_chromatin(
    const Config *cfg,
    double **R,           // chromatine
    double **R_rnap,      // sous-unités RNAP
    int mono,             // monomère p le plus proche
    double F_act,         // intensité de la force
    double Delta          // pas de temps
)
{
    const double eps = 1e-12;
    const int nsub = cfg->rnap_subunits;

    /* Sécurité bords */
    if (mono < 0 || mono >= cfg->N - 1) return;

    /* Tangente locale de la chromatine */
    double tx = R[mono + 1][0] - R[mono][0];
    double ty = R[mono + 1][1] - R[mono][1];
    double tz = R[mono + 1][2] - R[mono][2];

    double norm = sqrt(tx*tx + ty*ty + tz*tz);
    if (norm < eps) return;

    tx /= norm;
    ty /= norm;
    tz /= norm;

    /* Applique la force à chaque sous-unité RNAP */
    for (int sub = 0; sub < nsub; sub++) {
        R_rnap[sub][0] += Delta * F_act * tx;
        R_rnap[sub][1] += Delta * F_act * ty;
        R_rnap[sub][2] += Delta * F_act * tz;
    }

    // OPTIONNEL : réaction sur la chromatine (souvent NON souhaitée)
    R[mono][0]     -= Delta * F_act * tx / nsub;
    R[mono][1]     -= Delta * F_act * ty / nsub;
    R[mono][2]     -= Delta * F_act * tz / nsub;
    
}

void bond_rnap_chromatine(
    const Config *cfg, 
    double **R,
    double **R_rnap, 
    double a_transpt, 
    double K_transpt, 
    double Delta, 
    int mono_transcrpt, 
    double a, 
    double alpha
)
{
    const int NB_SUBUNITS = cfg->rnap_subunits;
    const double eps = 1e-12;   // évite division par 0

    int p = mono_transcrpt;
    int p1 = mono_transcrpt + 1;

    for(int sub = 0; sub < NB_SUBUNITS; sub++)
    {
        // Vecteur vers p et p+1
        double d0x = R[p][0]  - R_rnap[sub][0];
        double d0y = R[p][1]  - R_rnap[sub][1];
        double d0z = R[p][2]  - R_rnap[sub][2];

        double d1x = R[p1][0] - R_rnap[sub][0];
        double d1y = R[p1][1] - R_rnap[sub][1];
        double d1z = R[p1][2] - R_rnap[sub][2];

        double dist0 = distance(R[p], R_rnap[sub]); if (dist0 < eps) dist0 = eps;
        double dist1 = distance(R[p1], R_rnap[sub]); if (dist1 < eps) dist1 = eps;

        double f0 = K_transpt * (1.0 - a_transpt / dist0);
        double f1 = K_transpt * (1.0 - a_transpt / dist1);

        // Force exercée par p et p + 1 sur la sous-unité
        double Fx = d0x * f0 + d1x * f1;
        double Fy = d0y * f0 + d1y * f1;
        double Fz = d0z * f0 + d1z * f1;

        // Force exercée par p sur la sous-unité
        double Fx0 = d0x * f0;
        double Fy0 = d0y * f0;
        double Fz0 = d0z * f0;

        // Force exercée par p+1
        double Fx1 = d1x * f1;
        double Fy1 = d1y * f1;
        double Fz1 = d1z * f1;

        R_rnap[sub][0] += Delta * Fx;
        R_rnap[sub][1] += Delta * Fy;
        R_rnap[sub][2] += Delta * Fz;

        R[p][0]  -= Delta * Fx0;
        R[p][1]  -= Delta * Fy0;
        R[p][2]  -= Delta * Fz0;

        R[p1][0] -= Delta * Fx1;
        R[p1][1] -= Delta * Fy1;
        R[p1][2] -= Delta * Fz1;
    }
}


void bond_rnap_bead_progressive_mvt(
    const Config *cfg,
    double** R,                 // positions de la chromatine
    double** R_rnap,            // positions des sous-unités RNAP [8][3]
    double a_transpt,           // longueur à l'équilibre (si utile dans votre modèle)
    double K_transpt,           // raideur du couplage
    double Delta,               // pas de temps 
    int mono_transcrpt,         // index du monomère courant p
    double dx_avancement_rnap,  // progression dans [0,1] entre p et p+1
    double a,                   // (non utilisé ici, laissé pour compat)
    double alpha,               // paramètre de votre facteur (1 - (alpha+1)/(2*dist))
    FILE* fichier_force_lea,    // log optionnel
    int periode_enregistrement_force,
    int t
) {
    const int NB_SUBUNITS = cfg->rnap_subunits;
    const double eps = 1e-12;   // évite division par 0
    (void)a; // pas utilisé ici, mais conservé pour compat

    // Progression lissée pour éviter les discontinuités de vitesse
    // smoothstep: s -> s^2(3-2s)
    double s = dx_avancement_rnap;
    if (s < 0.0) s = 0.0;
    if (s > 1.0) s = 1.0;
    double s_eff = s * s * (3.0 - 2.0 * s);

    // Poids continus et lissés entre mono p et p+1
    double w0 = 1.0 - s_eff;   // poids appliqué au monomère p
    double w1 = s_eff;         // poids appliqué au monomère p+1

    // indices sécurité (assume mono_transcrpt+1 valide)
    int p   = mono_transcrpt;
    int p1  = mono_transcrpt + 1;

    if (t % periode_enregistrement_force == 0 && fichier_force_lea) {
        fprintf(fichier_force_lea, "t=%d RNAP-bond progressive (p=%d → p+1)\n", t, p);
        fprintf(fichier_force_lea, "weights: w0=%.6g w1=%.6g (s=%.6g s_eff=%.6g)\n", w0, w1, dx_avancement_rnap, s_eff);
    }

    for (int sub = 0; sub < NB_SUBUNITS; ++sub) {
        // Vecteurs vers p et p+1
        double d0x = R[p][0]  - R_rnap[sub][0];
        double d0y = R[p][1]  - R_rnap[sub][1];
        double d0z = R[p][2]  - R_rnap[sub][2];

        double d1x = R[p1][0] - R_rnap[sub][0];
        double d1y = R[p1][1] - R_rnap[sub][1];
        double d1z = R[p1][2] - R_rnap[sub][2];

        // Distances
        double dist0 = sqrt(d0x*d0x + d0y*d0y + d0z*d0z); if (dist0 < eps) dist0 = eps;
        double dist1 = sqrt(d1x*d1x + d1y*d1y + d1z*d1z); if (dist1 < eps) dist1 = eps;


        // Si vous préférez un ressort harmonique linéaire pur, remplacer f0/f1 par (dist - a_transpt)/dist
        double f0 = (1.0 - (alpha + 1.0) / (2.0 * dist0));
        double f1 = (1.0 - (alpha + 1.0) / (2.0 * dist1));

        // Force totale exercée par la chromatine sur la sous-unité (pondérée p/p+1)
        // F_subunit = K * [ w0 * d0 * f0 + w1 * d1 * f1 ]
        double Fx = K_transpt * ( w0 * d0x * f0 + w1 * d1x * f1 );
        double Fy = K_transpt * ( w0 * d0y * f0 + w1 * d1y * f1 );
        double Fz = K_transpt * ( w0 * d0z * f0 + w1 * d1z * f1 );

        // Mise à jour RNAP (r_new reçoit +F * Delta)
        R_rnap[sub][0] += Delta * Fx;
        R_rnap[sub][1] += Delta * Fy;
        R_rnap[sub][2] += Delta * Fz;

        // Réaction sur les monomères (opposée et pondérée)
        // p reçoit -w0 * F ; p+1 reçoit -w1 * F
        R[p][0]  -= Delta * w0 * Fx;
        R[p][1]  -= Delta * w0 * Fy;
        R[p][2]  -= Delta * w0 * Fz;

        R[p1][0] -= Delta * w1 * Fx;
        R[p1][1] -= Delta * w1 * Fy;
        R[p1][2] -= Delta * w1 * Fz;

        if (t % periode_enregistrement_force == 0 && fichier_force_lea) {
            fprintf(fichier_force_lea,
                " subunit=%d  dist0=%.4g dist1=%.4g  F=(%.4g,%.4g,%.4g)\n",
                sub, dist0, dist1, Fx, Fy, Fz
            );
        }
    }

    if (t % periode_enregistrement_force == 0 && fichier_force_lea) {
        fprintf(fichier_force_lea, "\n");
        fflush(fichier_force_lea);
    }
}


/*void lennard_jones_forces_rnap(const Config *cfg, double ***R_rnap, int nb_rnap, double **R, int N, NeighborList_rnap **neighbor_lists, double epsilon, double sigma6, double sigma12, double sigma6rnap2, double sigma12rnap2, double cut_rnap2, double Delta, int t, FILE* test, int T, FILE* fichier_force_rnap, int periode_enregistrement_force, int* l_rnap){
    if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap,"TIMESTEP: %d\n",t);}
    double deplacement = 0.0;
    for (int rnap = 0; rnap < nb_rnap; rnap++) 
    {
        if(l_rnap[rnap]<0){continue;}
        
        if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap,"RNAP: %d\n",rnap);}

        for (int subunit = 0; subunit < cfg->rnap_subunits ; subunit++)
        {
            if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap,"SUBUNIT: %d\n",subunit);}

            for (int k = 0; k < neighbor_lists[rnap][subunit].count; k++)
            {
                int j = neighbor_lists[rnap][subunit].neighbors[k];
                if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap,"NEIGHBOR: %d ",j);}
                double *Ri = R_rnap[rnap][subunit], *Rj = R[j];
                
                double dx = Ri[0] - Rj[0];
                double dy = Ri[1] - Rj[1];
                double dz = Ri[2] - Rj[2];

                double r2 = dx * dx + dy * dy + dz * dz;


                double r8 = r2 * r2 * r2 * r2; 
                double r14 = r8 * r2 * r2 * r2; 
                // double f = 4 * epsilon * (12 * sigma12/r14 - 6 * sigma6 / r8);
                double f = 4  * (12 * epsilon * sigma12/r14);

                double e = 300;
                if (f > e) {f = e;}


                if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap,"%f ", f * dx);}
                deplacement = - Ri[0];
                Ri[0] += Delta * f * dx;
                deplacement += Ri[0];
                if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap,"%f ", deplacement);}
                

                if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap,"%f ", f * dy);}
                deplacement = - Ri[1];
                Ri[1] += Delta * f * dy;
                deplacement += Ri[1];
                if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap,"%f ", deplacement);}
               
                if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap,"%f ", f * dz);}
                deplacement = - Ri[2];
                Ri[2] += Delta * f * dz;
                deplacement += Ri[2];
                if(t % periode_enregistrement_force == 0){fprintf(fichier_force_rnap,"%f\n", deplacement);}
                
                Rj[0] -= Delta * f * dx;
                Rj[1] -= Delta * f * dy;
                Rj[2] -= Delta * f * dz;
            }
            for (int subunit2 = 0; subunit2<cfg->rnap_subunits ; subunit2++)
            
            {
                double epsilon_test_att = epsilon; 
                double epsilon_test_rep = 1.5 * epsilon;
                if (subunit == subunit2) {continue;}
                double *Ri = R_rnap[rnap][subunit], *Rj = R_rnap[rnap][subunit2];
                double dx = Ri[0] - Rj[0];
                double dy = Ri[1] - Rj[1];
                double dz = Ri[2] - Rj[2];
                double r2 = dx * dx + dy * dy + dz * dz;
                if (sqrt(r2)>cut_rnap2) {continue;}
                

                double r8 = r2 * r2 * r2 * r2; 
                double r14 = r8 * r2 * r2 * r2; 
                double f = 4  * (12 * epsilon_test_rep * sigma12/r14 - 6 * epsilon_test_att * sigma6 / r8);
                

                double e = 300;
                if (f > e) {f = e;}

                Ri[0] += Delta * f * dx;
                Ri[1] += Delta * f * dy;
                Ri[2] += Delta * f * dz;
                Rj[0] -= Delta * f * dx;
                Rj[1] -= Delta * f * dy;
                Rj[2] -= Delta * f * dz;
            }
        } 
    }
}*/

/* void lennard_jones_forces_rnap(
    const Config *cfg,
    double ***R_rnap,
    int nb_rnap,
    double **R,
    int N,  // (optionnel)
    NeighborList_rnap **neighbor_lists,
    double epsilon,
    double sigma6,
    double sigma12,
    double sigma6rnap2,   // (optionnel)
    double sigma12rnap2,  // (optionnel)
    double cut_rnap2,     // cutoff^2
    double Delta,
    int t,
    FILE *test,           // (optionnel)
    int T,                // (optionnel)
    FILE *fichier_force_rnap,
    int periode_enregistrement_force,
    int *l_rnap
){
    const double eps2 = 1e-12;

    // Clamp sur le déplacement s = Delta * f
    const double s_max_rc = 0.02;   // RNAP-chromatine
    const double s_max_in = 0.02;   // intra-RNAP

    const int do_log = (fichier_force_rnap && periode_enregistrement_force > 0
                        && (t % periode_enregistrement_force == 0));

    if (do_log) fprintf(fichier_force_rnap, "TIMESTEP: %d\n", t);

    for (int rnap = 0; rnap < nb_rnap; rnap++) {
        if (l_rnap && l_rnap[rnap] < 0) continue;

        if (do_log) fprintf(fichier_force_rnap, "RNAP: %d\n", rnap);

        for (int subunit = 0; subunit < cfg->rnap_subunits; subunit++) {

            if (do_log) fprintf(fichier_force_rnap, "SUBUNIT: %d\n", subunit);

            for (int k = 0; k < neighbor_lists[rnap][subunit].count; k++) {
                int j = neighbor_lists[rnap][subunit].neighbors[k];
                if (j < 0 || j >= N) continue; // sécurité

                double *Ri = R_rnap[rnap][subunit];
                double *Rj = R[j];

                double dx = Ri[0] - Rj[0];
                double dy = Ri[1] - Rj[1];
                double dz = Ri[2] - Rj[2];

                double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 < eps2 || r2 > cut_rnap2) continue;

                // Calcul stable via inv_r2, inv_r6, inv_r12
                double inv_r2  = 1.0 / r2;
                double inv_r6  = inv_r2 * inv_r2 * inv_r2;
                double inv_r12 = inv_r6 * inv_r6;

                // Si tu veux LJ complet RNAP-chromatine, utilise ceci :
                double f = 24.0 * epsilon * (2.0 * sigma12 * inv_r12 - sigma6 * inv_r6) * inv_r2;

                // Si tu veux répulsion pure (ce que ton code faisait), remplace par :
                // double f = 48.0 * epsilon * sigma12 * inv_r12 * inv_r2; // ~ 4*12*epsilon*sigma12/r14 * rvec

                double s = Delta * f;
                if (s >  s_max_rc) s =  s_max_rc;
                if (s < -s_max_rc) s = -s_max_rc;

                if (do_log) {
                    // on logge "force*dx" et le déplacement effectué sur x/y/z
                    double dRx = s * dx, dRy = s * dy, dRz = s * dz;
                    fprintf(fichier_force_rnap, "NEIGHBOR: %d ", j);
                    fprintf(fichier_force_rnap, "%f %f %f %f %f %f\n",
                            f*dx, dRx, f*dy, dRy, f*dz, dRz);
                }

                Ri[0] += s * dx;  Ri[1] += s * dy;  Ri[2] += s * dz;
                Rj[0] -= s * dx;  Rj[1] -= s * dy;  Rj[2] -= s * dz;
            }

            for (int subunit2 = subunit + 1; subunit2 < cfg->rnap_subunits; subunit2++) {

                double epsilon_att = epsilon;
                double epsilon_rep = 1.5 * epsilon;

                double *Ri = R_rnap[rnap][subunit];
                double *Rj = R_rnap[rnap][subunit2];

                double dx = Ri[0] - Rj[0];
                double dy = Ri[1] - Rj[1];
                double dz = Ri[2] - Rj[2];

                double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 < eps2 || r2 > cut_rnap2) continue;

                double inv_r2  = 1.0 / r2;
                double inv_r6  = inv_r2 * inv_r2 * inv_r2;
                double inv_r12 = inv_r6 * inv_r6;

                double f = 24.0 * (2.0 * epsilon_rep * sigma12 * inv_r12
                                   -       epsilon_att * sigma6  * inv_r6) * inv_r2;

                double s = Delta * f;
                if (s >  s_max_in) s =  s_max_in;
                if (s < -s_max_in) s = -s_max_in;

                Ri[0] += s * dx;  Ri[1] += s * dy;  Ri[2] += s * dz;
                Rj[0] -= s * dx;  Rj[1] -= s * dy;  Rj[2] -= s * dz;
            }
        }
    }
}*/

void lennard_jones_forces_rnap(
    const Config *cfg,
    double ***R_rnap,
    int nb_rnap,
    double **R,
    int N,
    NeighborList_rnap **neighbor_lists,
    double epsilon,
    double sigma6,
    double sigma12,
    double sigma6rnap2,   // unused
    double sigma12rnap2,  // unused
    double cut_rnap2,     // cutoff^2
    double Delta,
    int t,
    FILE *test,           // unused
    int T,                // unused
    FILE *fichier_force_rnap,
    int periode_enregistrement_force,
    int *l_rnap
){
    const double eps2 = 1e-12;

    // --- Réglages : répulsion > attraction ---
    const double rep_factor_rc   = 5.0;   // RNAP–chromatine (essaye 3, 5, 10)
    const double rep_factor_intra= 2.0;   // intra-RNAP (anneau), à ajuster

    const double epsilon_rep_rc   = rep_factor_rc    * epsilon;
    const double epsilon_att_rc   = 1.0             * epsilon;

    const double epsilon_rep_intra= rep_factor_intra * epsilon;
    const double epsilon_att_intra= 1.0             * epsilon;

    // Clamp sur le déplacement s = Delta * f
    const double s_max_rc = 0.02;  // RNAP–chromatine
    const double s_max_in = 0.02;  // intra-RNAP

    const int do_log = (fichier_force_rnap && periode_enregistrement_force > 0
                        && (t % periode_enregistrement_force == 0));

    if (do_log) fprintf(fichier_force_rnap, "TIMESTEP: %d\n", t);

    for (int rnap = 0; rnap < nb_rnap; rnap++) {
        if (l_rnap && l_rnap[rnap] < 0) continue;

        if (do_log) fprintf(fichier_force_rnap, "RNAP: %d\n", rnap);

        for (int subunit = 0; subunit < cfg->rnap_subunits; subunit++) {

            if (do_log) fprintf(fichier_force_rnap, "SUBUNIT: %d\n", subunit);

            /* ===============================
               RNAP – CHROMATINE (neighbor list)
               répulsion > attraction
               =============================== */
            for (int k = 0; k < neighbor_lists[rnap][subunit].count; k++) {
                int j = neighbor_lists[rnap][subunit].neighbors[k];
                if (j < 0 || j >= N) continue;

                double *Ri = R_rnap[rnap][subunit];
                double *Rj = R[j];

                double dx = Ri[0] - Rj[0];
                double dy = Ri[1] - Rj[1];
                double dz = Ri[2] - Rj[2];

                double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 < eps2 || r2 > cut_rnap2) continue;

                double inv_r2  = 1.0 / r2;
                double inv_r6  = inv_r2 * inv_r2 * inv_r2;
                double inv_r12 = inv_r6 * inv_r6;

                // F_vec = f * r_vec
                // f = 24 [ 2 ε_rep σ^12/r^14  -  ε_att σ^6/r^8 ]
                double f = 24.0 * ( 2.0 * epsilon_rep_rc * sigma12 * inv_r12
                                   -      epsilon_att_rc * sigma6  * inv_r6 ) * inv_r2;

                double s = Delta * f;
                if (s >  s_max_rc) s =  s_max_rc;
                if (s < -s_max_rc) s = -s_max_rc;

                if (do_log) {
                    double dRx = s * dx, dRy = s * dy, dRz = s * dz;
                    fprintf(fichier_force_rnap, "NEIGHBOR: %d ", j);
                    fprintf(fichier_force_rnap, "%f %f %f %f %f %f\n",
                            f*dx, dRx, f*dy, dRy, f*dz, dRz);
                }

                Ri[0] += s * dx;  Ri[1] += s * dy;  Ri[2] += s * dz;
                Rj[0] -= s * dx;  Rj[1] -= s * dy;  Rj[2] -= s * dz;
            }

            /* ===============================
               Intra-RNAP (subunit2 > subunit)
               répulsion > attraction
               =============================== */
            for (int subunit2 = subunit + 1; subunit2 < cfg->rnap_subunits; subunit2++) {

                double *Ri = R_rnap[rnap][subunit];
                double *Rj = R_rnap[rnap][subunit2];

                double dx = Ri[0] - Rj[0];
                double dy = Ri[1] - Rj[1];
                double dz = Ri[2] - Rj[2];

                double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 < eps2 || r2 > cut_rnap2) continue;

                double inv_r2  = 1.0 / r2;
                double inv_r6  = inv_r2 * inv_r2 * inv_r2;
                double inv_r12 = inv_r6 * inv_r6;

                double f = 24.0 * ( 2.0 * epsilon_rep_intra * sigma12 * inv_r12
                                   -      epsilon_att_intra * sigma6  * inv_r6 ) * inv_r2;

                double s = Delta * f;
                if (s >  s_max_in) s =  s_max_in;
                if (s < -s_max_in) s = -s_max_in;

                Ri[0] += s * dx;  Ri[1] += s * dy;  Ri[2] += s * dz;
                Rj[0] -= s * dx;  Rj[1] -= s * dy;  Rj[2] -= s * dz;
            }
        }
    }
}


/*void lennard_jones_forces_rnap_rnap(
                                    const Config *cfg,
                                    double ***R_rnap,
                                    int nb_rnap,
                                    double epsilon,
                                    double sigma6,
                                    double sigma12,
                                    double cut_rnap2,
                                    double Delta,
                                    int* l_rnap)
{

    double epsilon_rep = 100 * epsilon; 
    double epsilon_att =  epsilon; 

    for (int rnap_i = 0; rnap_i < nb_rnap; rnap_i++)
    {
        if(l_rnap[rnap_i] < 0){continue;}
        for (int sub_i = 0; sub_i < cfg->rnap_subunits; sub_i++) {

            double *Ri = R_rnap[rnap_i][sub_i];

            // --- Interactions RNAP_i ↔ RNAP_j (j > i pour éviter les doublons) ---
            for (int rnap_j = rnap_i + 1; rnap_j < nb_rnap; rnap_j++) {
                if (l_rnap[rnap_j] < 0) continue;
                for (int sub_j = 0; sub_j < cfg->rnap_subunits; sub_j++) {

                    double *Rj = R_rnap[rnap_j][sub_j];

                    double dx = Ri[0] - Rj[0];
                    double dy = Ri[1] - Rj[1];
                    double dz = Ri[2] - Rj[2];

                    double r2 = dx * dx + dy * dy + dz * dz;
                    if (r2 > cut_rnap2 || r2 == 0.0)
                        continue;

                    // Pré-calculs pour le potentiel LJ
                    double r8  = r2 * r2 * r2 * r2;
                    double r14 = r8 * r2 * r2 * r2;

                    // Force de Lennard-Jones classique : F = -dU/dr
                    double f = 4.0 * (12.0 * epsilon_rep * sigma12 / r14
                                    - 6.0 * epsilon_att * sigma6 / r8);

                    // Saturation (pour éviter les explosions)
                    const double f_max = 1000.0;
                    if (f > f_max)
                        f = f_max;

                    // printf("f = %lf \n", f);

                    // Application symétrique de la force
                    Ri[0] += Delta * f * dx;
                    Ri[1] += Delta * f * dy;
                    Ri[2] += Delta * f * dz;

                    Rj[0] -= Delta * f * dx;
                    Rj[1] -= Delta * f * dy;
                    Rj[2] -= Delta * f * dz;
                }
            }

            // --- Auto-interactions entre sous-unités du même RNAP ---
            for (int sub_j = sub_i + 1; sub_j < cfg->rnap_subunits; sub_j++) {
                if (sub_i == sub_j) continue;

                double *Rj = R_rnap[rnap_i][sub_j];
                double dx = Ri[0] - Rj[0];
                double dy = Ri[1] - Rj[1];
                double dz = Ri[2] - Rj[2];
                double r2 = dx * dx + dy * dy + dz * dz;

                const double eps = 1e-12;
                if (r2 > cut_rnap2 || r2 < eps)
                    continue;

                
                double inv_r2 = 1.0/r2;
                double inv_r6 = inv_r2*inv_r2*inv_r2;
                double inv_r12 = inv_r6*inv_r6;

                // Force modifiée : plus répulsive entre sous-unités d’un même RNAP
                double f = 24.0 * (epsilon_rep*2.0*sigma12*inv_r12 - epsilon_att*sigma6*inv_r6) * inv_r2;

                const double f_max = 300.0;
                if (f > f_max)
                    f = f_max;

                Ri[0] += Delta * f * dx;
                Ri[1] += Delta * f * dy;
                Ri[2] += Delta * f * dz;

                Rj[0] -= Delta * f * dx;
                Rj[1] -= Delta * f * dy;
                Rj[2] -= Delta * f * dz;
            }
        }
    }
}*/

void lennard_jones_forces_rnap_rnap(
    const Config *cfg,
    double ***R_rnap,
    int nb_rnap,
    double epsilon,
    double sigma6,
    double sigma12,
    double cut_rnap2,   // ATTENTION: cutoff^2
    double Delta,
    int *l_rnap
)
{
    const double eps2 = 1e-12;

    // Tes choix
    const double epsilon_rep = 100.0 * epsilon;
    const double epsilon_att = 1.0   * epsilon;

    // Clamp sur le déplacement: s = Delta * f
    // s multiplie (dx,dy,dz) donc il a dimension ~ 1/r^2 (selon ton schéma),
    // mais en pratique ça borne |Δr| ~ |s| * r.
    const double s_max_inter = 0.02;  // inter-RNAP : à ajuster
    const double s_max_intra = 0.02;  // intra-RNAP : à ajuster

    for (int rnap_i = 0; rnap_i < nb_rnap; rnap_i++) {
        if (l_rnap && l_rnap[rnap_i] < 0) continue;

        for (int sub_i = 0; sub_i < cfg->rnap_subunits; sub_i++) {

            double *Ri = R_rnap[rnap_i][sub_i];

            /* =========================
               RNAP_i ↔ RNAP_j (j > i)
               ========================= */
            for (int rnap_j = rnap_i + 1; rnap_j < nb_rnap; rnap_j++) {
                if (l_rnap && l_rnap[rnap_j] < 0) continue;

                for (int sub_j = 0; sub_j < cfg->rnap_subunits; sub_j++) {

                    double *Rj = R_rnap[rnap_j][sub_j];

                    double dx = Ri[0] - Rj[0];
                    double dy = Ri[1] - Rj[1];
                    double dz = Ri[2] - Rj[2];

                    double r2 = dx*dx + dy*dy + dz*dz;
                    if (r2 < eps2 || r2 > cut_rnap2) continue;

                    // r^8 et r^14 à partir de r^2 (comme ton code)
                    double r4  = r2 * r2;
                    double r8  = r4 * r4;
                    double r14 = r8 * r4 * r2; // r^8 * r^4 * r^2 = r^14

                    // f est le scalaire tel que F_vec = f * r_vec
                    double f = 4.0 * (12.0 * epsilon_rep * sigma12 / r14
                                    - 6.0  * epsilon_att * sigma6  / r8);

                    // Clamp sur le déplacement (scalaire) s = Delta*f
                    double s = Delta * f;
                    if (s >  s_max_inter) s =  s_max_inter;
                    if (s < -s_max_inter) s = -s_max_inter;

                    Ri[0] += s * dx;  Ri[1] += s * dy;  Ri[2] += s * dz;
                    Rj[0] -= s * dx;  Rj[1] -= s * dy;  Rj[2] -= s * dz;
                }
            }

            /* =========================
               Intra-RNAP (sub_j > sub_i)
               ========================= */
            for (int sub_j = sub_i + 1; sub_j < cfg->rnap_subunits; sub_j++) {

                double *Rj = R_rnap[rnap_i][sub_j];

                double dx = Ri[0] - Rj[0];
                double dy = Ri[1] - Rj[1];
                double dz = Ri[2] - Rj[2];

                double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 < eps2 || r2 > cut_rnap2) continue;

                double r4  = r2 * r2;
                double r8  = r4 * r4;
                double r14 = r8 * r4 * r2;

                double f = 4.0 * (12.0 * epsilon_rep * sigma12 / r14
                                - 6.0  * epsilon_att * sigma6  / r8);

                double s = Delta * f;
                if (s >  s_max_intra) s =  s_max_intra;
                if (s < -s_max_intra) s = -s_max_intra;

                Ri[0] += s * dx;  Ri[1] += s * dy;  Ri[2] += s * dz;
                Rj[0] -= s * dx;  Rj[1] -= s * dy;  Rj[2] -= s * dz;
            }
        }
    }
}





/*******************************************************************************************************
 * 🧬 simu_LJ_RNAP_erwan
 * -----------------------------------------------------------------------------------------------------
 * • Description générale :
 *   Fonction principale de simulation de la dynamique brownienne d’une chaîne de chromatine 
 *   soumise à des interactions de Lennard-Jones (LJ) et à l’action de plusieurs ARN polymérases (RNAP).
 * 
 *   Elle constitue le moteur central de la simulation :
 *     - Initialise les structures de voisinage (NeighborList),
 *     - Met le système à l’équilibre (phase de relaxation),
 *     - Lance la boucle temporelle principale de la dynamique,
 *     - Gère les enregistrements (positions, RNAP, métriques),
 *     - Finalise la simulation et libère les ressources mémoire.
 *
 * • Entrées :
 *   - cfg  : paramètres physiques et numériques de la simulation (structure Config)
 *   - sv   : variables dynamiques du système (positions, forces, RNAP, etc.)
 *   - nbr_simu : index ou identifiant de la simulation courante
 *   - f    : structure regroupant les fichiers d’écriture
 *
 * • Sorties :
 *   - Écrit plusieurs fichiers de trajectoires et statistiques :
 *       → Trajectoire principale (.lammpstrj)
 *       → Positions des RNAP
 *       → Fichiers de métriques (forces, end-to-end, etc.)
 *   - Met à jour les structures de simulation
 *
 * • Étapes principales :
 *   1. Allocation des structures de voisinage
 *   2. Initialisation / enregistrement de l’état initial
 *   3. Mise à l’équilibre du système
 *   4. Boucle de simulation principale
 *   5. Calcul des grandeurs moyennes et libération mémoire
 *
 * • Auteur :
 *   Erwan Le Floch — LPT (Laboratoire de Physique Théorique, Toulouse)
 *   Projet TIRIS 2024–2028 — Modélisation multi-échelle de la transcription et de la dynamique chromatinienne
 *
 *******************************************************************************************************/

void simu_LJ_RNAP_erwan(const Config *cfg, SimVars *sv, const Files *f, int t_start)
{
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////   🕹️ Initialisation de la simulation   ///////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    clock_t start, end;
    double duree_tot = 0.0, temps_restant;
    int t = 0;

    double **R     = sv->R;
    double **R_new = sv->R_new;

    // --- Allocation des neighbor lists (particules + RNAP) ---
    NeighborList_rnap **neighbor_lists_rnap = allocate_neighbor_list_rnap(sv->nb_rnap, cfg->rnap_subunits);
    NeighborList *neighbor_lists = malloc(cfg->N * sizeof(NeighborList));

    if (!neighbor_lists) {
        fprintf(stderr, "❌ Erreur d'allocation : neighbor_lists\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < cfg->N; i++) {
        neighbor_lists[i].neighbors = malloc(10 * sizeof(int));
        neighbor_lists[i].capacity  = 10;
        if (!neighbor_lists[i].neighbors) {
            fprintf(stderr, "❌ Erreur d'allocation : neighbor_lists[%d].neighbors\n", i);
            exit(EXIT_FAILURE);
        }
    }

    if (cfg->quench == 1)
    {
        for ( ; sv->nb_rnap < cfg->nb_rnap_initial; sv->nb_rnap++) {

            sv->l_rnap[sv->nb_rnap] = 1;

            int old_nb = sv->nb_rnap;
            int new_nb = sv->nb_rnap + 1;

            // 🔧 Mise à jour dynamique de la liste des voisins RNAP
            resize_neighbor_list_rnap(
                &neighbor_lists_rnap,
                old_nb,
                new_nb,
                cfg->N,      // nombre de lignes (typiquement nombre de monomères)
                t            // pas de temps (pour debug)
            );

            // 💡 Position de la nouvelle RNAP
            sv->positions_bille_rnap[sv->nb_rnap] = cfg->debut_segment + 3 * sv->nb_rnap;

            // 🧬 Création physique de la nouvelle RNAP
            creation_1_rnap_erwan(
                cfg, sv->R, sv->nb_rnap, sv->positions_bille_rnap, sv->R_rnap,
                cfg->N, cfg->rnap_subunits,
                cfg->debut_segment, cfg->fin_segment, cfg->nb_rnap_initial
            );
        }
    }






    for(int i = 0; i < cfg->nb_rnap_initial; i++)
    {
        for(int j = 0; j < cfg->rnap_subunits; j++)
        {
            for(int k = 0; k < 3; k++)
            {
                printf("R_rnap[%d][%d][%d] = %lf ", i, j, k, sv->R_rnap[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////   💾 État initial et enregistrements initiaux   ////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    if (cfg->nb_rnap_initial > 0 && cfg->resume_from_checkpoint == 0)
    {
        enregistrement_RNAP(
            f->fichier,                 
            sv->R,                     
            cfg->N,                     
            sv->R_rnap,                 
            sv->nb_rnap,                
            t,                          
            sv->positions_bille_rnap,   
            cfg->rnap_subunits,         
            cfg->nb_rnap_initial,      
            sv->l_rnap,
            sv,
            cfg
        );
    }

//     //////////////////////////////////////////////////////////////////////////////////////////////////////
//     //////////////////////////////////////   ⚖️ Mise à l'équilibre   //////////////////////////////////////
//     //////////////////////////////////////////////////////////////////////////////////////////////////////
    if(cfg->resume_from_checkpoint == 0)
    {
        f_equilibriate(sv, cfg, f, neighbor_lists, neighbor_lists_rnap);
    }

    
//     //////////////////////////////////////////////////////////////////////////////////////////////////////
//     ///////////////////////////////////////   🧮 Boucle principale   //////////////////////////////////////
//     //////////////////////////////////////////////////////////////////////////////////////////////////////
calcul(sv, cfg, f, neighbor_lists, neighbor_lists_rnap, t_start);

//     endtoend /= cfg->T_endtoend;

//     //////////////////////////////////////////////////////////////////////////////////////////////////////
//     //////////////////////////////////////   🧾 Finalisation et sortie   //////////////////////////////////
//     //////////////////////////////////////////////////////////////////////////////////////////////////////
//     finalize_simulationaaaaa(sv, cfg, f, endtoend, nbr_simu);

//     for (int i = 0; i < cfg->N; i++)
//         free(neighbor_lists[i].neighbors);
//     free(neighbor_lists);

//     printf("🏁 Simulation terminée avec succès.\n\n");
}