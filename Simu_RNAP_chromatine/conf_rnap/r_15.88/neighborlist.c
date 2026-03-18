#include "simulation.h"
#include "config.h"
#include "basic_functions.h"
#include "transcription.h"
#include "movement.h"   
#include "neighborlist.h" 
#include <assert.h>


/**
 * @brief Builds the neighbor list for each particle based on a cutoff radius.
 *
 * For every particle i, this function finds all other particles j such that
 * the squared distance r_ij² < (RCUT + SKIN)², and stores their indices
 * in neighbor_lists[i]. The list expands dynamically as needed.
 *
 * @param R                Positions of all particles [N][3]
 * @param neighbor_lists   Array of NeighborList structures, one per particle
 * @param N                Total number of particles
 * @param RCUT             Cutoff radius for neighbor search
 * @param SKIN             Additional buffer distance for Verlet list updates
 *
 * @note This function reallocates memory dynamically if capacity is exceeded.
 *       Make sure neighbor_lists[i].neighbors and capacity are properly initialized.
 */
void build_neighbor_list(double **R, NeighborList *neighbor_lists, int N, double RCUT, double SKIN) {

    double rcut_sq = (RCUT + SKIN) * (RCUT + SKIN);

    for (int i = 0; i < N; i++) {
        neighbor_lists[i].count = 0;

        for (int j = 0; j < N; j++) {
            if (i == j) continue;

            double dx = R[i][0] - R[j][0];
            double dy = R[i][1] - R[j][1];
            double dz = R[i][2] - R[j][2];
            double dist_sq = dx * dx + dy * dy + dz * dz;

            if (dist_sq < rcut_sq) {
                if (neighbor_lists[i].count >= neighbor_lists[i].capacity) {
                    neighbor_lists[i].capacity *= 2;
                    int *tmp = realloc(neighbor_lists[i].neighbors, neighbor_lists[i].capacity * sizeof(int));
                    if (!tmp) {
                        fprintf(stderr, "❌ Error: realloc failed for neighbor list of particle %d\n", i);
                        exit(EXIT_FAILURE);
                    }
                    neighbor_lists[i].neighbors = tmp;
                }
                neighbor_lists[i].neighbors[neighbor_lists[i].count++] = j;
            }
        }
    }
}


NeighborList_rnap** allocate_neighbor_list_rnap(int nb_rnap, int rows) {
    NeighborList_rnap** neighbor_lists = (NeighborList_rnap**)malloc(nb_rnap * sizeof(NeighborList_rnap*));
    if (neighbor_lists == NULL) {
        perror("Erreur d'allocation mémoire pour neighbor_lists");
        exit(EXIT_FAILURE);
    }

    for (int rnap = 0; rnap < nb_rnap; rnap++) {
        neighbor_lists[rnap] = (NeighborList_rnap*)malloc(rows * sizeof(NeighborList_rnap));
        if (neighbor_lists[rnap] == NULL) {
            perror("Erreur d'allocation mémoire pour neighbor_lists[rnap]");
            exit(EXIT_FAILURE);
        }

        for (int subunit = 0; subunit < rows; subunit++) {
            neighbor_lists[rnap][subunit].neighbors = (int*)malloc(10 * sizeof(int));
            if (neighbor_lists[rnap][subunit].neighbors == NULL) {
                perror("Erreur d'allocation mémoire pour neighbor_lists[rnap][subunit].neighbors");
                exit(EXIT_FAILURE);
            }
            neighbor_lists[rnap][subunit].capacity = 10;
            neighbor_lists[rnap][subunit].count = 0;
        }
    }

    return neighbor_lists;
}


// Libération de la mémoire pour NeighborList_rnap
void free_neighbor_list_rnap(NeighborList_rnap** neighbor_lists, int nb_rnap, int rows) {
    for (int rnap = 0; rnap < nb_rnap; rnap++) {
        for (int subunit = 0; subunit < rows; subunit++) {
            free(neighbor_lists[rnap][subunit].neighbors);
        }
        free(neighbor_lists[rnap]);
    }
    free(neighbor_lists);
}


/**
 * @brief Builds the neighbor list between RNAP subunits and chromatin monomers.
 *
 * For each RNAP (rnap) and each of its subunits, this function finds
 * all chromatin beads (monomers) within a given cutoff distance
 * (rayon_ecrantage_LJ_rnap). The indices of these monomers are stored
 * in neighbor_lists[rnap][subunit].neighbors.
 *
 * Memory is reallocated dynamically when the current capacity is exceeded.
 *
 * @param R_rnap                  3D positions of RNAP subunits [nb_rnap][rnap_subunits][3]
 * @param nb_rnap                 Number of active RNAPs
 * @param R                       3D positions of chromatin monomers [N][3]
 * @param N                       Total number of chromatin monomers
 * @param neighbor_lists           Array of neighbor lists for each RNAP subunit
 * @param rayon_ecrantage_LJ_rnap Cutoff radius for interaction (Lennard-Jones screening)
 * @param t                       Current timestep (for debug or logging)
 */
void build_neighbor_list_rnap_chrom(const Config *cfg, 
                                    double ***R_rnap,
                                    int nb_rnap,
                                    double **R,
                                    int N,
                                    NeighborList_rnap **neighbor_lists,
                                    double rayon_ecrantage_LJ_rnap,
                                    int t)
{
    const double rcut_sq = rayon_ecrantage_LJ_rnap * rayon_ecrantage_LJ_rnap;
    const int rows = cfg->rnap_subunits;  

    for (int rnap = 0; rnap < nb_rnap; rnap++) {
        for (int subunit = 0; subunit < rows; subunit++) {

            // Réinitialisation du compteur de voisins
            neighbor_lists[rnap][subunit].count = 0;

            // Parcours de tous les monomères de la chromatine
            for (int j = 0; j < N; j++) {
                double d2 = distance(R_rnap[rnap][subunit], R[j]);

                if (d2 < rcut_sq) {

                    // Réallocation si capacité dépassée
                    if (neighbor_lists[rnap][subunit].count >= neighbor_lists[rnap][subunit].capacity) {
                        int old_cap = neighbor_lists[rnap][subunit].capacity;
                        int new_cap = old_cap * 2;

                        int *tmp = realloc(
                            neighbor_lists[rnap][subunit].neighbors,
                            new_cap * sizeof(int)
                        );

                        if (!tmp) {
                            fprintf(stderr,
                                    "❌ [t=%d] realloc failed for RNAP %d subunit %d (old_cap=%d)\n",
                                    t, rnap, subunit, old_cap);
                            exit(EXIT_FAILURE);
                        }

                        neighbor_lists[rnap][subunit].neighbors = tmp;
                        neighbor_lists[rnap][subunit].capacity = new_cap;
                    }

                    // Ajout du voisin j à la liste
                    neighbor_lists[rnap][subunit]
                        .neighbors[
                            neighbor_lists[rnap][subunit].count++
                        ] = j;
                }
            }
        }
    }
}

void build_neighbor_list_rnap(double **R, NeighborList *neighbor_lists, int N) {
    int RCUT = 60 ; int SKIN =5 ;
    double rcut_sq = (RCUT + SKIN); //* (RCUT + SKIN);  // Rayon de coupure au carré
    for (int i = 0; i < N; i++) {
        neighbor_lists[i].count = 0;  // Réinitialise le nombre de voisins
        for (int j = 0; j < N; j++) {
            if (i == j) continue;  // Ignore la même particule
            double dx = R[i][0] - R[j][0];
            double dy = R[i][1] - R[j][1];
            double dz = R[i][2] - R[j][2];
            //double dist_sq = distance(R[i], R[j]);
            double dist_sq = dx * dx + dy * dy + dz * dz;
            if (dist_sq < rcut_sq) {  // Si la distance est inférieure à rcut_sq
                if (neighbor_lists[i].count >= neighbor_lists[i].capacity) {
                    neighbor_lists[i].capacity *= 2;
                    neighbor_lists[i].neighbors = realloc(neighbor_lists[i].neighbors, neighbor_lists[i].capacity * sizeof(int));
                }
                neighbor_lists[i].neighbors[neighbor_lists[i].count++] = j;  // Ajoute j à la liste des voisins de i
            }
        }
    }
}



#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


void resize_neighbor_list_rnap(NeighborList_rnap ***plists,
                               int old_nb,
                               int new_nb,
                               int rows,
                               int t)
{
    assert(plists != NULL);
    if (old_nb < 0 || new_nb < 0) {
        fprintf(stderr, "[t=%d] ERROR: old_nb=%d, new_nb=%d invalid\n",
                t, old_nb, new_nb);
        exit(EXIT_FAILURE);
    }

    NeighborList_rnap **lists = *plists;

    // 1) Si on supprime des RNAP
    if (new_nb < old_nb) {
        for (int i = new_nb; i < old_nb; i++) {
            for (int j = 0; j < rows; j++) {
                free(lists[i][j].neighbors);
                lists[i][j].neighbors = NULL;
            }
            free(lists[i]);
            lists[i] = NULL;
        }
    }

    // 2) Realloc du tableau principal
    lists = realloc(lists, new_nb * sizeof *lists);
    if (!lists && new_nb > 0) {
        perror("ERROR realloc neighbor_lists_rnap");
        exit(EXIT_FAILURE);
    }

    // 3) Si on ajoute des RNAP
    if (new_nb > old_nb) {
        for (int i = old_nb; i < new_nb; i++) {
            lists[i] = malloc(rows * sizeof *lists[i]);
            if (!lists[i]) {
                perror("ERROR malloc new neighbor_lists_rnap[i]");
                exit(EXIT_FAILURE);
            }
            for (int j = 0; j < rows; j++) {
                lists[i][j].capacity  = 10;
                lists[i][j].count     = 0;
                lists[i][j].neighbors = malloc(
                    lists[i][j].capacity * sizeof *lists[i][j].neighbors);
                if (!lists[i][j].neighbors) {
                    perror("ERROR malloc new neighbors");
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    // 4) Mise à jour du pointeur de sortie
    *plists = lists;

    // 5) Assertion finale
    assert( (new_nb > 0 && *plists != NULL) || (new_nb == 0) );

}