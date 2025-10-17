#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H

// Définition de vos deux structs :
typedef struct {
    int *neighbors;
    int  count;
    int  capacity;
} NeighborList;

typedef struct {
    int *neighbors;
    int    count;
    int    capacity;
} NeighborList_rnap;



void build_neighbor_list(double **R, NeighborList *neighbor_lists, int N, int RCUT, int SKIN);
NeighborList_rnap** allocate_neighbor_list_rnap(int nb_rnap, int rows);
void free_neighbor_list_rnap(NeighborList_rnap** neighbor_lists, int nb_rnap, int rows);
void build_neighbor_list_rnap_chrom(double ***R_rnap, int nb_rnap, int rnap_subunits, double **R, int N, NeighborList_rnap **neighbor_lists, double rayon_ecrantage_LJ_rnap, int t);
void build_neighbor_list_rnap(double **R, NeighborList *neighbor_lists, int N);
void resize_neighbor_list_rnap(NeighborList_rnap ***plists, int old_nb, int new_nb, int rows, int t);


#endif // NEIGHBORLIST_H