#include "simulation.h"
#include "config.h"
#include "basic_functions.h"
#include "transcription_erwan.h"
#include "movement.h"   
#include "neighborlist.h"   
#include "time.h"      


#include <stdlib.h>  // pour malloc, calloc
#include <math.h>    // pour log10()


void init_sim_vars(SimVars *sv, const Config *cfg) {

    // --- Periode 
    sv->T_enregistrement = cfg->T_enregistrement;
    sv->T_msd = cfg->T_msd;
    sv->T_correlation = cfg->T_correlation;
    sv->T_endtoend = cfg->T_endtoend;
    sv->T_centre_de_masse = cfg->T_centre_de_masse;
    sv->T_force = cfg->T_force;
    sv->T_voisin = cfg->T_voisin; 

    // Initialisation tableaux   
    sv->bending_forces = allocate_matrix(cfg->N, 3);     
    sv->R = allocate_matrix(cfg->N, 3);
    sv->R_new = allocate_matrix(cfg->N, 3);
    sv->t_link = allocate_matrix(cfg->N-1,3);
    sv->R_matrix = (double**)malloc(sizeof(double*)); 

    sv->compteur_mono_rnap = calloc(cfg->N, sizeof(int));

    sv->list_monomere = (int*)malloc(cfg->Nm * sizeof(int));
    for(int i = 0; i < cfg->Nm; i++) {
        sv->list_monomere[i] = i;
    }


    // Corrélations et temps
    sv->stock_correlation         = calloc(sv->T_correlation, sizeof(double));
    sv->stock_correlation_segment = calloc(sv->T_correlation, sizeof(double));
    sv->time                      = malloc(sv->T_correlation * sizeof(double));
    sv->log10_time                = malloc(sv->T_correlation * sizeof(double));

    // Stock end-to-end
    sv->stock         = allocate_matrix(sv->T_endtoend, 2);
    sv->Rbb_segment   = allocate_matrix(sv->T_endtoend, 3);
    sv->Rbb_avant     = allocate_matrix(sv->T_endtoend, 3);
    sv->Rbb_apres     = allocate_matrix(sv->T_endtoend, 3);
    sv->nombre_voisins= malloc(cfg->N * sizeof(int));

    // Corrélations spaciales
    sv->Rbb       = allocate_matrix(sv->T_correlation, 3);
    sv->R_segment = allocate_matrix(sv->T_correlation, 3);

    // Centre de masse & rayon de giration
    sv->R_centre_de_masse = allocate_matrix(sv->T_centre_de_masse, 3);
    sv->stock_cdm         = allocate_matrix(sv->T_centre_de_masse, 3);
    sv->gyration_radius   = malloc(sv->T_centre_de_masse * sizeof(double));

    // MSD
    sv->stock_msd = malloc(cfg->nbr_total_simu * sizeof(double**));
    for (int i = 0; i < cfg->nbr_total_simu; i++) {
        sv->stock_msd[i] = malloc(sv->T_msd * sizeof(double*));
        for (int j = 0; j < sv->T_msd; j++)
            sv->stock_msd[i][j] = calloc(cfg->Nm, sizeof(double));
    }

    // Trajectoires des monomères
    sv->R_monomere_arrays = malloc(cfg->nbr_total_simu * sizeof(double***));
    for (int k = 0; k < cfg->nbr_total_simu; k++) {
        sv->R_monomere_arrays[k] = malloc(sv->T_msd * sizeof(double**));
        for (int i = 0; i < cfg->T_msd; i++) {
            sv->R_monomere_arrays[k][i] = malloc(cfg->Nm * sizeof(double*));
            for (int j = 0; j < cfg->Nm; j++)
                sv->R_monomere_arrays[k][i][j] = calloc(3, sizeof(double));
        }
    }

    // Initialisation de time/log10_time
    for (int i = 0; i < sv->T_correlation; i++) {
        sv->time[i]       = (double)i;
        sv->log10_time[i] = log10(sv->time[i] + 1e-12); // éviter log10(0)
    }

    // ---- RNAP


     if(cfg->nb_rnap_initial == 0){
        sv->nb_rnap = 0;
    } else {
        sv->nb_rnap = 1;
    }
    //sv->nb_rnap = sv->nb_rnap_initial;

    sv->sortie = 0;


    sv->positions_bille_rnap = calloc(sv->nb_rnap, sizeof(int));

    sv->R_rnap = matrix_rnap(sv->nb_rnap, cfg->rnap_subunits); // Matrice cube (rnap, sub_units, positons sub units) 
    sv->R_rnap_new = matrix_rnap(sv->nb_rnap, cfg->rnap_subunits); // Matrice de stockage

    sv->avancement_transcription = calloc(sv->nb_rnap, sizeof(double));

    sv->compteur_mono_rnap = calloc(cfg->N, sizeof(int));

    sv->num_rnap = calloc(cfg->nb_rnap_initial, sizeof(int)); 

    for (int i = 0; i < cfg->nb_rnap_initial; i++){
        sv->num_rnap[i] = -1;
    }

    // Compteur 

    sv->premier_passage = 1;
    sv->compteur_grand_deplacement = 0;
    sv->compteur_train = 0; 
    sv->variable_test_2 = 0;
    sv->variable_test_3 = 0;
}

void cleanup_sim_vars(SimVars *sv, Config *cfg){

    // --- Chromatine
    free_matrix_cube_if_allocated(sv->stock_msd, cfg->nbr_total_simu, cfg->T_msd);
    free_matrix_if_allocated(&sv->Rbb, cfg->T_correlation);
    free_if_allocated((void **)&sv->list_monomere);
    free_matrix_if_allocated(&sv->R, cfg->N);
    free_matrix_if_allocated(&sv->t_link, cfg->N-1);
    free_matrix_if_allocated(&sv->bending_forces,cfg->N);
    free_matrix_rnap(sv->R_rnap, sv->nb_rnap, cfg->rnap_subunits);
    free_matrix_rnap(sv->R_rnap_new, sv->nb_rnap, cfg->rnap_subunits);
}

void ajouter_rnap(SimVars *sv, const Config *cfg, NeighborList *neighbor_list, NeighborList_rnap ***p_neighbor_lists_rnap, int t){

    NeighborList_rnap **neighbor_lists_rnap = *p_neighbor_lists_rnap;
     
    if(sv->nb_rnap == 0){

        sv->num_rnap[sv->nb_rnap] = 0;
        // printf("num rnap = %d \n ",sv->num_rnap[sv->nb_rnap]);

        sv->nb_rnap +=1;
        
        sv->R_rnap = matrix_rnap(sv->nb_rnap, cfg->rnap_subunits); // Matrice cube (rnap, sub_units, positons sub units) 
        sv->R_rnap_new = matrix_rnap(sv->nb_rnap, cfg->rnap_subunits); // Matrice de stockage


        sv->positions_bille_rnap = calloc(sv->nb_rnap, sizeof(int));
        sv->avancement_transcription = calloc(sv->nb_rnap, sizeof(double));

        creation_1_rnap_erwan(sv->R, sv->nb_rnap, sv->positions_bille_rnap, sv->R_rnap, cfg->N, cfg->rnap_subunits, cfg->debut_segment, cfg->fin_segment, cfg->nb_rnap_initial);

    }

    if( (sv->nb_rnap > 0) && (sv->positions_bille_rnap[0] > cfg->debut_segment + 1) && (sv->nb_rnap < cfg->nb_rnap_initial) ){
                sv->nb_rnap += 1;

                for(int i = 0; i < sv->nb_rnap; i++){
                    sv->num_rnap[i] = sv->nb_rnap - i -1;
                }

                if(sv->nb_rnap == cfg->nb_rnap_initial){
                    sv->premier_passage = 0;
                    printf("A t = %d premier passage = %d \n \n \n", t, sv->premier_passage);
                    sv->variable_test_3 = 0;
                }

                double ***tmp_rnap_new = realloc(sv->R_rnap_new, sv->nb_rnap * sizeof(double **));
                if (tmp_rnap_new == NULL) {
                    perror("Erreur realloc R_rnap_new");
                    exit(EXIT_FAILURE);
                }
                sv->R_rnap_new = tmp_rnap_new;

                sv->R_rnap_new[sv->nb_rnap-1] = malloc(cfg->rnap_subunits * sizeof(double *));
                if(sv->R_rnap_new[sv->nb_rnap-1] == NULL){
                    perror("Erreur malloc R_rnap_new");
                    exit(EXIT_FAILURE);
                }
                
                for(int i = 0; i < cfg->rnap_subunits; i++ ){
                    sv->R_rnap_new[sv->nb_rnap - 1][i] = malloc(3 * sizeof(double));
                    if(sv->R_rnap_new[sv->nb_rnap - 1][i] == NULL){
                        perror("Erreur malloc R_rnap_new");
                        exit(EXIT_FAILURE);
                    }
                }

                for(int i = 0; i < cfg->rnap_subunits; i++){
                    for(int j = 0; j < 3; j++){
                        sv->R_rnap_new[sv->nb_rnap - 1][i][j] = 0;
                    }
                }

                double ***tmp_rnap = realloc(sv->R_rnap, sv->nb_rnap * sizeof(double **));
                if (tmp_rnap == NULL) {
                    perror("Erreur realloc R_rnap");
                    exit(EXIT_FAILURE);
                }
                sv->R_rnap = tmp_rnap;

                sv->R_rnap[sv->nb_rnap-1] = malloc(cfg->rnap_subunits * sizeof(double *));
                if(sv->R_rnap[sv->nb_rnap - 1] == NULL){
                    perror("Erreur malloc R_rnap");
                    exit(EXIT_FAILURE);
                }

                for(int i = 0; i < cfg->rnap_subunits; i++ ){
                    sv->R_rnap[sv->nb_rnap - 1][i] = malloc(3 * sizeof(double));
                    if(sv->R_rnap[sv->nb_rnap - 1][i] == NULL){
                        perror("Erreur malloc R_rnap"); 
                        exit(EXIT_FAILURE);
                    }
                }
                for(int i = 0; i < cfg->rnap_subunits; i++){
                    for(int j = 0; j < 3; j++){
                        sv->R_rnap[sv->nb_rnap - 1][i][j] = 0;
                    }
                }

                int* tmp_positions = realloc(sv->positions_bille_rnap, sv->nb_rnap * sizeof(int));
                if (tmp_positions == NULL) {
                    perror("Erreur realloc positions_bille_rnap");
                    exit(EXIT_FAILURE);
                }

                

                double* tmp_avancement = realloc(sv->avancement_transcription, sv->nb_rnap * sizeof(double));
                if (tmp_avancement == NULL) {
                    perror("Erreur realloc avancement_transcription");
                    exit(EXIT_FAILURE);
                }

                

                
                // Mise à jour des pointeurs après réussite
                sv->R_rnap_new = tmp_rnap_new;
                sv->R_rnap = tmp_rnap;
                sv->positions_bille_rnap = tmp_positions;
                sv->avancement_transcription = tmp_avancement;
                sv->avancement_transcription[sv->nb_rnap - 1] = 0;
                sv->positions_bille_rnap[sv->nb_rnap - 1] = 0;
                

                for( int i = sv->nb_rnap - 1; i > 0; i--){
                    sv->avancement_transcription[i] = sv->avancement_transcription[i-1];
                    sv->positions_bille_rnap[i] = sv->positions_bille_rnap[i-1];
                    for(int j = 0; j < cfg->rnap_subunits; j++){
                        for (int k = 0; k < 3; k++){
                            sv->R_rnap[i][j][k] = sv->R_rnap[i-1][j][k];
                            sv->R_rnap_new[i][j][k] = sv->R_rnap_new[i-1][j][k];
                        }
                    }
                }
                sv->avancement_transcription[0] = 0;
               
                creation_1_rnap_erwan(sv->R, sv->nb_rnap, sv->positions_bille_rnap, sv->R_rnap, cfg->N, cfg->rnap_subunits, cfg->debut_segment, cfg->fin_segment, cfg->nb_rnap_initial);
                
    }
    build_neighbor_list(sv->R, neighbor_list, cfg->N, 2, 0);
    if(cfg->nb_rnap_initial > 0) {
        resize_neighbor_list_rnap(&neighbor_lists_rnap, sv->nb_rnap-1, sv->nb_rnap, cfg->rnap_subunits, t);
        *p_neighbor_lists_rnap = neighbor_lists_rnap;
        build_neighbor_list_rnap_chrom(sv->R_rnap, sv->nb_rnap, sv->R, cfg->N, neighbor_lists_rnap, cfg->rayon_ecrantage_LJ_chrom, t);}
}

void retirer_rnap(SimVars *sv, const Config *cfg, NeighborList *neighbor_list, NeighborList_rnap** neighbor_lists_rnap, int rnap, int prout, int t){

    sv->positions_bille_rnap[rnap] += 1; 

    sv->compteur_mono_rnap[sv->positions_bille_rnap[rnap]]+=1;
    sv->avancement_transcription[rnap] = 0;



                    if(sv->positions_bille_rnap[rnap] > cfg->fin_segment){

                    if (sv->nb_rnap > 1){
                        prout = 1;
                        
                        int last_index = sv->nb_rnap - 1;
                        if(last_index == 0){
                            sv->variable_test_2 = 1;
                            sv->premier_passage = 0; 
                        }
                        sv->nb_rnap --;
                        
                        double ***tmp_rnap_new = realloc(sv->R_rnap_new, sv->nb_rnap * sizeof(double **));
                      

                        if(tmp_rnap_new != NULL){
                            sv->R_rnap_new = tmp_rnap_new;
                        }
                        else{
                            perror("Erreur lors du realloc de R_rnap_new 1 ");
                            exit(EXIT_FAILURE);
                        }

                        double ***tmp_rnap = realloc(sv->R_rnap, sv->nb_rnap * sizeof(double **));

                        if(tmp_rnap != NULL){
                            sv->R_rnap = tmp_rnap;
                        }
                        else{
                            perror("Erreur lors du realloc de R_rnap 1 ");
                            exit(EXIT_FAILURE);
                        }
                    }else {sv->nb_rnap--;}
                    
                    if(sv->nb_rnap == 0){
                        free_matrix_rnap(sv->R_rnap, sv->nb_rnap, cfg->rnap_subunits);
                        sv->R_rnap = NULL;
                        free_matrix_rnap(sv->R_rnap_new, sv->nb_rnap, cfg->rnap_subunits);
                        sv->R_rnap_new = NULL;
                        free(sv->positions_bille_rnap);
                        sv->positions_bille_rnap = NULL; 
                        free(sv->avancement_transcription);
                        sv->avancement_transcription = NULL; 
                        sv->premier_passage = 0;
                    }
                    build_neighbor_list(sv->R, neighbor_list, cfg->N, 2, 0);
                    resize_neighbor_list_rnap(&neighbor_lists_rnap, sv->nb_rnap-1, sv->nb_rnap, cfg->rnap_subunits, t );
                    build_neighbor_list_rnap_chrom(sv->R_rnap, sv->nb_rnap, sv->R, cfg->N, neighbor_lists_rnap, cfg->rayon_ecrantage_LJ_chrom, t);
                }

}

void retirer_rnap_2(SimVars *sv, const Config *cfg,
                    NeighborList *neighborlist,
                    NeighborList_rnap*** p_lists,
                    int rnap, int prout, int t)
{
    // 1) Mise à jour des compteurs de position/avancement
    sv->positions_bille_rnap[rnap] += 1;
    sv->compteur_mono_rnap[sv->positions_bille_rnap[rnap]]++;
    sv->avancement_transcription[rnap] = 0;
    if (sv->positions_bille_rnap[rnap] < cfg->fin_segment)
        return;  // pas encore sorti → rien à faire

    // printf("nb_rnap = %d \n ", sv->nb_rnap);
    // printf("num rnap = %d \n ", sv->num_rnap[cfg->nb_rnap_initial - sv->nb_rnap]);
    sv->num_rnap[sv->nb_rnap - 1] = - 1;
    sv->sortie = 1;

    // printf("Apres num rnap = %d \n ", sv->num_rnap[cfg->nb_rnap_initial - sv->nb_rnap]);
    int old_nb = sv->nb_rnap;
    if (old_nb <= 0) return;
    int new_nb = old_nb - 1;
    sv->nb_rnap = new_nb;

    // 2) On récupère le tableau de listes
    NeighborList_rnap **lists = *p_lists;
    int rows = cfg->rnap_subunits;

    // 3) Libération de la sous-liste à l’indice 'rnap'
    for (int j = 0; j < rows; ++j)
        free(lists[rnap][j].neighbors);
    free(lists[rnap]);

    // 4) Décalage des pointeurs vers la gauche pour combler le trou
    for (int i = rnap; i < old_nb - 1; ++i)
        lists[i] = lists[i+1];

    // 5) Realloc du tableau principal à new_nb
    if (new_nb > 0) {
      lists = realloc(lists, new_nb * sizeof *lists);
      if (!lists) { perror("realloc neighbor_lists_rnap"); exit(EXIT_FAILURE); }
    } else {
      free(lists);
      lists = NULL;
    }
    *p_lists = lists;

    // 6) Même chose pour positions_bille_rnap et avancement_transcription
    sv->positions_bille_rnap = realloc(
      sv->positions_bille_rnap,
      new_nb * sizeof *sv->positions_bille_rnap
    );
    sv->avancement_transcription = realloc(
      sv->avancement_transcription,
      new_nb * sizeof *sv->avancement_transcription
    );

    // 7) Et pour R_rnap et R_rnap_new
    //    - on libère la dernière matrice (celle qui n’existe plus)
    for (int sub = 0; sub < rows; ++sub) {
        free(sv->R_rnap[old_nb-1][sub]);
        free(sv->R_rnap_new[old_nb-1][sub]);
    }
    free(sv->R_rnap[old_nb-1]);
    free(sv->R_rnap_new[old_nb-1]);
    //    - on décale les pointeurs 
    for (int i = rnap; i < old_nb - 1; ++i) {
      sv->R_rnap[i]     = sv->R_rnap[i+1];
      sv->R_rnap_new[i] = sv->R_rnap_new[i+1];
    }
    //    - puis realloc
    if (new_nb > 0) {
      sv->R_rnap     = realloc(sv->R_rnap,     new_nb * sizeof *sv->R_rnap);
      sv->R_rnap_new = realloc(sv->R_rnap_new, new_nb * sizeof *sv->R_rnap_new);
      if (!sv->R_rnap || !sv->R_rnap_new) {
        perror("realloc R_rnap"); exit(EXIT_FAILURE);
      }
    } else {
      free(sv->R_rnap);     sv->R_rnap     = NULL;
      free(sv->R_rnap_new); sv->R_rnap_new = NULL;
    }

    // 8) On reconstruit les neighbor-lists pour la chromatine et les RNAPs restants
    build_neighbor_list(sv->R, neighborlist, cfg->N, 2, 0);
    if (new_nb > 0) {
      build_neighbor_list_rnap_chrom(
         sv->R_rnap, new_nb, sv->R,
         cfg->N, lists, cfg->rayon_ecrantage_LJ_chrom, t
      );
    }
}


void enregistrement_data(SimVars *sv, const Config *cfg, const Files *f, int t, double endtoend, int nbr_simu){

    if(t%cfg->periode_msd == 0){
            for(int i = 0; i < cfg->Nm; i++){
                for(int j = 0; j < 3; j++){
                    sv->R_monomere_arrays[nbr_simu][(int)t/cfg->periode_msd][i][j] = sv->R[sv->list_monomere[i]][j];
                }
            }
        }

    if(t%cfg->periode_correlation == 0){
            for ( int i = 0; i < 3; i++){
                sv->Rbb[t/cfg->periode_correlation][i] = sv->R[cfg->N-2][i] - sv->R[0][i];
                sv->R_segment[t/cfg->periode_correlation][i] = sv->R[cfg->fin_segment][i] - sv->R[cfg->debut_segment][i];
                
            }
            fprintf(f->fichier_endtoend, "%f %d \n", distance(sv->R[0], sv->R[cfg->N-1]), t/cfg->periode_correlation);
            fprintf(f->fichier_correl_segment, "%f %d \n", distance(sv->R[cfg->debut_segment], sv->R[cfg->fin_segment]), t/cfg->periode_correlation);


        }
    
        if(t%cfg->periode_centre_de_masse == 0){
            for(int i = 0; i < cfg->N; i ++){
                
                sv->R_centre_de_masse[(int)t/cfg->periode_centre_de_masse][0] += sv->R[i][0];
                sv->R_centre_de_masse[(int)t/cfg->periode_centre_de_masse][1] += sv->R[i][1];
                sv->R_centre_de_masse[(int)t/cfg->periode_centre_de_masse][2] += sv->R[i][2];

                sv->stock_cdm[(int)t/cfg->periode_centre_de_masse][0] = sv->R[i][0];
                sv->stock_cdm[(int)t/cfg->periode_centre_de_masse][1] = sv->R[i][1];
                sv->stock_cdm[(int)t/cfg->periode_centre_de_masse][2] = sv->R[i][2];
                 
            }
        }

        if(t%cfg->periode_endtoend == 0){

            sv->stock[t/cfg->periode_endtoend][0] = distance(sv->R[0], sv->R[cfg->N-1]);
            sv->stock[t/cfg->periode_endtoend][1] = t/cfg->periode_endtoend; 

            sv->Rbb_segment[t/cfg->periode_endtoend][0] = distance(sv->R[(3*cfg->N/10)], sv->R[(4*cfg->N)/10]);
            sv->Rbb_segment[t/cfg->periode_endtoend][1] = t/cfg->periode_endtoend;
            fprintf(f->fichier_endtoend_segment, "%f %d \n", distance(sv->R[(3*cfg->N)/10], sv->R[(4*cfg->N)/10]), t/cfg->periode_endtoend);

            sv->Rbb_avant[t/cfg->periode_endtoend][0] = distance(sv->R[(2*cfg->N)/10], sv->R[(3*cfg->N)/10]);
            sv->Rbb_avant[t/cfg->periode_endtoend][1] = t/cfg->periode_endtoend;
            fprintf(f->fichier_endtoend_avant, "%f %d \n", distance(sv->R[(2*cfg->N)/10], sv->R[(3*cfg->N)/10]), t/cfg->periode_endtoend);

            sv->Rbb_apres[t/cfg->periode_endtoend][0] = distance(sv->R[(4*cfg->N)/10], sv->R[(5*cfg->N)/10]);
            sv->Rbb_apres[t/cfg->periode_endtoend][1] = t/cfg->periode_endtoend;
            fprintf(f->fichier_endtoend_apres, "%f %d \n", distance(sv->R[(4*cfg->N)/10], sv->R[(5*cfg->N)/10]), t/cfg->periode_endtoend);

            endtoend += distance(sv->R[0], sv->R[cfg->N-1]);
        }

        if(t%cfg->periode_enregistrement == 0){
            if(sv->nb_rnap > 0){
                enregistrement_RNAP(f->fichier, sv->R, cfg->N, sv->R_rnap, sv->nb_rnap, cfg->T_eq + t, cfg->mono_transcrpt, sv->positions_bille_rnap, cfg->rnap_subunits, cfg->nb_rnap_initial, sv->num_rnap, sv->sortie); 
                enregistrement_RNAP_position(f->fichier_rnap, sv->nb_rnap, cfg->T_eq + t, sv->positions_bille_rnap, sv->avancement_transcription[0]);
            }
            else{
                enregistrement(f->fichier, sv->R, cfg->N, cfg->T_eq + t);
            }
        }

}

void finalize_simulationaaaaa(SimVars *sv, const Config *cfg, const Files *f, double endtoend, int nbr_simu){

    printf("End to end distance Gaussian : %lf \n", endtoend);

    printf("Simulation %d Nombre de deplacements totals : %d Nombre de deplacement plus grand que 0.1a : %d \n", nbr_simu, cfg->T, sv->compteur_grand_deplacement);

    calculate_autocorrelation(sv->Rbb, sv->stock_correlation, cfg->T_correlation, cfg->N);
    calculate_autocorrelation(sv->R_segment, sv->stock_correlation_segment, cfg->T_correlation, cfg->N);

    calculate_msd(sv->R_monomere_arrays, sv->stock_msd, cfg->T_msd, cfg->Nm, sv->list_monomere, nbr_simu, cfg->nbr_total_simu);

    for(int i = 0; i < cfg->T_centre_de_masse; i++){
        for(int j = 0; j < 3; j++){
            sv->R_centre_de_masse[i][j] /= cfg->N; 
        }
    }

    save_centre_de_masse(sv->R_centre_de_masse, cfg->N, cfg->T_centre_de_masse, nbr_simu, cfg->periode_centre_de_masse);
    calculate_gyration_radius(sv->R_centre_de_masse, sv->gyration_radius, sv->R, cfg->N, cfg->T_centre_de_masse, sv->stock_cdm);
    save_gyration(sv->gyration_radius, nbr_simu, cfg->T_centre_de_masse, cfg->periode_centre_de_masse);

    for(int i = 0; i < cfg->N; i++){
        printf("monomere %d a vu %d rnap \n",i, sv->compteur_mono_rnap[i]);
    }
}

void calcul(SimVars *sv, const Config *cfg, const Files *f, NeighborList *neighbor_lists, NeighborList_rnap **neighbor_lists_rnap, int nbr_simu, double endtoend){

    clock_t start_2, end_2; 
    clock_t start, end;double duree_boucle, duree_tot = 0, temps_restant;

    struct timespec last, now;
    double interval = 60; // en secondes

    clock_gettime(CLOCK_MONOTONIC, &last);


    for (int t = 0; t < cfg->T; t++){
        // printf("nb_rnap = %d premier_passage = %d \n",sv->nb_rnap, sv->premier_passage);

        clock_gettime(CLOCK_MONOTONIC, &now); 
        double elapsed = now.tv_sec - last.tv_sec + (now.tv_nsec - last.tv_nsec)*1e-9; 

        if (elapsed >= interval){
            printf("Itération %d  (après %.1f s)\n", t, elapsed);
            last = now;

        }

        if( cfg->nb_rnap_initial > 0){
            if( /*(sv->variable_test_3 == 1 && sv->nb_rnap < cfg->nb_rnap_initial) ||*/ (sv->premier_passage == 1 && sv->nb_rnap < cfg->nb_rnap_initial) ){
                ajouter_rnap(sv, cfg, neighbor_lists, &neighbor_lists_rnap, t);
            }
        }
        start = clock();

        if (t%1000==0){
            build_neighbor_list(sv->R, neighbor_lists, cfg->N, 2, 0);

            if(cfg->nb_rnap_initial > 0) {
                build_neighbor_list_rnap_chrom(sv->R_rnap, sv->nb_rnap, sv->R, cfg->N, neighbor_lists_rnap, cfg->rayon_ecrantage_LJ_chrom, t);
            }

            for(int i = 0; i < cfg->N; i++){
                fprintf(f->fichier_voisin,"%d ",neighbor_lists[i].count);
            }
            fprintf(f->fichier_voisin, "\n");
        }
        sv->R_new = polymere_brownian_motion(sv->R, cfg->K, cfg->Delta, cfg->N, sv->R_new, cfg->K_bend, sv->bending_forces, cfg->attache, cfg->plan, t, f->test, cfg->bending, sv->truc, cfg->T, f->fichier_force, cfg->periode_force, f->fichier_force_thermique, cfg->temperature);

        //update_link_vectors(R, t_link, N);
        //f_bending_forces(R, t_link, bending_forces, K_bend, N, t);

        if(cfg->nb_rnap_initial > 0) {matrix_rnap_0(sv->R_rnap_new, sv->nb_rnap, cfg->rnap_subunits);}// Remplis une matrice 3D de 0

        /////////// Boucle sur les RNAPS /////////////
        if(cfg->nb_rnap_initial > 0){
            for(int rnap = 0; rnap < sv->nb_rnap; rnap++){

                int prout = 0;
                
                polymere_brownian_motion_ring_force(sv->R_rnap[rnap], cfg->alpha, 1 * cfg->K_rnap, cfg->Delta, cfg->rnap_subunits, rnap, sv->R_rnap_new[rnap], f->test, t, cfg->periode_force, f->fichier_force_rnap, f->fichier_force_thermique, cfg->temperature);
                
                liaison_sup(3*cfg->a_transpt, 2*cfg->K_rnap, cfg->Delta, sv->R_rnap_new[rnap], 0, 4, f->fichier_force_rnap_2, t, cfg->periode_force); // Liaison entre les monomères opposés
                liaison_sup(3*cfg->a_transpt, 2*cfg->K_rnap, cfg->Delta, sv->R_rnap_new[rnap], 2, 6, f->fichier_force_rnap_2, t, cfg->periode_force); // Liaison entre les monomères opposés
                liaison_sup(3*cfg->a_transpt, 2*cfg->K_rnap, cfg->Delta, sv->R_rnap_new[rnap], 1, 5, f->fichier_force_rnap_2, t, cfg->periode_force); // Liaison entre les monomères opposés
                liaison_sup(3*cfg->a_transpt, 2*cfg->K_rnap, cfg->Delta, sv->R_rnap_new[rnap], 3, 7, f->fichier_force_rnap_2, t, cfg->periode_force); // Liaison entre les monomères opposés
                
                sv->avancement_transcription[rnap] += cfg->dx_avancement_rnap;
                sv->R_rnap_new[rnap] = bond_rnap_bead_progressive_mvt(sv->R, sv->R_rnap[rnap], cfg->a_transpt, cfg->K_transpt, cfg->Delta, sv->R_rnap_new[rnap], sv->positions_bille_rnap[rnap], sv->avancement_transcription[rnap], cfg->a, cfg->alpha, f->fichier_force_lea, cfg->periode_force, t);
                copy_matrix(sv->R_rnap[rnap], sv->R_rnap_new[rnap],cfg->rnap_subunits);

                if(1-sv->avancement_transcription[rnap]<0.0000001){

                    retirer_rnap_2(sv, cfg, neighbor_lists, &neighbor_lists_rnap, rnap, prout, t);
                }  
            }
        }
        
        if(sv->variable_test_2 == 1){
            sv->compteur_train +=1; 
            if(sv->compteur_train > cfg->attente_train ){
                sv->variable_test_3 = 1; 
                sv->variable_test_2 = 0; 
                sv->compteur_train = 0;
            }
        }

        sv->R = sv->R_new;
        
        lennard_jones_forces(sv->R, neighbor_lists, cfg->N, cfg->epsilon, cfg->sigma6, cfg->sigma12, cfg->Delta, cfg->attache, cfg->periode_force, f->fichier_force_LJ, t);
        lennard_jones_forces_rnap(sv->R_rnap, sv->nb_rnap, sv->R, cfg->N, neighbor_lists_rnap, cfg->epsilon_rnap, cfg->sigma6_rnap, cfg->sigma12_rnap, cfg->sigma6_rnap2, cfg->sigma12_rnap2, cfg->rayon_ecrantage_LJ_rnap, cfg->Delta, t, f->test, cfg->T, f->fichier_force_rnap_LJ, cfg->periode_force);
        lennard_jones_forces_rnap_rnap(sv->R_rnap, sv->nb_rnap, cfg->epsilon, cfg->sigma6_rnap, cfg->sigma12_rnap, cfg->rayon_ecrantage_LJ_rnap, cfg->Delta);
        compteur_grands_deplacements(cfg->N, cfg->T, sv->R, sv->R_new, sv->compteur_grand_deplacement);


        if(cfg->confinement == 1){
            confinement_sphere(sv->R, cfg->N, cfg->r_sphere);
        }

        enregistrement_data(sv, cfg, f, t, endtoend, nbr_simu);

        end = clock();  
        duree_boucle = (double)(end - start)/CLOCKS_PER_SEC;
        duree_tot += duree_boucle ;

        if (t%(cfg->T/10) == 0){
            Mesures mesures = calcul_mesures (sv->R, cfg->N);
            double duree_min = (int)(duree_tot / 60);
            double duree_sec = duree_tot - (duree_min * 60);
            temps_restant = duree_boucle * (cfg->T-t-1) / 60;
            printf("%d/%.d %.f:%.f %.2fmin std %.10f moy %.10f      %d/%d \n",t, cfg->T, duree_min, duree_sec, temps_restant, mesures.std, mesures.moyenne, nbr_simu, cfg->nbr_total_simu);   
        }
    }
}


void f_equilibriate(SimVars *sv, const Config *cfg, const Files *f, NeighborList *neighbor_lists, NeighborList_rnap **neighbor_lists_rnap, int nbr_simu){
    printf("//// Avant la simulation //// \n");

    clock_t start_2, end_2; 
    clock_t start, end;double duree_boucle, duree_tot = 0, temps_restant;

    for (int t = 0; t < cfg->T_eq; t++){

        start = clock();

        if (t%1000==0){
            build_neighbor_list(sv->R, neighbor_lists, cfg->N, 2, 0);

            if(cfg->nb_rnap_initial > 0) {build_neighbor_list_rnap_chrom(sv->R_rnap, sv->nb_rnap, sv->R, cfg->N, neighbor_lists_rnap, cfg->rayon_ecrantage_LJ_chrom, t);}

            for(int i = 0; i < cfg->N; i++){
                fprintf(f->fichier_voisin,"%d ",neighbor_lists[i].count);
            }
            fprintf(f->fichier_voisin, "\n");
        }

        sv->R_new = polymere_brownian_motion(sv->R, cfg->K, cfg->Delta, cfg->N, sv->R_new, cfg->K_bend, sv->bending_forces, cfg->attache, cfg->plan, t, f->test, cfg->bending, sv->truc, cfg->T, f->fichier_force, cfg->periode_force, f->fichier_force_thermique, cfg->temperature);

        //update_link_vectors(R, t_link, N);
        //f_bending_forces(R, t_link, bending_forces, K_bend, N, t);


        /////////// Boucle sur les RNAPS /////////////

        sv->R = sv->R_new;
        
        lennard_jones_forces(sv->R, neighbor_lists, cfg->N, cfg->epsilon, cfg->sigma6, cfg->sigma12, cfg->Delta, cfg->attache, cfg->periode_force, f->fichier_force_LJ, t);
        lennard_jones_forces_rnap(sv->R_rnap, sv->nb_rnap, sv->R, cfg->N, neighbor_lists_rnap, cfg->epsilon_rnap, cfg->sigma6_rnap, cfg->sigma12_rnap, cfg->sigma6_rnap2, cfg->sigma12_rnap2, cfg->rayon_ecrantage_LJ_rnap, cfg->Delta, t, f->test, cfg->T, f->fichier_force_rnap_LJ, cfg->periode_force);
        compteur_grands_deplacements(cfg->N, cfg->T, sv->R, sv->R_new, sv->compteur_grand_deplacement);

        if(cfg->confinement == 1){
            confinement_sphere(sv->R, cfg->N, cfg->r_sphere);
        }


        end = clock();  
        duree_boucle = (double)(end - start)/CLOCKS_PER_SEC;
        duree_tot += duree_boucle ;

        if (t%(cfg->T_eq/10) == 0){
            Mesures mesures = calcul_mesures (sv->R, cfg->N);
            double duree_min = (int)(duree_tot / 60);
            double duree_sec = duree_tot - (duree_min * 60);
            temps_restant = duree_boucle * (cfg->T-t-1) / 60;
            printf("%d/%.d %.f:%.f %.2fmin std %.10f moy %.10f      %d/%d \n",t, cfg->T_eq, duree_min, duree_sec, temps_restant, mesures.std, mesures.moyenne, nbr_simu, cfg->nbr_total_simu);   
        }

        if(t%cfg->periode_enregistrement == 0){
            if(sv->nb_rnap > 0){
                enregistrement_RNAP(f->fichier, sv->R, cfg->N, sv->R_rnap, sv->nb_rnap, t, cfg->mono_transcrpt, sv->positions_bille_rnap, cfg->rnap_subunits, cfg->nb_rnap_initial, sv->num_rnap, sv->sortie);
                enregistrement_RNAP_position(f->fichier_rnap, sv->nb_rnap, t, sv->positions_bille_rnap, sv->avancement_transcription[0]);
            }
            else{
                enregistrement(f->fichier, sv->R, cfg->N, t);
            }
        }
    }
}