#include "config.h"
#include "time.h" 
#include "simulation.h" 
#include "basic_functions.h"
#include "transcription.h"
#include "file.h"
#include "colors.h"
#include "mt19937ar.h"


#include <stdlib.h>  // pour malloc, calloc
#include <math.h>    // pour log10()

#define MAX_RNAP 50
#define MAX_RNAP_SUBUNITS 8 
#define DIM 3


void init_sim_vars(SimVars *sv, Config *cfg) {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////                 Période                    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    sv->T_enregistrement = cfg->T_enregistrement;
    sv->T_msd = cfg->T_msd;
    sv->T_correlation = cfg->T_correlation;
    sv->T_endtoend = cfg->T_endtoend;
    sv->T_centre_de_masse = cfg->T_centre_de_masse;
    sv->T_force = cfg->T_force;
    sv->T_voisin = cfg->T_voisin; 

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////          Initialisation Tableaux           ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


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


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////           Corrélations et temps            ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    sv->stock_correlation         = calloc(sv->T_correlation, sizeof(double));
    sv->stock_correlation_segment = calloc(sv->T_correlation, sizeof(double));
    sv->time                      = malloc(sv->T_correlation * sizeof(double));
    sv->log10_time                = malloc(sv->T_correlation * sizeof(double));

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////              Stock end-to-end              ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    sv->stock         = allocate_matrix(sv->T_endtoend, 2);
    sv->Rbb_segment   = allocate_matrix(sv->T_endtoend, 3);
    sv->Rbb_avant     = allocate_matrix(sv->T_endtoend, 3);
    sv->Rbb_apres     = allocate_matrix(sv->T_endtoend, 3);
    sv->nombre_voisins= malloc(cfg->N * sizeof(int));
    sv->Rbb           = allocate_matrix(sv->T_correlation, 3);
    sv->R_segment     = allocate_matrix(sv->T_correlation, 3);


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////     Centre de masse & rayon de giration    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    sv->R_centre_de_masse = allocate_matrix(sv->T_centre_de_masse, 3);
    sv->stock_cdm         = allocate_matrix(sv->T_centre_de_masse, 3);
    sv->gyration_radius   = malloc(sv->T_centre_de_masse * sizeof(double));

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////                     MSD                    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    sv->stock_msd = malloc(cfg->nbr_total_simu * sizeof(double**));
    for (int i = 0; i < cfg->nbr_total_simu; i++) {
        sv->stock_msd[i] = malloc(sv->T_msd * sizeof(double*));
        for (int j = 0; j < sv->T_msd; j++)
            sv->stock_msd[i][j] = calloc(cfg->Nm, sizeof(double));
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////         Trajectoires des monomères         ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    sv->R_monomere_arrays = malloc(sv->T_msd * sizeof(double**));
    for (int i = 0; i < cfg->T_msd; i++) {
        sv->R_monomere_arrays[i] = malloc(cfg->Nm * sizeof(double*));
        for (int j = 0; j < cfg->Nm; j++)
            sv->R_monomere_arrays[i][j] = calloc(3, sizeof(double));
    }
    

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////      Initialisation de time/log10_time     ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (int i = 0; i < sv->T_correlation; i++) {
        sv->time[i]       = (double)i;
        sv->log10_time[i] = log10(sv->time[i] + 1e-12); // éviter log10(0)
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////      RNAP     ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    sv->nb_rnap = 0;

    sv->l_rnap = malloc(MAX_RNAP * sizeof(int)); // -1 non active // 1 active // -2 sortie + ne revient plus // 2 active deuxième passage // 3 active troisième passage ....
    if (!sv->l_rnap)
    {
        fprintf(stderr, "Erreur malloc l_rnap\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < MAX_RNAP; i ++ )
    {
        sv->l_rnap[i] = -1; 
    }

    sv->positions_bille_rnap = calloc(MAX_RNAP, sizeof(int)); 
    sv->R_rnap = allocate_matrix_3D(MAX_RNAP, cfg->rnap_subunits); 
    sv->R_rnap_new = allocate_matrix_3D(MAX_RNAP, cfg->rnap_subunits);
    sv->avancement_transcription = calloc(MAX_RNAP, sizeof(double));
    sv->compteur_mono_rnap = calloc(cfg->N, sizeof(int));
    sv->is_rnap = calloc(cfg->N, sizeof(int));
    sv->last_rnap_add_time = -1000000000;


    // Compteur 

    sv->compteur_grand_deplacement = 0;
    sv->compteur_train = 0; 
    sv->variable_test_2 = 0;
    sv->variable_test_3 = 0;


    sv->cdm_conf = malloc(3 * sizeof(double));
    sv->nb_move_large = 0;


}

void cleanup_sim_vars(SimVars *sv, Config *cfg)
{
    // --- Matrices simples ---
    free_matrix_if_allocated(&sv->bending_forces, cfg->N);
    free_matrix_if_allocated(&sv->R, cfg->N);
    free_matrix_if_allocated(&sv->R_new, cfg->N);
    free_matrix_if_allocated(&sv->t_link, cfg->N - 1);
    free_matrix_if_allocated(&sv->Rbb, sv->T_correlation);
    free_matrix_if_allocated(&sv->R_segment, sv->T_correlation);
    free_matrix_if_allocated(&sv->Rbb_segment, sv->T_endtoend);
    free_matrix_if_allocated(&sv->Rbb_avant, sv->T_endtoend);
    free_matrix_if_allocated(&sv->Rbb_apres, sv->T_endtoend);
    free_matrix_if_allocated(&sv->R_centre_de_masse, sv->T_centre_de_masse);
    free_matrix_if_allocated(&sv->stock_cdm, sv->T_centre_de_masse);
    free_matrix_if_allocated(&sv->stock, sv->T_endtoend);

    // --- Tableaux simples ---
    free_if_allocated((void**)&sv->nombre_voisins);
    free_if_allocated((void**)&sv->gyration_radius);
    free_if_allocated((void**)&sv->stock_correlation);
    free_if_allocated((void**)&sv->stock_correlation_segment);
    free_if_allocated((void**)&sv->time);
    free_if_allocated((void**)&sv->log10_time);
    free_if_allocated((void**)&sv->list_monomere);
    free_if_allocated((void**)&sv->compteur_mono_rnap);
    free_if_allocated((void**)&sv->R_matrix);

    // --- MSD ---
    free_matrix_cube_if_allocated(sv->stock_msd, cfg->nbr_total_simu, sv->T_msd);

    // --- RNAP ---
    free_matrix_3D(sv->R_rnap, MAX_RNAP, cfg->rnap_subunits);
    free_matrix_3D(sv->R_rnap_new, MAX_RNAP, cfg->rnap_subunits);
    free_if_allocated((void**)&sv->l_rnap);
    free_if_allocated((void**)&sv->positions_bille_rnap);
    free_if_allocated((void**)&sv->avancement_transcription);
    free_if_allocated((void**)&sv->is_rnap);

    // --- R_monomere_arrays ---
    if (sv->R_monomere_arrays)
    {
            if (sv->R_monomere_arrays)
            {
                for (int i = 0; i < sv->T_msd; i++)
                {
                    if (sv->R_monomere_arrays[i])
                    {
                        for (int j = 0; j < cfg->Nm; j++)
                        {
                            free_if_allocated((void**)&sv->R_monomere_arrays[i][j]);
                        }
                        free(sv->R_monomere_arrays[i]);
                    }
                }
            }
        free(sv->R_monomere_arrays);
    }
}

void ajouter_rnap(SimVars *sv, const Config *cfg,
                  NeighborList *neighbor_list,
                  NeighborList_rnap ***p_neighbor_lists_rnap,
                  int t) {

    NeighborList_rnap **neighbor_lists_rnap = *p_neighbor_lists_rnap;

    // --- Trouver un slot libre ---
    int new_id = -1;
    for (int i = 0; i < MAX_RNAP; i++) {
        if (sv->l_rnap[i] == -1) { // Libre et jamais sortie
            new_id = i;
            break;
        }
    }

    if (new_id == -1) {
        fprintf(stderr, C_RED "⚠️  Aucune RNAP disponible (MAX_RNAP=%d)\n" C_RESET, MAX_RNAP);
        return;
    }

    // --- Trouver la dernière RNAP active ---
    int last_active = -1;
    for (int i = MAX_RNAP - 1; i >= 0; i--) {
        if (sv->l_rnap[i] == 1) {
            last_active = i;
            break;
        }
    }

    // --- Initialisation de base ---
    sv->l_rnap[new_id] = 1;
    sv->positions_bille_rnap[new_id] = cfg->debut_segment;
    sv->avancement_transcription[new_id] = 0.0;
    sv->nb_rnap++;

    printf(C_CYAN "🧬 [t=%d] RNAP ajoutée (id=%d, nb_rnap=%d)\n" C_RESET,
           t, new_id, sv->nb_rnap);

    // --- Copier la position de la RNAP précédente pour continuité ---
    if (last_active >= 0) {
        for (int sub = 0; sub < cfg->rnap_subunits; sub++) {
            for (int d = 0; d < 3; d++) {
                sv->R_rnap[new_id][sub][d] = sv->R_rnap[last_active][sub][d];
                sv->R_rnap_new[new_id][sub][d] = sv->R_rnap_new[last_active][sub][d];
            }
        }
    }

    // --- Génération géométrique de la nouvelle RNAP ---
    creation_1_rnap_erwan(
        cfg, sv->R, new_id, sv->positions_bille_rnap, sv->R_rnap,
        cfg->N, cfg->rnap_subunits,
        cfg->debut_segment, cfg->fin_segment, cfg->nb_rnap_initial
    );

    // --- Mise à jour des neighbors ---
    build_neighbor_list(sv->R, neighbor_list, cfg->N, 2, 0);
    if (cfg->nb_rnap_initial > 0) {
        resize_neighbor_list_rnap(&neighbor_lists_rnap,
                                  sv->nb_rnap - 1,
                                  sv->nb_rnap,
                                  cfg->rnap_subunits,
                                  t);
        *p_neighbor_lists_rnap = neighbor_lists_rnap;
        build_neighbor_list_rnap_chrom(
            cfg, sv->R_rnap, sv->nb_rnap, sv->R,
            cfg->N, neighbor_lists_rnap, cfg->rayon_ecrantage_LJ_chrom, t
        );
    }
}



void retirer_rnap(SimVars *sv, const Config *cfg,
                  NeighborList *neighborlist,
                  NeighborList_rnap ***p_lists,
                  int rnap, int prout, int t)
{
    if (rnap < 0 || rnap >= MAX_RNAP) {
        fprintf(stderr, C_RED "⚠️ Retrait RNAP invalide : index %d hors limites (MAX_RNAP=%d)\n" C_RESET,
                rnap, MAX_RNAP);
        return;
    }

    // --- Mise à jour logique ---
    sv->positions_bille_rnap[rnap] += 1;
    if (sv->positions_bille_rnap[rnap] < cfg->N)
        sv->compteur_mono_rnap[sv->positions_bille_rnap[rnap]]++;

    sv->avancement_transcription[rnap] = 0.0;

    // --- Si la RNAP est encore dans le segment actif, on ne la retire pas encore
    if (sv->positions_bille_rnap[rnap] < cfg->fin_segment)
        return;
    for (int i = 0; i < cfg->nb_rnap_initial; i++)
        {
            printf("Retire : l_rnap[%d]=%d\n",i,sv->l_rnap[i]);
        }

    // --- Marque la RNAP comme inactive ---
    sv->l_rnap[rnap] = -2;
    sv->sortie = 1;

    if (sv->nb_rnap > 0)
        sv->nb_rnap--;

    printf(C_MAGENTA "\n=== [t=%d] RNAP %d retirée ===\n" C_RESET, t, rnap);
    printf(C_MAGENTA "   → nb_rnap restant = %d\n" C_RESET, sv->nb_rnap);

    // --- Nettoyage complet des coordonnées ---
    // ⚠️ Mettre à (0,0,0) *et* loin du domaine pour éviter les interactions
    const double COORD_INACTIVE = -1e6; // sécurité pour sortir de la simulation

    for (int sub = 0; sub < cfg->rnap_subunits; ++sub) {
        for (int d = 0; d < 3; ++d) {
            sv->R_rnap[rnap][sub][d]     = COORD_INACTIVE;
            sv->R_rnap_new[rnap][sub][d] = COORD_INACTIVE;
        }
    }

    sv->positions_bille_rnap[rnap] = 0;
    sv->avancement_transcription[rnap] = 0.0;

    printf("   🧼 Données RNAP[%d] réinitialisées à (%.1e,%.1e,%.1e)\n",
           rnap, COORD_INACTIVE, COORD_INACTIVE, COORD_INACTIVE);

    // --- Reconstruction des neighbor lists ---
    build_neighbor_list(sv->R, neighborlist, cfg->N, 2, 0);

    // Comptage des RNAP encore actives
    int nb_actives = 0;
    for (int i = 0; i < MAX_RNAP; ++i)
        if (sv->l_rnap[i] == 1)
            nb_actives++;

    NeighborList_rnap **lists = *p_lists;
    if (nb_actives > 0) {
        build_neighbor_list_rnap_chrom(
            cfg, sv->R_rnap, nb_actives, sv->R,
            cfg->N, lists, cfg->rayon_ecrantage_LJ_chrom, t
        );
        printf(C_GREEN "   🔁 Reconstruction des neighbor lists (%d actives)\n" C_RESET, nb_actives);
    } else {
        printf(C_YELLOW "   ⚠️ Aucune RNAP active restante\n" C_RESET);
    }

    // --- Vérif de cohérence ---
    if (sv->nb_rnap != nb_actives) {
        fprintf(stderr, C_RED "⚠️ Incohérence détectée : sv->nb_rnap=%d mais %d actives selon l_rnap[]\n" C_RESET,
                sv->nb_rnap, nb_actives);
    }

    // --- Affichage d’état ---
    printf("   État RNAP : ");
    for (int i = 0; i < MAX_RNAP; i++)
        printf("%d", sv->l_rnap[i]);
    printf("\n\n");
}




void creation_1_rnap_erwan(const Config *cfg, double **R, int id_rnap, int *positions_bille_rnap,
                           double ***R_rnap, int N, int rnap_subunits,
                           int debut_segment, int fin_segment, int nb_rnap_initial) {
    
    if (id_rnap < 0 || id_rnap >= MAX_RNAP) return;

    int promoteur;
    if(cfg->quench == 0)
    {
        promoteur = debut_segment;
    }
    else
    {
        promoteur = positions_bille_rnap[id_rnap];
        printf("promoteur = %d \n ", promoteur);
    }



    double taille_pol = 2.0;
    double threshold = 1.0;

    double x0 = R[promoteur][0], y0 = R[promoteur][1], z0 = R[promoteur][2];

    // Axe local autour du promoteur
    double axe[3] = {
        R[promoteur + 2][0] - R[promoteur - 2][0],
        R[promoteur + 2][1] - R[promoteur - 2][1],
        R[promoteur + 2][2] - R[promoteur - 2][2]
    };
    double norme_axe = sqrt(axe[0]*axe[0] + axe[1]*axe[1] + axe[2]*axe[2]);
    for (int i = 0; i < 3; i++) axe[i] /= norme_axe;

    // Vecteur normal (choisi orthogonal à axe)
    double ref[3] = {0, 0, 1};
    if (fabs(axe[2]) > 0.9) { ref[0] = 1; ref[2] = 0; } // éviter la colinéarité
    double normal[3] = {
        axe[1]*ref[2] - axe[2]*ref[1],
        axe[2]*ref[0] - axe[0]*ref[2],
        axe[0]*ref[1] - axe[1]*ref[0]
    };
    double norme_normal = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
    for (int i = 0; i < 3; i++) normal[i] /= norme_normal;

    // Tangente orthogonale
    double tangent[3] = {
        axe[1]*normal[2] - axe[2]*normal[1],
        axe[2]*normal[0] - axe[0]*normal[2],
        axe[0]*normal[1] - axe[1]*normal[0]
    };

    // Génération circulaire
    for (int p = 0; p < rnap_subunits; p++) {
        double angle = (p * 2.0 * M_PI) / rnap_subunits;
        double xi = cos(angle) * normal[0] * taille_pol/2 + sin(angle) * tangent[0] * taille_pol/2;
        double yi = cos(angle) * normal[1] * taille_pol/2 + sin(angle) * tangent[1] * taille_pol/2;
        double zi = cos(angle) * normal[2] * taille_pol/2 + sin(angle) * tangent[2] * taille_pol/2;
        R_rnap[id_rnap][p][0] = x0 + xi;
        R_rnap[id_rnap][p][1] = y0 + yi;
        R_rnap[id_rnap][p][2] = z0 + zi;
    }

    positions_bille_rnap[id_rnap] = promoteur;
}





void enregistrement_data(SimVars *sv, const Config *cfg, const Files *f, int t){

    if(t%cfg->periode_msd == 0){
            for(int i = 0; i < cfg->Nm; i++){
                for(int j = 0; j < 3; j++){
                    sv->R_monomere_arrays[(int)t/cfg->periode_msd][i][j] = sv->R[sv->list_monomere[i]][j];
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
        }

        if(t%cfg->periode_enregistrement == 0){

            if(sv->nb_rnap > 0){
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
                // enregistrement_RNAP_position(f->fichier_rnap, sv->nb_rnap, cfg->T_eq + t, sv->positions_bille_rnap, sv->avancement_transcription[0]);
            }
            else{
                enregistrement(f->fichier, sv->R, cfg->N, cfg->T_eq + t);
            }
        }
}





void calcul(SimVars *sv, const Config *cfg, const Files *f, NeighborList *neighbor_lists, NeighborList_rnap **neighbor_lists_rnap, int t_start)
{
    clock_t start_2, end_2;
    clock_t start, end;
    double duree_boucle, duree_tot = 0, temps_restant;

    struct timespec chrono_start, last, now;
    double interval = 60.0; // affichage de progression toutes les 60 s
    double checkpoint_limit_h = 47; // 47h30min
    double checkpoint_limit_s = checkpoint_limit_h * 3600.0;

    // --- Initialisation du timer global ---
    clock_gettime(CLOCK_MONOTONIC, &chrono_start);
    clock_gettime(CLOCK_MONOTONIC, &last);

    int time_depart = (t_start != 0) ? t_start : 0;
    printf("▶️  Simulation lancée à t=%d avec limite de %.1f h (%.0f s)\n",
           time_depart, checkpoint_limit_h, checkpoint_limit_s);

    
    
   

    for (int t = time_depart; t < cfg->T; t++) {

        
        clock_gettime(CLOCK_MONOTONIC, &now);
        double elapsed = (now.tv_sec - last.tv_sec) +
                         (now.tv_nsec - last.tv_nsec) * 1e-9;

        // Affichage de progression périodique
        if (elapsed >= interval) {
            double total_elapsed =
                (now.tv_sec - chrono_start.tv_sec) +
                (now.tv_nsec - chrono_start.tv_nsec) * 1e-9;
            printf("Itération %d (%.1f s depuis début)\n", t, total_elapsed);
            fflush(stdout);
            last = now;
        }

        // -----------------------------------------------------------------
        // ⏱️ Vérification de la durée totale → sauvegarde avant limite SLURM
        // -----------------------------------------------------------------
        double total_elapsed =
            (now.tv_sec - chrono_start.tv_sec) +
            (now.tv_nsec - chrono_start.tv_nsec) * 1e-9;

        // if (total_elapsed >= checkpoint_limit_s) {
        //     printf("\n💾 Limite de %.1f h atteinte (%.2f h écoulées)\n",
        //            checkpoint_limit_h, total_elapsed / 3600.0);
        //     printf("   → Sauvegarde automatique du checkpoint et arrêt propre.\n");

        //     char checkpoint_name[256];
        //     snprintf(checkpoint_name, sizeof(checkpoint_name),
        //              "checkpoint_t%d.dat", t);

        //     save_checkpoint(sv, cfg, t);

        //     printf("✅ Checkpoint enregistré sous %s\n", checkpoint_name);
        //     fflush(stdout);
        //     exit(0);
        // }
        

        // -----------------------------------------------------------------
        //  🧩 Code de simulation existant
        // -----------------------------------------------------------------
        
        // Ajout conditionnel des RNAP
        if (cfg->nb_rnap_initial > 0) {
            if (t - sv->last_rnap_add_time < cfg->rnap_refract_time) continue;
            int last_active = -1;
            for (int i = MAX_RNAP - 1; i >= 0; i--) {
                if (sv->l_rnap[i] == 1) {
                    last_active = i;
                    break;
                }
            }

            if (last_active == -1 && sv->nb_rnap < cfg->nb_rnap_initial) {
                ajouter_rnap(sv, cfg, neighbor_lists, &neighbor_lists_rnap, t);
            } else if (last_active >= 0) {
                int pos = sv->positions_bille_rnap[last_active];
                if (pos >= cfg->debut_segment + 3 &&
                    sv->nb_rnap < cfg->nb_rnap_initial) {
                    int has_available = 0;
                    for (int i = 0; i < cfg->nb_rnap_initial; i++) {
                        if (sv->l_rnap[i] == -1) {
                            has_available = 1;
                            break;
                        }
                    }
                    if (has_available) {
                        ajouter_rnap(sv, cfg, neighbor_lists,
                                     &neighbor_lists_rnap, t);
                        sv-> last_rnap_add_time = t;
                    }
                }
            }
        }

        start = clock();


        if (t % 1000 == 0) {
            build_neighbor_list(sv->R, neighbor_lists, cfg->N, 2, 0);
            if (cfg->nb_rnap_initial > 0) {
                build_neighbor_list_rnap_chrom(
                    cfg, sv->R_rnap, sv->nb_rnap, sv->R, cfg->N, neighbor_lists_rnap,
                    cfg->rayon_ecrantage_LJ_chrom, t);
            }
            
            for (int i = 0; i < cfg->N; i++) {
                fprintf(f->fichier_voisin, "%d ", neighbor_lists[i].count);
            }
            fprintf(f->fichier_voisin, "\n");
        }

        if (cfg->confinement == 1)
        {
             confinement_sphere(cfg, sv, t);
        }
        
        polymere_brownian_motion(
            sv->R, cfg->K, cfg->Delta, cfg->N, cfg->K_bend, sv->bending_forces,
            cfg->attache, cfg->plan, t, f->test, cfg->bending, sv->truc,
            cfg->T, f->fichier_force, cfg->periode_force,
            f->fichier_force_thermique, cfg->temperature);

        
        if (cfg->nb_rnap_initial > 0)
        {
            for (int rnap = 0; rnap < cfg->nb_rnap_initial; rnap++)
            {
                if (sv->l_rnap[rnap] < 0)
                    continue;
                int prout = 0;
                
                
                if(cfg->rnap_subunits == 8)
                {
                    polymere_brownian_motion_ring_force(
                        sv->R_rnap[rnap], cfg->alpha, cfg->K_rnap, cfg->Delta,
                        cfg->rnap_subunits, rnap, f->test, t, cfg->periode_force,
                        f->fichier_force_rnap, f->fichier_force_thermique,
                        cfg->temperature);

                    liaison_sup(3 * cfg->a_transpt, 2 * cfg->K_rnap, cfg->Delta,
                                sv->R_rnap[rnap], 0, 4, f->fichier_force_rnap_2, t,
                                cfg->periode_force);
                    liaison_sup(3 * cfg->a_transpt, 2 * cfg->K_rnap, cfg->Delta,
                                sv->R_rnap[rnap], 2, 6, f->fichier_force_rnap_2, t,
                                cfg->periode_force);
                    liaison_sup(3 * cfg->a_transpt, 2 * cfg->K_rnap, cfg->Delta,
                                sv->R_rnap[rnap], 1, 5, f->fichier_force_rnap_2, t,
                                cfg->periode_force);
                    liaison_sup(3 * cfg->a_transpt, 2 * cfg->K_rnap, cfg->Delta,
                                sv->R_rnap[rnap], 3, 7, f->fichier_force_rnap_2, t,
                                cfg->periode_force);
                }
               
                // if(cfg->quench == 0)
                // {
                //     sv->avancement_transcription[rnap] += cfg->dx_avancement_rnap;
                // }

                double com_rnap[3];
                int mono_proche;

                centre_masse_rnap(
                    sv->R_rnap[rnap],
                    cfg->rnap_subunits,
                    com_rnap
                );

                int p_old = sv->positions_bille_rnap[rnap];
                int W = 5;   // ex: 5–20

                mono_proche = monomere_plus_proche(
                    sv->R,
                    cfg->N,
                    com_rnap,
                    p_old,
                    W
                );

                sv->positions_bille_rnap[rnap] = mono_proche;

                force_rnap_along_chromatin(
                    cfg,
                    sv->R,
                    sv->R_rnap[rnap],
                    mono_proche,
                    3.52e-4,   // nouveau paramètre
                    cfg->Delta
                );

                bond_rnap_chromatine(
                    cfg, 
                    sv->R,
                    sv->R_rnap[rnap], 
                    cfg->a_transpt, 
                    cfg->K_transpt, 
                    cfg->Delta, 
                    sv->positions_bille_rnap[rnap], 
                    cfg->a, 
                    cfg->alpha
                );

                if (sv->positions_bille_rnap[rnap] == 400){
                    retirer_rnap(sv, cfg, neighbor_lists, &neighbor_lists_rnap,
                                 rnap, prout, t);

                }

                if (1 - sv->avancement_transcription[rnap] < 1e-7) {
                    sv->is_rnap[sv->positions_bille_rnap[rnap]] = 0;
                    retirer_rnap(sv, cfg, neighbor_lists, &neighbor_lists_rnap,
                                 rnap, prout, t);
                }

                if (sv->positions_bille_rnap[rnap] >= 0 && sv->positions_bille_rnap[rnap] < cfg->N)
                {
                    sv->is_rnap[sv->positions_bille_rnap[rnap]] = 1; 
                }
            }
        }

        lennard_jones_forces(sv->R, neighbor_lists, cfg->N, cfg->epsilon,
                             cfg->sigma6, cfg->sigma12, cfg->Delta,
                             cfg->attache, cfg->periode_force,
                             f->fichier_force_LJ, t);
        lennard_jones_forces_rnap(
            cfg, sv->R_rnap, sv->nb_rnap, sv->R, cfg->N, neighbor_lists_rnap,
            cfg->epsilon_rnap, cfg->sigma6_rnap, cfg->sigma12_rnap,
            cfg->sigma6_rnap2, cfg->sigma12_rnap2,
            cfg->rayon_ecrantage_LJ_rnap, cfg->Delta, t, f->test, cfg->T,
            f->fichier_force_rnap_LJ, cfg->periode_force, sv->l_rnap);

        lennard_jones_forces_rnap_rnap(
            cfg, sv->R_rnap, sv->nb_rnap, cfg->epsilon, cfg->sigma6_rnap,
            cfg->sigma12_rnap, cfg->rayon_ecrantage_LJ_rnap, cfg->Delta,
            sv->l_rnap);

        compteur_grands_deplacements(cfg->N, cfg->T, sv->R, sv->R_new,
                                     sv->compteur_grand_deplacement);

        enregistrement_data(sv, cfg, f, t);

        end = clock();
        duree_boucle = (double)(end - start) / CLOCKS_PER_SEC;
        duree_tot += duree_boucle;

        if (t % (cfg->T / 10) == 0) {
            Mesures mesures = calcul_mesures(sv->R, cfg->N);
            double duree_min = (int)(duree_tot / 60);
            double duree_sec = duree_tot - (duree_min * 60);
            temps_restant = duree_boucle * (cfg->T - t - 1) / 60;
            printf("%d/%d %.f:%.f %.2fmin std %.10f moy %.10f\n", t, cfg->T,
                   duree_min, duree_sec, temps_restant, mesures.std,
                   mesures.moyenne);
        }
    }
}

void f_equilibriate(SimVars *sv, const Config *cfg, const Files *f, NeighborList *neighbor_lists, NeighborList_rnap **neighbor_lists_rnap)
{
    printf("//// Avant la simulation //// \n");

    clock_t start_2, end_2; 
    clock_t start, end;double duree_boucle, duree_tot = 0, temps_restant;

    for (int t = 0; t < cfg->T_eq; t++){

        start = clock();

        if (t%1000==0){
            build_neighbor_list(sv->R, neighbor_lists, cfg->N, 2, 0);

            if(cfg->nb_rnap_initial > 0) {build_neighbor_list_rnap_chrom(cfg, sv->R_rnap, sv->nb_rnap, sv->R, cfg->N, neighbor_lists_rnap, cfg->rayon_ecrantage_LJ_chrom, t);}

            for(int i = 0; i < cfg->N; i++){
                fprintf(f->fichier_voisin,"%d ",neighbor_lists[i].count);
            }
            fprintf(f->fichier_voisin, "\n");
        }

        
        polymere_brownian_motion(
            sv->R, cfg->K, cfg->Delta, cfg->N, cfg->K_bend, sv->bending_forces,
            cfg->attache, cfg->plan, t, f->test, cfg->bending, sv->truc,
            cfg->T, f->fichier_force, cfg->periode_force,
            f->fichier_force_thermique, cfg->temperature);

        //update_link_vectors(R, t_link, N);
        //f_bending_forces(R, t_link, bending_forces, K_bend, N, t);


        /////////// Boucle sur les RNAPS /////////////


        
        lennard_jones_forces(sv->R, neighbor_lists, cfg->N, cfg->epsilon, cfg->sigma6, cfg->sigma12, cfg->Delta, cfg->attache, cfg->periode_force, f->fichier_force_LJ, t);
        lennard_jones_forces_rnap(
            cfg, sv->R_rnap, sv->nb_rnap, sv->R, cfg->N, neighbor_lists_rnap,
            cfg->epsilon_rnap, cfg->sigma6_rnap, cfg->sigma12_rnap,
            cfg->sigma6_rnap2, cfg->sigma12_rnap2,
            cfg->rayon_ecrantage_LJ_rnap, cfg->Delta, t, f->test, cfg->T,
            f->fichier_force_rnap_LJ, cfg->periode_force, sv->l_rnap);
        compteur_grands_deplacements(cfg->N, cfg->T, sv->R, sv->R_new, sv->compteur_grand_deplacement);


        end = clock();  
        duree_boucle = (double)(end - start)/CLOCKS_PER_SEC;
        duree_tot += duree_boucle ;

        if (t%(cfg->T_eq/10) == 0){
            Mesures mesures = calcul_mesures (sv->R, cfg->N);
            double duree_min = (int)(duree_tot / 60);
            double duree_sec = duree_tot - (duree_min * 60);
            temps_restant = duree_boucle * (cfg->T-t-1) / 60;
            printf("%d/%.d %.f:%.f %.2fmin std %.10f moy %.10f\n",t, cfg->T_eq, duree_min, duree_sec, temps_restant, mesures.std, mesures.moyenne);   
        }

        if(t%cfg->periode_enregistrement == 0){
            if(sv->nb_rnap > 0){
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
                // enregistrement_RNAP_position(f->fichier_rnap, sv->nb_rnap, t, sv->positions_bille_rnap, sv->avancement_transcription[0]);
            }
            else{
                enregistrement(f->fichier, sv->R, cfg->N, t);
            }
        }
    }
}



// void save_checkpoint(SimVars *sv, const Config *cfg, int t)
// {
//     unsigned long state[624];
//     int index;
//     get_mt_state(state, &index);

//     char current[256], backup[256];
//     snprintf(current, sizeof(current), "./Resultats/checkpoint_last.dat");
//     snprintf(backup, sizeof(backup), "./Resultats/checkpoint_prev.dat");

//     rename(current, backup);

//     FILE *f = fopen(current, "wb");
//     if (!f) {
//         perror("❌ fopen checkpoint");
//         return;
//     }

//     // --- META DATA ---
//     fwrite(&t,          sizeof(int), 1, f);
//     fwrite(&sv->nb_rnap, sizeof(int), 1, f);

//     // --- POLYMERE ---
//     for (int i = 0; i < cfg->N; i++)
//         fwrite(sv->R[i], sizeof(double), 3, f);

//     // --- RNAP ---

//     for (int i = 0; i < MAX_RNAP; i++)
//         for (int j = 0; j < cfg->rnap_subunits; j++)
//             for (int k = 0; k < 3; k++)
//                 fwrite(&sv->R_rnap[i][j][k], sizeof(double), 1, f);


//     // --- ARRAYS RNAP 1D ---
   
//     for (int i = 0; i < MAX_RNAP; i++)
//         fwrite(&sv->l_rnap[i], sizeof(int), 1, f);

//     for (int i = 0; i < MAX_RNAP; i++)
//         fwrite(&sv->positions_bille_rnap[i], sizeof(int), 1, f);
    
//     for (int i = 0; i < MAX_RNAP; i++)
//         fwrite(&sv->avancement_transcription[i], sizeof(double), 1, f);

//     for (int i = 0; i < MAX_RNAP; i++)
//         fwrite(&sv->compteur_mono_rnap[i], sizeof(int), 1, f);

//     // --- RNG ---
//     fwrite(state, sizeof(state), 1, f);
//     fwrite(&index, sizeof(int), 1, f);

//     fclose(f);

//     printf("💾 Checkpoint sauvegardé à t = %d\n", t);
// }


// #define CHECKPOINT_FILE "./Resultats/checkpoint_last.dat"

// int load_checkpoint_metadata(Config *cfg, int *t_start)
// {
//     FILE *f = fopen("./Resultats/checkpoint_last.dat", "rb");
//     if (!f) return 0;

//     // L'ordre doit matcher EXACTEMENT save_checkpoint()

//     fread(t_start, sizeof(int), 1, f);     // ✅ lis t
//     // fread(&cfg->nb_rnap_initial, sizeof(int), 1, f);  // ✅ lis nb_rnap (ne l'ignore pas)

//     fclose(f);
//     return 1;
// }


// int load_checkpoint_data(SimVars *sv, Config *cfg, int *t_start)
// {
//     FILE *f = fopen(CHECKPOINT_FILE, "rb");
//     if (!f)
//         return 0;

//     // Lire t_start et cfg sauvegardée
//     fread(t_start, sizeof(int), 1, f);
//     fread(&sv->nb_rnap, sizeof(int), 1, f);

//     // POLYMERE

//     for (int i = 0; i < cfg->N; i++)
//         fread(sv->R[i], sizeof(double), 3, f);
    

//     // RNAP matrices
//     for (int i = 0; i < MAX_RNAP; i++)
//         for (int j = 0; j < cfg->rnap_subunits; j++)
//             for (int k = 0; k < 3; k++)
//                 fread(&sv->R_rnap[i][j][k], sizeof(double), 1, f);

//     // RNAP tableaux 1D
//     for (int i = 0; i < MAX_RNAP; i++)
//         fread(&sv->l_rnap[i], sizeof(int), 1, f);

//     for (int i = 0; i < MAX_RNAP; i++)
//         fread(&sv->positions_bille_rnap[i], sizeof(int), 1, f);
    
//     for (int i = 0; i < MAX_RNAP; i++)
//         fread(&sv->avancement_transcription[i], sizeof(double), 1, f);

//     for (int i = 0; i < MAX_RNAP; i++)
//         fread(&sv->compteur_mono_rnap[i], sizeof(int), 1, f);

//     // RNG
//     unsigned long state[624];
//     int index;
//     fread(state, sizeof(state), 1, f);
//     fread(&index, sizeof(int), 1, f);
//     set_mt_state(state, index);

//     fclose(f);
//     return 1;
// }

void confinement_sphere(const Config *cfg, SimVars *sv, int t)
{
    // --- INITIALISATION DU CENTRE DE MASSE À t = 0 ---
    if (t == 0)
    {
        sv->cdm_conf[0] = 0.0;
        sv->cdm_conf[1] = 0.0;
        sv->cdm_conf[2] = 0.0;

        for (int i = 0; i < cfg->N; i++)
        {
            sv->cdm_conf[0] += sv->R[i][0];
            sv->cdm_conf[1] += sv->R[i][1];
            sv->cdm_conf[2] += sv->R[i][2];
        }
        sv->cdm_conf[0] /= cfg->N;
        sv->cdm_conf[1] /= cfg->N;
        sv->cdm_conf[2] /= cfg->N;
    }

    // --- PARAMÈTRES ---
    double R = cfg->r_conf;
    double sigma = cfg->sigma_conf;
    double eps = cfg->epsilon_conf;

    double sigma6 = pow(sigma, 6);
    double sigma12 = sigma6 * sigma6;

    // seuil de déplacement (important)
    double seuil = 20 * cfg->a;

    // compteur local
    int deplacement_trop_grand = 0;


    // --- BOUCLE SUR LES MONOMÈRES ---
    for (int i = 0; i < cfg->N; i++)
    {
        double dx = sv->R[i][0] - sv->cdm_conf[0];
        double dy = sv->R[i][1] - sv->cdm_conf[1];
        double dz = sv->R[i][2] - sv->cdm_conf[2];

        double d = sqrt(dx*dx + dy*dy + dz*dz);
        double wall_dist = R - d;

        if (wall_dist < sigma)  
        {
            double inv = sigma / (sigma - wall_dist);
            double inv6 = pow(inv, 6);
            double inv12 = inv6 * inv6;

            double delta = 0.02 * (2.0 * inv12 - inv6);

            // clamp pour éviter les explosions
            double delta_max = 1e-4 * cfg->a;    // 5% de la taille d'une bille
            if (delta > delta_max)
                delta = delta_max;

            // direction vers le centre
            double nx = -dx / d;
            double ny = -dy / d;
            double nz = -dz / d;

            // --- On calcule le déplacement ---
            double dx_move = delta * nx;
            double dy_move = delta * ny;
            double dz_move = delta * nz;

            double move_norm = sqrt(dx_move*dx_move + dy_move*dy_move + dz_move*dz_move);

            // --- Vérification du déplacement ---
            if (move_norm > seuil)
                deplacement_trop_grand++;

            // --- On applique le déplacement ---
            sv->R[i][0] += dx_move;
            sv->R[i][1] += dy_move;
            sv->R[i][2] += dz_move;

            // Sécurité : clamp sur la sphère
            double d2 = sqrt(
                (sv->R[i][0] - sv->cdm_conf[0]) * (sv->R[i][0] - sv->cdm_conf[0]) +
                (sv->R[i][1] - sv->cdm_conf[1]) * (sv->R[i][1] - sv->cdm_conf[1]) +
                (sv->R[i][2] - sv->cdm_conf[2]) * (sv->R[i][2] - sv->cdm_conf[2])
            );

            if (d2 > R)
            {
                double k = R / d2;
                sv->R[i][0] = sv->cdm_conf[0] + k * (sv->R[i][0] - sv->cdm_conf[0]);
                sv->R[i][1] = sv->cdm_conf[1] + k * (sv->R[i][1] - sv->cdm_conf[1]);
                sv->R[i][2] = sv->cdm_conf[2] + k * (sv->R[i][2] - sv->cdm_conf[2]);
            }
        }
    }

    // Boucle sur les RNAP

    for (int i = 0; i < MAX_RNAP; i++)
    {
        if(sv->l_rnap < 0)
            continue;
        
        for (int j = 0; j < cfg->rnap_subunits; j++ )
        {
            double dx = sv->R_rnap[i][j][0] - sv->cdm_conf[0];
            double dy = sv->R_rnap[i][j][1] - sv->cdm_conf[1];
            double dz = sv->R_rnap[i][j][2] - sv->cdm_conf[2];

            double d = sqrt(dx*dx + dy*dy + dz*dz);
            double wall_dist = R - d;

            if (wall_dist < sigma)
            {

                double inv = sigma / (sigma - wall_dist);
                double inv6 = pow(inv, 6);
                double inv12 = inv6 * inv6;

                double delta = 0.02 * (2.0 * inv12 - inv6);

                // clamp pour éviter les explosions
                double delta_max = 1e-4 * cfg->a;    // 5% de la taille d'une bille
                if (delta > delta_max)
                    delta = delta_max;

                // direction vers le centre
                double nx = -dx / d;
                double ny = -dy / d;
                double nz = -dz / d;
            
                // --- On calcule le déplacement ---
                double dx_move = delta * nx;
                double dy_move = delta * ny;
                double dz_move = delta * nz;

                double move_norm = sqrt(dx_move*dx_move + dy_move*dy_move + dz_move*dz_move);

                // --- Vérification du déplacement ---
                if (move_norm > seuil)
                    deplacement_trop_grand++;

                // --- On applique le déplacement ---
                sv->R_rnap[i][j][0] += dx_move;
                sv->R_rnap[i][j][1] += dy_move;
                sv->R_rnap[i][j][2] += dz_move;

                double d2 = sqrt(
                    (sv->R_rnap[i][j][0] - sv->cdm_conf[0]) * (sv->R_rnap[i][j][0] - sv->cdm_conf[0]) +
                    (sv->R_rnap[i][j][1] - sv->cdm_conf[1]) * (sv->R_rnap[i][j][1] - sv->cdm_conf[1]) +
                    (sv->R_rnap[i][j][2] - sv->cdm_conf[2]) * (sv->R_rnap[i][j][2] - sv->cdm_conf[2])
                );
                
                if (d2 > R)
                {
                    double k = R / d2;
                    sv->R_rnap[i][j][0] = sv->cdm_conf[0] + k * (sv->R_rnap[i][j][0] - sv->cdm_conf[0]);
                    sv->R_rnap[i][j][1] = sv->cdm_conf[1] + k * (sv->R_rnap[i][j][1] - sv->cdm_conf[1]);
                    sv->R_rnap[i][j][2] = sv->cdm_conf[2] + k * (sv->R_rnap[i][j][2] - sv->cdm_conf[2]);
                }
            }
        }
    }

    // --- On stocke le compteur dans sv ---
    sv->nb_move_large += deplacement_trop_grand;

    // // Tu peux afficher périodiquement
    // if (t % 1000 == 0 && t > 0)
    //     printf("[CONF] Grandes corrections cumulées = %d\n", sv->nb_move_large);
}
