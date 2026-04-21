#ifndef TRAJ_BINARY_H
#define TRAJ_BINARY_H

#include <stdio.h>
#include "config.h"
#include "simulation.h"

/* ============================================================
 * Format du fichier binaire de trajectoire
 * ------------------------------------------------------------
 *
 * [ HEADER ]  (écrit une seule fois à l'ouverture)
 *   int32   magic         = 0x54524A42  ("TRJB") — signature
 *   int32   version       = 1
 *   int32   N             — nombre total de monomères
 *   int32   N_segment     — nb monomères du segment (fin - debut)
 *   int32   debut_segment
 *   int32   fin_segment
 *   int32   nb_rnap       — nb_rnap_initial
 *   int32   rnap_subunits
 *   float64 Delta         — pas de temps
 *   float64 a             — taille d'un monomère
 *   int32   periode       — période d'enregistrement (en pas)
 *
 * [ FRAME ]   (répété à chaque enregistrement)
 *   int32   timestep
 *   float32 R[debut..fin][3]        — positions segment chromatine
 *   float32 R_rnap[nb_rnap][nsub][3] — positions sous-unités RNAP
 *                                      (0.0 si RNAP inactive)
 * ============================================================ */

#define TRAJ_MAGIC   0x54524A42
#define TRAJ_VERSION 1

/* Ouvre le fichier et écrit le header.
 * Retourne le FILE* (NULL si erreur). */
FILE *traj_open(const char *path, const Config *cfg);

/* Écrit une frame (appelé à chaque période d'enregistrement). */
void traj_write_frame(FILE *f, int timestep,
                      double **R,
                      double ***R_rnap,
                      const int *l_rnap,
                      const Config *cfg);

/* Ferme proprement le fichier. */
void traj_close(FILE *f);

#endif /* TRAJ_BINARY_H */
