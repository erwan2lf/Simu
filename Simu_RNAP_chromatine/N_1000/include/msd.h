#ifndef MSD_H
#define MSD_H

/* ============================================================
 * msd.h — Calcul du MSD par autocorrélation via FFT
 * ------------------------------------------------------------
 *
 * Workflow :
 *   1. msd_compute_from_file()
 *      Lit le fichier binaire produit par traj_binary.c,
 *      calcule le MSD par monomère (300→400) via FFT,
 *      moyenne sur les monomères, écrit le résultat.
 *
 * Format du fichier de sortie (texte, lisible en Python) :
 *   # lag_index   lag_time   msd_mean   msd_std
 *   0             0.0        0.0        0.0
 *   1             0.011      0.00432    0.00021
 *   ...
 *
 * Dépendances : FFTW3  (-lfftw3)
 * ============================================================ */

/* Point d'entrée principal.
 *
 *  bin_path   : chemin vers le fichier trajectoire.bin
 *  out_path   : chemin du fichier MSD de sortie (texte)
 *  n_lags     : nombre de lags voulus (≤ n_frames, typiquement n_frames)
 *
 * Retourne 0 si succès, -1 si erreur. */
int msd_compute_from_file(const char *bin_path,
                          const char *out_path,
                          int         n_lags);

#endif /* MSD_H */
