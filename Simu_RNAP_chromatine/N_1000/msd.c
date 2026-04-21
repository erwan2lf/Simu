#include "msd.h"
#include "traj_binary.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <fftw3.h>

/* ============================================================
 * Structures internes
 * ============================================================ */

/* Header lu depuis le fichier .bin */
typedef struct {
    int32_t magic;
    int32_t version;
    int32_t N;
    int32_t N_segment;
    int32_t debut_segment;
    int32_t fin_segment;
    int32_t nb_rnap;
    int32_t rnap_subunits;
    double  Delta;
    double  a;
    int32_t periode;
} BinHeader;

/* ============================================================
 * Lecture du header
 * ============================================================ */
static int read_header(FILE *f, BinHeader *h)
{
    if (fread(&h->magic,         sizeof(int32_t), 1, f) != 1) return -1;
    if (fread(&h->version,       sizeof(int32_t), 1, f) != 1) return -1;
    if (fread(&h->N,             sizeof(int32_t), 1, f) != 1) return -1;
    if (fread(&h->N_segment,     sizeof(int32_t), 1, f) != 1) return -1;
    if (fread(&h->debut_segment, sizeof(int32_t), 1, f) != 1) return -1;
    if (fread(&h->fin_segment,   sizeof(int32_t), 1, f) != 1) return -1;
    if (fread(&h->nb_rnap,       sizeof(int32_t), 1, f) != 1) return -1;
    if (fread(&h->rnap_subunits, sizeof(int32_t), 1, f) != 1) return -1;
    if (fread(&h->Delta,         sizeof(double),  1, f) != 1) return -1;
    if (fread(&h->a,             sizeof(double),  1, f) != 1) return -1;
    if (fread(&h->periode,       sizeof(int32_t), 1, f) != 1) return -1;

    if (h->magic != (int32_t)TRAJ_MAGIC) {
        fprintf(stderr, "❌ msd: magic incorrect (0x%X attendu, 0x%X lu)\n",
                TRAJ_MAGIC, h->magic);
        return -1;
    }
    return 0;
}

/* ============================================================
 * Comptage des frames
 * ============================================================ */
static int count_frames(FILE *f, const BinHeader *h)
{
    /* Taille d'une frame en octets */
    long frame_bytes =
        sizeof(int32_t)                                              /* timestep */
        + (long)h->N_segment * 3 * sizeof(float)                    /* chromatine */
        + (long)h->nb_rnap * h->rnap_subunits * 3 * sizeof(float);  /* RNAP */

    long pos_start = ftell(f);
    fseek(f, 0, SEEK_END);
    long pos_end = ftell(f);
    fseek(f, pos_start, SEEK_SET);

    long data_bytes = pos_end - pos_start;
    int n_frames = (int)(data_bytes / frame_bytes);

    printf("📊 Fichier binaire : %d frames détectées (%.1f Mo)\n",
           n_frames, (double)(pos_end) / (1024.0 * 1024.0));

    return n_frames;
}

/* ============================================================
 * Lecture de toutes les trajectoires d'un monomère donné
 *
 * Retourne un tableau [n_frames][3] alloué par la fonction.
 * Le appelant doit le libérer avec free().
 * ============================================================ */
static float (*read_monomer_traj(FILE *f, const BinHeader *h,
                                 int n_frames, int local_idx,
                                 long pos_data_start))[3]
{
    long frame_bytes =
        sizeof(int32_t)
        + (long)h->N_segment * 3 * sizeof(float)
        + (long)h->nb_rnap * h->rnap_subunits * 3 * sizeof(float);

    long mono_offset =
        sizeof(int32_t)
        + (long)local_idx * 3 * sizeof(float);

    float (*traj)[3] = malloc((size_t)n_frames * sizeof(*traj));
    if (!traj) { perror("malloc traj"); return NULL; }

    for (int fr = 0; fr < n_frames; fr++) {
        long pos = pos_data_start
                 + (long)fr * frame_bytes
                 + mono_offset;

        if (fseek(f, pos, SEEK_SET) != 0) {
            fprintf(stderr, "❌ fseek échoué frame %d mono %d\n", fr, local_idx);
            free(traj);
            return NULL;
        }

        size_t n_read = fread(traj[fr], sizeof(float), 3, f);
        if (n_read != 3) {
            fprintf(stderr, "❌ fread frame %d mono %d : lu %zu/3 floats "
                            "(pos=%ld frame_bytes=%ld)\n",
                    fr, local_idx, n_read, pos, frame_bytes);
            free(traj);
            return NULL;
        }
    }

    return traj;
}

/* ============================================================
 * Calcul du MSD d'une trajectoire 1D via FFT
 * (algorithme de Calandrini / Kneller, O(T log T))
 *
 *  x       : signal de longueur n
 *  msd_out : tableau de longueur n_lags, déjà alloué
 *  count   : tableau de longueur n_lags (accumulation du nb de termes)
 *
 * MSD(m) = S1(m) - 2 * S2(m)
 *   S2(m) = Re[ IFFT( |FFT(x)|² ) ] / n   (autocorrélation via FFT)
 *   S1(m) = somme cumulative de x²
 * ============================================================ */
static void msd_1d_fft(const double *x, int n, int n_lags,
                       double *msd_out, long *count)
{
    int nfft = 1;
    while (nfft < 2 * n) nfft <<= 1;   /* zero-padding pour autocorrélation */

    fftw_complex *in  = fftw_alloc_complex(nfft);
    fftw_complex *out = fftw_alloc_complex(nfft);

    fftw_plan plan_fwd = fftw_plan_dft_1d(nfft, in, out, FFTW_FORWARD,  FFTW_ESTIMATE);
    fftw_plan plan_inv = fftw_plan_dft_1d(nfft, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    /* Remplissage avec zero-padding */
    for (int i = 0;  i < n;    i++) { in[i][0] = x[i]; in[i][1] = 0.0; }
    for (int i = n; i < nfft; i++) { in[i][0] = 0.0;  in[i][1] = 0.0; }

    /* FFT → |FFT|² → IFFT = autocorrélation */
    fftw_execute(plan_fwd);
    for (int i = 0; i < nfft; i++) {
        double re = out[i][0], im = out[i][1];
        out[i][0] = re*re + im*im;
        out[i][1] = 0.0;
    }
    fftw_execute(plan_inv);

    /* S2(m) = autocorrélation normalisée */
    double *S2 = malloc((size_t)n * sizeof(double));
    for (int m = 0; m < n; m++)
        S2[m] = in[m][0] / (double)nfft;

    /* S1 : somme cumulative de x[t]² + x[t+m]²
     * Calcul récursif : S1(m) = S1(m-1) - x[m-1]² - x[n-m]² */
    double S1 = 0.0;
    for (int i = 0; i < n; i++) S1 += 2.0 * x[i] * x[i];

    for (int m = 0; m < n_lags; m++) {
        if (m > 0) S1 -= x[m-1]*x[m-1] + x[n-m]*x[n-m];

        int denom = n - m;
        if (denom <= 0) break;

        double msd_m = (S1 - 2.0 * S2[m]) / (double)denom;
        msd_out[m] += msd_m;
        count[m]   += 1;
    }

    free(S2);
    fftw_destroy_plan(plan_fwd);
    fftw_destroy_plan(plan_inv);
    fftw_free(in);
    fftw_free(out);
}

/* ============================================================
 * msd_compute_from_file  — point d'entrée public
 * ============================================================ */
int msd_compute_from_file(const char *bin_path,
                          const char *out_path,
                          int         n_lags)
{
    /* --- Ouverture du fichier binaire --- */
    FILE *f = fopen(bin_path, "rb");
    if (!f) {
        fprintf(stderr, "❌ msd: impossible d'ouvrir '%s'\n", bin_path);
        return -1;
    }

    /* --- Lecture du header --- */
    BinHeader h;
    if (read_header(f, &h) != 0) { fclose(f); return -1; }

    printf("📂 Header lu :\n");
    printf("   N=%d  segment=[%d,%d)  N_segment=%d\n",
           h.N, h.debut_segment, h.fin_segment, h.N_segment);
    printf("   nb_rnap=%d  nsub=%d  Delta=%.2e  periode=%d\n",
           h.nb_rnap, h.rnap_subunits, h.Delta, h.periode);

    /* --- Comptage des frames --- */
    int n_frames = count_frames(f, &h);
    if (n_frames <= 0) {
        fprintf(stderr, "❌ msd: aucune frame dans le fichier\n");
        fclose(f);
        return -1;
    }

    /* Clamp n_lags */
    if (n_lags <= 0 || n_lags > n_frames)
        n_lags = n_frames;

    printf("   n_frames=%d  n_lags=%d\n\n", n_frames, n_lags);

    /* Temps réel entre deux frames */
    // double dt_frame = h.Delta * (double)h.periode; // en τ₀
    double dt_frame = (double)h.periode; // en timestep

    /* --- Allocation des accumulateurs ---
     * msd_per_mono[i][m] : MSD du monomère i au lag m
     * msd_mean[m]        : moyenne sur les monomères
     * msd_std[m]         : écart-type sur les monomères */
    int  nseg = h.N_segment;

    double **msd_per_mono = malloc((size_t)nseg * sizeof(double *));
    long   **count_mono   = malloc((size_t)nseg * sizeof(long *));
    for (int i = 0; i < nseg; i++) {
        msd_per_mono[i] = calloc((size_t)n_lags, sizeof(double));
        count_mono[i]   = calloc((size_t)n_lags, sizeof(long));
        if (!msd_per_mono[i] || !count_mono[i]) {
            perror("calloc msd_per_mono");
            fclose(f);
            return -1;
        }
    }

    long pos_data_start = ftell(f);
    printf("   pos_data_start = %ld\n", pos_data_start); // debug


    /* --- Boucle sur les monomères du segment --- */
    printf("⏳ Calcul MSD en cours...\n");

    double *x = malloc((size_t)n_frames * sizeof(double));
    if (!x) { perror("malloc x"); fclose(f); return -1; }

    for (int i = 0; i < nseg; i++) {

        if (i % 10 == 0)
            printf("   monomère %d/%d\n", h.debut_segment + i, h.fin_segment);

        /* Lecture de la trajectoire du monomère i */
        float (*traj)[3] = read_monomer_traj(f, &h, n_frames, i, pos_data_start);
        if (!traj) { free(x); fclose(f); return -1; }

        /* MSD sur chaque dimension puis somme */
        for (int dim = 0; dim < 3; dim++) {
            for (int fr = 0; fr < n_frames; fr++)
                x[fr] = (double)traj[fr][dim];

            msd_1d_fft(x, n_frames, n_lags,
                       msd_per_mono[i], count_mono[i]);
        }

        free(traj);
    }
    free(x);
    fclose(f);

    /* --- Calcul moyenne et écart-type sur les monomères --- */
    double *msd_mean = calloc((size_t)n_lags, sizeof(double));
    double *msd_std  = calloc((size_t)n_lags, sizeof(double));

    for (int m = 0; m < n_lags; m++) {
        /* Normalise chaque monomère (3 dimensions déjà sommées) */
        for (int i = 0; i < nseg; i++) {
            if (count_mono[i][m] > 0)
                msd_per_mono[i][m] /= (double)count_mono[i][m];
        }

        /* Moyenne */
        double sum = 0.0;
        for (int i = 0; i < nseg; i++) sum += msd_per_mono[i][m];
        msd_mean[m] = sum / (double)nseg;

        /* Écart-type */
        double sum2 = 0.0;
        for (int i = 0; i < nseg; i++) {
            double d = msd_per_mono[i][m] - msd_mean[m];
            sum2 += d * d;
        }
        msd_std[m] = sqrt(sum2 / (double)nseg);
    }

    /* --- Écriture du fichier de sortie --- */
    FILE *fout = fopen(out_path, "w");
    if (!fout) {
        fprintf(stderr, "❌ msd: impossible d'ouvrir '%s' en écriture\n", out_path);
        return -1;
    }

    fprintf(fout, "# MSD — segment [%d, %d)\n", h.debut_segment, h.fin_segment);
    fprintf(fout, "# n_monomeres=%d  n_frames=%d  Delta=%.6e  periode=%d\n",
            nseg, n_frames, h.Delta, h.periode);
    fprintf(fout, "# colonnes: lag_index  lag_time  msd_mean  msd_std\n");
    fprintf(fout, "# [adim]    [timestep]    [a^2]     [a^2]\n");

    for (int m = 0; m < n_lags; m++) {
        double lag_time = (double)m * dt_frame;
        fprintf(fout, "%d\t%.6e\t%.6e\t%.6e\n",
                m, lag_time, msd_mean[m], msd_std[m]);
    }

    fclose(fout);
    printf("✅ MSD écrit dans '%s'\n", out_path);

    /* --- Libération mémoire --- */
    for (int i = 0; i < nseg; i++) {
        free(msd_per_mono[i]);
        free(count_mono[i]);
    }
    free(msd_per_mono);
    free(count_mono);
    free(msd_mean);
    free(msd_std);

    return 0;
}
