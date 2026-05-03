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
    long frame_bytes =
        sizeof(int32_t)
        + (long)h->N_segment * 3 * sizeof(float)
        + (long)h->nb_rnap * h->rnap_subunits * 3 * sizeof(float);

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
 * Lecture SÉQUENTIELLE de toutes les trajectoires
 *
 * Lit le fichier UNE SEULE FOIS de début à fin — pas de fseek.
 * Retourne traj[mono][frame][3] alloué par la fonction.
 *
 * RAM : nseg × n_frames × 3 × sizeof(float)
 *       Pour nseg=100, n_frames=100000 : ~115 Mo
 * ============================================================ */
static float ***read_all_trajs_sequential(FILE *f, const BinHeader *h,
                                          int n_frames)
{
    int  nseg      = h->N_segment;
    long rnap_bytes = (long)h->nb_rnap * h->rnap_subunits * 3 * sizeof(float);

    /* Allocation traj[mono][frame][3] */
    float ***traj = malloc((size_t)nseg * sizeof(float **));
    if (!traj) { perror("malloc traj"); return NULL; }

    for (int i = 0; i < nseg; i++) {
        traj[i] = malloc((size_t)n_frames * sizeof(float *));
        if (!traj[i]) { perror("malloc traj[i]"); goto error; }
        for (int fr = 0; fr < n_frames; fr++) {
            traj[i][fr] = malloc(3 * sizeof(float));
            if (!traj[i][fr]) { perror("malloc traj[i][fr]"); goto error; }
        }
    }

    /* Buffer pour absorber les données RNAP qu'on ignore */
    float *rnap_buf = (rnap_bytes > 0) ? malloc(rnap_bytes) : NULL;

    /* Lecture séquentielle — UNE frame à la fois, sans fseek */
    for (int fr = 0; fr < n_frames; fr++) {

        /* Timestep — ignoré */
        int32_t ts;
        if (fread(&ts, sizeof(int32_t), 1, f) != 1) {
            fprintf(stderr, "❌ fread timestep frame %d\n", fr);
            free(rnap_buf);
            goto error;
        }

        /* Positions chromatine — stockées */
        for (int i = 0; i < nseg; i++) {
            if (fread(traj[i][fr], sizeof(float), 3, f) != 3) {
                fprintf(stderr, "❌ fread chrom frame %d mono %d\n", fr, i);
                free(rnap_buf);
                goto error;
            }
        }

        /* Positions RNAP — ignorées */
        if (rnap_bytes > 0) {
            if (fread(rnap_buf, 1, rnap_bytes, f) != (size_t)rnap_bytes) {
                fprintf(stderr, "❌ fread rnap frame %d\n", fr);
                free(rnap_buf);
                goto error;
            }
        }

        if (fr % 10000 == 0)
            printf("   lecture frame %d/%d\n", fr, n_frames);
    }

    free(rnap_buf);
    return traj;

error:
    for (int i = 0; i < nseg; i++) {
        if (!traj[i]) break;
        for (int fr = 0; fr < n_frames; fr++)
            free(traj[i][fr]);
        free(traj[i]);
    }
    free(traj);
    return NULL;
}

/* ============================================================
 * Libération des trajectoires
 * ============================================================ */
static void free_all_trajs(float ***traj, int nseg, int n_frames)
{
    if (!traj) return;
    for (int i = 0; i < nseg; i++) {
        if (!traj[i]) continue;
        for (int fr = 0; fr < n_frames; fr++)
            free(traj[i][fr]);
        free(traj[i]);
    }
    free(traj);
}

/* ============================================================
 * Calcul du MSD d'une trajectoire 1D via FFT
 * (algorithme de Calandrini / Kneller, O(T log T))
 * ============================================================ */
static void msd_1d_fft(const double *x, int n, int n_lags,
                       double *msd_out, long *count)
{
    int nfft = 1;
    while (nfft < 2 * n) nfft <<= 1;

    fftw_complex *in  = fftw_alloc_complex(nfft);
    fftw_complex *out = fftw_alloc_complex(nfft);

    fftw_plan plan_fwd = fftw_plan_dft_1d(nfft, in, out, FFTW_FORWARD,  FFTW_ESTIMATE);
    fftw_plan plan_inv = fftw_plan_dft_1d(nfft, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    for (int i = 0;  i < n;    i++) { in[i][0] = x[i]; in[i][1] = 0.0; }
    for (int i = n; i < nfft; i++) { in[i][0] = 0.0;  in[i][1] = 0.0; }

    fftw_execute(plan_fwd);
    for (int i = 0; i < nfft; i++) {
        double re = out[i][0], im = out[i][1];
        out[i][0] = re*re + im*im;
        out[i][1] = 0.0;
    }
    fftw_execute(plan_inv);

    double *S2 = malloc((size_t)n * sizeof(double));
    for (int m = 0; m < n; m++)
        S2[m] = in[m][0] / (double)nfft;

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
 * msd_compute_from_file — point d'entrée public
 * ============================================================ */
int msd_compute_from_file(const char *bin_path,
                          const char *out_path,
                          int         n_lags)
{
    FILE *f = fopen(bin_path, "rb");
    if (!f) {
        fprintf(stderr, "❌ msd: impossible d'ouvrir '%s'\n", bin_path);
        return -1;
    }

    BinHeader h;
    if (read_header(f, &h) != 0) { fclose(f); return -1; }

    printf("📂 Header lu :\n");
    printf("   N=%d  segment=[%d,%d)  N_segment=%d\n",
           h.N, h.debut_segment, h.fin_segment, h.N_segment);
    printf("   nb_rnap=%d  nsub=%d  Delta=%.2e  periode=%d\n",
           h.nb_rnap, h.rnap_subunits, h.Delta, h.periode);

    int n_frames = count_frames(f, &h);
    if (n_frames <= 0) {
        fprintf(stderr, "❌ msd: aucune frame dans le fichier\n");
        fclose(f);
        return -1;
    }

    if (n_lags <= 0 || n_lags > n_frames)
        n_lags = n_frames;

    printf("   n_frames=%d  n_lags=%d\n", n_frames, n_lags);

    double dt_frame = (double)h.periode;
    int    nseg     = h.N_segment;

    double ram_mb = (double)nseg * n_frames * 3 * sizeof(float) / (1024.0 * 1024.0);
    printf("   RAM trajectoires : %.1f Mo\n\n", ram_mb);

    /* ── Lecture séquentielle — UN SEUL PASSAGE dans le fichier ─── */
    printf("📖 Lecture séquentielle...\n");
    float ***traj = read_all_trajs_sequential(f, &h, n_frames);
    fclose(f);

    if (!traj) {
        fprintf(stderr, "❌ msd: échec lecture\n");
        return -1;
    }
    printf("✅ Lecture terminée.\n\n");

    /* ── Allocation accumulateurs ────────────────────────────────── */
    double **msd_per_mono = malloc((size_t)nseg * sizeof(double *));
    long   **count_mono   = malloc((size_t)nseg * sizeof(long *));
    for (int i = 0; i < nseg; i++) {
        msd_per_mono[i] = calloc((size_t)n_lags, sizeof(double));
        count_mono[i]   = calloc((size_t)n_lags, sizeof(long));
        if (!msd_per_mono[i] || !count_mono[i]) {
            perror("calloc");
            free_all_trajs(traj, nseg, n_frames);
            return -1;
        }
    }

    /* ── Calcul MSD ──────────────────────────────────────────────── */
    printf("⏳ Calcul MSD...\n");
    double *x = malloc((size_t)n_frames * sizeof(double));
    if (!x) { perror("malloc x"); free_all_trajs(traj, nseg, n_frames); return -1; }

    for (int i = 0; i < nseg; i++) {
        if (i % 10 == 0)
            printf("   monomère %d/%d\n", h.debut_segment + i, h.fin_segment);

        for (int dim = 0; dim < 3; dim++) {
            for (int fr = 0; fr < n_frames; fr++)
                x[fr] = (double)traj[i][fr][dim];
            msd_1d_fft(x, n_frames, n_lags, msd_per_mono[i], count_mono[i]);
        }
    }
    free(x);
    free_all_trajs(traj, nseg, n_frames);

    /* ── Moyenne et écart-type ───────────────────────────────────── */
    double *msd_mean = calloc((size_t)n_lags, sizeof(double));
    double *msd_std  = calloc((size_t)n_lags, sizeof(double));

    for (int m = 0; m < n_lags; m++) {
        for (int i = 0; i < nseg; i++)
            if (count_mono[i][m] > 0)
                msd_per_mono[i][m] /= (double)count_mono[i][m];

        double sum = 0.0;
        for (int i = 0; i < nseg; i++) sum += msd_per_mono[i][m];
        msd_mean[m] = sum / (double)nseg;

        double sum2 = 0.0;
        for (int i = 0; i < nseg; i++) {
            double d = msd_per_mono[i][m] - msd_mean[m];
            sum2 += d * d;
        }
        msd_std[m] = sqrt(sum2 / (double)nseg);
    }

    /* ── Écriture ────────────────────────────────────────────────── */
    FILE *fout = fopen(out_path, "w");
    if (!fout) {
        fprintf(stderr, "❌ msd: impossible d'ouvrir '%s'\n", out_path);
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

    /* ── Libération ──────────────────────────────────────────────── */
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