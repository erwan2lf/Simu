#include "traj_binary.h"
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

/* ============================================================
 * Helpers internes
 * ============================================================ */

/* Écrit un int32 en little-endian (portable). */
static void write_i32(FILE *f, int32_t v)
{
    fwrite(&v, sizeof(int32_t), 1, f);
}

/* Écrit un float64 (double). */
static void write_f64(FILE *f, double v)
{
    uint8_t buf[8];
    memcpy(buf, &v, 8);
    fwrite(buf, 1, 8, f);
}

/* Écrit un float32.
 * Les positions n'ont pas besoin de double précision pour le MSD. */
static void write_f32(FILE *f, double v)
{
    float fv = (float)v;
    fwrite(&fv, sizeof(float), 1, f);
}


/* ============================================================
 * traj_open
 * ============================================================ */
FILE *traj_open(const char *path, const Config *cfg)
{
    FILE *f = fopen(path, "wb");
    if (!f) {
        fprintf(stderr, "❌ traj_open : impossible d'ouvrir '%s'\n", path);
        return NULL;
    }

    int N_segment = cfg->fin_segment - cfg->debut_segment;

    /* --- Header --- */
    write_i32(f, (int32_t)TRAJ_MAGIC);
    write_i32(f, (int32_t)TRAJ_VERSION);
    write_i32(f, (int32_t)cfg->N);
    write_i32(f, (int32_t)N_segment);
    write_i32(f, (int32_t)cfg->debut_segment);
    write_i32(f, (int32_t)cfg->fin_segment);
    write_i32(f, (int32_t)cfg->nb_rnap_initial);
    write_i32(f, (int32_t)cfg->rnap_subunits);
    write_f64(f, cfg->Delta);
    printf("\n");
    printf("DEBUG a = %f\n", cfg->a);
    printf("DEBUG Delta = %e\n", cfg->Delta);
    printf("\n");
    write_f64(f, cfg->a);
    write_i32(f, (int32_t)cfg->periode_enregistrement);

    fflush(f);

    printf("DEBUG traj_open : pos après header = %ld\n", ftell(f));

    printf("📂 traj_binary ouvert : %s\n", path);
    printf("   N_segment=%d  [%d → %d]  nb_rnap=%d  nsub=%d\n",
           N_segment,
           cfg->debut_segment, cfg->fin_segment,
           cfg->nb_rnap_initial,
           cfg->rnap_subunits);

    /* Estimation du poids final */
    long bytes_per_frame =
        sizeof(int32_t)                                          /* timestep  */
        + (long)N_segment * 3 * sizeof(float)                   /* chromatine */
        + (long)cfg->nb_rnap_initial * cfg->rnap_subunits * 3 * sizeof(float); /* RNAP */

    long n_frames = cfg->T_enregistrement;
    double size_mb = (double)(bytes_per_frame * n_frames) / (1024.0 * 1024.0);

    printf("   ~%.1f Mo estimés (%ld frames × %ld octets/frame)\n\n",
           size_mb, n_frames, bytes_per_frame);

    return f;
}


/* ============================================================
 * traj_write_frame
 * ============================================================ */
void traj_write_frame(FILE *f, int timestep,
                      double **R,
                      double ***R_rnap,
                      const int *l_rnap,
                      const Config *cfg)
{
    if (!f) return;

    /* --- Timestep --- */
    write_i32(f, (int32_t)timestep);

    /* --- Segment chromatine [debut_segment, fin_segment) --- */
    for (int i = cfg->debut_segment; i < cfg->fin_segment; i++) {
        write_f32(f, R[i][0]);
        write_f32(f, R[i][1]);
        write_f32(f, R[i][2]);
    }

    /* --- Sous-unités RNAP ---
     * On écrit toujours nb_rnap_initial blocs pour garder
     * une taille de frame constante (plus simple à lire ensuite).
     * Si la RNAP est inactive (l_rnap < 0), on écrit (0, 0, 0). */
    for (int rnap = 0; rnap < cfg->nb_rnap_initial; rnap++) {
        int active = (l_rnap != NULL && l_rnap[rnap] >= 1);

        for (int s = 0; s < cfg->rnap_subunits; s++) {
            if (active) {
                write_f32(f, R_rnap[rnap][s][0]);
                write_f32(f, R_rnap[rnap][s][1]);
                write_f32(f, R_rnap[rnap][s][2]);
            } else {
                write_f32(f, 0.0);
                write_f32(f, 0.0);
                write_f32(f, 0.0);
            }
        }
    }
}


/* ============================================================
 * traj_close
 * ============================================================ */
void traj_close(FILE *f)
{
    if (f) {
        fflush(f);
        fclose(f);
        printf("✅ traj_binary fermé.\n");
    }
}
