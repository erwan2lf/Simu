"""
Lecteur et visualiseur pour le format binaire de trajectoire.

Calcule le vecteur bout à bout entre deux billes i et j :
    r(t) = R[j](t) - R[i](t)
et trace |r(t)|² ainsi que ses composantes x, y, z.
"""

import sys
import struct
import numpy as np
import matplotlib.pyplot as plt


# ──────────────────────────────────────────────────────────────
# Lecture
# ──────────────────────────────────────────────────────────────

def read_header(f):
    def ri32(): return struct.unpack("<i", f.read(4))[0]
    def rf64(): return struct.unpack("<d", f.read(8))[0]

    hdr = {}
    hdr["magic"]                  = ri32()
    hdr["version"]                = ri32()
    hdr["N"]                      = ri32()
    hdr["N_segment"]              = ri32()
    hdr["debut_segment"]          = ri32()
    hdr["fin_segment"]            = ri32()
    hdr["nb_rnap_initial"]        = ri32()
    hdr["rnap_subunits"]          = ri32()
    hdr["Delta"]                  = rf64()
    hdr["a"]                      = rf64()
    hdr["periode_enregistrement"] = ri32()
    return hdr


def load_trajectory(filepath: str):
    with open(filepath, "rb") as f:
        hdr = read_header(f)

        N_seg  = hdr["N_segment"]
        n_rnap = hdr["nb_rnap_initial"]
        n_sub  = hdr["rnap_subunits"]

        frame_bytes = 4 + N_seg * 3 * 4 + n_rnap * n_sub * 3 * 4
        raw = f.read()

    n_frames, remainder = divmod(len(raw), frame_bytes)
    if remainder != 0:
        print(f"⚠️  {remainder} octets résiduels ignorés (frame incomplète).")

    print(f"[load] Header : N={hdr['N']}, N_segment={N_seg}, "
          f"debut={hdr['debut_segment']}, fin={hdr['fin_segment']}, "
          f"nb_rnap={n_rnap}, rnap_subunits={n_sub}, "
          f"Delta={hdr['Delta']}, a={hdr['a']}")
    print(f"[load] {n_frames} frames trouvées dans '{filepath}'")

    timesteps = np.empty(n_frames, dtype=np.int32)
    chrom     = np.empty((n_frames, N_seg, 3), dtype=np.float32)
    rnap      = np.empty((n_frames, n_rnap, n_sub, 3), dtype=np.float32)

    offset = 0
    for k in range(n_frames):
        timesteps[k] = struct.unpack_from("<i", raw, offset)[0]
        offset += 4

        n_f = N_seg * 3
        chrom[k] = np.frombuffer(raw, dtype="<f4", count=n_f, offset=offset).reshape(N_seg, 3)
        offset += n_f * 4

        n_f_rnap = n_rnap * n_sub * 3
        rnap[k] = np.frombuffer(raw, dtype="<f4", count=n_f_rnap, offset=offset).reshape(n_rnap, n_sub, 3)
        offset += n_f_rnap * 4

    return hdr, timesteps, chrom, rnap


# ──────────────────────────────────────────────────────────────
# Visualisation
# ──────────────────────────────────────────────────────────────

def plot_end_to_end(filepath: str, bead_i: int = 300, bead_j: int = 400):
    """
    Trace le vecteur bout à bout r(t) = R[bead_j](t) - R[bead_i](t) :
      1. |r(t)|²
      2. Composantes rx(t), ry(t), rz(t)

    Les indices bead_i et bead_j sont des indices GLOBAUX dans la chaîne
    (le script soustrait debut_segment pour accéder au bon indice dans le tableau).
    """
    hdr, timesteps, chrom, rnap = load_trajectory(filepath)

    debut = hdr["debut_segment"]
    N_seg = hdr["N_segment"]

    # Conversion indices globaux → indices locaux dans chrom
    li = bead_i - debut
    lj = bead_j - debut

    if not (0 <= li < N_seg and 0 <= lj < N_seg):
        raise ValueError(
            f"Les billes {bead_i} et {bead_j} doivent être dans le segment "
            f"[{debut}, {debut + N_seg - 1}] enregistré."
        )

    # Vecteur bout à bout
    r  = chrom[:, lj, :] - chrom[:, li, :]   # (n_frames, 3)
    rx, ry, rz = r[:, 0], r[:, 1], r[:, 2]
    r2 = rx**2 + ry**2 + rz**2
    t  = timesteps.astype(float)

    fig, axes = plt.subplots(2, 1, figsize=(11, 8))
    label = f"billes {bead_i}→{bead_j}"

    # --- |r(t)|² ---
    axes[0].plot(t, r2, color="steelblue", linewidth=1.2)
    axes[0].set_xlabel("Timestep")
    axes[0].set_ylabel(r"$|\mathbf{r}_{" + str(bead_i) + r"\to" + str(bead_j) + r"}(t)|^2$")
    axes[0].set_title(f"Carré du vecteur bout à bout — {label}")
    axes[0].grid(True, alpha=0.3)

    # --- Composantes ---
    axes[1].plot(t, rx, label="x", color="tomato",       linewidth=1.0)
    axes[1].plot(t, ry, label="y", color="seagreen",     linewidth=1.0)
    axes[1].plot(t, rz, label="z", color="mediumpurple", linewidth=1.0)
    axes[1].set_xlabel("Timestep")
    axes[1].set_ylabel("Composante")
    axes[1].set_title(f"Composantes du vecteur bout à bout — {label}")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    out = filepath.replace(".bin", f"_r2_{bead_i}_{bead_j}.png")
    plt.savefig(out, dpi=150)
    plt.show()
    print(f"[plot] Graphe sauvegardé dans '{out}'")

    return hdr, timesteps, chrom, rnap, r2


# ──────────────────────────────────────────────────────────────
# Analyse statistique gaussienne
# ──────────────────────────────────────────────────────────────

def plot_gaussian_analysis(filepath: str, bead_i: int = None, bead_j: int = None):
    """
    Analyse statistique du vecteur bout à bout r = R[j] - R[i] :

      1. Histogrammes de rx, ry, rz avec ajustement gaussien
         → vérifie que chaque composante ~ N(0, σ²)
      2. Histogramme de |r|² avec ajustement chi² à 3 degrés de liberté
         → pour une chaîne gaussienne : |r|²/σ² ~ χ²(3)
      3. Q-Q plot de chaque composante
         → déviation visuelle par rapport à la droite gaussienne
      4. Impression des statistiques : <|r|²>, σ², skewness, kurtosis

    Paramètres
    ----------
    filepath : chemin vers le fichier binaire
    bead_i, bead_j : indices globaux (défaut : debut et fin du segment)
    """
    from scipy import stats

    hdr, timesteps, chrom, rnap = load_trajectory(filepath)

    debut = hdr["debut_segment"]
    N_seg = hdr["N_segment"]

    if bead_i is None: bead_i = debut
    if bead_j is None: bead_j = debut + N_seg - 1

    li = bead_i - debut
    lj = bead_j - debut

    if not (0 <= li < N_seg and 0 <= lj < N_seg):
        raise ValueError(
            f"Les billes {bead_i} et {bead_j} doivent être dans [{debut}, {debut + N_seg - 1}]."
        )

    r  = chrom[:, lj, :] - chrom[:, li, :]
    rx, ry, rz = r[:, 0], r[:, 1], r[:, 2]
    r2 = rx**2 + ry**2 + rz**2

    # ── Statistiques texte ──────────────────────────────────────
    print(f"\n{'─'*50}")
    print(f"  Analyse gaussienne : billes {bead_i} → {bead_j}")
    print(f"{'─'*50}")
    print(f"  <|r|²>      = {r2.mean():.4f}")
    print(f"  σ²  (via rx) = {rx.var():.4f}  (attendu : <|r|²>/3 = {r2.mean()/3:.4f})")
    for name, comp in [("rx", rx), ("ry", ry), ("rz", rz)]:
        sk = stats.skew(comp)
        ku = stats.kurtosis(comp)   # excès de kurtosis (0 = gaussien)
        _, p = stats.normaltest(comp)
        print(f"  {name} : skew={sk:+.3f}  kurt_exc={ku:+.3f}  p_normaltest={p:.3e}"
              + ("  ✅" if p > 0.05 else "  ❌"))
    print(f"{'─'*50}\n")

    # ── Figure ──────────────────────────────────────────────────
    fig = plt.figure(figsize=(14, 10))
    fig.suptitle(f"Analyse gaussienne — billes {bead_i}→{bead_j}", fontsize=13)

    components = [("rx", rx, "tomato"), ("ry", ry, "seagreen"), ("rz", rz, "mediumpurple")]

    # --- Ligne 1 : histogrammes des composantes avec fit gaussien ---
    for col, (name, comp, color) in enumerate(components):
        ax = fig.add_subplot(3, 3, col + 1)
        mu, sigma = comp.mean(), comp.std()
        ax.hist(comp, bins=60, density=True, color=color, alpha=0.6, label="données")
        xs = np.linspace(comp.min(), comp.max(), 300)
        ax.plot(xs, stats.norm.pdf(xs, mu, sigma), "k--", linewidth=1.5, label=f"N({mu:.2f}, {sigma:.2f}²)")
        ax.set_title(f"Histogramme {name}")
        ax.set_xlabel(name)
        ax.set_ylabel("Densité")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    # --- Ligne 2 : Q-Q plots ---
    for col, (name, comp, color) in enumerate(components):
        ax = fig.add_subplot(3, 3, col + 4)
        (osm, osr), (slope, intercept, _) = stats.probplot(comp, dist="norm")
        ax.scatter(osm, osr, s=1, color=color, alpha=0.4)
        xs = np.array([osm[0], osm[-1]])
        ax.plot(xs, slope * xs + intercept, "k--", linewidth=1.5)
        ax.set_title(f"Q-Q plot {name}")
        ax.set_xlabel("Quantiles théoriques")
        ax.set_ylabel("Quantiles observés")
        ax.grid(True, alpha=0.3)

    # --- Ligne 3 : histogramme de |r|² avec fit χ²(3)·σ² ---
    ax_r2 = fig.add_subplot(3, 3, 7)
    sigma2 = rx.var()   # estimé sur rx (isotropie supposée)
    r2_norm = r2 / sigma2  # devrait suivre χ²(3)
    ax_r2.hist(r2_norm, bins=60, density=True, color="steelblue", alpha=0.6, label=r"$|r|^2/\sigma^2$")
    xs = np.linspace(0, r2_norm.max(), 300)
    ax_r2.plot(xs, stats.chi2.pdf(xs, df=3), "k--", linewidth=1.5, label=r"$\chi^2(3)$")
    ax_r2.set_title(r"Histogramme $|r|^2/\sigma^2$")
    ax_r2.set_xlabel(r"$|r|^2 / \sigma^2$")
    ax_r2.set_ylabel("Densité")
    ax_r2.legend(fontsize=8)
    ax_r2.grid(True, alpha=0.3)

    # --- Évolution temporelle de |r|² ---
    ax_t = fig.add_subplot(3, 3, (8, 9))
    t = timesteps.astype(float)
    ax_t.plot(t, r2, color="steelblue", linewidth=0.6, alpha=0.7)
    ax_t.axhline(r2.mean(), color="k", linestyle="--", linewidth=1.2, label=f"<|r|²> = {r2.mean():.2f}")
    ax_t.set_xlabel("Timestep")
    ax_t.set_ylabel(r"$|r|^2$")
    ax_t.set_title(r"$|r(t)|^2$ au cours du temps")
    ax_t.legend(fontsize=8)
    ax_t.grid(True, alpha=0.3)

    plt.tight_layout()
    out = filepath.replace(".bin", f"_gauss_{bead_i}_{bead_j}.png")
    plt.savefig(out, dpi=150)
    plt.show()
    print(f"[plot] Analyse sauvegardée dans '{out}'")

    return r, r2


# ──────────────────────────────────────────────────────────────
# Point d'entrée
# ──────────────────────────────────────────────────────────────

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage : python3 plot_trajectoire.py <trajectoire.bin> [bead_i] [bead_j]")
        sys.exit(1)

    path = sys.argv[1]

    with open(path, "rb") as f:
        hdr_preview = read_header(f)
    debut = hdr_preview["debut_segment"]
    fin   = hdr_preview["fin_segment"] - 1

    bead_i = int(sys.argv[2]) if len(sys.argv) > 2 else debut
    bead_j = int(sys.argv[3]) if len(sys.argv) > 3 else fin

    print(f"[info] Vecteur bout à bout : bille {bead_i} → bille {bead_j}")

    # Trajectoire + composantes
    plot_end_to_end(path, bead_i=bead_i, bead_j=bead_j)

    # Analyse statistique gaussienne
    plot_gaussian_analysis(path, bead_i=bead_i, bead_j=bead_j)