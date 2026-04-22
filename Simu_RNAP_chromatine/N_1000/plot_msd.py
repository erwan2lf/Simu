"""
plot_msd.py
-----------
Lit le fichier MSD produit par msd.c et trace MSD(τ) en log-log.
Le temps est en unités de pas de temps (timestep).
Usage : python3 plot_msd.py Resultats/msd.txt
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# ─────────────────────────────────────────────────────────────
# Paramètres physiques — à adapter si tu changes la config
# ─────────────────────────────────────────────────────────────
K      = 10.0    # raideur du ressort (cfg->K)
DELTA  = 1e-4    # pas de temps (cfg->Delta)
N_POLY = 1000    # nombre de monomères (cfg->N)
A      = 1.0     # taille d'un monomère (cfg->a)

# Diffusion libre : D = kT/γ = 1 en unités réduites (kT=1, γ=1)
D_FREE = 1.0

# ── Temps caractéristiques ────────────────────────────────────
# τ₀ : temps pour diffuser d'une distance a (taille d'un monomère)
#   MSD = 6 D Δ τ₀ = a²  →  τ₀ = a² / (6 D Δ)
TAU_0 = A**2 / (6.0 * D_FREE * DELTA)   # en pas de temps

# τ_R : temps de Rouse — relaxation de toute la chaîne
#   τ_R = N² τ₀ / π²
TAU_R = N_POLY**2 * TAU_0 / (np.pi**2)  # en pas de temps


def load_msd(path):
    data = np.loadtxt(path, comments='#')
    lag_idx  = data[:, 0].astype(int)
    lag_time = data[:, 1]
    msd_mean = data[:, 2]
    msd_std  = data[:, 3]
    return lag_idx, lag_time, msd_mean, msd_std


def fit_powerlaw(t, msd, t_min=None, t_max=None):
    """Ajustement log-log sur [t_min, t_max]."""
    mask = (msd > 0) & (t > 0)
    if t_min is not None: mask &= (t >= t_min)
    if t_max is not None: mask &= (t <= t_max)
    if mask.sum() < 3:
        return None, None
    slope, intercept, r, p, se = stats.linregress(
        np.log10(t[mask]), np.log10(msd[mask])
    )
    return slope, 10**intercept


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "Resultats/msd.txt"

    lag_idx, lag_time, msd_mean, msd_std = load_msd(path)

    # --- Supprimer lag 0 ---
    mask = lag_idx > 0
    t    = lag_time[mask]
    msd  = msd_mean[mask]
    err  = msd_std[mask]

    # --- Tronquer à T/4 ---
    n_valid = len(t) // 4
    t   = t[:n_valid]
    msd = msd[:n_valid]
    err = err[:n_valid]

    print(f"Fichier : {path}")
    print(f"Points affichés : {len(t)} (tronqué à T/4)")
    print(f"Lag min : {t[0]:.1f} pas  |  Lag max : {t[-1]:.1f} pas")
    print(f"\nTemps caractéristiques :")
    print(f"  τ₀  = a²/(6DΔ) = {TAU_0:.1f} pas  "
          f"(D={D_FREE:.1f}, a={A:.1f}, Δ={DELTA:.0e})")
    print(f"  τ_R = N²τ₀/π²  = {TAU_R:.2e} pas  (N={N_POLY})")
    print(f"  Fenêtre simulée : {t[-1]:.1f} pas  "
          f"→ {t[-1]/TAU_0:.1f} τ₀  /  {t[-1]/TAU_R:.2e} τ_R")

    # --- Fits ---
    t_early_max = t[len(t) // 10]
    alpha_early, A_early = fit_powerlaw(t, msd, t_max=t_early_max)

    t_late_min = t[len(t) // 10]
    t_late_max = t[len(t) * 7 // 10]
    alpha_late, A_late = fit_powerlaw(t, msd, t_min=t_late_min, t_max=t_late_max)

    if alpha_early is not None:
        print(f"\nRégime précoce  (τ < {t_early_max:.0f} pas) : α = {alpha_early:.3f}")
    if alpha_late is not None:
        print(f"Régime tardif   ({t_late_min:.0f} < τ < {t_late_max:.0f} pas) : α = {alpha_late:.3f}")

    # =========================================================================
    # Figure
    # =========================================================================
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # ── Panneau gauche : log-log ──────────────────────────────────────────────
    ax = axes[0]

    ax.loglog(t, msd, color='steelblue', lw=1.5, zorder=3, label='MSD moyen')
    ax.fill_between(t, np.maximum(msd - err, 1e-15), msd + err,
                    alpha=0.2, color='steelblue', label='± std monomères')

    # Fits loi de puissance
    if alpha_early is not None:
        t_plot = np.array([t[1], t[len(t) // 10]])
        ax.loglog(t_plot, A_early * t_plot**alpha_early,
                  'r--', lw=1.5, label=f'α = {alpha_early:.2f} (précoce)')
    if alpha_late is not None:
        t_plot_late = t[(t >= t_late_min) & (t <= t_late_max)]
        ax.loglog(t_plot_late, A_late * t_plot_late**alpha_late,
                  'g--', lw=1.5, label=f'α = {alpha_late:.2f} (tardif)')

    # Référence diffusion libre : MSD = 6 D Δ τ
    t_ref = np.logspace(np.log10(t[0]), np.log10(t[-1]), 300)
    ax.loglog(t_ref, 6.0 * D_FREE * DELTA * t_ref,
              color='orange', lw=1.5, ls='-.',
              label=f'MSD = 6DΔτ  (D={D_FREE:.0f}, libre)')

    # Références théoriques Rouse / diffusif
    ax.loglog(t, msd[0] * (t / t[0])**0.5, 'k:',  lw=0.8, alpha=0.4,
              label='∝ τ^0.5 (Rouse)')
    ax.loglog(t, msd[0] * (t / t[0])**1.0, 'k-.', lw=0.8, alpha=0.4,
              label='∝ τ^1.0 (diffusif)')

    # ── Annotation τ₀ ────────────────────────────────────────────────────────
    msd_at_tau0 = 6.0 * D_FREE * DELTA * TAU_0  # = a² par définition
    ax.axvline(TAU_0, color='green', lw=1.0, ls='--', alpha=0.7)
    ax.annotate(f'τ₀ ≈ {TAU_0:.0f} pas\n(MSD ~ a²)',
                xy=(TAU_0, msd_at_tau0),
                xytext=(TAU_0 * 0.15, msd_at_tau0 * 4),
                fontsize=8, color='darkgreen',
                arrowprops=dict(arrowstyle='->', color='darkgreen', lw=0.8))

    # ── Annotation τ_R ───────────────────────────────────────────────────────
    # τ_R est probablement hors de la fenêtre — on indique juste sa valeur
    ax.annotate(f'τ_R ≈ {TAU_R:.1e} pas\n(hors fenêtre)',
                xy=(t[-1], msd[-1]),
                xytext=(t[-1] * 0.1, msd[-1] * 0.15),
                fontsize=8, color='purple',
                arrowprops=dict(arrowstyle='->', color='purple', lw=0.8))

    # Annotation premier lag
    ax.axvline(t[0], color='gray', lw=0.8, ls='--', alpha=0.5)
    ax.annotate(f'1er lag\n= {t[0]:.0f} pas',
                xy=(t[0], msd[0]),
                xytext=(t[0] * 3, msd[0] * 0.3),
                fontsize=8, color='gray',
                arrowprops=dict(arrowstyle='->', color='gray', lw=0.8))

    ax.set_xlabel('Lag τ [pas de temps]', fontsize=12)
    ax.set_ylabel('MSD [a²]', fontsize=12)
    ax.set_title('MSD segment transcrit (log-log)', fontsize=12)
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(True, which='both', alpha=0.3)

    # ── Panneau droit : régime précoce linéaire ───────────────────────────────
    ax2 = axes[1]
    n_show = max(len(t) // 4, 2)

    ax2.plot(t[:n_show], msd[:n_show], color='steelblue', lw=1.5,
             label='MSD moyen')
    ax2.fill_between(t[:n_show],
                     msd[:n_show] - err[:n_show],
                     msd[:n_show] + err[:n_show],
                     alpha=0.2, color='steelblue', label='± std')

    # Annotation τ₀ sur le panneau linéaire si dans la fenêtre
    if TAU_0 <= t[n_show - 1]:
        ax2.axvline(TAU_0, color='green', lw=1.0, ls='--', alpha=0.7,
                    label=f'τ₀ = {TAU_0:.0f} pas')

    ax2.set_xlabel('Lag τ [pas de temps]', fontsize=12)
    ax2.set_ylabel('MSD [a²]', fontsize=12)
    ax2.set_title('MSD — régime précoce (linéaire)', fontsize=12)
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    out_fig = path.replace('.txt', '.png')
    plt.savefig(out_fig, dpi=150, bbox_inches='tight')
    print(f"\n📊 Figure sauvegardée : {out_fig}")
    plt.show()


if __name__ == '__main__':
    main()