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
D_FREE = DELTA   # diffusion libre : D = kT/γ = Delta en unités réduites

# Temps de crossover en pas de temps : τ_c = 1 / (K * Delta)
TAU_C = 1.0 / (K * DELTA)   # = 1000 pas


def load_msd(path):
    data = np.loadtxt(path, comments='#')
    lag_idx  = data[:, 0].astype(int)
    lag_time = data[:, 1]   # en pas de temps (après modif msd.c)
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


def estimate_D_eff(t, msd, n_pts=5):
    """D_eff estimé depuis les premiers lags : D = MSD / (6τ*Delta)."""
    mask = (t > 0) & (msd > 0)
    t_e   = t[mask][:n_pts]
    msd_e = msd[mask][:n_pts]
    # t est en pas de temps, MSD = 6 D_eff Delta t
    D_vals = msd_e / (6.0 * D_FREE * t_e)
    return np.mean(D_vals)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "Resultats/msd.txt"

    lag_idx, lag_time, msd_mean, msd_std = load_msd(path)

    # --- Supprimer lag 0 ---
    mask = lag_idx > 0
    t    = lag_time[mask]   # en pas de temps
    msd  = msd_mean[mask]
    err  = msd_std[mask]

    # --- Tronquer à T/4 ---
    n_valid = len(t)
    t   = t[:n_valid]
    msd = msd[:n_valid]
    err = err[:n_valid]

    print(f"Fichier : {path}")
    print(f"Points affichés : {len(t)} (tronqué à T/4)")
    print(f"Lag min : {t[0]:.1f} pas  |  Lag max : {t[-1]:.1f} pas")
    print(f"τ_c théorique : {TAU_C:.0f} pas  (K={K}, Δ={DELTA})")
    print(f"D libre théorique : {D_FREE:.2e} a²/pas")

    # --- Fits ---
    t_early_max = t[len(t) // 10]
    alpha_early, A_early = fit_powerlaw(t, msd, t_max=t_early_max)

    t_late_min = t[len(t) // 10]
    t_late_max = t[len(t) * 7 // 10]
    alpha_late, A_late = fit_powerlaw(t, msd, t_min=t_late_min, t_max=t_late_max)

    D_eff = estimate_D_eff(t, msd)

    if alpha_early is not None:
        print(f"\nRégime précoce  (τ < {t_early_max:.0f} pas) : α = {alpha_early:.3f}")
    if alpha_late is not None:
        print(f"Régime tardif   ({t_late_min:.0f} < τ < {t_late_max:.0f} pas) : α = {alpha_late:.3f}")
    print(f"D_eff mesuré : {D_eff:.4f}  (ratio D_eff/D_free = {D_eff/D_FREE:.2f})")

    # =========================================================================
    # Figure
    # =========================================================================
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # ── Panneau gauche : log-log ──────────────────────────────────────────────
    ax = axes[0]

    ax.loglog(t, msd, color='steelblue', lw=1.5, zorder=3, label='MSD moyen')
    ax.fill_between(t, np.maximum(msd - err, 1e-15), msd + err,
                    alpha=0.2, color='steelblue', label='± std monomères')

    # Fits
    if alpha_early is not None:
        t_plot = np.array([t[1], t[len(t) // 10]])
        ax.loglog(t_plot, A_early * t_plot**alpha_early,
                  'r--', lw=1.5, label=f'α = {alpha_early:.2f} (précoce)')

    if alpha_late is not None:
        t_plot_late = t[(t >= t_late_min) & (t <= t_late_max)]
        ax.loglog(t_plot_late, A_late * t_plot_late**alpha_late,
                  'g--', lw=1.5, label=f'α = {alpha_late:.2f} (tardif)')

    # Référence diffusion libre : MSD = 6 D_free Delta t = 6 Delta^2 t
    t_free = np.logspace(np.log10(max(TAU_C * 0.01, t[0] * 0.1)),
                         np.log10(t[-1]), 300)
    ax.loglog(t_free, 6.0 * D_FREE * t_free,
              color='orange', lw=1.5, ls='-.',
              label=f'MSD = 6DΔτ  (D={D_FREE:.0e}, libre)')

    # Annotation τ_c
    msd_at_tc = 6.0 * D_FREE * TAU_C
    ax.axvline(TAU_C, color='orange', lw=0.8, ls=':', alpha=0.7)
    ax.annotate(f'τ_c ≈ {TAU_C:.0f} pas',
                xy=(TAU_C, msd_at_tc),
                xytext=(TAU_C * 5, msd_at_tc * 0.2),
                fontsize=8, color='darkorange',
                arrowprops=dict(arrowstyle='->', color='darkorange', lw=0.8))

    # Annotation premier lag
    ax.axvline(t[0], color='gray', lw=0.8, ls='--', alpha=0.6)
    ax.annotate(f'1er lag\n= {t[0]:.0f} pas',
                xy=(t[0], msd[0]),
                xytext=(t[0] * 3, msd[0] * 0.3),
                fontsize=8, color='gray',
                arrowprops=dict(arrowstyle='->', color='gray', lw=0.8))

    # Références théoriques
    ax.loglog(t, msd[0] * (t / t[0])**0.5, 'k:',  lw=0.8, alpha=0.4,
              label='∝ τ^0.5 (Rouse)')
    ax.loglog(t, msd[0] * (t / t[0])**1.0, 'k-.', lw=0.8, alpha=0.4,
              label='∝ τ^1.0 (diffusif)')

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

    # # Référence libre
    # t_lin = np.linspace(0, t[n_show - 1], 300)
    # ax2.plot(t_lin, 6.0 * D_FREE * t_lin,
    #          color='orange', lw=1.5, ls='-.',
    #          label=f'MSD = 6DΔτ (libre, D={D_FREE:.0e})')
    # ax2.plot(t_lin, 6.0 * D_eff * D_FREE * t_lin,
    #          color='red', lw=1.2, ls=':',
    #          label=f'MSD = 6D_eff·Δτ (D_eff={D_eff:.3f})')

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
