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
from scipy.optimize import curve_fit

# ─────────────────────────────────────────────────────────────
# Paramètres physiques
# ─────────────────────────────────────────────────────────────
K      = 10.0
DELTA  = 1e-4
N_POLY = 1000
A      = 1.0
D_FREE = 1.0

TAU_0 = 1.0 / DELTA
TAU_R = N_POLY**2 * TAU_0 / (np.pi**2)


def load_msd(path):
    data = np.loadtxt(path, comments='#')
    return data[:, 0].astype(int), data[:, 1], data[:, 2], data[:, 3]


def fit_powerlaw(t, msd, t_min=None, t_max=None):
    mask = (msd > 0) & (t > 0)
    if t_min is not None: mask &= (t >= t_min)
    if t_max is not None: mask &= (t <= t_max)
    if mask.sum() < 3:
        return None, None
    slope, intercept, *_ = stats.linregress(
        np.log10(t[mask]), np.log10(msd[mask])
    )
    return slope, 10**intercept


def fit_powerlaw_v2(t, msd, t_min=None, t_max=None):
    """Fit MSD ~ A*t^alpha + v²*t²"""
    mask = (msd > 0) & (t > 0)
    if t_min is not None: mask &= (t >= t_min)
    if t_max is not None: mask &= (t <= t_max)
    if mask.sum() < 6:
        return None

    t_fit   = t[mask]
    msd_fit = msd[mask]

    def model(t, A, alpha, v2):
        return A * t**alpha + v2 * t**2

    try:
        p0 = [msd_fit[0] / t_fit[0]**0.75, 0.75, 0.0]
        popt, pcov = curve_fit(
            model, t_fit, msd_fit, p0=p0,
            bounds=([0, 0.1, 0], [np.inf, 1.0, np.inf]),
            maxfev=5000
        )
        perr = np.sqrt(np.diag(pcov))
        return {
            'A'        : popt[0],
            'alpha'    : popt[1],
            'v2'       : popt[2],
            'v'        : np.sqrt(max(popt[2], 0.0)),
            'A_err'    : perr[0],
            'alpha_err': perr[1],
            'v_err'    : perr[2] / (2*np.sqrt(popt[2])) if popt[2] > 0 else 0.0
        }
    except Exception as e:
        print(f"  ⚠ fit_powerlaw_v2 échoué : {e}")
        return None


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "Resultats/msd.txt"

    lag_idx, lag_time, msd_mean, msd_std = load_msd(path)

    # Supprimer lag 0
    mask = lag_idx > 0
    t    = lag_time[mask]
    msd  = msd_mean[mask]
    err  = msd_std[mask]

    # Tronquer à T/4 pour le log-log
    n_valid = len(t) // 4
    t   = t[:n_valid]
    msd = msd[:n_valid]
    err = err[:n_valid]

    # 20% du T total pour le fit balistique
    t_max_lin = t[-1] * 0.2 * 4   # 20% de T (pas de T/4)
    mask_lin  = t <= t_max_lin
    t_plot    = t[mask_lin]
    msd_plot  = msd[mask_lin]
    err_plot  = err[mask_lin]

    print(f"Fichier : {path}")
    print(f"τ₀  = {TAU_0:.1f} pas")
    print(f"τ_R = {TAU_R:.2e} pas")
    print(f"T/4 = {t[-1]:.1e} pas  |  T_20% = {t_max_lin:.1e} pas")

    # ── Fit log-log (régime intermédiaire) ───────────────────────────
    t_fit_min = t[len(t) // 10]
    t_fit_max = t[len(t) * 7 // 10]
    alpha_ll, A_ll = fit_powerlaw(t, msd, t_min=t_fit_min, t_max=t_fit_max)

    if alpha_ll is not None:
        print(f"\nFit log-log : α = {alpha_ll:.3f}")

    # ── Fit balistique sur les 20% ────────────────────────────────────
    fr = fit_powerlaw_v2(t_plot, msd_plot)

    if fr is not None:
        print(f"Fit balistique (20%) :")
        print(f"  A     = {fr['A']:.3e} ± {fr['A_err']:.1e}")
        print(f"  alpha = {fr['alpha']:.3f} ± {fr['alpha_err']:.3f}")
        print(f"  v     = {fr['v']:.3e} ± {fr['v_err']:.1e} a/pas")
        print(f"  v²    = {fr['v2']:.3e}")
        SNR = fr['v'] / fr['v_err'] if fr['v_err'] > 0 else np.inf
        print(f"  SNR   = {SNR:.1f}")

    # =========================================================================
    # Figure
    # =========================================================================
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(f'MSD — {path}', fontsize=12, fontweight='bold')

    # ── Panneau gauche : log-log ──────────────────────────────────────
    ax = axes[0]

    ax.loglog(t, msd, color='steelblue', lw=1.5, zorder=3, label='MSD moyen')
    ax.fill_between(t, np.maximum(msd - err, 1e-15), msd + err,
                    alpha=0.2, color='steelblue', label='± std')

    # Fit log-log
    if alpha_ll is not None:
        t_arr = t[(t >= t_fit_min) & (t <= t_fit_max)]
        ax.loglog(t_arr, A_ll * t_arr**alpha_ll,
                  'r--', lw=1.5, label=f'α = {alpha_ll:.2f}')

    # Références Rouse / diffusif
    ax.loglog(t, msd[0] * (t/t[0])**0.5, 'k:',  lw=0.8, alpha=0.4,
              label='∝ τ^0.5 (Rouse)')
    ax.loglog(t, msd[0] * (t/t[0])**1.0, 'k-.', lw=0.8, alpha=0.4,
              label='∝ τ^1.0 (diffusif)')

    # τ₀
    ax.axvline(TAU_0, color='green', lw=1.0, ls='--', alpha=0.7,
               label=f'τ₀ = {TAU_0:.0f} pas')

    # Limite T_20%
    ax.axvline(t_max_lin, color='orange', lw=1.0, ls='--', alpha=0.7,
               label=f'T_20% = {t_max_lin:.1e} pas')

    ax.set_xlabel('Lag τ [pas de temps]', fontsize=12)
    ax.set_ylabel('MSD [a²]', fontsize=12)
    ax.set_title('Log-log', fontsize=11)
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(True, which='both', alpha=0.3)

    # ── Panneau droit : régime précoce + fit balistique ───────────────
    ax2 = axes[1]

    ax2.plot(t_plot, msd_plot, color='steelblue', lw=1.5, label='MSD moyen')
    ax2.fill_between(t_plot, msd_plot - err_plot, msd_plot + err_plot,
                     alpha=0.2, color='steelblue', label='± std')

    # Fit balistique
    if fr is not None:
        t_dense = np.linspace(t_plot[0], t_plot[-1], 500)
        msd_fit_curve = fr['A'] * t_dense**fr['alpha'] + fr['v2'] * t_dense**2
        ax2.plot(t_dense, msd_fit_curve,
                 color='black', lw=2.0, alpha=0.9,
                 linestyle=(0, (3, 1)),
                 label=f'fit: α={fr["alpha"]:.2f}, v={fr["v"]:.2e}')

    # τ₀ si dans la fenêtre
    if TAU_0 <= t_plot[-1]:
        ax2.axvline(TAU_0, color='green', lw=1.0, ls='--', alpha=0.7,
                    label=f'τ₀ = {TAU_0:.0f} pas')

    ax2.set_xlabel('Lag τ [pas de temps]', fontsize=12)
    ax2.set_ylabel('MSD [a²]', fontsize=12)
    ax2.set_title(f'Régime précoce — 20% de T ({t_max_lin:.1e} pas)', fontsize=11)
    ax2.set_xlim(left=0)
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    out_fig = path.replace('.txt', '.png')
    plt.savefig(out_fig, dpi=150, bbox_inches='tight')
    print(f"\n📊 Figure sauvegardée : {out_fig}")
    plt.show()


if __name__ == '__main__':
    main()