"""Redraw Figures 1 and 2 of the paper from results/summary.csv (n = 3,000).

Figure 1: median paired absolute-loss difference with 95% percentile-bootstrap intervals, k = 3, 9, 25.
Figure 2: empirical coverage of the Student-t interval with binomial Monte Carlo error bars, k = 9.
"""
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402  (backend must be set first)

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / 'results'
REPS = 1200


def main():
    s = pd.read_csv(RESULTS / 'summary.csv')
    s = s[s.n == 3000]

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    for ax, alpha in zip(axes, (1.5, 2.5)):
        for k in (3, 9, 25):
            g = s[(s.alpha == alpha) & (s.k == k)]
            ax.errorbar(g.p, g.median_D, yerr=[g.median_D - g.median_D_low, g.median_D_high - g.median_D],
                        marker='o', capsize=3, label=f'k={k}')
        ax.axhline(0, color='black', lw=0.8, ls=':')
        ax.set(title=f'alpha={alpha}', xlabel='Positive sign probability', ylabel='Median paired absolute-loss difference')
        ax.grid(alpha=0.2)
        ax.legend()
    fig.tight_layout()
    fig.savefig(RESULTS / 'figure1_paired_loss.png', dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    for ax, alpha in zip(axes, (1.5, 2.5)):
        g = s[(s.alpha == alpha) & (s.k == 9)]
        err = 1.96 * np.sqrt(g.coverage * (1 - g.coverage) / REPS)
        ax.errorbar(g.p, g.coverage, yerr=err, marker='o', capsize=3, label='k=9')
        ax.axhline(0.95, color='black', lw=0.8, ls=':')
        ax.set(title=f'alpha={alpha}', xlabel='Positive sign probability', ylabel='Studentized 95% interval coverage')
        ax.grid(alpha=0.2)
        ax.legend()
    fig.tight_layout()
    fig.savefig(RESULTS / 'figure2_coverage.png', dpi=200)
    plt.close(fig)
    print('Saved figure1_paired_loss.png and figure2_coverage.png in', RESULTS)


if __name__ == '__main__':
    main()
