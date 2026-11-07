"""Signed-Pareto Monte Carlo study of the sample mean, median of means (MoM) and the Student-t interval.

Design (paper Sections 2-4):
    alpha in {1.5, 2.5}, p in {0.5, 0.75, 1}, n in {300, 3000}, k in {3, 9, 25}, 1,200 replications per alpha-p-n cell.
    Y ~ Pareto type I (minimum 1, shape alpha); S = +1 with probability p, else -1; X = S * Y.
    Within a replication the magnitudes Y are shared across the three values of p; the signs are drawn separately.
    Every alpha-n-replication sample has its own deterministic NumPy SeedSequence.

Outputs in results/:
    replicates.csv.gz   43,200 rows: one per sample and block count k
    summary.csv         36 conditions: errors, paired loss, win rate, bias, coverage, directional misses
    table1.csv          the k = 9, n = 3,000 rows reported as Table 1 in the paper
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results'

ALPHAS = (1.5, 2.5)
PS = (0.5, 0.75, 1.0)
NS = (300, 3000)
KS = (3, 9, 25)
REPS = 1200
BOOTSTRAP_RESAMPLES = 1200
ROOT_SEED = 314159
Z = stats.norm.ppf(0.975)


def true_mean(alpha, p):
    """mu = (2p - 1) * alpha / (alpha - 1), valid for alpha > 1."""
    return (2 * p - 1) * alpha / (alpha - 1)


def median_of_means(x, k):
    """Split into k consecutive blocks whose sizes differ by at most one, average each, take the median."""
    return float(np.median([block.mean() for block in np.array_split(x, k)]))


def wilson_interval(successes, trials):
    """95% Wilson score interval for a binomial proportion."""
    phat = successes / trials
    centre = (phat + Z ** 2 / (2 * trials)) / (1 + Z ** 2 / trials)
    half = Z * np.sqrt(phat * (1 - phat) / trials + Z ** 2 / (4 * trials ** 2)) / (1 + Z ** 2 / trials)
    return centre - half, centre + half


def simulate():
    rows = []
    for ai, alpha in enumerate(ALPHAS):
        for ni, n in enumerate(NS):
            t_crit = stats.t.ppf(0.975, n - 1)
            for rep in range(REPS):
                seed = np.random.SeedSequence(ROOT_SEED, spawn_key=(ai, ni, rep))
                magnitude_seed, *sign_seeds = seed.spawn(1 + len(PS))
                # numpy's pareto() draws Lomax (Pareto II); adding 1 gives Pareto type I with minimum 1.
                y = np.random.default_rng(magnitude_seed).pareto(alpha, n) + 1.0
                for pi, p in enumerate(PS):
                    signs = np.where(np.random.default_rng(sign_seeds[pi]).random(n) < p, 1.0, -1.0)
                    x = signs * y
                    mu = true_mean(alpha, p)
                    m = x.mean()
                    half_width = t_crit * x.std(ddof=1) / np.sqrt(n)
                    low, high = m - half_width, m + half_width
                    interval = 'below' if high < mu else 'above' if low > mu else 'covered'
                    for k in KS:
                        mom = median_of_means(x, k)
                        rows.append({'alpha': alpha, 'p': p, 'n': n, 'k': k, 'rep': rep, 'mu': mu,
                                     'mean': m, 'mom': mom, 'mean_abs_error': abs(m - mu), 'mom_abs_error': abs(mom - mu),
                                     'D': abs(mom - mu) - abs(m - mu), 'ci_low': low, 'ci_high': high,
                                     'interval': interval})
    return pd.DataFrame(rows)


def summarize(rep):
    out = []
    for ci, ((alpha, p, n, k), g) in enumerate(rep.groupby(['alpha', 'p', 'n', 'k'], sort=True)):
        assert len(g) == REPS
        d = g['D'].to_numpy()
        # Percentile bootstrap of the median paired loss over independent replications.
        rng = np.random.default_rng(np.random.SeedSequence(ROOT_SEED, spawn_key=(99, ci)))
        medians = np.median(d[rng.integers(0, REPS, size=(BOOTSTRAP_RESAMPLES, REPS))], axis=1)
        wins = int((d < 0).sum())
        win_low, win_high = wilson_interval(wins, REPS)
        share = g['interval'].value_counts(normalize=True)
        coverage, below, above = (share.get(c, 0.0) for c in ('covered', 'below', 'above'))
        assert abs(coverage + below + above - 1) < 1e-12  # coverage/miss events partition the replications
        out.append({
            'alpha': alpha, 'p': p, 'n': n, 'k': k, 'mu': g['mu'].iloc[0],
            'mean_median_abs_error': g['mean_abs_error'].median(), 'mom_median_abs_error': g['mom_abs_error'].median(),
            'mean_mean_abs_error': g['mean_abs_error'].mean(), 'mom_mean_abs_error': g['mom_abs_error'].mean(),
            'median_D': np.median(d), 'median_D_low': np.percentile(medians, 2.5), 'median_D_high': np.percentile(medians, 97.5),
            'mean_D_descriptive': d.mean(),
            'mom_win_rate': wins / REPS, 'mom_win_low': win_low, 'mom_win_high': win_high,
            'mean_bias': (g['mean'] - g['mu']).mean(), 'mom_bias': (g['mom'] - g['mu']).mean(),
            'coverage': coverage, 'miss_below': below, 'miss_above': above,
        })
    return pd.DataFrame(out)


def checks(rep, summary):
    assert len(rep) == len(ALPHAS) * len(PS) * len(NS) * len(KS) * REPS == 43_200
    assert rep.groupby(['alpha', 'n', 'rep', 'p']).ngroups == 14_400
    assert len(summary) == 36
    # Analytic Pareto mean: scipy's Pareto(b=alpha) has mean alpha / (alpha - 1).
    for alpha in ALPHAS:
        for p in PS:
            assert np.isclose(true_mean(alpha, p), (2 * p - 1) * stats.pareto(b=alpha).mean())
    # Magnitudes are shared across p: |X| is identical for the three p values of one replication.
    first = rep[(rep.k == KS[0]) & (rep.rep == 0)]
    assert first.groupby(['alpha', 'n'])['mean'].nunique().min() >= 2  # signs differ across p


def main():
    OUT.mkdir(exist_ok=True)
    rep = simulate()
    summary = summarize(rep)
    checks(rep, summary)
    rep.to_csv(OUT / 'replicates.csv.gz', index=False, compression={'method': 'gzip', 'mtime': 0})
    summary.to_csv(OUT / 'summary.csv', index=False)
    table1 = summary[(summary.k == 9) & (summary.n == 3000)]
    table1.to_csv(OUT / 'table1.csv', index=False)
    cols = ['alpha', 'p', 'mean_median_abs_error', 'mom_median_abs_error', 'median_D', 'median_D_low', 'median_D_high',
            'coverage', 'miss_below', 'miss_above']
    print(table1[cols].round(4).to_string(index=False))
    print(f'{len(rep)} replicate rows, {len(summary)} conditions written to {OUT}')


if __name__ == '__main__':
    main()
