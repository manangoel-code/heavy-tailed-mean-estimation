"""Compare results/table1.csv from simulate.py with Table 1 as printed in the paper.

The paper's run used different random seeds, so the two are independent Monte Carlo replications of the same design.
For coverage the script reports z = (new - paper) / sqrt(SE_new^2 + SE_paper^2) with binomial standard errors;
for the median paired loss it checks whether the two 95% bootstrap intervals overlap.
Writes results/table1_comparison.csv.
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
REPS = 1200

# Table 1 of the paper: k = 9, n = 3,000.
PAPER = pd.DataFrame([
    # alpha, p, mean med., MoM med., median D, D low, D high, coverage, below, above
    (1.5, 0.50, 0.125, 0.102, -0.019, -0.026, -0.012, 0.957, 0.024, 0.019),
    (1.5, 0.75, 0.132, 0.121, -0.009, -0.012, -0.002, 0.910, 0.085, 0.005),
    (1.5, 1.00, 0.165, 0.205, +0.041, +0.034, +0.049, 0.794, 0.206, 0.000),
    (2.5, 0.50, 0.027, 0.032, +0.004, +0.002, +0.006, 0.960, 0.019, 0.021),
    (2.5, 0.75, 0.024, 0.029, +0.003, +0.001, +0.004, 0.964, 0.020, 0.016),
    (2.5, 1.00, 0.016, 0.020, +0.003, +0.003, +0.004, 0.944, 0.049, 0.007),
], columns=['alpha', 'p', 'mean_med', 'mom_med', 'D', 'D_low', 'D_high', 'coverage', 'below', 'above'])


def main():
    new = pd.read_csv(ROOT / 'results' / 'table1.csv')
    df = PAPER.merge(new.rename(columns={'coverage': 'coverage_new'}), on=['alpha', 'p'])
    se = np.sqrt(df.coverage * (1 - df.coverage) / REPS + df.coverage_new * (1 - df.coverage_new) / REPS)
    out = pd.DataFrame({
        'alpha': df.alpha, 'p': df.p,
        'mean_med_paper': df.mean_med, 'mean_med_new': df.mean_median_abs_error.round(3),
        'mom_med_paper': df.mom_med, 'mom_med_new': df.mom_median_abs_error.round(3),
        'D_paper': df.D, 'D_new': df.median_D.round(3),
        'D_intervals_overlap': (df.median_D_low <= df.D_high) & (df.D_low <= df.median_D_high),
        'coverage_paper': df.coverage, 'coverage_new': df.coverage_new.round(3),
        'coverage_z': ((df.coverage_new - df.coverage) / se).round(2),
        'below_above_paper': [f'{b:.1%}/{a:.1%}' for b, a in zip(df.below, df.above)],
        'below_above_new': [f'{b:.1%}/{a:.1%}' for b, a in zip(df.miss_below, df.miss_above)],
    })
    out.to_csv(ROOT / 'results' / 'table1_comparison.csv', index=False)
    print(out.to_string(index=False))
    print(f"\nD intervals overlap in {out.D_intervals_overlap.sum()}/6 rows; "
          f"largest |coverage z| = {out.coverage_z.abs().max():.2f} (|z| > 2 in {(out.coverage_z.abs() > 2).sum()}/6 rows)")


if __name__ == '__main__':
    main()
