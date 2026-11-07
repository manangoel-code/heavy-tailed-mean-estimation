# When Do Monte Carlo Mean Estimates and Confidence Intervals Fail Under Heavy Tails?

Simulation code and results for the paper

> **When Do Monte Carlo Mean Estimates and Confidence Intervals Fail Under Heavy Tails?**
> *A controlled simulation of tail index, sign asymmetry, and median-of-means block choice*
> Manan Goel. [Read the paper (PDF)](paper/When_Do_Monte_Carlo_Mean_Estimates.pdf)

The study compares the ordinary sample mean, median of means (MoM) and the conventional Student-t interval on a
signed-Pareto family. The tail index α sets how heavy the tails are (α = 1.5 has infinite variance, α = 2.5 does
not), and the sign probability p sets how one-sided the data are. Every condition has an exact known mean, so errors
and interval coverage can be measured directly.

Main findings, from Table 1 of the paper (k = 9, n = 3,000):

- MoM has lower typical error than the mean at α = 1.5 when the data are symmetric (p = 0.5), but higher error when
  they are one-sided (p = 1), because the median of block means drops the rare large values that the true mean includes.
- With finite variance (α = 2.5) the ordinary mean does slightly better in every case.
- The Student-t interval covers the true mean about 95% of the time at α = 2.5 and at α = 1.5 with symmetric data.
  At α = 1.5 coverage falls to 91.0% when p = 0.75 and 79.4% when p = 1, and almost every miss lies below the true
  mean.

## Repository layout

```
paper/                       the paper (PDF)
simulate.py                  runs the full design and writes results/
compare_table1.py            compares results/table1.csv with Table 1 of the paper
make_figures.py              redraws Figures 1 and 2 from results/summary.csv
results/
    replicates.csv.gz        43,200 rows: one per sample (14,400) and block count k (3)
    summary.csv              36 conditions: errors, paired loss, win rate, bias, coverage, misses
    table1.csv               the k = 9, n = 3,000 rows (Table 1)
    table1_comparison.csv    this run next to the paper's Table 1
    figure1_paired_loss.png  median paired loss difference with bootstrap intervals
    figure2_coverage.png     interval coverage with binomial error bars
```

## Design

| Setting | Values |
|---|---|
| Tail index α | 1.5, 2.5 |
| Sign probability p | 0.5, 0.75, 1 |
| Sample size n | 300, 3,000 |
| MoM blocks k | 3, 9, 25 |
| Replications | 1,200 per α-p-n cell |

Y follows a Pareto type I distribution with minimum 1 and shape α, S is +1 with probability p and −1 otherwise, and
X = SY, so the true mean is μ = (2p − 1)α/(α − 1). Within a replication the magnitudes are shared across the three
values of p and only the signs are redrawn, so results for different p are paired. Each α-n-replication sample has its
own NumPy `SeedSequence`, so every run gives the same numbers.

For each sample the code records the mean, MoM for each k, the paired loss D = |MoM − μ| − |m − μ| and whether the
95% Student-t interval covers μ or misses it from below or above. The summary adds the median of D with a percentile
bootstrap interval (1,200 resamples), the MoM win rate with a Wilson interval, and bias.

## Run it

You need Python 3.12 and the pinned versions in `requirements.txt`.

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate    macOS/Linux: source .venv/bin/activate
python -m pip install -r requirements.txt

python simulate.py         # about 30 seconds; rewrites results/
python compare_table1.py   # compares the run with Table 1
python make_figures.py     # redraws the two figures
```

`simulate.py` checks the row counts (43,200 rows from 14,400 samples), that coverage and the two kinds of miss add up
to 1 in every condition, and the analytic Pareto mean against SciPy.

## Agreement with the paper

The random seeds of the run reported in the paper are not part of this repository, so `results/` is an independent
replication of the same design and its numbers differ from Table 1 by Monte Carlo error. `compare_table1.py` shows:

- the 95% intervals for the median paired loss overlap with the paper's in all 6 rows
- coverage differs by at most 1.8 percentage points; the largest standardized difference is z = −2.08 (α = 2.5,
  p = 0.5), with |z| ≤ 1.11 in the other five rows
- the direction of interval misses matches: at α = 1.5, p = 1, all misses fall below the true mean (21.8% here,
  20.6% in the paper)

## Scope

This is a computational working paper and has not been peer reviewed. The signed-Pareto family is a synthetic stress
test, not a model of financial returns, and the study proposes no new estimator.

## Citation

See [`CITATION.cff`](CITATION.cff).
