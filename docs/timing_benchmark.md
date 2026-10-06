# codeml run times and the automatic time limit

These measurements set the default time limit for each codeml run
(`src/backend/timeouts.py`), in the window and on the command line.

## How it was measured

codeml 4.9j fitted one model per run, as EasyPAML does, with the default settings
(CodonFreq = 2, ncatG = 10, cleandata = 1, initial ω 0.5). The data were simulated
with `evolverNSsites` from PAML 4.10.10: random trees of 10, 30 and 60 taxa, 150,
500 and 1500 codons, an M3 model with three classes (proportions 0.6, 0.3, 0.1;
ω = 0.05, 1, 4) and κ = 2. Branch and Branch-site runs labelled the first pair of
sister tips as `#1`. The machine was an Intel Core i7-13620H (16 threads, 16 GB,
Linux) running 8 jobs at once under `nice 5`, so a single run is somewhat faster.
Of 81 runs (9 sizes × 9 models), 79 finished; M8 and M8a on 60 × 1500 were stopped
at the 8.3-hour cap.

Measured codeml time (minutes):

| Taxa × codons | M0 | M1a | M2a | M7 | M8 | M8a | Branch | Branch-site | BS null |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 × 150 | 0.1 | 0.2 | 0.2 | 1.6 | 0.7 | 1.1 | 0.1 | 0.3 | 0.2 |
| 10 × 500 | 0.1 | 0.2 | 0.3 | 1.4 | 1.4 | 2.6 | 0.1 | 0.8 | 1.0 |
| 10 × 1500 | 0.2 | 0.6 | 1.0 | 3.9 | 3.2 | 4.8 | 0.3 | 1.9 | 1.2 |
| 30 × 150 | 0.8 | 2.5 | 3.3 | 11 | 14 | 18 | 1.3 | 5.3 | 5.5 |
| 30 × 500 | 2.1 | 6.5 | 8.8 | 24 | 33 | 33 | 2.2 | 13 | 10 |
| 30 × 1500 | 5.4 | 12 | 22 | 67 | 74 | 82 | 6.6 | 35 | 40 |
| 60 × 150 | 7.8 | 17 | 27 | 57 | 99 | 73 | 6.7 | 61 | 36 |
| 60 × 500 | 20 | 33 | 67 | 203 | 237 | 266 | 16 | 74 | 88 |
| 60 × 1500 | 43 | 84 | 189 | 469 | > 500¹ | > 500¹ | 40 | 240 | 230 |

¹ Stopped at the 30,000 s cap of the benchmark.

Time grows steeply with the number of taxa and slowly with length: doubling the
taxa multiplies it by about 6.6, tripling the codons by about 2.2. M7, M8 and M8a
take 3.5 to 7 times as long as M1a and M2a.

## Formula

A least-squares fit on the log scale, with exponents shared by all models and one
constant per model:

    time ≈ T_ref[model] × (taxa / 30)^2.73 × (codons / 500)^0.71

| Model | M0 | M1a | M2a | M7 | M8 | M8a | Branch | Branch-site | BS null |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| T_ref (s, 30 taxa × 500 codons) | 140 | 340 | 530 | 1890 | 1950 | 2400 | 160 | 890 | 810 |

The time limit is 5 times the estimate, and at least 30 minutes. Models not in the
table use the Branch-site constant if their name starts with "Branch", otherwise
the largest constant.

The largest measured time in the benchmark was 2.4 times the estimate (M7, 10 × 150).
In 125 runs on 25 real *Cereus* genes (M1a, M2a, M7, M8, M8a; 11 to 24 taxa, 268 to
1997 codons) the largest was also 2.4 times the estimate, which is 0.47 of the limit;
the median was 0.38 times the estimate.

Resulting automatic time limit:

| Taxa × codons | M0 | M1a | M2a | M7 | M8 | M8a | Branch | Branch-site | BS null |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 × 150 | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min |
| 10 × 500 | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min |
| 10 × 1500 | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min |
| 30 × 150 | 30 min | 30 min | 30 min | 1.1 h | 1.2 h | 1.4 h | 30 min | 32 min | 30 min |
| 30 × 500 | 30 min | 30 min | 44 min | 2.6 h | 2.7 h | 3.3 h | 30 min | 1.2 h | 1.1 h |
| 30 × 1500 | 30 min | 1.0 h | 1.6 h | 5.7 h | 5.9 h | 7.3 h | 30 min | 2.7 h | 2.5 h |
| 60 × 150 | 33 min | 1.3 h | 2.1 h | 7.4 h | 7.6 h | 9.4 h | 38 min | 3.5 h | 3.2 h |
| 60 × 500 | 1.3 h | 3.1 h | 4.9 h | 17.4 h | 18.0 h | 22.1 h | 1.5 h | 8.2 h | 7.5 h |
| 60 × 1500 | 2.8 h | 6.8 h | 10.7 h | 38.0 h | 39.2 h | 48.2 h | 3.2 h | 17.9 h | 16.3 h |

## Changing the limit

In the window: Advanced settings › "Time limit per model (min)"; empty means
automatic. On the command line: `--timeout SECONDS` (0, the default, means
automatic), or the `"timeout"` key with `--config`. The limit used for each gene
and model is written to `batch_analysis_log.txt`.

A stuck codeml is caught earlier by the idle check: no CPU use for 5 minutes
(`--idle-timeout`, 0 turns it off).
