# Phase 1 baselines

Median absolute relative error (%), median across folds. Memory length 20 x delta99. W1 excluded (not a well-behaved TBL).

> **Read with the limits in the README.** The `+history` columns use a
> kernel that is unresolvable on the 50-station wing grid, and the model
> comparison is scored on the same held-out folds used to select it. See
> *Identifiability* in the README.

## leave-one-Re-out

| target   |   M0 constant |   M1 ZPG correlation |   M2 linear, local |   M3 linear, +history |   M4 trees, local |   M5 trees, +history |   Mp linear, +placebo |
|:---------|--------------:|---------------------:|-------------------:|----------------------:|------------------:|---------------------:|----------------------:|
| Cf       |          30.7 |                 21.4 |                3.7 |                   2.3 |              12.1 |                  3.6 |                   4   |
| H        |           5.8 |                nan   |                1.2 |                   0.5 |               4.5 |                  4.5 |                   1   |
| Pk_outer |          59.2 |                nan   |                9.8 |                   8.2 |              47.3 |                 48.5 |                   9.6 |
| Ue_plus  |          19.1 |                  9.3 |                1.9 |                   1.2 |               6.4 |                  1.8 |                   2.1 |
| uu_outer |          46.8 |                nan   |                8.6 |                   7.5 |              14.4 |                  9.2 |                   8.1 |
| uv_outer |          38.8 |                nan   |                4.5 |                   3.7 |               4.8 |                  3.1 |                   5   |


## blocked streamwise

| target   |   M0 constant |   M1 ZPG correlation |   M2 linear, local |   M3 linear, +history |   M4 trees, local |   M5 trees, +history |   Mp linear, +placebo |
|:---------|--------------:|---------------------:|-------------------:|----------------------:|------------------:|---------------------:|----------------------:|
| Cf       |          40.5 |                 24.2 |                4.5 |                   2.2 |               9.8 |                  8   |                   4.9 |
| H        |           8.9 |                nan   |                1.1 |                   0.8 |               2.2 |                  2.3 |                   1   |
| Pk_outer |          85.9 |                nan   |               12.6 |                  12.4 |              15.2 |                 18.3 |                  10   |
| Ue_plus  |          29.4 |                 10.1 |                2.2 |                   1.1 |               5   |                  4.1 |                   2.4 |
| uu_outer |          65.6 |                nan   |               11.2 |                   8.6 |               9.7 |                 10   |                  11.1 |
| uv_outer |          57.6 |                nan   |                4.2 |                   1.6 |              10.3 |                 12.4 |                   3.8 |


## Error reduction from adding history (leave-one-Re-out)

Percent reduction in median relative error when the memory-kernel
features are added, within each model family. Positive = history helps.

| target   |   linear |   trees |   placebo |
|:---------|---------:|--------:|----------:|
| Cf       |     36.3 |    70.6 |      -9.6 |
| H        |     58.8 |     0.1 |      15.8 |
| Pk_outer |     15.6 |    -2.5 |       1.4 |
| Ue_plus  |     37   |    71.5 |     -10   |
| uu_outer |     13.4 |    36.5 |       6.2 |
| uv_outer |     17.8 |    35.5 |     -10.6 |


## Per-fold error, M3 linear + history (%)

Column = the Reynolds-number case held out.

| target   |   W10 |   W2 |   W4 |
|:---------|------:|-----:|-----:|
| Cf       |   5.8 |  2.3 |  0.8 |
| H        |   0.5 |  1.2 |  0.3 |
| Pk_outer |   8.2 |  8.7 |  4.3 |
| Ue_plus  |   2.8 |  1.2 |  0.4 |
| uu_outer |   7.5 |  7.9 |  4.2 |
| uv_outer |   5.1 |  3.7 |  2   |
