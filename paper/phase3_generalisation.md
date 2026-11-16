# Phase 3a: generalisation

Median relative error (%), by test.

| target   |   A interpolate in Re |   B extrapolate in beta |   C change the geometry |
|:---------|----------------------:|------------------------:|------------------------:|
| H        |                  0.57 |                    2.99 |                    2.54 |
| Cf       |                  2.03 |                   19.96 |                    7.42 |
| Ue_plus  |                  1.01 |                    8.7  |                    3.93 |
| uu_outer |                  4.72 |                   20.96 |                    9.81 |
| uv_outer |                  2.59 |                   13.54 |                    4.24 |
| Pk_outer |                  5.11 |                   34.98 |                   17.72 |

- **A** holds out Re_c = 400k and trains on 200k and 1M.
- **B** trains on `x/c <= 0.70` (beta <= 2.5) and predicts
  `x/c > 0.70` (beta up to 49.5).
- **C** trains on the NACA4412 at 5 deg and predicts the NACA0012 at
  0 deg, a different section with a far milder pressure gradient.

The pressure side cannot be tested: it is a favourable pressure
gradient with beta < 0 at 88-95% of stations, and a law containing
beta^a has no real value there. The fitted form is an
adverse-pressure-gradient law and does not continue through beta = 0.
