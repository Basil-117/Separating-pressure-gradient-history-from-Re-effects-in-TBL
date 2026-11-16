# Phase 3b: symbolic search

Exhaustive over all 255 expressions of up to 4 terms from a 9-term library, at each of 4 memory lengths.

## Best expression per target

| target   |   k |   ell |   searched_loo_% |   searched_0012_% |   powerlaw_loo_% |   powerlaw_0012_% |   optimism_x |
|:---------|----:|------:|-----------------:|------------------:|-----------------:|------------------:|-------------:|
| H        |   4 |  1    |             0.31 |              1.77 |             0.74 |              2.18 |          5.6 |
| Cf       |   4 |  0.45 |             0.97 |              4.24 |             2.32 |              7.92 |          4.4 |
| Ue_plus  |   4 |  0.45 |             0.49 |              2.19 |             1.14 |              4.21 |          4.5 |
| uu_outer |   4 |  3    |             3.44 |             11.09 |             7.57 |              9.74 |          3.2 |
| uv_outer |   4 |  3    |             1.5  |              0.86 |             3.91 |              4.26 |          0.6 |
| Pk_outer |   4 |  1    |             3.67 |             21.77 |             6.93 |             18.47 |          5.9 |

`searched_loo_%` was used for selection and is therefore
optimistic. `searched_0012_%` and `powerlaw_0012_%` are both fitted
on the NACA4412 and scored once on the NACA0012, which took no part
in selection -- the only fair comparison of the two forms.
`optimism_x` is how much worse the honest number is than the one
selection reported.

Terms chosen:

- `H`: log R + log b log R + log m log R + log b log m
- `Cf`: log^2 b + log^2 R + log b log R + log b log m
- `Ue_plus`: log^2 b + log^2 R + log b log R + log b log m
- `uu_outer`: log m + log R + log^2 b + log b log m
- `uv_outer`: log b + log^2 b + log b log R + log b log m
- `Pk_outer`: log R + log^2 b + log m log R + log b log m

## Term frequency across the six targets

|             |   times chosen (of 6) |
|:------------|----------------------:|
| log b log m |                     6 |
| log^2 b     |                     5 |
| log b log R |                     4 |
| log R       |                     3 |
| log m log R |                     2 |
| log^2 R     |                     2 |
| log m       |                     1 |
| log b       |                     1 |

## Pareto front

| target   |   k |        err |   ell | terms                                           |
|:---------|----:|-----------:|------:|:------------------------------------------------|
| H        |   1 | 0.0423878  |  1    | m - 1                                           |
| H        |   2 | 0.0183133  |  1    | log^2 R + log m log R                           |
| H        |   3 | 0.00478829 |  0.2  | log^2 b + log^2 R + log b log R                 |
| H        |   4 | 0.00313752 |  1    | log R + log b log R + log m log R + log b log m |
| Cf       |   1 | 0.143642   |  3    | log m log R                                     |
| Cf       |   2 | 0.0296558  |  0.45 | log m + log R                                   |
| Cf       |   3 | 0.0119715  |  0.2  | log^2 R + log b log R + log m log R             |
| Cf       |   4 | 0.0097014  |  0.45 | log^2 b + log^2 R + log b log R + log b log m   |
| Ue_plus  |   1 | 0.0736049  |  0.2  | log^2 b                                         |
| Ue_plus  |   2 | 0.0151661  |  0.45 | log m + log R                                   |
| Ue_plus  |   3 | 0.00598569 |  0.2  | log^2 R + log b log R + log m log R             |
| Ue_plus  |   4 | 0.00485239 |  0.45 | log^2 b + log^2 R + log b log R + log b log m   |
| uu_outer |   1 | 0.0811994  |  1    | log m                                           |
| uu_outer |   2 | 0.0514792  |  1    | log m + log^2 R                                 |
| uu_outer |   3 | 0.0392459  |  1    | log m + log^2 b + log^2 R                       |
| uu_outer |   4 | 0.0343604  |  3    | log m + log R + log^2 b + log b log m           |
| uv_outer |   1 | 0.0774109  |  0.45 | m - 1                                           |
| uv_outer |   2 | 0.0299401  |  0.45 | log b + log m log R                             |
| uv_outer |   3 | 0.0185954  |  0.2  | log b + log m + m - 1                           |
| uv_outer |   4 | 0.0150329  |  3    | log b + log^2 b + log b log R + log b log m     |
| Pk_outer |   1 | 0.471547   |  1    | m - 1                                           |
| Pk_outer |   2 | 0.0989221  |  0.45 | log R + log m log R                             |
| Pk_outer |   3 | 0.0491121  |  3    | log R + log b log m + m - 1                     |
| Pk_outer |   4 | 0.0367055  |  1    | log R + log^2 b + log m log R + log b log m     |
