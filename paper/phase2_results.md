# Phase 2: fitted memory and the Re / pressure-gradient split

> **Read with the limits in the README.** The pressure-gradient / Reynolds
> split below is **withdrawn**. The design matrix conditions at 18-57, so
> the coefficient split is not identified, and the pressure-gradient share
> pools the `beta` and `Lambda` terms rather than isolating the `beta`
> mismatch. The fitted lengths are also below the station spacing. Kept as
> the record of what was computed. See *Identifiability* in the README.
## Held-out accuracy, memory length refitted inside each fold

| target   |   ell_median |   ell_min |   ell_max |   a_beta |   a_hist |   a_Re |   err_% |
|:---------|-------------:|----------:|----------:|---------:|---------:|-------:|--------:|
| H        |        0.587 |     0.574 |     0.7   |    0.028 |   -0.875 | -0.07  |    1.15 |
| Cf       |        0.425 |     0.405 |     0.481 |   -0.039 |    8.143 | -0.345 |    2.36 |
| Ue_plus  |        0.425 |     0.405 |     0.481 |    0.019 |   -4.072 |  0.172 |    1.2  |
| uu_outer |        0.463 |     0.461 |     0.604 |    0.09  |   -7.452 | -0.139 |    9.24 |
| uv_outer |        0.449 |     0.414 |     0.521 |    0.145 |   -6.666 |  0.062 |    3.82 |
| Pk_outer |        0.469 |     0.461 |     0.579 |    0.168 |  -11.364 | -0.965 |    6.08 |

`ell` is in multiples of the local `delta99`; `err_%` is the median
relative error on the held-out Reynolds number.

## Decomposition of the W2 -> W10 change

| target   |   observed |   model |   pressure_gradient |   reynolds |   PG_share_% |
|:---------|-----------:|--------:|--------------------:|-----------:|-------------:|
| H        |    -0.1198 | -0.1172 |             -0.0399 |    -0.0773 |      34.0305 |
| Cf       |    -0.1828 | -0.1794 |              0.2066 |    -0.386  |      34.8668 |
| Ue_plus  |     0.0914 |  0.0897 |             -0.1033 |     0.193  |      34.8668 |
| uu_outer |    -0.3837 | -0.3754 |             -0.2388 |    -0.1366 |      63.6209 |
| uv_outer |    -0.1627 | -0.1454 |             -0.2083 |     0.0628 |      76.8245 |
| Pk_outer |    -1.4315 | -1.4168 |             -0.3756 |    -1.0412 |      26.5108 |

Mean over `0.30 <= x/c <= 0.90` of the change in `log y`, split
exactly by the fitted law into a pressure-gradient part and a
Reynolds-number part.

## Robustness: the same split for a second pair of cases

| target   |   W2->W10 |   W4->W10 |
|:---------|----------:|----------:|
| H        |      34   |      31.2 |
| Cf       |      34.9 |      31.1 |
| Ue_plus  |      34.9 |      31.1 |
| uu_outer |      63.6 |      60.4 |
| uv_outer |      76.8 |      73.5 |
| Pk_outer |      26.5 |      23.7 |

## Sensitivity to the memory length

The fitted length is identifiable only to tens of percent, because a
shorter memory and a larger history exponent trade off. The split is
built from their product, so it barely moves when the length is
forced to half or double the fitted value.

| target   |   fitted |   ell x 0.5 |   ell x 2 |
|:---------|---------:|------------:|----------:|
| H        |     34   |        34.3 |      33.4 |
| Cf       |     34.9 |        33.3 |      36.1 |
| Ue_plus  |     34.9 |        33.3 |      36.1 |
| uu_outer |     63.6 |        61.1 |      59.3 |
| uv_outer |     76.8 |        83.7 |      80.2 |
| Pk_outer |     26.5 |        24.7 |      25.5 |
