# Phase 4: the flat-plate test

## Was the memory length ever resolved?

Ratio of a 0.5 `delta99` memory length to the station spacing:
wing **0.47**, flat plates **4.34**. Below 1 the kernel decays
inside a single quadrature panel, so the wing fit reported a length
its own grid could not represent.

## Does the third term earn its place?

| dataset    | target   |   local_only |   memory |   dbeta_dx |
|:-----------|:---------|-------------:|---------:|-----------:|
| wing       | H        |         3.67 |     1.33 |       1    |
| wing       | Cf       |        26.37 |     3    |       3.93 |
| wing       | Ue_plus  |        16.4  |     1.46 |       1.95 |
| wing       | uu_outer |        30.98 |     8.93 |       6.32 |
| wing       | uv_outer |        20.68 |     3.34 |       3.09 |
| wing       | Pk_outer |        41.88 |     6.05 |       7.74 |
| flat plate | H        |         0.96 |     1.07 |       0.93 |
| flat plate | Cf       |         3.44 |     4.02 |       3.69 |
| flat plate | Ue_plus  |         1.67 |     1.95 |       1.8  |
| flat plate | uu_outer |         5.04 |     4.73 |       4.82 |
| flat plate | uv_outer |         3.68 |     2.61 |       2.92 |
| flat plate | Pk_outer |         4.1  |     4.08 |       4.02 |

On the wing it removes most of the error; on the near-equilibrium
flat plates it does almost nothing, and the kernel and the plain
`d beta/dx` gradient are interchangeable throughout.

## Does the headline split survive?

|          |   memory |   dbeta_dx |
|:---------|---------:|-----------:|
| H        |     33.3 |       25.9 |
| Cf       |     34.8 |       30.5 |
| Ue_plus  |     34.8 |       30.5 |
| uu_outer |     62.1 |       46.9 |
| uv_outer |     76.4 |       91.6 |
| Pk_outer |     25.9 |       19.5 |

The ordering is unchanged, but `uu_outer` crosses 50%, so that
share is parameterisation dependent and must be quoted as a range.
