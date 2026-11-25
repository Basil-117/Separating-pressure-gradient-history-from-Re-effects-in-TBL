# Memory-length scan

> **Read with the limits in the README.** The fitted memory lengths below
> are near or under one station spacing on the 50-station wing grid, where
> the kernel decays inside a single quadrature panel and cannot be
> resolved. They are a property of the grid, not a physical length. See
> *Identifiability* in the README.
Leave-one-Reynolds-number-out median relative error, as a function of
the kernel memory length in units of the local `delta99`.

## Best length per target

| target   |   length |   err |   local_only |   reduction_% | region   |   within10%_lo |   within10%_hi |
|:---------|---------:|------:|-------------:|--------------:|:---------|---------------:|---------------:|
| Cf       |    0.202 | 0.009 |        0.037 |        76.786 | wall     |          0.202 |          0.202 |
| Ue_plus  |    0.202 | 0.004 |        0.019 |        77.439 | integral |          0.202 |          0.202 |
| H        |    0.267 | 0.003 |        0.012 |        75.79  | integral |          0.267 |          0.616 |
| uv_outer |    5.749 | 0.027 |        0.045 |        39.565 | outer    |          4.349 |          5.749 |
| uu_outer |    7.6   | 0.057 |        0.086 |        33.581 | outer    |          7.6   |         10.046 |
| Pk_outer |   10.046 | 0.057 |        0.098 |        41.685 | outer    |         10.046 |         10.046 |

## Median best length by region

| region   |   length_over_delta99 |
|:---------|----------------------:|
| integral |                   0.2 |
| outer    |                   7.6 |
| wall     |                   0.2 |
