# wing-memory

**Separating local pressure gradient, flow history and Reynolds number in
turbulent wing boundary layers.**

Built on the VinuesaLab WING database: well-resolved LES of the turbulent
boundary layers around a NACA4412 wing section at `Re_c` = 100k, 200k, 400k
and 1M (5° angle of attack), plus a NACA0012 at `Re_c` = 400k.

---

## Status

| phase | state |
|---|---|
| 0. Loader, boundary-layer quantities, reproduction of the paper | **done** |
| 1. Baselines (classical correlations, local-feature GBM) | next |
| 2. History-aware model with a learned memory kernel + SHAP | planned |
| 3. Symbolic regression to a compact law | planned |
| 4. Generalisation: leave-one-Re-out, cross-geometry, APG → FPG | planned |

## Phase 0 result: the published results come back out of the published files

![Reproduction of Fig. 5](paper/figures/fig05_bl_development.png)

Every quantity tabulated in Vinuesa et al. (2018), recomputed from the
distributed data and interpolated to the exact chordwise stations the paper
quotes:

| | |
|---|---|
| quantities compared | **48** |
| median absolute error | **0.75%** |
| mean absolute error | **1.24%** |
| within 2% | **38 / 48** |
| within 5% | **47 / 48** |

The residual is dominated by the paper quoting `beta` to two significant
figures. Full table: [`paper/reproduction_report.md`](paper/reproduction_report.md).
Pinned as regression tests in [`tests/test_reproduction.py`](tests/test_reproduction.py).

```bash
python scripts/reproduce_paper.py      # figures + report
python -m pytest tests/ -q             # 31 tests
```

## Two findings from the audit

### 1. The velocity scale is `max(U)`, not the `Ue` field

Every shipped integral quantity — `Cf`, `Re_delta*`, `Re_theta`, `delta*`,
`theta` — is built on the peak tangential velocity of each profile, not on
the `Ue` field stored alongside it. Verified exactly, across all nine cases,
by inverting each quantity for the velocity it implies.

Using `Ue` instead, which is the obvious reading of the files, puts the
thicknesses out by **9–13%** and `Cf` by about **6%**. `wingtbl` exposes the
correct scale as `U_ref`.

### 2. The matched-`beta` premise holds only to tens of percent

![The beta confound](paper/figures/fig_beta_confound.png)

The paper's design holds the pressure-gradient history fixed and varies the
Reynolds number, describing the `Re_c >= 200k` cases as seeing
"approximately the same" `beta(x)` with differences "on the order of 10%"
only beyond `x/c > 0.9`. Over the analysis window `0.25 <= x/c <= 0.95` the
data gives:

| comparison | mean ratio | range |
|---|---|---|
| `beta(200k) / beta(1M)` | **1.42** | 1.17 – 1.74 |
| `beta(400k) / beta(1M)` | **1.19** | 1.14 – 1.43 |

with `beta` decreasing monotonically with `Re_c` at every station. The
paper's own Table 4 shows the same (2.6 / 2.4 / 2.0 at `x/c` = 0.7). The
apparent collapse is an artefact of a four-decade log axis.

This matters because the paper's headline result — that the wing-to-ZPG
ratios all decrease with `Re_c`, so low-Re boundary layers are more sensitive
to a given pressure gradient — rests on `beta` being matched. It is not, and
`beta` moves in the *same* direction as the effect being measured. The two
causes are confounded.

**That confound is what this project sets out to resolve**, by fitting a
model whose terms are the local pressure gradient, its upstream history and
the Reynolds number, and reporting how the measured sensitivity divides
between them.

Both findings are documented in full, with the evidence, in
[`docs/DATA_NOTES.md`](docs/DATA_NOTES.md).

## Install

```bash
pip install -e ".[dev]"
```

## Data

The download links on `vinuesalab.com/wing/` point at a personal KTH OneDrive
account that is no longer served. The files were supplied directly by the
author on request (November 2026); contact rvinuesa@umich.edu. Place
`naca4412.mat` and `naca0012.mat` in `data/raw/`. The data is not
redistributed here.

```python
from wingtbl.loader import load_wing, concat_cases

cases = load_wing("data/raw/naca4412.mat")
suction = concat_cases([cases[k] for k in ("top1n", "top2n", "top4n", "top10n")])
suction.Re_tau.sel(case="W10").max()
```

Each case is an `xarray.Dataset` over `(station, j)` — 50 chordwise stations,
52 wall-normal points — carrying the mean flow, the Reynolds-stress tensor,
the TKE budget terms and the shipped integral quantities, all in the local
tangential / wall-normal frame.

## Layout

```
wingtbl/        loader, boundary-layer quantities, plotting style
scripts/        reproduce_paper.py, beta_confound.py
tests/          regression tests against the published values
docs/           DATA_NOTES.md -- provenance audit
paper/          figures and the reproduction report
```

## Caveats carried forward

- **W1** (`Re_c` = 100k) is not a well-behaved TBL and its last station has
  separated. Out-of-distribution test case only, never training data.
- **W10** ran at coarser resolution; the paper notes a slightly
  over-predicted near-wall `u_t^2` peak.
- **`beta` is used as shipped.** It cannot be recomputed from the 50 stored
  stations to better than ~7%, because they under-resolve `dPe/dx`.
- The pressure-side TKE budgets are empty in the distributed file.

## References

1. R. Vinuesa, P. S. Negi, M. Atzori, A. Hanifi, D. S. Henningson and
   P. Schlatter. Turbulent boundary layers around wing sections up to
   `Re_c` = 1,000,000. *Int. J. Heat Fluid Flow* **72**, 86–99 (2018).
2. A. Tanarro, R. Vinuesa and P. Schlatter. Effect of adverse pressure
   gradients on turbulent wing boundary layers. *J. Fluid Mech.* **883**,
   A8 (2020).
3. A. Bobke, R. Vinuesa, R. Örlü and P. Schlatter. History effects and near
   equilibrium in adverse-pressure-gradient turbulent boundary layers.
   *J. Fluid Mech.* **820**, 667–692 (2017).
4. R. Vinuesa, A. Bobke, R. Örlü and P. Schlatter. On determining
   characteristic length scales in pressure-gradient turbulent boundary
   layers. *Phys. Fluids* **28**, 055101 (2016).
5. P. Schlatter and R. Örlü. Assessment of direct numerical simulation data
   of turbulent boundary layers. *J. Fluid Mech.* **659**, 116–126 (2010).

## Licence

MIT for the code. The WING database is the property of its authors and is
not redistributed here; cite references 1 and 2 when using it.
