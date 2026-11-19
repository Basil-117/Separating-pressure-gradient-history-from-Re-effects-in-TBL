# WING database: provenance notes

Findings from auditing the published `.mat` files against Vinuesa et al.
(2018). Anyone building on this database should read the first section before
computing anything, because the obvious reading of the files is wrong by
several percent.

## 1. The velocity scale is `max(U)`, not the `Ue` field

Each station ships both a `Ue` field and a velocity profile `U`. The shipped
integral quantities are **not** built on `Ue`. They are built on the maximum
tangential velocity in the profile:

```
U_ref = max(U)        over the stored wall-normal range
```

Verified by inverting each shipped quantity for the velocity it implies:

| implied from | agreement with `max(U)` |
|---|---|
| `Cf`  = 2 (u_tau/U_ref)^2        | exact to machine precision |
| `Reds` = U_ref delta* / nu       | exact to machine precision |
| `Reth` = U_ref theta / nu        | exact to machine precision |

across **all nine cases** in both files (`top`/`bot` x four Reynolds numbers
for the NACA4412, plus the NACA0012), maximum deviation 0.000%.

The thicknesses follow the same convention, integrated to `delta99`:

```
delta* = int_0^delta99 (1 - U/U_ref)        dyn
theta  = int_0^delta99 (U/U_ref)(1 - U/U_ref) dyn
```

Recomputing them this way reproduces the shipped values to a median of 0.21%
(`delta*`) and 0.34% (`theta`). Using the `Ue` field instead gives errors of
**9-13%**, and `Cf` is then wrong by about 6%.

`Ue` sits around 0.96-0.98 of `max(U)`, with a ratio that varies from station
to station, and the sequence is visibly noisy (it is a fitted quantity, not a
sampled one). It is best read as the diagnostic-method edge velocity that
accompanies `delta99`, not as a normalising scale.

**Consequence.** `wingtbl` computes `U_ref` from the profile and uses it for
every derived quantity, and keeps `Ue` as a separate recorded field.

## 2. `delta99` does not follow from `U_ref` algebraically

`delta99` is not the height where `U = 0.99 U_ref`: that ratio is 0.984-0.988
on the suction side and varies by up to 0.056 within a single case. Nor is it
where `U = 0.99 Ue`. Both come from the diagnostic-scaling method of Vinuesa,
Bobke, Orlu & Schlatter (Phys. Fluids 28, 055101, 2016), which deliberately
does not rely on the mean profile reaching a plateau -- the whole point, since
under a pressure gradient it never does.

Take `delta99` and `Ue` as given. Do not try to re-derive them from the 52
stored points.

## 3. `beta` cannot be recomputed from the stored stations

The shipped `beta` is the published quantity and should be used as given.
Recomputing it needs `dPe/dxt`, and the files store only 50 chordwise
stations, which under-resolves that derivative:

- over `0.25 < x/c < 0.95`, recomputation agrees to a **median of ~7%**
  (best variant: pressure interpolated to `delta99`, differentiated along
  surface arc length);
- over the first few stations the recomputed value disagrees **in sign**. The
  stored pressure rises monotonically from the first station, implying an
  adverse gradient, while the shipped `beta` is negative there.

Inverting the shipped `beta` for the pressure gradient it implies gives a
ratio to the recomputed gradient of 1.06 +/- 0.07 over the analysis window --
a bias and a scatter, not a clean constant. The most likely explanation is
that the authors computed `beta` from the full simulation field, at far
higher chordwise resolution than the 50 stations distributed here.

**Consequence.** `wingtbl` uses the shipped `beta` (exposed as `beta_db`).
The history kernel integrates an interpolant of that field.

## 3b. The matched-`beta` premise holds only to tens of percent

This is the most consequential finding, and it is the reason the project
exists in its current form.

The paper's design is to hold the pressure-gradient history fixed and vary
only the Reynolds number, so that the remaining differences can be read as
Reynolds-number effects. Section 3.1 describes the `Re_c >= 200k` cases as
subjected to "approximately the same" `beta(x)`, with relative differences
"on the order of 10%" arising only beyond `x/c > 0.9`.

In the distributed data, over the analysis window `0.25 <= x/c <= 0.95`:

| comparison | mean ratio | range |
|---|---|---|
| `beta(Re_c=200k) / beta(Re_c=1M)` | **1.42** | 1.17 - 1.74 |
| `beta(Re_c=400k) / beta(Re_c=1M)` | **1.19** | 1.14 - 1.43 |

and `beta` decreases monotonically with Reynolds number at **every** station.
The paper's own Table 4 shows the same thing: at `x/c = 0.7` it lists
`beta` = 2.6, 2.4 and 2.0 for W2, W4 and W10 -- a ratio of 1.30, not 1.10.

The apparent collapse in Fig. 5 is an artefact of plotting `beta` on a
four-decade log axis, where a factor of 1.4 is a hairline offset. Plotted as
a ratio it is unmistakable (`paper/figures/fig_beta_confound.png`).

**Why it matters.** The paper's headline result is that the wing-to-ZPG
ratios (`Phi_Ue+`, `Phi_H`, `Phi_u_t^2`, `Phi_uv`, `Phi_Pk`) all *decrease*
with `Re_c`, concluding that low-Re boundary layers are more sensitive to a
given pressure gradient. But `beta` itself also decreases with `Re_c`, by
roughly 20-40%, and every one of those ratios increases with `beta`. The two
effects are **confounded and act in the same direction**, so an unknown share
of the reported Reynolds-number sensitivity is residual pressure-gradient
sensitivity.

Nothing here says the paper's conclusion is wrong -- the mechanism it
proposes is physically well motivated and supported by the wall-normal
velocity data. What is missing is the decomposition that would say how much
of the measured effect belongs to each cause. Supplying that decomposition,
with history included as its own term, is the project's contribution.

## 4. Analysis window

The paper restricts itself to `x_ss/c < 0.95`, and its quoted maxima confirm
this: the published "maximum" `Re_theta` values of 1050 / 1730 / 2930 / 6000
match the database at `x/c = 0.95` (1029 / 1743 / 2927 / 6203), not at the
last stored station `x/c = 0.983` (1112 / 1900 / 3187 / 6900).

Lower bound: the boundary layers are tripped at `x/c = 0.1`, and `beta` is
negative over roughly the first five stations. This project uses
**`0.25 <= x/c <= 0.95`**, matching the paper and the later XAI study on the
same wing.

## 5. Per-case caveats

- **W1 (`Re_c` = 100k)** is not a well-behaved TBL. The paper demonstrates
  this via the modified diagnostic plot, and it reaches `Re_delta*` = 450 only
  at `x/c` ~ 0.35. Its **last station (`x/c` = 0.983) has separated**: the
  shipped `beta` there is -7506, i.e. `u_tau` has passed through zero. Treat
  W1 as an out-of-distribution test case, never as training data.
- **W10 (`Re_c` = 1M)** was run at coarser resolution than the other three.
  The paper reports a slight over-prediction of the near-wall `u_t^2` peak and
  mild attenuation of the wall-normal and spanwise fluctuations.
- **Pressure side** (`bot*`): the TKE budget fields (`Pk`, `Dk`, `Tk`, `VDk`,
  `VPGk`) are **empty** in the distributed file. The loader fills them with
  NaN. Mean flow and Reynolds stresses are present.
- **NACA0012** ships a single case, `top4n12` (suction side, `Re_c` = 400k).
  Its statistics exploit the symmetry of the flow.

## 6. Units and frame

Lengths are normalised by chord, velocities by freestream, so `nu = 1/Re_c`
exactly -- the loader asserts this as a check on the Reynolds-number code in
the variable names. All turbulence statistics are in the local tangential /
wall-normal frame with the tensors rotated accordingly.

## 7. Obtaining the data

The download links on `vinuesalab.com/wing/` point at a personal KTH OneDrive
account that is no longer served. The files used here were supplied directly
by the author on request (November 2026). Contact rvinuesa@umich.edu.

Files, after extraction:

| file | md5 | contents |
|---|---|---|
| `naca4412.mat` | `a019ac85b30bbd428288f3456d41d864` | 8 cases, 50 stations each |
| `naca0012.mat` | `803540516df1503eb5b041e1c94fa95b` | 1 case, 50 stations |

## 8. Modelling notes from Phase 2

Two findings about the *model* rather than the data, recorded here because
they govern how its numbers should be read.

**The memory length is weakly identified.** The fitted law
`y = C beta^a1 (beta_bar/beta)^a2 Re_tau^a3` recovers all four parameters to
four significant figures from noiseless synthetic data. Adding 0.2% noise in
log space -- well below the scatter in the LES statistics -- moves the fitted
`ell` by 15-27% and `a2` by up to 65%, because a shorter memory shrinks
`log(beta_bar/beta)` and a larger `|a2|` compensates almost exactly. Quote
`ell` to a factor of ~1.3, never to two digits.

The decomposition uses the product `a2 * log(beta_bar/beta)`, in which the
trade-off cancels: forcing `ell` to half or double its fitted value moves the
pressure-gradient share by at most 7 points, usually 1-3.

**Do not parameterise with `beta` and `beta_bar` as separate powers.** They
are nearly collinear -- at short memory `beta_bar -> beta` -- and the fit
returns equal-and-opposite exponents around +/- 7 that are unstable and
meaningless. Since `beta^a1 beta_bar^a2 = beta^(a1+a2) (beta_bar/beta)^a2`,
the ratio form is the same model class with separable terms.

**Correction to an earlier Phase 1 claim.** A scan of `ell` against held-out
error appeared to show wall quantities preferring `ell ~ 0.2 delta99` and
outer-region quantities `5.7-10 delta99`. That scan selected `ell` on the
test set and used a collinear feature set. With the length fitted inside each
training fold and the ratio parameterisation, every target lands at
`ell ~ 0.4-0.6 delta99`. There is no order-of-magnitude separation by region.

## 9. Phase 3: what generalises and what does not

**The law transfers across geometry.** Fitted on the NACA4412 at 5 deg and
applied unchanged to the NACA0012 at 0 deg -- a different section, a
different angle of attack, and a pressure gradient reaching beta = 3.3
rather than 46 -- the median error is 2.5-9.8% on five of six targets
(17.7% on TKE production, which is the hardest quantity everywhere).

**It does not extrapolate far in beta.** Trained on `x/c <= 0.70`
(beta <= 2.5) and asked about `x/c > 0.70` (beta up to 49), the error rises
to 20-35% on Cf, uu and Pk. Use the law inside its calibration range.

**It has no favourable-pressure-gradient limit.** The pressure side has
beta < 0 at 88-95% of stations in the analysis window, and `beta^a` is not
real there. Any extension to FPG needs a different functional form, not a
refit.

**Local pressure gradient and history multiply.** An exhaustive search over
255 expressions from a nine-term library chose the interaction term
`log(beta) * log(beta_bar/beta)` for all six targets -- out of nine
candidates, in six independent selections. The exponent on beta therefore
depends on how far the upstream history departs from the local value, so
the two are not separable factors as the Phase 2 power law assumes. Next
most frequent: `log^2(beta)` (5/6) and `log(beta) * log(Re_tau)` (4/6).

**Selection optimism is large.** Leave-one-Reynolds-number-out error, used
to pick an expression out of ~1000 candidates, understates the true held-out
error by a median factor of **4.4** (range 3.2-5.9). Scored honestly on the
unseen NACA0012, the searched four-term expressions beat the three-term
power law on only four of six targets. Any symbolic-regression result on a
dataset this size should be reported with a genuinely held-out score, not
with the criterion it was selected on.

## 10. The APG flat-plate database, and what it revised

Obtained on request from R. Vinuesa, November 2026 (`APG.mat`, 306 MB). Five
cases from Bobke et al. (2017): `b1`, `b2` at nominally constant `beta` of 1
and 2, and `m13`, `m16`, `m18` following a power-law freestream
`U_e ~ x^m`. 1570 chordwise stations per case, 301-401 wall-normal points.

### The velocity scale is the opposite of WING's

Every shipped integral quantity here is built on the **`Ue` field**. In WING
they are built on **`max(U)`** and `Ue` is something else (section 1).
Verified by inverting `Cf`, `Reds` and `Reth` for the velocity each implies:
all three agree with `Ue` to six decimal places across all five cases.
Two databases, one lab, one file format, opposite conventions. Reusing the
wing loader here would be silently wrong.

Units are inlet-based: `x` runs 200-2500, `nu` = 1/450, lengths in inlet
displacement thicknesses.

### Analysis window

`700 <= x <= 1900`. Below 700 the pressure gradient is still being
established; above 1900 `beta` collapses towards zero and turns negative as
the outflow boundary is approached, which is a domain-end artefact.

### The memory length was never resolved on the wing

A memory length of `0.5 delta99` corresponds to `ell/dx` of **0.47** on the
wing grid and **4.3** on the flat plates (at stride 4; 20.5 at full
resolution). Below 1 the exponential kernel decays inside a single
quadrature panel, `beta_bar` collapses onto the local `beta`, and the
feature degenerates into a finite-difference gradient. **The wing fit of
`ell` = 0.43-0.59 reported a length its own grid could not represent**, and
the Phase 1 and Phase 2 statements about memory length should be read as
statements about a local gradient.

### What the third term actually measures

The local rate of change of `beta`. On the wing, replacing the kernel with
a plain `(d beta/dx) * delta99 / beta` term gives interchangeable accuracy
(`H` 1.33 vs 1.00%, `Cf` 3.00 vs 3.93%, `uv` 3.34 vs 3.09%). The two
parameterisations are the same feature to first order, since the kernel's
`ell -> 0` limit is exactly that gradient.

### The term tracks departure from equilibrium

| | wing | flat plate |
|---|---|---|
| local features only | 3.7-41.9% | 1.0-5.0% |
| best third term | 1.0-6.3% | 0.9-4.0% |
| error removed | **73-91%** | **-8 to +29%** |

On the wing, local `beta` and `Re_tau` alone are badly insufficient and the
third term removes most of the error. On the flat plates they are already
sufficient to 1-5%, and the third term adds nothing (it makes `Cf` and
`Ue_plus` slightly worse).

This is coherent rather than contradictory. The wing's `beta` climbs from
0.3 to 46 across the window; the flat-plate cases were constructed to sit
near equilibrium and their `beta` is nearly flat. The term measures how far
the layer is from equilibrium, and correctly reports "not far" where that is
true. Bobke et al. titled their paper "History effects and *near
equilibrium*" for the same reason.

### Consequence for the headline split

Pressure-gradient share of the `W2 -> W10` change, under both
parameterisations:

| target | kernel | `d beta/dx` |
|---|---|---|
| `uv_outer` | 76.4% | 91.6% |
| `uu_outer` | 62.1% | 46.9% |
| `Cf` | 34.8% | 30.5% |
| `Ue_plus` | 34.8% | 30.5% |
| `H` | 33.3% | 25.9% |
| `Pk_outer` | 25.9% | 19.5% |

The ordering is unchanged and the conclusion stands for `uv_outer`
(76-92%, pressure-gradient dominated) and for the wall, integral and
production quantities (70-80% Reynolds number). But `uu_outer` crosses 50%,
so it must be quoted as a range, 47-62%, and not as "64% pressure
gradient".
