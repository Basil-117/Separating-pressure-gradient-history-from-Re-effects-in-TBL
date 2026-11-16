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
