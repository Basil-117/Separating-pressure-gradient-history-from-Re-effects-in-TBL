# Reproduction of Vinuesa et al. (2018)

Values are interpolated to the exact chordwise station quoted in the
paper. Residual differences in `beta` are dominated by the paper
quoting it to two significant figures.

## Summary

- quantities compared: **48**
- median absolute error: **0.75%**
- mean absolute error: **1.24%**
- within 2%: **38/48**, within 5%: **47/48**

## Tables

### x_ss/c = 0.4 (paper Table 3)

| case | quantity | paper | reproduced | diff |
|---|---|---|---|---|
| W1 | beta | 1.6 | 1.62 | +1.0% |
| W1 | Re_tau | 76 | 76 | -0.5% |
| W1 | Re_theta | 276 | 289 | +4.8% |
| W1 | Cf | 2.800e-03 | 2.752e-03 | -1.7% |
| W1 | H | 2.040 | 2.023 | -0.8% |
| W2 | beta | 0.67 | 0.67 | +0.0% |
| W2 | Re_tau | 160 | 157 | -2.1% |
| W2 | Re_theta | 452 | 461 | +1.9% |
| W2 | Cf | 4.400e-03 | 4.430e-03 | +0.7% |
| W2 | H | 1.680 | 1.689 | +0.6% |
| W4 | beta | 0.67 | 0.641 | -4.4% |
| W4 | Re_tau | 241 | 239 | -1.0% |
| W4 | Re_theta | 759 | 753 | -0.8% |
| W4 | Cf | 3.900e-03 | 3.944e-03 | +1.1% |
| W4 | H | 1.610 | 1.609 | -0.1% |
| W10 | beta | 0.58 | 0.559 | -3.6% |
| W10 | Re_tau | 449 | 446 | -0.7% |
| W10 | Re_theta | 1465 | 1487 | +1.5% |
| W10 | Cf | 3.400e-03 | 3.401e-03 | +0.0% |
| W10 | H | 1.500 | 1.498 | -0.1% |

### x_ss/c = 0.7 (paper Table 4)

| case | quantity | paper | reproduced | diff |
|---|---|---|---|---|
| W1 | beta | 3.7 | 3.65 | -1.5% |
| W1 | Re_tau | 136 | 135 | -0.6% |
| W1 | Re_theta | 566 | 561 | -0.9% |
| W1 | Cf | 3.000e-03 | 3.021e-03 | +0.7% |
| W1 | H | 1.920 | 1.924 | +0.2% |
| W2 | beta | 2.6 | 2.75 | +5.7% |
| W2 | Re_tau | 228 | 223 | -2.3% |
| W2 | Re_theta | 901 | 901 | +0.1% |
| W2 | Cf | 3.000e-03 | 3.024e-03 | +0.8% |
| W2 | H | 1.730 | 1.766 | +2.1% |
| W4 | beta | 2.4 | 2.39 | -0.4% |
| W4 | Re_tau | 360 | 353 | -2.0% |
| W4 | Re_theta | 1506 | 1474 | -2.1% |
| W4 | Cf | 2.800e-03 | 2.801e-03 | +0.0% |
| W4 | H | 1.680 | 1.685 | +0.3% |
| W10 | beta | 2 | 1.96 | -1.9% |
| W10 | Re_tau | 671 | 669 | -0.3% |
| W10 | Re_theta | 2877 | 2896 | +0.7% |
| W10 | Cf | 2.500e-03 | 2.483e-03 | -0.7% |
| W10 | H | 1.580 | 1.578 | -0.2% |

### Maxima quoted in section 3.1

| quantity | case | paper | reproduced | diff |
|---|---|---|---|---|
| Re_theta at x/c=0.95 | W1 | 1050 | 1029 | -2.0% |
| Re_theta at x/c=0.95 | W2 | 1730 | 1743 | +0.7% |
| Re_theta at x/c=0.95 | W4 | 2930 | 2927 | -0.1% |
| Re_theta at x/c=0.95 | W10 | 6000 | 6203 | +3.4% |
| max Re_tau | W1 | 139 | 138 | -0.4% |
| max Re_tau | W2 | 231 | 228 | -1.5% |
| max Re_tau | W4 | 366 | 367 | +0.3% |
| max Re_tau | W10 | 707 | 706 | -0.2% |
