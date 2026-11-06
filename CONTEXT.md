# RealiseD simulation

A simulation study comparing methods for testing a treatment-by-time interaction in longitudinal trials with missing outcomes.

## LSPIM inference

**Deviation contrast**:
The difference between one visit's treatment effect and the mean of all visits' treatment effects; there are K of them for K visits, and they sum to zero.
_Avoid_: mean contrast, centered contrast

**Treatment hypothesis**:
The null hypothesis that all deviation contrasts are zero, i.e. the treatment effect is the same at every visit. It is the LSPIM interaction hypothesis.
_Avoid_: interaction family

**Time contrast**:
The difference between the within-arm time trends of the treated and the control arm (trend_treat − trend_ctrl).
_Avoid_: trend interaction

**Time hypothesis**:
The null hypothesis that the time contrast is zero. It is a secondary hypothesis, tested on its own and never part of the interaction decision.

**Single-step test**:
The max-|t| test over the deviation contrasts, with its p-value from the joint (multivariate t) distribution of the contrast statistics. It is the primary LSPIM interaction test.
_Avoid_: max-t test, maximum test, single_step

**Holm test**:
The Holm-adjusted per-contrast t-tests over the deviation contrasts. It is a secondary LSPIM test.

**Reference distribution**:
The null distribution used for an LSPIM test statistic: t with df N − 1 (N = number of subjects) for a single dataset, Barnard-Rubin t after multiple imputation.
_Avoid_: finite/asymptotic reference
