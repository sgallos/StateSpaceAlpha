# 2 State-Space Model

This folder contains the model code.

- `SS_age_diff.m`: the 4-D Kalman filter and RTS smoother for the baseline
  and CP-Control difference trajectories.
- `SS_age_diff_em.m`: EM estimation for the model hyperparameters.
- `SS_age_diff_bootstrap_simultaneous.m`: optional bootstrap sensitivity
  wrapper.
- `SS_age_pooled_iwp.m`: earlier 2-D pooled IWP smoother used as a scaffold.
- `estimate_sigmaBio_from_subepochs.m`: estimates empirical between-subject
  biological variance from the sub-epoch table.

The primary no-resampling inference uses `SS_age_diff.m` and reads the
posterior credible band for `delta(a)` directly from the smoother covariance.
