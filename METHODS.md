# Methods

## Data

Each subject contributes frontal EEG from F3 during stable sevoflurane
maintenance. The original MFDB pipeline extracts alpha power from a 120 s
epoch. The sub-epoch workflow splits that epoch into four 30 s windows, runs
the multitaper frequency-domain bootstrap (MFDB) for each window, and then
collapses the four alpha estimates into a 19-row subject table.

For each subject, the SSM table contains:

- `subjectID`
- `groupLabel`: CP or Control
- `ageYears`
- `alpha_dB`: alpha-band power in dB
- `mfdb_var`: observation variance for alpha_dB

The alpha band is 8 to 13 Hz. Power is integrated on the linear scale and
then converted to dB:

```matlab
alpha_linear = trapz(freq(alphaMask), 10.^(spectrum_dB(alphaMask) / 10));
alpha_dB = 10 * log10(max(alpha_linear, eps));
```

## State-Space Model

The latent state at age `a` is:

```text
x(a) = [f0(a); f0_slope(a); delta(a); delta_slope(a)]
```

where:

- `f0(a)` is the baseline trajectory. With the current group coding, this is
  the Control trajectory.
- `delta(a)` is the CP minus Control alpha-power difference.
- The two slope terms allow local-linear evolution across age.

The observation equation for subject `k` is:

```text
y_k = [1, 0, g_k, 0] * x(a_k) + e_k
```

where `g_k = 0` for Control and `g_k = 1` for CP. Therefore:

- Control: `y_k = f0(a_k) + e_k`
- CP: `y_k = f0(a_k) + delta(a_k) + e_k`

The observation variance is:

```text
V_k = mfdb_var_k + sigmaBio
```

`mfdb_var_k` captures uncertainty in the subject alpha estimate. `sigmaBio`
captures residual between-subject biological variability.

## Smoothness Prior

Both `f0(a)` and `delta(a)` follow an integrated Wiener process (IWP) prior.
For an age gap `h`, the transition block is:

```text
A(h) = [1 h
        0 1]
```

and the process covariance block is:

```text
Q_iwp(h) = [h^3/3  h^2/2
            h^2/2  h]
```

The full 4-D process covariance is:

```text
Q(h) = blockdiag(q_f0 * Q_iwp(h), q_delta * Q_iwp(h))
```

The two process-noise values `q_f0` and `q_delta` control how much the
baseline and group-difference trajectories are allowed to bend with age.

## Inferential Target

The inferential target is the smoother posterior for `delta(a)`.

The posterior mean is:

```matlab
delta_mean = smoothedState(3, :).';
```

The posterior standard deviation is:

```matlab
delta_sd = sqrt(smoothedCovariance(3, 3, :));
```

The pointwise 95% credible band is:

```matlab
delta_mean +/- 1.96 * delta_sd
```

No subject resampling is used for the primary no-resampling result.

## Hyperparameter Sensitivity

Because the cohort is small, joint variance-component estimation is unstable
in some settings. We therefore report a no-subject-resampling matrix of six
hyperparameter strategies:

- A: fixed `qInit`, fixed empirical `sigmaBio`
- B: EM-estimated `q_f0` and `q_delta`, fixed empirical `sigmaBio`
- C: joint EM for `q_f0`, `q_delta`, and `sigmaBio`
- D: joint EM on the direct 76-row sub-epoch table
- E: fixed `qInit`, EM-estimated `sigmaBio`
- F: profile-likelihood `q_f0` and `q_delta`, fixed empirical `sigmaBio`

Runs A, E, and F are the cleanest no-resampling sensitivity checks. Runs B
and C are useful but did not converge cleanly. Run D is diagnostic only
because it hit the smoothness ceiling.

## What This Analysis Does Not Claim

This analysis does not claim that the exact age interval where zero is
excluded is invariant. It varies across reasonable hyperparameter choices.
The stable finding is qualitative: CP alpha power is lower than Control alpha
power, by roughly 5 dB, with the difference most clearly separated in the
older part of the 10 to 15.5 year range.
