# MFDB Implementation Checklist

## Locked Parameters

- `electrodeLabel = 'F3'`
- `epochDurationSec = 120`
- `TW = 6`
- `K = 11`
- `Bsubj = 2000`
- `Bgroup = 2000`
- `fpass = [1 40]`
- Subject-level MFDB in linear power
- Convert subject spectra to dB after MFDB
- Group analysis in dB
- Primary significance: pointwise CI exclusion of zero
- No contiguous-frequency criterion in the primary pipeline

## Subject-Level Pipeline

1. Finalize inclusion list and F3 QC.
- If `F3` is bad or missing, exclude the subject.
- Record exclusion reason in the registry.

2. Run `run_subject_mfdb_bootstrap.m` for each included subject.
- Input is one 2-minute continuous F3 signal.
- Use the full 2-minute epoch, not short windows.
- Save linear and dB spectra plus metadata.

3. Validate subject-level outputs.
- Same `freq`
- Same `Fs`
- Same `T`
- Same `N`
- Same `TW`
- Same `K`
- Same `Bsubj`
- Same `fpass`

## Group-Level Pipeline

4. Build a subject registry.
- `subjectID`
- `groupLabel`
- `filePath`
- `F3QualityPass`
- `excluded`
- `exclusionReason`

5. Compute observed group medians in dB.
- CP median spectrum
- Control median spectrum
- Group difference

6. Run group-level bootstrap with `Bgroup = 2000`.
- Draw one bootstrap spectrum per subject per iteration
- Compute within-group medians
- Store CP minus Control difference

7. Compute pointwise 95% CIs and `sigMask`.
- `sigMask = (CI_low > 0) | (CI_high < 0)`

8. Plot group-level results.
- Original group spectra
- Group difference with 95% CI
- Significance mask

9. Save final outputs.
- Group spectra
- Bootstrap differences
- CIs
- `sigMask`
- Metadata
- Excluded subjects and reasons
