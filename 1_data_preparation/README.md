# 1 Data Preparation

This folder contains scripts that turn the existing MFDB spectral outputs into
the alpha-power tables used by the state-space model.

Run order:

```matlab
cd('/Users/gallo/projects/StateSpaceAlpha')
run('1_data_preparation/smoke_test_alpha_extraction.m')
run('1_data_preparation/check_subepoch_consistency.m')
run('1_data_preparation/collapse_subepoch_alpha_table_for_ssm.m')
```

`load_one_subject_data.m` is a minimal inspection script. Edit `subjectID` at
the top, run it, and it loads that subject's raw preprocessed file into the
workspace as `data` and `HDR`.
`smoke_test_alpha_extraction.m` extracts alpha power from the saved MFDB
spectra and checks that the values match the existing pipeline.  
`check_subepoch_consistency.m` summarizes and plots the four 30 s alpha
values per subject so unstable epochs are easy to spot.
`collapse_subepoch_alpha_table_for_ssm.m` collapses the 4 x 30 s sub-epoch
alpha table into a 19-row subject table with a more realistic observation
variance.
