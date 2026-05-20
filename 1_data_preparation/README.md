# 1 Data Preparation

This folder contains scripts that turn the existing MFDB spectral outputs into
the alpha-power tables used by the state-space model.

Run order:

```matlab
cd('/Users/gallo/projects/StateSpaceAlpha')
run('1_data_preparation/smoke_test_alpha_extraction.m')
run('1_data_preparation/collapse_subepoch_alpha_table_for_ssm.m')
```

`smoke_test_alpha_extraction.m` extracts alpha power from the saved MFDB
spectra and checks that the values match the existing pipeline.  
`collapse_subepoch_alpha_table_for_ssm.m` collapses the 4 x 30 s sub-epoch
alpha table into a 19-row subject table with a more realistic observation
variance.
