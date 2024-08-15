Test app
  $ python -m tcr_antigen_prediction.apps.select_tcr_pmhc_structures \
  > --log-level error \
  > --seed 123 \
  > --tcr-types abTCR \
  > --mhc-types MH1 MH2 \
  > --antigen-types peptide \
  > -o output \
  > $TESTDIR/data/stcrdab-mock/

  $ diff output/stcrdab_split.csv $TESTDIR/reference/stcrdab_split.csv

  $ diff output/3qiw_CDEAB.pdb $TESTDIR/reference/3qiw_CDEAB.pdb
  $ diff output/7q9b_DECAB.pdb $TESTDIR/reference/7q9b_DECAB.pdb
  $ diff output/7q9b_IJHFG.pdb $TESTDIR/reference/7q9b_IJHFG.pdb
