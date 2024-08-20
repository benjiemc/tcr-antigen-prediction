Test app (standard)
  $ python -m tcr_antigen_prediction.apps.select_tcr_pmhc_structures \
  > --log-level error \
  > --seed 123 \
  > --tcr-types abTCR \
  > --mhc-types MH1 MH2 \
  > --antigen-types peptide \
  > -o output \
  > $TESTDIR/data/stcrdab-mock/

  $ diff output/stcrdab_split.csv $TESTDIR/reference/standard/stcrdab_split.csv

  $ diff output/3qiw_CDEAB.pdb $TESTDIR/reference/standard/3qiw_CDEAB.pdb
  $ diff output/7q9b_DECAB.pdb $TESTDIR/reference/standard/7q9b_DECAB.pdb
  $ diff output/7q9b_IJHFG.pdb $TESTDIR/reference/standard/7q9b_IJHFG.pdb

Test excluding PDB IDs
  $  python -m tcr_antigen_prediction.apps.select_tcr_pmhc_structures \
  > --log-level error \
  > --seed 123 \
  > --tcr-types abTCR \
  > --mhc-types MH1 MH2 \
  > --antigen-types peptide \
  > --pdb-ids-to-exclude 3qiw \
  > -o output-exclude-pdb-ids \
  > $TESTDIR/data/stcrdab-mock/

  $ diff output-exclude-pdb-ids/stcrdab_split.csv $TESTDIR/reference/exclude-pdb-ids/stcrdab_split.csv

Test removing structures missing residues
  $ python -m tcr_antigen_prediction.apps.select_tcr_pmhc_structures \
  > --log-level error \
  > --seed 123 \
  > --tcr-types abTCR \
  > --mhc-types MH1 MH2 \
  > --antigen-types peptide \
  > --remove-structures-missing-residues \
  > -o output-no-missing-residues \
  > $TESTDIR/data/stcrdab-mock/

  $ diff output-no-missing-residues/stcrdab_split.csv $TESTDIR/reference/output-no-missing-residues/stcrdab_split.csv