Test app (standard)
  $ python -m tcr_antigen_prediction.apps.select_stcrdab_tcr_pmhc_structures \
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

Test app adding MHC-TCR Pseudo Sequences
  $  python -m tcr_antigen_prediction.apps.select_stcrdab_tcr_pmhc_structures \
  > --log-level error \
  > --seed 123 \
  > --tcr-types abTCR \
  > --mhc-types MH1 MH2 \
  > --antigen-types peptide \
  > --mhc-class-I-tcr-contact-residues 1058 1061A 1062 1063 1065 1066 1068 1069 1070 1073 1076 1077 58 62 65 66 68 69 72 73 75 76 79 \
  > --mhc-class-II-alpha-chain-tcr-contact-residues 63 65 66 68 69 70 72 73 76 \
  > --mhc-class-II-beta-chain-tcr-contact-residues 58 61A 61B 62 63 65 66 69 72 72A 76 \
  > -o output-mhc-tcr-pseudo-seqs \
  > $TESTDIR/data/stcrdab-mock/

  $ diff output-mhc-tcr-pseudo-seqs/stcrdab_split.csv $TESTDIR/reference/mhc-tcr-pseudo-seqs/stcrdab_split.csv

Test excluding PDB IDs
  $  python -m tcr_antigen_prediction.apps.select_stcrdab_tcr_pmhc_structures \
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
  $ python -m tcr_antigen_prediction.apps.select_stcrdab_tcr_pmhc_structures \
  > --log-level error \
  > --seed 123 \
  > --tcr-types abTCR \
  > --mhc-types MH1 MH2 \
  > --antigen-types peptide \
  > --remove-structures-missing-residues \
  > -o output-no-missing-residues \
  > $TESTDIR/data/stcrdab-mock/

  $ diff output-no-missing-residues/stcrdab_split.csv $TESTDIR/reference/output-no-missing-residues/stcrdab_split.csv

Test structural similarity cutoff
  $ python -m tcr_antigen_prediction.apps.select_stcrdab_tcr_pmhc_structures \
  > --log-level error \
  > --seed 123 \
  > --tcr-types abTCR \
  > --mhc-types MH1 MH2 \
  > --antigen-types peptide \
  > --structural-similarity-cutoff 2.0 \
  > -o output-structural-similarity-cutoff \
  > $TESTDIR/data/stcrdab-mock/

  $ diff output-structural-similarity-cutoff/stcrdab_split.csv $TESTDIR/reference/structural-similarity-cutoff/stcrdab_split.csv
