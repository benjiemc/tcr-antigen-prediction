Test app (standard)
  $ python -m tcr_antigen_prediction.data.apps.select_stcrdab_tcr_pmhc_structures \
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
  $  python -m tcr_antigen_prediction.data.apps.select_stcrdab_tcr_pmhc_structures \
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

Test removing structures missing residues
  $ python -m tcr_antigen_prediction.data.apps.select_stcrdab_tcr_pmhc_structures \
  > --log-level error \
  > --seed 123 \
  > --tcr-types abTCR \
  > --mhc-types MH1 MH2 \
  > --antigen-types peptide \
  > --remove-structures-missing-residues \
  > -o output-no-missing-residues \
  > $TESTDIR/data/stcrdab-mock/

  $ diff output-no-missing-residues/stcrdab_split.csv $TESTDIR/reference/output-no-missing-residues/stcrdab_split.csv

Test fixing structures missing reisdues
  $ python -m tcr_antigen_prediction.data.apps.select_stcrdab_tcr_pmhc_structures \
  > --log-level error \
  > --seed 123 \
  > --tcr-types abTCR \
  > --mhc-types MH1 MH2 \
  > --antigen-types peptide \
  > --remove-structures-missing-residues \
  > --fix-structures-missing-residues \
  > -o output-fixed-missing-residues \
  > $TESTDIR/data/stcrdab-mock/

  $ diff output-fixed-missing-residues/stcrdab_split.csv $TESTDIR/reference/fixed-missing-residues/stcrdab_split.csv
  $ diff output-fixed-missing-residues/3qiw_CDEAB.pdb $TESTDIR/reference/fixed-missing-residues/3qiw_CDEAB.pdb
  $ diff output-fixed-missing-residues/6v19_DECAB.pdb $TESTDIR/reference/fixed-missing-residues/6v19_DECAB.pdb

Test structural similarity cutoff
  $ python -m tcr_antigen_prediction.data.apps.select_stcrdab_tcr_pmhc_structures \
  > --log-level error \
  > --seed 123 \
  > --tcr-types abTCR \
  > --mhc-types MH1 MH2 \
  > --antigen-types peptide \
  > --structural-similarity-cutoff 2.0 \
  > -o output-structural-similarity-cutoff \
  > $TESTDIR/data/stcrdab-mock/

  $ diff output-structural-similarity-cutoff/stcrdab_split.csv $TESTDIR/reference/structural-similarity-cutoff/stcrdab_split.csv

Test crop structures
  $ python -m tcr_antigen_prediction.data.apps.select_stcrdab_tcr_pmhc_structures \
  > --log-level error \
  > --seed 123 \
  > --tcr-types abTCR \
  > --mhc-types MH1 MH2 \
  > --antigen-types peptide \
  > --crop-structures \
  > -o output-crop \
  > $TESTDIR/data/stcrdab-mock/

  $ diff output-crop/stcrdab_split.csv $TESTDIR/reference/crop/stcrdab_split.csv

  $ diff output-crop/3qiw_CDEAB.pdb $TESTDIR/reference/crop/3qiw_CDEAB.pdb
  $ diff output-crop/7q9b_DECAB.pdb $TESTDIR/reference/crop/7q9b_DECAB.pdb
  $ diff output-crop/7q9b_IJHFG.pdb $TESTDIR/reference/crop/7q9b_IJHFG.pdb

Test remove HETATMs
  $ python -m tcr_antigen_prediction.data.apps.select_stcrdab_tcr_pmhc_structures \
  > --log-level error \
  > --seed 123 \
  > --tcr-types abTCR \
  > --mhc-types MH1 MH2 \
  > --antigen-types peptide \
  > --remove-het-atoms \
  > -o output-no-het-atoms \
  > $TESTDIR/data/stcrdab-mock/

  $ diff output-no-het-atoms/stcrdab_split.csv $TESTDIR/reference/no-het-atoms/stcrdab_split.csv

  $ diff output-no-het-atoms/3qiw_CDEAB.pdb $TESTDIR/reference/no-het-atoms/3qiw_CDEAB.pdb
  $ diff output-no-het-atoms/7q9b_DECAB.pdb $TESTDIR/reference/no-het-atoms/7q9b_DECAB.pdb
  $ diff output-no-het-atoms/7q9b_IJHFG.pdb $TESTDIR/reference/no-het-atoms/7q9b_IJHFG.pdb
