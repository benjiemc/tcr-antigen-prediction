Test 1oa7 chains D E C A B
  $ python -m tcr_antigen_prediction.data.apps.extract_chains_from_structure \
  > --chains D E C A B \
  > -o 1ao7_DECAB.pdb \
  > $TESTDIR/data/1ao7.pdb

  $ diff 1ao7_DECAB.pdb $TESTDIR/reference/1ao7_DECAB.pdb
