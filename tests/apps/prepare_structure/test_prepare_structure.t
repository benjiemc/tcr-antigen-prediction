Creating a mesh for a TCR
  $ python -m tcr_antigen_prediction.apps.prepare_structure \
  > -o ./ \
  > $TESTDIR/data/1ao7_DECAB.pdb \
  > --chains D E

  $ diff 1ao7_DECAB_DE.ply $TESTDIR/reference/1ao7_DECAB_DE.ply

Creating a mesh for a pMHC
  $ python -m tcr_antigen_prediction.apps.prepare_structure \
  > -o ./ \
  > $TESTDIR/data/1ao7_DECAB.pdb \
  > --chains C A

  $ diff 1ao7_DECAB_CA.ply $TESTDIR/reference/1ao7_DECAB_CA.ply