Test 1ao7 - abTCR with MH1
  $ python -m tcr_antigen_prediction.apps.renumber_tcr_pmhc_structure -o 1ao7_renumbered.pdb $TESTDIR/data/1ao7.pdb
  $ diff 1ao7_renumbered.pdb $TESTDIR/reference/1ao7_renumbered.pdb

Test 1j8h - abTCR with MH2
  $ python -m tcr_antigen_prediction.apps.renumber_tcr_pmhc_structure -o 1j8h_renumbered.pdb $TESTDIR/data/1j8h.pdb
  $ diff 1j8h_renumbered.pdb $TESTDIR/reference/1j8h_renumbered.pdb