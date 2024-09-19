Test multiple TCR:pMHC-Is
  $ python -m tcr_antigen_prediction.apps.identify_tcr_pmhc_interactions "$TESTDIR/data/7pbe.pdb"
  Achain\tBchain\tantigen_chain\tmhc_chain1\tmhc_chain2\tmhc_type (esc)
  D\tE\tC\tA\tB\tMH1 (esc)
  I\tJ\tH\tF\tG\tMH1 (esc)

  $ python -m tcr_antigen_prediction.apps.identify_tcr_pmhc_interactions "$TESTDIR/data/4e41.pdb"
  Achain\tBchain\tantigen_chain\tmhc_chain1\tmhc_chain2\tmhc_type (esc)
  D\tE\tC\tA\tB\tMH2 (esc)
  I\tJ\tH\tF\tG\tMH2 (esc)
