Test app.
  $ python -m tcr_antigen_prediction.data.apps.annotate_mhc_pseudo_sequences \
  > --log-level warning \
  > -o summary.csv \
  > --summary-csv "$TESTDIR/data/structure_summary.csv" \
  > --mhc-pseudo-sequence-imgt-numbers "$TESTDIR/data/mhc_pseudo_seq_imgt_positions.json" \
  > $TESTDIR/data/*.pdb
  * - WARNING: No path given for 3tf7_cCBA.pdb, defaulting to None (glob)

  $ diff summary.csv "$TESTDIR/reference/summary.csv"
