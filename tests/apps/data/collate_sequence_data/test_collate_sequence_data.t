Test app.
  $ python -m tcr_antigen_prediction.data.apps.collate_sequence_data \
  > --iedb-path $TESTDIR/data/iedb.csv \
  > --vdjdb-path $TESTDIR/data/vdjdb.tsv \
  > --itrap-path $TESTDIR/data/itrap.csv \
  > --mcpas-tcr-path $TESTDIR/data/mcpas_tcr.csv \
  > -o sequences.csv

  $ diff sequences.csv $TESTDIR/reference/sequences.csv
