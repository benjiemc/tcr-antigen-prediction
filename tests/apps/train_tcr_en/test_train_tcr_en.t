Test training
  $ python -m tcr_antigen_prediction.apps.train_tcr_en \
  > -o TCRen_probabilities.csv \
  > --summary-csv "$TESTDIR/data/stcrdab_split.csv" \
  > $(cat "$TESTDIR/data/stcrdab_split.csv" | grep "train" | awk -F, -v dir="$TESTDIR/data" '{ printf "%s/%s_%s%s%s%s%s.pdb ", dir, $1, $2, $3, $4, $5, $6 }')

  $ diff TCRen_probabilities.csv $TESTDIR/reference/TCRen_probabilities.csv
