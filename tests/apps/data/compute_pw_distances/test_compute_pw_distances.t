  $ python -m tcr_antigen_prediction.data.apps.compute_pw_distances -o distances.txt $(cat $TESTDIR/data/sequences.txt | tr '\n' ' ')
  $ diff distances.txt $TESTDIR/reference/distances.txt