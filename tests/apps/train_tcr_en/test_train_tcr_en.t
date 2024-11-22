Test training
  $ python -m tcr_antigen_prediction.apps.train_tcr_en \
  > -o TCRen_probabilities.csv \
  > --summary-csv "$TESTDIR/data/stcrdab_split.csv" \
  > $(cat "$TESTDIR/data/stcrdab_split.csv" | grep "train" | awk -F, -v dir="$TESTDIR/data" '{ printf "%s/%s_%s%s%s%s%s.pdb ", dir, $1, $2, $3, $4, $5, $6 }')

  $ cut -d, -f1-2 TCRen_probabilities.csv > test_entries
  $ cut -d, -f1-2 $TESTDIR/reference/TCRen_probabilities.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f3 TCRen_probabilities.csv | sed 1d > test_values
  $ cut -d, -f3 $TESTDIR/reference/TCRen_probabilities.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

Test training with validation
  $ python -m tcr_antigen_prediction.apps.train_tcr_en \
  > -o TCRen_probabilities.csv \
  > --summary-csv "$TESTDIR/data/stcrdab_split.csv" \
  > $(cat "$TESTDIR/data/stcrdab_split.csv" | grep "train" | awk -F, -v dir="$TESTDIR/data" '{ printf "%s/%s_%s%s%s%s%s.pdb ", dir, $1, $2, $3, $4, $5, $6 }') \
  > --validation-data $(cat "$TESTDIR/data/stcrdab_split.csv" | grep "validation" | awk -F, -v dir="$TESTDIR/data" '{ printf "%s/%s_%s%s%s%s%s.pdb ", dir, $1, $2, $3, $4, $5, $6 }')

  $ cut -d, -f1-2 TCRen_probabilities.csv > test_entries
  $ cut -d, -f1-2 $TESTDIR/reference/TCRen_probabilities.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f3 TCRen_probabilities.csv | sed 1d > test_values
  $ cut -d, -f3 $TESTDIR/reference/TCRen_probabilities.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

Test LOO training
  $ python -m tcr_antigen_prediction.apps.train_tcr_en \
  > --strategy leave-one-out \
  > -o TCRen_probabilities.csv \
  > --summary-csv "$TESTDIR/data/stcrdab_split.csv" \
  > $(cat "$TESTDIR/data/stcrdab_split.csv" | grep "train" | awk -F, -v dir="$TESTDIR/data" '{ printf "%s/%s_%s%s%s%s%s.pdb ", dir, $1, $2, $3, $4, $5, $6 }')

  $ cut -d, -f1-2 TCRen_probabilities_LOO_1.csv > test_entries
  $ cut -d, -f1-2 $TESTDIR/reference/TCRen_probabilities_LOO_1.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f3 TCRen_probabilities_LOO_1.csv | sed 1d > test_values
  $ cut -d, -f3 $TESTDIR/reference/TCRen_probabilities_LOO_1.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

  $ cut -d, -f1-2 TCRen_probabilities_LOO_2.csv > test_entries
  $ cut -d, -f1-2 $TESTDIR/reference/TCRen_probabilities_LOO_2.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f3 TCRen_probabilities_LOO_2.csv | sed 1d > test_values
  $ cut -d, -f3 $TESTDIR/reference/TCRen_probabilities_LOO_2.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"
