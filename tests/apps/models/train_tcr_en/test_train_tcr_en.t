Test training
  $ python -m tcr_antigen_prediction.models.apps.train_tcr_en \
  > -o TCRen_probabilities.csv \
  > --summary-csv "$TESTDIR/data/stcrdab_split.csv" \
  > $TESTDIR/data/3qdj_DECAB.pdb $TESTDIR/data/6eqa_DECAB.pdb $TESTDIR/data/6eqb_DECAB.pdb

  $ cut -d, -f1-2 TCRen_probabilities.csv > test_entries
  $ cut -d, -f1-2 $TESTDIR/reference/TCRen_probabilities.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f3 TCRen_probabilities.csv | sed 1d > test_values
  $ cut -d, -f3 $TESTDIR/reference/TCRen_probabilities.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

Test training with validation
  $ python -m tcr_antigen_prediction.models.apps.train_tcr_en \
  > -o TCRen_probabilities.csv \
  > --summary-csv "$TESTDIR/data/stcrdab_split.csv" \
  > $TESTDIR/data/3qdj_DECAB.pdb $TESTDIR/data/6eqa_DECAB.pdb $TESTDIR/data/6eqb_DECAB.pdb \
  > --validation-data $TESTDIR/data/2pxy_ABPCD.pdb

  $ cut -d, -f1-2 TCRen_probabilities.csv > test_entries
  $ cut -d, -f1-2 $TESTDIR/reference/TCRen_probabilities.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f3 TCRen_probabilities.csv | sed 1d > test_values
  $ cut -d, -f3 $TESTDIR/reference/TCRen_probabilities.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

Test LOO training
  $ python -m tcr_antigen_prediction.models.apps.train_tcr_en \
  > --strategy leave-one-out \
  > -o TCRen_probabilities.csv \
  > --summary-csv "$TESTDIR/data/stcrdab_split.csv" \
  > $TESTDIR/data/3qdj_DECAB.pdb $TESTDIR/data/6eqa_DECAB.pdb $TESTDIR/data/6eqb_DECAB.pdb $TESTDIR/data/2pxy_ABPCD.pdb

  $ cut -d, -f1-2 TCRen_probabilities_LOO_AAGIGILTV.csv > test_entries
  $ cut -d, -f1-2 $TESTDIR/reference/TCRen_probabilities_LOO_AAGIGILTV.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f3 TCRen_probabilities_LOO_AAGIGILTV.csv | sed 1d > test_values
  $ cut -d, -f3 $TESTDIR/reference/TCRen_probabilities_LOO_AAGIGILTV.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

  $ cut -d, -f1-2 TCRen_probabilities_LOO_RGGASQYRPSQ.csv > test_entries
  $ cut -d, -f1-2 $TESTDIR/reference/TCRen_probabilities_LOO_RGGASQYRPSQ.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f3 TCRen_probabilities_LOO_RGGASQYRPSQ.csv | sed 1d > test_values
  $ cut -d, -f3 $TESTDIR/reference/TCRen_probabilities_LOO_RGGASQYRPSQ.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"
