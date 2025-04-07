Test counts.
  $ python -m tcr_antigen_prediction.data.apps.count_ppi_interactions \
  >  -o counts.csv \
  > $TESTDIR/data/test_ppi3d.csv

  $ diff counts.csv $TESTDIR/reference/counts.csv

Test normalising counts.
  $ python -m tcr_antigen_prediction.data.apps.count_ppi_interactions \
  > --normalise \
  > -o normalised.csv \
  > $TESTDIR/data/test_ppi3d.csv

  $ cut -d, -f2 normalised.csv > test_entries
  $ cut -d, -f2 $TESTDIR/reference/normalised.csv > ref_entries
  $ diff test_entries ref_entries

  $ cut -d, -f3 normalised.csv | sed 1d > test_vals
  $ cut -d, -f3 $TESTDIR/reference/normalised.csv | sed 1d > ref_vals
  $ python -c "import numpy as np; test = np.loadtxt('test_vals'); ref = np.loadtxt('ref_vals'); np.testing.assert_array_almost_equal(test, ref)"

Test sampling input.
  $ python -m tcr_antigen_prediction.data.apps.count_ppi_interactions \
  > --seed 123 \
  > --sample-size 3 \
  > -o sampled.csv \
  > $TESTDIR/data/test_ppi3d.csv

  $ diff sampled.csv $TESTDIR/reference/sampled.csv

Test separating interaction types in output.
  $ python -m tcr_antigen_prediction.data.apps.count_ppi_interactions \
  > --separate-interaction-types \
  > -o interaction_types.csv \
  > $TESTDIR/data/test_ppi3d.csv

  $ diff interaction_types.csv $TESTDIR/reference/interaction_types.csv
