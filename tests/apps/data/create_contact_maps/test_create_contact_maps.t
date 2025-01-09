Test relative numbering
  $ python -m tcr_antigen_prediction.data.apps.create_contact_maps \
  > --log-level warning \
  > --tcr-norm relative_pos_centre \
  > --peptide-norm relative_pos_centre \
  > --mhc-norm imgt_number \
  > -o relative_numbering.csv \
  > --summary-csv $TESTDIR/data/stcrdab_split.csv \
  > $TESTDIR/data/*.pdb

  $ cut -d, -f 1-6 relative_numbering.csv > test_entries
  $ cut -d, -f 1-6 $TESTDIR/reference/test_relative_numbering.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f 7 relative_numbering.csv | sed 1d > test_values
  $ cut -d, -f 7 $TESTDIR/reference/test_relative_numbering.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

Test right offset on relative numbering
  $ python -m tcr_antigen_prediction.data.apps.create_contact_maps \
  > --log-level warning \
  > --tcr-norm relative_pos_centre \
  > --peptide-norm relative_pos_centre \
  > --mhc-norm imgt_number \
  > --even-offset-side right \
  > -o relative_numbering_right_offset.csv \
  > --summary-csv $TESTDIR/data/stcrdab_split.csv \
  > $TESTDIR/data/*.pdb

  $ cut -d, -f 1-6 relative_numbering_right_offset.csv > test_entries
  $ cut -d, -f 1-6 $TESTDIR/reference/test_relative_numbering_right_offset.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f 7 relative_numbering_right_offset.csv | sed 1d > test_values
  $ cut -d, -f 7 $TESTDIR/reference/test_relative_numbering_right_offset.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

Test filtering
  $ python -m tcr_antigen_prediction.data.apps.create_contact_maps \
  > --log-level warning \
  > --tcr-norm \
  > --peptide-norm \
  > --mhc-norm imgt_number \
  > --probability-filter 0.01 \
  > -o mhc_resi_filtered.csv \
  > --summary-csv $TESTDIR/data/stcrdab_split.csv \
  > $TESTDIR/data/*.pdb

  $ cut -d, -f 1-3 mhc_resi_filtered.csv > test_entries
  $ cut -d, -f 1-3 $TESTDIR/reference/test_mhc_resi_filtered.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f 4 mhc_resi_filtered.csv | sed 1d > test_values
  $ cut -d, -f 4 $TESTDIR/reference/test_mhc_resi_filtered.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"


Test blank
  $ python -m tcr_antigen_prediction.data.apps.create_contact_maps \
  > --log-level warning \
  > --tcr-norm \
  > --peptide-norm \
  > --mhc-norm \
  > -o blank.csv \
  > --summary-csv $TESTDIR/data/stcrdab_split.csv \
  > $TESTDIR/data/*.pdb

  $ diff blank.csv $TESTDIR/reference/test_blank.csv