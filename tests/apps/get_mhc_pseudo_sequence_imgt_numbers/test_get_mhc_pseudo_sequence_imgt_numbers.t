Test app.
  $ python -m tcr_antigen_prediction.apps.get_mhc_pseudo_sequence_imgt_numbers \
  > --probability-filter 0.05 \
  > -o mhc_pseudo_sequence_imgt_numbers.json \
  > $TESTDIR/data/tcr_pmhc_contact_map.csv

  $ diff mhc_pseudo_sequence_imgt_numbers.json $TESTDIR/reference/mhc_pseudo_sequence_imgt_numbers.json
