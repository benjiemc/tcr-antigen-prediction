Test MH1 contact maps with padding.
  $ python -m tcr_antigen_prediction.data.apps.process_contact_maps \
  > --mhc-types MH1 \
  > --pad \
  > --mhc-pseudo-sequence-imgt-numbers $TESTDIR/data/mhc_pseudo_seq_imgt_positions.json \
  > -o tcr_mhc1_pad_contact_maps.h5 \
  > $TESTDIR/data/tcr_pmhc_contacts.csv

  $ python <<EOF
  > import os
  > import h5py
  > import numpy as np
  > fh = h5py.File('tcr_mhc1_pad_contact_maps.h5')
  > ref_peptide_cdr1_alpha = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'mhc1_pad', 'peptide_cdr1_alpha.txt'))
  > ref_peptide_cdr2_alpha = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'mhc1_pad', 'peptide_cdr2_alpha.txt'))
  > ref_peptide_cdr3_alpha = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'mhc1_pad', 'peptide_cdr3_alpha.txt'))
  > ref_peptide_cdr1_beta = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'mhc1_pad', 'peptide_cdr1_beta.txt'))
  > ref_peptide_cdr2_beta = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'mhc1_pad', 'peptide_cdr2_beta.txt'))
  > ref_peptide_cdr3_beta = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'mhc1_pad', 'peptide_cdr3_beta.txt'))
  > ref_mhc_pseudo_cdr1_alpha = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'mhc1_pad', 'mhc_pseudo_cdr1_alpha.txt'))
  > ref_mhc_pseudo_cdr2_alpha = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'mhc1_pad', 'mhc_pseudo_cdr2_alpha.txt'))
  > ref_mhc_pseudo_cdr3_alpha = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'mhc1_pad', 'mhc_pseudo_cdr3_alpha.txt'))
  > ref_mhc_pseudo_cdr1_beta = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'mhc1_pad', 'mhc_pseudo_cdr1_beta.txt'))
  > ref_mhc_pseudo_cdr2_beta = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'mhc1_pad', 'mhc_pseudo_cdr2_beta.txt'))
  > ref_mhc_pseudo_cdr3_beta = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'mhc1_pad', 'mhc_pseudo_cdr3_beta.txt'))
  > np.testing.assert_array_almost_equal(ref_peptide_cdr1_alpha, fh['peptide']['cdr1_alpha'][:])
  > np.testing.assert_array_almost_equal(ref_peptide_cdr2_alpha, fh['peptide']['cdr2_alpha'][:])
  > np.testing.assert_array_almost_equal(ref_peptide_cdr3_alpha, fh['peptide']['cdr3_alpha'][:])
  > np.testing.assert_array_almost_equal(ref_peptide_cdr1_beta, fh['peptide']['cdr1_beta'][:])
  > np.testing.assert_array_almost_equal(ref_peptide_cdr2_beta, fh['peptide']['cdr2_beta'][:])
  > np.testing.assert_array_almost_equal(ref_peptide_cdr3_beta, fh['peptide']['cdr3_beta'][:])
  > np.testing.assert_array_almost_equal(ref_mhc_pseudo_cdr1_alpha, fh['mhc_pseudo']['cdr1_alpha'][:])
  > np.testing.assert_array_almost_equal(ref_mhc_pseudo_cdr2_alpha, fh['mhc_pseudo']['cdr2_alpha'][:])
  > np.testing.assert_array_almost_equal(ref_mhc_pseudo_cdr3_alpha, fh['mhc_pseudo']['cdr3_alpha'][:])
  > np.testing.assert_array_almost_equal(ref_mhc_pseudo_cdr1_beta, fh['mhc_pseudo']['cdr1_beta'][:])
  > np.testing.assert_array_almost_equal(ref_mhc_pseudo_cdr2_beta, fh['mhc_pseudo']['cdr2_beta'][:])
  > np.testing.assert_array_almost_equal(ref_mhc_pseudo_cdr3_beta, fh['mhc_pseudo']['cdr3_beta'][:])
  > fh.close()
  > EOF

Test combined MHC with padding.
  $ python -m tcr_antigen_prediction.data.apps.process_contact_maps \
  > --mhc-types MH1 MH2 \
  > --pad \
  > --mhc-pseudo-sequence-imgt-numbers $TESTDIR/data/mhc_pseudo_seq_imgt_positions.json \
  > -o combined_pad_contact_maps.h5 \
  > $TESTDIR/data/tcr_pmhc_contacts.csv

  $ python <<EOF
  > import os
  > import h5py
  > import numpy as np
  > fh = h5py.File('combined_pad_contact_maps.h5')
  > ref_peptide_cdr1_alpha = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined_pad', 'peptide_cdr1_alpha.txt'))
  > ref_peptide_cdr2_alpha = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined_pad', 'peptide_cdr2_alpha.txt'))
  > ref_peptide_cdr3_alpha = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined_pad', 'peptide_cdr3_alpha.txt'))
  > ref_peptide_cdr1_beta = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined_pad', 'peptide_cdr1_beta.txt'))
  > ref_peptide_cdr2_beta = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined_pad', 'peptide_cdr2_beta.txt'))
  > ref_peptide_cdr3_beta = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined_pad', 'peptide_cdr3_beta.txt'))
  > ref_mhc_pseudo_cdr1_alpha = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined_pad', 'mhc_pseudo_cdr1_alpha.txt'))
  > ref_mhc_pseudo_cdr2_alpha = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined_pad', 'mhc_pseudo_cdr2_alpha.txt'))
  > ref_mhc_pseudo_cdr3_alpha = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined_pad', 'mhc_pseudo_cdr3_alpha.txt'))
  > ref_mhc_pseudo_cdr1_beta = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined_pad', 'mhc_pseudo_cdr1_beta.txt'))
  > ref_mhc_pseudo_cdr2_beta = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined_pad', 'mhc_pseudo_cdr2_beta.txt'))
  > ref_mhc_pseudo_cdr3_beta = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined_pad', 'mhc_pseudo_cdr3_beta.txt'))
  > np.testing.assert_array_almost_equal(ref_peptide_cdr1_alpha, fh['peptide']['cdr1_alpha'][:])
  > np.testing.assert_array_almost_equal(ref_peptide_cdr2_alpha, fh['peptide']['cdr2_alpha'][:])
  > np.testing.assert_array_almost_equal(ref_peptide_cdr3_alpha, fh['peptide']['cdr3_alpha'][:])
  > np.testing.assert_array_almost_equal(ref_peptide_cdr1_beta, fh['peptide']['cdr1_beta'][:])
  > np.testing.assert_array_almost_equal(ref_peptide_cdr2_beta, fh['peptide']['cdr2_beta'][:])
  > np.testing.assert_array_almost_equal(ref_peptide_cdr3_beta, fh['peptide']['cdr3_beta'][:])
  > np.testing.assert_array_almost_equal(ref_mhc_pseudo_cdr1_alpha, fh['mhc_pseudo']['cdr1_alpha'][:])
  > np.testing.assert_array_almost_equal(ref_mhc_pseudo_cdr2_alpha, fh['mhc_pseudo']['cdr2_alpha'][:])
  > np.testing.assert_array_almost_equal(ref_mhc_pseudo_cdr3_alpha, fh['mhc_pseudo']['cdr3_alpha'][:])
  > np.testing.assert_array_almost_equal(ref_mhc_pseudo_cdr1_beta, fh['mhc_pseudo']['cdr1_beta'][:])
  > np.testing.assert_array_almost_equal(ref_mhc_pseudo_cdr2_beta, fh['mhc_pseudo']['cdr2_beta'][:])
  > np.testing.assert_array_almost_equal(ref_mhc_pseudo_cdr3_beta, fh['mhc_pseudo']['cdr3_beta'][:])
  > fh.close()
  > EOF

Test combined without padding.
  $ python -m tcr_antigen_prediction.data.apps.process_contact_maps \
  > --mhc-types MH1 MH2 \
  > --mhc-pseudo-sequence-imgt-numbers $TESTDIR/data/mhc_pseudo_seq_imgt_positions.json \
  > -o combined_contact_maps.h5 \
  > $TESTDIR/data/tcr_pmhc_contacts.csv

  $ python <<EOF
  > import os
  > import h5py
  > import numpy as np
  > fh = h5py.File('combined_contact_maps.h5')
  > ref_peptide_cdr1_alpha = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined', 'peptide_cdr1_alpha.txt'))
  > ref_peptide_cdr2_alpha = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined', 'peptide_cdr2_alpha.txt'))
  > ref_peptide_cdr3_alpha = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined', 'peptide_cdr3_alpha.txt'))
  > ref_peptide_cdr1_beta = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined', 'peptide_cdr1_beta.txt'))
  > ref_peptide_cdr2_beta = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined', 'peptide_cdr2_beta.txt'))
  > ref_peptide_cdr3_beta = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined', 'peptide_cdr3_beta.txt'))
  > ref_mhc_pseudo_cdr1_alpha = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined', 'mhc_pseudo_cdr1_alpha.txt'))
  > ref_mhc_pseudo_cdr2_alpha = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined', 'mhc_pseudo_cdr2_alpha.txt'))
  > ref_mhc_pseudo_cdr3_alpha = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined', 'mhc_pseudo_cdr3_alpha.txt'))
  > ref_mhc_pseudo_cdr1_beta = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined', 'mhc_pseudo_cdr1_beta.txt'))
  > ref_mhc_pseudo_cdr2_beta = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined', 'mhc_pseudo_cdr2_beta.txt'))
  > ref_mhc_pseudo_cdr3_beta = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'combined', 'mhc_pseudo_cdr3_beta.txt'))
  > np.testing.assert_array_almost_equal(ref_peptide_cdr1_alpha, fh['peptide']['cdr1_alpha'][:])
  > np.testing.assert_array_almost_equal(ref_peptide_cdr2_alpha, fh['peptide']['cdr2_alpha'][:])
  > np.testing.assert_array_almost_equal(ref_peptide_cdr3_alpha, fh['peptide']['cdr3_alpha'][:])
  > np.testing.assert_array_almost_equal(ref_peptide_cdr1_beta, fh['peptide']['cdr1_beta'][:])
  > np.testing.assert_array_almost_equal(ref_peptide_cdr2_beta, fh['peptide']['cdr2_beta'][:])
  > np.testing.assert_array_almost_equal(ref_peptide_cdr3_beta, fh['peptide']['cdr3_beta'][:])
  > np.testing.assert_array_almost_equal(ref_mhc_pseudo_cdr1_alpha, fh['mhc_pseudo']['cdr1_alpha'][:])
  > np.testing.assert_array_almost_equal(ref_mhc_pseudo_cdr2_alpha, fh['mhc_pseudo']['cdr2_alpha'][:])
  > np.testing.assert_array_almost_equal(ref_mhc_pseudo_cdr3_alpha, fh['mhc_pseudo']['cdr3_alpha'][:])
  > np.testing.assert_array_almost_equal(ref_mhc_pseudo_cdr1_beta, fh['mhc_pseudo']['cdr1_beta'][:])
  > np.testing.assert_array_almost_equal(ref_mhc_pseudo_cdr2_beta, fh['mhc_pseudo']['cdr2_beta'][:])
  > np.testing.assert_array_almost_equal(ref_mhc_pseudo_cdr3_beta, fh['mhc_pseudo']['cdr3_beta'][:])
  > fh.close()
  > EOF