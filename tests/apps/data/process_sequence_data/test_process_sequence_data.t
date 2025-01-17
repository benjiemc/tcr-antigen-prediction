  $ python -m tcr_antigen_prediction.data.apps.process_sequence_data \
  > --seed 123 \
  > --mhc-sequences $TESTDIR/data/h2_d.json $TESTDIR/data/hla_a.json $TESTDIR/data/hla_b.json \
  > --mhc-pseudo-sequence-imgt-numbers $TESTDIR/data/mhc_pseudo_seq_imgt_positions.json \
  > -o processed_sequences.h5 \
  > $TESTDIR/data/sequences.csv

  $ python <<EOF
  > import os
  > import h5py
  > import numpy as np
  > fh = h5py.File('processed_sequences.h5')
  > ref_cdr1a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'cdr_1a.txt'))
  > ref_cdr2a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'cdr_2a.txt'))
  > ref_cdr3a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'cdr_3a.txt'))
  > ref_cdr1b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'cdr_1b.txt'))
  > ref_cdr2b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'cdr_2b.txt'))
  > ref_cdr3b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'cdr_3b.txt'))
  > ref_peptide = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'peptide.txt'))
  > ref_mhc = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'mhc_pseudo_sequence.txt'))
  > ref_label = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'label.txt'))
  > ref_fold = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'fold.txt'))
  > np.testing.assert_array_equal(ref_cdr1a, fh['cdr_1a'][:].reshape(fh['cdr_1a'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr2a, fh['cdr_2a'][:].reshape(fh['cdr_2a'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr3a, fh['cdr_3a'][:].reshape(fh['cdr_3a'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr1b, fh['cdr_1b'][:].reshape(fh['cdr_1b'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr2b, fh['cdr_2b'][:].reshape(fh['cdr_2b'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr3b, fh['cdr_3b'][:].reshape(fh['cdr_3b'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_peptide, fh['peptide'][:].reshape(fh['peptide'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_mhc, fh['mhc_pseudo_sequence'][:].reshape(fh['mhc_pseudo_sequence'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_label, fh['label'][:])
  > np.testing.assert_array_equal(ref_fold, fh['fold'][:])
  > fh.close()
  > EOF
