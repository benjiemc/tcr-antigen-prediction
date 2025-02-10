Test app centre-pad one-hot.
  $ python -m tcr_antigen_prediction.data.apps.process_sequence_data \
  > --seed 123 \
  > -o processed_sequences_centre_pad_one_hot.h5 \
  > $TESTDIR/data/sequences.csv

  $ python <<EOF
  > import os
  > import h5py
  > import numpy as np
  > fh = h5py.File('processed_sequences_centre_pad_one_hot.h5')
  > ref_cdr1a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'centre-pad-one-hot', 'cdr_1a.txt'))
  > ref_cdr2a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'centre-pad-one-hot', 'cdr_2a.txt'))
  > ref_cdr3a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'centre-pad-one-hot', 'cdr_3a.txt'))
  > ref_cdr1b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'centre-pad-one-hot', 'cdr_1b.txt'))
  > ref_cdr2b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'centre-pad-one-hot', 'cdr_2b.txt'))
  > ref_cdr3b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'centre-pad-one-hot', 'cdr_3b.txt'))
  > ref_peptide = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'centre-pad-one-hot', 'peptide.txt'))
  > ref_mhc = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'centre-pad-one-hot', 'mhc_pseudo_sequence.txt'))
  > ref_label = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'centre-pad-one-hot', 'label.txt'))
  > ref_fold = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'centre-pad-one-hot', 'fold.txt'))
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

Test right-pad blosum-50.
  $ python -m tcr_antigen_prediction.data.apps.process_sequence_data \
  > --seed 123 \
  > --cdr-1a-length 7 \
  > --cdr-2a-length 8 \
  > --cdr-3a-length 22 \
  > --cdr-1b-length 6 \
  > --cdr-2b-length 7 \
  > --cdr-3b-length 23 \
  > --peptide-length 12 \
  > --pad-direction right \
  > --encoding blosum50 \
  > -o processed_sequences_right_pad_blosum.h5 \
  > $TESTDIR/data/sequences.csv

  $ python <<EOF
  > import os
  > import h5py
  > import numpy as np
  > fh = h5py.File('processed_sequences_right_pad_blosum.h5')
  > ref_cdr1a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50', 'cdr_1a.txt'))
  > ref_cdr2a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50', 'cdr_2a.txt'))
  > ref_cdr3a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50', 'cdr_3a.txt'))
  > ref_cdr1b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50', 'cdr_1b.txt'))
  > ref_cdr2b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50', 'cdr_2b.txt'))
  > ref_cdr3b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50', 'cdr_3b.txt'))
  > ref_peptide = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50', 'peptide.txt'))
  > ref_mhc = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50', 'mhc_pseudo_sequence.txt'))
  > ref_label = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50', 'label.txt'))
  > ref_fold = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50', 'fold.txt'))
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

Test right-pad blosum-50 normalised.
  $ python -m tcr_antigen_prediction.data.apps.process_sequence_data \
  > --seed 123 \
  > --cdr-1a-length 7 \
  > --cdr-2a-length 8 \
  > --cdr-3a-length 22 \
  > --cdr-1b-length 6 \
  > --cdr-2b-length 7 \
  > --cdr-3b-length 23 \
  > --peptide-length 12 \
  > --pad-direction right \
  > --encoding blosum50 \
  > --normalisation-factor 5.0 \
  > -o processed_sequences_right_pad_blosum_norm.h5 \
  > $TESTDIR/data/sequences.csv

  $ python <<EOF
  > import os
  > import h5py
  > import numpy as np
  > fh = h5py.File('processed_sequences_right_pad_blosum_norm.h5')
  > ref_cdr1a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50-norm', 'cdr_1a.txt'))
  > ref_cdr2a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50-norm', 'cdr_2a.txt'))
  > ref_cdr3a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50-norm', 'cdr_3a.txt'))
  > ref_cdr1b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50-norm', 'cdr_1b.txt'))
  > ref_cdr2b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50-norm', 'cdr_2b.txt'))
  > ref_cdr3b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50-norm', 'cdr_3b.txt'))
  > ref_peptide = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50-norm', 'peptide.txt'))
  > ref_mhc = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50-norm', 'mhc_pseudo_sequence.txt'))
  > ref_label = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50-norm', 'label.txt'))
  > ref_fold = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50-norm', 'fold.txt'))
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
