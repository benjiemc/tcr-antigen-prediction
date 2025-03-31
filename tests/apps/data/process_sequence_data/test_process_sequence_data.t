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
  > ref_cdr1a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'centre-pad-one-hot', 'cdr1_alpha.txt'))
  > ref_cdr2a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'centre-pad-one-hot', 'cdr2_alpha.txt'))
  > ref_cdr3a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'centre-pad-one-hot', 'cdr3_alpha.txt'))
  > ref_cdr1b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'centre-pad-one-hot', 'cdr1_beta.txt'))
  > ref_cdr2b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'centre-pad-one-hot', 'cdr2_beta.txt'))
  > ref_cdr3b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'centre-pad-one-hot', 'cdr3_beta.txt'))
  > ref_peptide = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'centre-pad-one-hot', 'peptide.txt'))
  > ref_mhc = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'centre-pad-one-hot', 'mhc_pseudo.txt'))
  > ref_label = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'centre-pad-one-hot', 'label.txt'))
  > ref_fold = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'centre-pad-one-hot', 'fold.txt'))
  > np.testing.assert_array_equal(ref_cdr1a, fh['cdr1_alpha'][:].reshape(fh['cdr1_alpha'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr2a, fh['cdr2_alpha'][:].reshape(fh['cdr2_alpha'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr3a, fh['cdr3_alpha'][:].reshape(fh['cdr3_alpha'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr1b, fh['cdr1_beta'][:].reshape(fh['cdr1_beta'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr2b, fh['cdr2_beta'][:].reshape(fh['cdr2_beta'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr3b, fh['cdr3_beta'][:].reshape(fh['cdr3_beta'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_peptide, fh['peptide'][:].reshape(fh['peptide'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_mhc, fh['mhc_pseudo'][:].reshape(fh['mhc_pseudo'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_label, fh['label'][:])
  > np.testing.assert_array_equal(ref_fold, fh['fold'][:])
  > fh.close()
  > EOF

Test right-pad blosum-50.
  $ python -m tcr_antigen_prediction.data.apps.process_sequence_data \
  > --seed 123 \
  > --cdr1-alpha-length 7 \
  > --cdr2-alpha-length 8 \
  > --cdr3-alpha-length 22 \
  > --cdr1-beta-length 6 \
  > --cdr2-beta-length 7 \
  > --cdr3-beta-length 23 \
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
  > ref_cdr1a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50', 'cdr1_alpha.txt'))
  > ref_cdr2a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50', 'cdr2_alpha.txt'))
  > ref_cdr3a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50', 'cdr3_alpha.txt'))
  > ref_cdr1b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50', 'cdr1_beta.txt'))
  > ref_cdr2b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50', 'cdr2_beta.txt'))
  > ref_cdr3b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50', 'cdr3_beta.txt'))
  > ref_peptide = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50', 'peptide.txt'))
  > ref_mhc = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50', 'mhc_pseudo.txt'))
  > ref_label = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50', 'label.txt'))
  > ref_fold = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50', 'fold.txt'))
  > np.testing.assert_array_equal(ref_cdr1a, fh['cdr1_alpha'][:].reshape(fh['cdr1_alpha'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr2a, fh['cdr2_alpha'][:].reshape(fh['cdr2_alpha'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr3a, fh['cdr3_alpha'][:].reshape(fh['cdr3_alpha'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr1b, fh['cdr1_beta'][:].reshape(fh['cdr1_beta'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr2b, fh['cdr2_beta'][:].reshape(fh['cdr2_beta'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr3b, fh['cdr3_beta'][:].reshape(fh['cdr3_beta'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_peptide, fh['peptide'][:].reshape(fh['peptide'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_mhc, fh['mhc_pseudo'][:].reshape(fh['mhc_pseudo'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_label, fh['label'][:])
  > np.testing.assert_array_equal(ref_fold, fh['fold'][:])
  > fh.close()
  > EOF

Test right-pad blosum-50 normalised.
  $ python -m tcr_antigen_prediction.data.apps.process_sequence_data \
  > --seed 123 \
  > --cdr1-alpha-length 7 \
  > --cdr2-alpha-length 8 \
  > --cdr3-alpha-length 22 \
  > --cdr1-beta-length 6 \
  > --cdr2-beta-length 7 \
  > --cdr3-beta-length 23 \
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
  > ref_cdr1a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50-norm', 'cdr1_alpha.txt'))
  > ref_cdr2a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50-norm', 'cdr2_alpha.txt'))
  > ref_cdr3a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50-norm', 'cdr3_alpha.txt'))
  > ref_cdr1b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50-norm', 'cdr1_beta.txt'))
  > ref_cdr2b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50-norm', 'cdr2_beta.txt'))
  > ref_cdr3b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50-norm', 'cdr3_beta.txt'))
  > ref_peptide = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50-norm', 'peptide.txt'))
  > ref_mhc = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50-norm', 'mhc_pseudo.txt'))
  > ref_label = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50-norm', 'label.txt'))
  > ref_fold = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'right-pad-blosum-50-norm', 'fold.txt'))
  > np.testing.assert_array_equal(ref_cdr1a, fh['cdr1_alpha'][:].reshape(fh['cdr1_alpha'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr2a, fh['cdr2_alpha'][:].reshape(fh['cdr2_alpha'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr3a, fh['cdr3_alpha'][:].reshape(fh['cdr3_alpha'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr1b, fh['cdr1_beta'][:].reshape(fh['cdr1_beta'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr2b, fh['cdr2_beta'][:].reshape(fh['cdr2_beta'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr3b, fh['cdr3_beta'][:].reshape(fh['cdr3_beta'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_peptide, fh['peptide'][:].reshape(fh['peptide'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_mhc, fh['mhc_pseudo'][:].reshape(fh['mhc_pseudo'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_label, fh['label'][:])
  > np.testing.assert_array_equal(ref_fold, fh['fold'][:])
  > fh.close()
  > EOF

Test splitting on pMHCs
  $ python -m tcr_antigen_prediction.data.apps.process_sequence_data \
  > --seed 123 \
  > --split-type pMHC \
  > -o processed_sequences_pmhc_split.h5 \
  > $TESTDIR/data/sequences.csv
  * - WARNING: Insufficient data to create 5 folds. Only 3 fold(s) will be created (glob)

  $ python <<EOF
  > import os
  > import h5py
  > import numpy as np
  > fh = h5py.File('processed_sequences_pmhc_split.h5')
  > ref_cdr1a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'pmhc-split', 'cdr1_alpha.txt'))
  > ref_cdr2a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'pmhc-split', 'cdr2_alpha.txt'))
  > ref_cdr3a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'pmhc-split', 'cdr3_alpha.txt'))
  > ref_cdr1b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'pmhc-split', 'cdr1_beta.txt'))
  > ref_cdr2b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'pmhc-split', 'cdr2_beta.txt'))
  > ref_cdr3b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'pmhc-split', 'cdr3_beta.txt'))
  > ref_peptide = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'pmhc-split', 'peptide.txt'))
  > ref_mhc = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'pmhc-split', 'mhc_pseudo.txt'))
  > ref_label = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'pmhc-split', 'label.txt'))
  > ref_fold = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'pmhc-split', 'fold.txt'))
  > np.testing.assert_array_equal(ref_cdr1a, fh['cdr1_alpha'][:].reshape(fh['cdr1_alpha'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr2a, fh['cdr2_alpha'][:].reshape(fh['cdr2_alpha'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr3a, fh['cdr3_alpha'][:].reshape(fh['cdr3_alpha'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr1b, fh['cdr1_beta'][:].reshape(fh['cdr1_beta'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr2b, fh['cdr2_beta'][:].reshape(fh['cdr2_beta'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr3b, fh['cdr3_beta'][:].reshape(fh['cdr3_beta'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_peptide, fh['peptide'][:].reshape(fh['peptide'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_mhc, fh['mhc_pseudo'][:].reshape(fh['mhc_pseudo'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_label, fh['label'][:])
  > np.testing.assert_array_equal(ref_fold, fh['fold'][:])
  > fh.close()
  > EOF

Test splitting on TCRs
  $ python -m tcr_antigen_prediction.data.apps.process_sequence_data \
  > --seed 123 \
  > --split-type tcr \
  > -o processed_sequences_tcr_split.h5 \
  > $TESTDIR/data/sequences.csv

  $ python <<EOF
  > import os
  > import h5py
  > import numpy as np
  > fh = h5py.File('processed_sequences_tcr_split.h5')
  > ref_cdr1a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'tcr-split', 'cdr1_alpha.txt'))
  > ref_cdr2a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'tcr-split', 'cdr2_alpha.txt'))
  > ref_cdr3a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'tcr-split', 'cdr3_alpha.txt'))
  > ref_cdr1b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'tcr-split', 'cdr1_beta.txt'))
  > ref_cdr2b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'tcr-split', 'cdr2_beta.txt'))
  > ref_cdr3b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'tcr-split', 'cdr3_beta.txt'))
  > ref_peptide = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'tcr-split', 'peptide.txt'))
  > ref_mhc = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'tcr-split', 'mhc_pseudo.txt'))
  > ref_label = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'tcr-split', 'label.txt'))
  > ref_fold = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'tcr-split', 'fold.txt'))
  > np.testing.assert_array_equal(ref_cdr1a, fh['cdr1_alpha'][:].reshape(fh['cdr1_alpha'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr2a, fh['cdr2_alpha'][:].reshape(fh['cdr2_alpha'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr3a, fh['cdr3_alpha'][:].reshape(fh['cdr3_alpha'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr1b, fh['cdr1_beta'][:].reshape(fh['cdr1_beta'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr2b, fh['cdr2_beta'][:].reshape(fh['cdr2_beta'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr3b, fh['cdr3_beta'][:].reshape(fh['cdr3_beta'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_peptide, fh['peptide'][:].reshape(fh['peptide'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_mhc, fh['mhc_pseudo'][:].reshape(fh['mhc_pseudo'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_label, fh['label'][:])
  > np.testing.assert_array_equal(ref_fold, fh['fold'][:])
  > fh.close()
  > EOF

Test split on levenshtein distance
  $ python <<EOF
  > import os
  > import h5py
  > import numpy as np
  > with h5py.File('distances.h5', 'w') as hdf5:
  >     for entity in [
  >         'cdr1_alpha', 'cdr2_alpha', 'cdr3_alpha', 'cdr1_beta', 'cdr2_beta', 'cdr3_beta', 'peptide', 'mhc_pseudo'
  >     ]:
  >         group = hdf5.create_group(entity)
  >         with open(os.path.join(os.environ['TESTDIR'], 'data', f'{entity}s.txt')) as fh:
  >             group['names'] = [line.strip() for line in fh.readlines() if line]
  >         group['distance_matrix'] = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'data', f'{entity}_distances.txt'), dtype=int)
  > EOF

  $ python -m tcr_antigen_prediction.data.apps.process_sequence_data \
  > --seed 123 \
  > --split-type levenshtein \
  > --split-entities cdr1_alpha cdr3_alpha \
  > --split-distance 5 \
  > --distances distances.h5 \
  > -o processed_sequences_levenshtein_split.h5 \
  > $TESTDIR/data/sequences.csv

  $ python <<EOF
  > import os
  > import h5py
  > import numpy as np
  > fh = h5py.File('processed_sequences_levenshtein_split.h5')
  > ref_cdr1a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'levenshtein-split', 'cdr1_alpha.txt'))
  > ref_cdr2a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'levenshtein-split', 'cdr2_alpha.txt'))
  > ref_cdr3a = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'levenshtein-split', 'cdr3_alpha.txt'))
  > ref_cdr1b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'levenshtein-split', 'cdr1_beta.txt'))
  > ref_cdr2b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'levenshtein-split', 'cdr2_beta.txt'))
  > ref_cdr3b = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'levenshtein-split', 'cdr3_beta.txt'))
  > ref_peptide = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'levenshtein-split', 'peptide.txt'))
  > ref_mhc = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'levenshtein-split', 'mhc_pseudo.txt'))
  > ref_label = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'levenshtein-split', 'label.txt'))
  > ref_fold = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'levenshtein-split', 'fold.txt'))
  > np.testing.assert_array_equal(ref_cdr1a, fh['cdr1_alpha'][:].reshape(fh['cdr1_alpha'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr2a, fh['cdr2_alpha'][:].reshape(fh['cdr2_alpha'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr3a, fh['cdr3_alpha'][:].reshape(fh['cdr3_alpha'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr1b, fh['cdr1_beta'][:].reshape(fh['cdr1_beta'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr2b, fh['cdr2_beta'][:].reshape(fh['cdr2_beta'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_cdr3b, fh['cdr3_beta'][:].reshape(fh['cdr3_beta'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_peptide, fh['peptide'][:].reshape(fh['peptide'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_mhc, fh['mhc_pseudo'][:].reshape(fh['mhc_pseudo'][:].shape[0], -1))
  > np.testing.assert_array_equal(ref_label, fh['label'][:])
  > np.testing.assert_array_equal(ref_fold, fh['fold'][:])
  > fh.close()
  > EOF
