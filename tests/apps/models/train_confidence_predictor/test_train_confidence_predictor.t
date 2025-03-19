Test training model.
  $ python -m tcr_antigen_prediction.models.apps.train_confidence_predictor \
  > --seed 123 \
  > --test-size 0.5 \
  > -o test.onnx \
  > $TESTDIR/data/test.h5

  $ python <<EOF
  > import os
  > import h5py
  > import numpy as np
  > import onnxruntime as rt
  > sess_options = rt.SessionOptions()
  > sess_options.intra_op_num_threads = 1
  > sess_options.inter_op_num_threads = 1
  > with h5py.File(os.path.join(os.environ['TESTDIR'], 'data', 'test.h5')) as fh:
  >     cdr1_alphas = fh['cdr1_alpha'][:]
  >     cdr2_alphas = fh['cdr2_alpha'][:]
  >     cdr3_alphas = fh['cdr3_alpha'][:]
  >     cdr1_betas = fh['cdr1_beta'][:]
  >     cdr2_betas = fh['cdr2_beta'][:]
  >     cdr3_betas = fh['cdr3_beta'][:]
  >     peptides = fh['peptide'][:]
  >     mhcs = fh['mhc_pseudo'][:]
  > tcr_pmhcs = np.concatenate(
  >     [cdr1_alphas, cdr2_alphas, cdr3_alphas, cdr1_betas, cdr2_betas, cdr3_betas, peptides, mhcs], axis=1
  > )
  > tcr_pmhcs_flat = tcr_pmhcs.reshape(tcr_pmhcs.shape[0], -1)
  > ref_sess1 = rt.InferenceSession(os.path.join(os.environ['TESTDIR'], 'reference', 'test_1.onnx'), sess_options=sess_options)
  > ref_confidence1, *_ = ref_sess1.run(None, {'X': tcr_pmhcs_flat.astype(np.float32)})
  > test_sess1 = rt.InferenceSession('test_1.onnx', sess_options=sess_options)
  > test_confidence1, *_ = test_sess1.run(None, {'X': tcr_pmhcs_flat.astype(np.float32)})
  > np.testing.assert_array_almost_equal(test_confidence1, ref_confidence1)
  > ref_sess2 = rt.InferenceSession(os.path.join(os.environ['TESTDIR'], 'reference', 'test_2.onnx'), sess_options=sess_options)
  > ref_confidence2, *_ = ref_sess2.run(None, {'X': tcr_pmhcs_flat.astype(np.float32)})
  > test_sess2 = rt.InferenceSession('test_2.onnx', sess_options=sess_options)
  > test_confidence2, *_ = test_sess2.run(None, {'X': tcr_pmhcs_flat.astype(np.float32)})
  > np.testing.assert_array_almost_equal(test_confidence2, ref_confidence2)
  > EOF
