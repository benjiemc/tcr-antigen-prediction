Test training model.
  $ python -m tcr_antigen_prediction.models.apps.train_tcr_struct_map \
  > --seed 123 \
  > --batch-size 1 \
  > --num-epochs 1 \
  > --eval-interval 1 \
  > -o model \
  > --contact-maps $TESTDIR/data/contact_maps.h5 \
  > $TESTDIR/data/test.h5
  *: UndefinedMetricWarning: Only one class is present in y_true. ROC AUC score is not defined in that case. (glob)
    warnings.warn(
  *: UndefinedMetricWarning: Only one class is present in y_true. ROC AUC score is not defined in that case. (glob)
    warnings.warn(
  *: UndefinedMetricWarning: Only one class is present in y_true. ROC AUC score is not defined in that case. (glob)
    warnings.warn(
  *: UndefinedMetricWarning: Only one class is present in y_true. ROC AUC score is not defined in that case. (glob)
    warnings.warn(
  *: UndefinedMetricWarning: Only one class is present in y_true. ROC AUC score is not defined in that case. (glob)
    warnings.warn(
  *: UndefinedMetricWarning: Only one class is present in y_true. ROC AUC score is not defined in that case. (glob)
    warnings.warn(
  *: UndefinedMetricWarning: Only one class is present in y_true. ROC AUC score is not defined in that case. (glob)
    warnings.warn(
  *: UndefinedMetricWarning: Only one class is present in y_true. ROC AUC score is not defined in that case. (glob)
    warnings.warn(

  $ python <<EOF
  > import os
  > import h5py
  > import numpy as np
  > import torch
  > from tcr_antigen_prediction.models import TCRStructMap
  > with h5py.File(os.path.join(os.environ['TESTDIR'], 'data', 'test.h5')) as fh:
  >     cdr1_alphas = fh['cdr1_alpha'][:]
  >     cdr2_alphas = fh['cdr2_alpha'][:]
  >     cdr3_alphas = fh['cdr3_alpha'][:]
  >     cdr1_betas = fh['cdr1_beta'][:]
  >     cdr2_betas = fh['cdr2_beta'][:]
  >     cdr3_betas = fh['cdr3_beta'][:]
  >     peptides = fh['peptide'][:]
  >     mhcs = fh['mhc_pseudo'][:]
  > ref_model1 = TCRStructMap()
  > ref_model1.load_state_dict(torch.load(
  >     os.path.join(os.environ['TESTDIR'], 'reference', 'model_1.pt'),
  >     map_location=torch.device('cpu'),
  >     weights_only=True,
  > ))
  > ref_model1.eval()
  > ref_predictions1 = ref_model1(
  >     torch.tensor(cdr1_alphas, dtype=torch.float32),
  >     torch.tensor(cdr2_alphas, dtype=torch.float32),
  >     torch.tensor(cdr3_alphas, dtype=torch.float32),
  >     torch.tensor(cdr1_betas, dtype=torch.float32),
  >     torch.tensor(cdr2_betas, dtype=torch.float32),
  >     torch.tensor(cdr3_betas, dtype=torch.float32),
  >     torch.tensor(peptides, dtype=torch.float32),
  >     torch.tensor(mhcs, dtype=torch.float32),
  > )
  > test_model1 = TCRStructMap()
  > test_model1.load_state_dict(torch.load(
  >     os.path.join('model', 'model_1.pt'),
  >     map_location=torch.device('cpu'),
  >     weights_only=True,
  > ))
  > test_model1.eval()
  > test_predictions1 = test_model1(
  >     torch.tensor(cdr1_alphas, dtype=torch.float32),
  >     torch.tensor(cdr2_alphas, dtype=torch.float32),
  >     torch.tensor(cdr3_alphas, dtype=torch.float32),
  >     torch.tensor(cdr1_betas, dtype=torch.float32),
  >     torch.tensor(cdr2_betas, dtype=torch.float32),
  >     torch.tensor(cdr3_betas, dtype=torch.float32),
  >     torch.tensor(peptides, dtype=torch.float32),
  >     torch.tensor(mhcs, dtype=torch.float32),
  > )
  > np.testing.assert_array_almost_equal(test_predictions1.detach().numpy(), ref_predictions1.detach().numpy())
  > ref_model2 = TCRStructMap()
  > ref_model2.load_state_dict(torch.load(
  >     os.path.join(os.environ['TESTDIR'], 'reference', 'model_2.pt'),
  >     map_location=torch.device('cpu'),
  >     weights_only=True,
  > ))
  > ref_model2.eval()
  > ref_predictions2 = ref_model2(
  >     torch.tensor(cdr1_alphas, dtype=torch.float32),
  >     torch.tensor(cdr2_alphas, dtype=torch.float32),
  >     torch.tensor(cdr3_alphas, dtype=torch.float32),
  >     torch.tensor(cdr1_betas, dtype=torch.float32),
  >     torch.tensor(cdr2_betas, dtype=torch.float32),
  >     torch.tensor(cdr3_betas, dtype=torch.float32),
  >     torch.tensor(peptides, dtype=torch.float32),
  >     torch.tensor(mhcs, dtype=torch.float32),
  > )
  > test_model2 = TCRStructMap()
  > test_model2.load_state_dict(torch.load(
  >     os.path.join('model', 'model_2.pt'),
  >     map_location=torch.device('cpu'),
  >     weights_only=True,
  > ))
  > test_model2.eval()
  > test_predictions2 = test_model2(
  >     torch.tensor(cdr1_alphas, dtype=torch.float32),
  >     torch.tensor(cdr2_alphas, dtype=torch.float32),
  >     torch.tensor(cdr3_alphas, dtype=torch.float32),
  >     torch.tensor(cdr1_betas, dtype=torch.float32),
  >     torch.tensor(cdr2_betas, dtype=torch.float32),
  >     torch.tensor(cdr3_betas, dtype=torch.float32),
  >     torch.tensor(peptides, dtype=torch.float32),
  >     torch.tensor(mhcs, dtype=torch.float32),
  > )
  > np.testing.assert_array_almost_equal(test_predictions2.detach().numpy(), ref_predictions2.detach().numpy())
  > EOF
