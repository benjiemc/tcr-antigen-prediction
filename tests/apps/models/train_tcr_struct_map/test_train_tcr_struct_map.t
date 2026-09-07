Test training model.
  $ python -m tcr_antigen_prediction.models.apps.train_tcr_struct_map \
  > --log-level error \
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
  >     os.path.join(os.environ['TESTDIR'], 'reference', 'full', 'model_1.pt'),
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
  >     os.path.join(os.environ['TESTDIR'], 'reference', 'full', 'model_2.pt'),
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

Test training model without MHC.
  $ python -m tcr_antigen_prediction.models.apps.train_tcr_struct_map \
  > --log-level error \
  > --seed 123 \
  > --features-to-include cdr1_alpha cdr2_alpha cdr3_alpha cdr1_beta cdr2_beta cdr3_beta peptide \
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
  > with h5py.File(os.path.join(os.environ['TESTDIR'], 'data', 'contact_maps.h5')) as fh:
  >     cdr_peptide_contact_maps = (
  >         torch.tensor(fh['peptide']['cdr1_alpha'][:], dtype=torch.float32),
  >         torch.tensor(fh['peptide']['cdr2_alpha'][:], dtype=torch.float32),
  >         torch.tensor(fh['peptide']['cdr3_alpha'][:], dtype=torch.float32),
  >         torch.tensor(fh['peptide']['cdr1_beta'][:], dtype=torch.float32),
  >         torch.tensor(fh['peptide']['cdr2_beta'][:], dtype=torch.float32),
  >         torch.tensor(fh['peptide']['cdr3_beta'][:], dtype=torch.float32),
  >     )
  >     cdr_mhc_contact_maps = (
  >         torch.tensor(fh['mhc_pseudo']['cdr1_alpha'][:], dtype=torch.float32),
  >         torch.tensor(fh['mhc_pseudo']['cdr2_alpha'][:], dtype=torch.float32),
  >         torch.tensor(fh['mhc_pseudo']['cdr3_alpha'][:], dtype=torch.float32),
  >         torch.tensor(fh['mhc_pseudo']['cdr1_beta'][:], dtype=torch.float32),
  >         torch.tensor(fh['mhc_pseudo']['cdr2_beta'][:], dtype=torch.float32),
  >         torch.tensor(fh['mhc_pseudo']['cdr3_beta'][:], dtype=torch.float32),
  >     )
  > with h5py.File(os.path.join(os.environ['TESTDIR'], 'data', 'test.h5')) as fh:
  >     cdr1_alphas = fh['cdr1_alpha'][:]
  >     cdr2_alphas = fh['cdr2_alpha'][:]
  >     cdr3_alphas = fh['cdr3_alpha'][:]
  >     cdr1_betas = fh['cdr1_beta'][:]
  >     cdr2_betas = fh['cdr2_beta'][:]
  >     cdr3_betas = fh['cdr3_beta'][:]
  >     peptides = fh['peptide'][:]
  > ref_model1 = TCRStructMap(cdr_peptide_contact_maps, cdr_mhc_contact_maps, mhc_length=0)
  > ref_model1.load_state_dict(torch.load(
  >     os.path.join(os.environ['TESTDIR'], 'reference', 'no-mhc', 'model_1.pt'),
  >     map_location=torch.device('cpu'),
  >     weights_only=True,
  > ))
  > ref_model1.eval()
  > ref_predictions1 = ref_model1(**{
  >     'cdr1_alpha': torch.tensor(cdr1_alphas, dtype=torch.float32),
  >     'cdr2_alpha': torch.tensor(cdr2_alphas, dtype=torch.float32),
  >     'cdr3_alpha': torch.tensor(cdr3_alphas, dtype=torch.float32),
  >     'cdr1_beta': torch.tensor(cdr1_betas, dtype=torch.float32),
  >     'cdr2_beta': torch.tensor(cdr2_betas, dtype=torch.float32),
  >     'cdr3_beta': torch.tensor(cdr3_betas, dtype=torch.float32),
  >     'peptide': torch.tensor(peptides, dtype=torch.float32),
  > })
  > test_model1 = TCRStructMap(cdr_peptide_contact_maps, cdr_mhc_contact_maps, mhc_length=0)
  > test_model1.load_state_dict(torch.load(
  >     os.path.join('model', 'model_1.pt'),
  >     map_location=torch.device('cpu'),
  >     weights_only=True,
  > ))
  > test_model1.eval()
  > test_predictions1 = test_model1(**{
  >     'cdr1_alpha': torch.tensor(cdr1_alphas, dtype=torch.float32),
  >     'cdr2_alpha': torch.tensor(cdr2_alphas, dtype=torch.float32),
  >     'cdr3_alpha': torch.tensor(cdr3_alphas, dtype=torch.float32),
  >     'cdr1_beta': torch.tensor(cdr1_betas, dtype=torch.float32),
  >     'cdr2_beta': torch.tensor(cdr2_betas, dtype=torch.float32),
  >     'cdr3_beta': torch.tensor(cdr3_betas, dtype=torch.float32),
  >     'peptide': torch.tensor(peptides, dtype=torch.float32),
  > })
  > np.testing.assert_array_almost_equal(test_predictions1.detach().numpy(), ref_predictions1.detach().numpy())
  > ref_model2 = TCRStructMap(cdr_peptide_contact_maps, cdr_mhc_contact_maps, mhc_length=0)
  > ref_model2.load_state_dict(torch.load(
  >     os.path.join(os.environ['TESTDIR'], 'reference', 'no-mhc', 'model_2.pt'),
  >     map_location=torch.device('cpu'),
  >     weights_only=True,
  > ))
  > ref_model2.eval()
  > ref_predictions2 = ref_model2(**{
  >     'cdr1_alpha': torch.tensor(cdr1_alphas, dtype=torch.float32),
  >     'cdr2_alpha': torch.tensor(cdr2_alphas, dtype=torch.float32),
  >     'cdr3_alpha': torch.tensor(cdr3_alphas, dtype=torch.float32),
  >     'cdr1_beta': torch.tensor(cdr1_betas, dtype=torch.float32),
  >     'cdr2_beta': torch.tensor(cdr2_betas, dtype=torch.float32),
  >     'cdr3_beta': torch.tensor(cdr3_betas, dtype=torch.float32),
  >     'peptide': torch.tensor(peptides, dtype=torch.float32),
  > })
  > test_model2 = TCRStructMap(cdr_peptide_contact_maps, cdr_mhc_contact_maps, mhc_length=0)
  > test_model2.load_state_dict(torch.load(
  >     os.path.join('model', 'model_2.pt'),
  >     map_location=torch.device('cpu'),
  >     weights_only=True,
  > ))
  > test_model2.eval()
  > test_predictions2 = test_model2(**{
  >     'cdr1_alpha': torch.tensor(cdr1_alphas, dtype=torch.float32),
  >     'cdr2_alpha': torch.tensor(cdr2_alphas, dtype=torch.float32),
  >     'cdr3_alpha': torch.tensor(cdr3_alphas, dtype=torch.float32),
  >     'cdr1_beta': torch.tensor(cdr1_betas, dtype=torch.float32),
  >     'cdr2_beta': torch.tensor(cdr2_betas, dtype=torch.float32),
  >     'cdr3_beta': torch.tensor(cdr3_betas, dtype=torch.float32),
  >     'peptide': torch.tensor(peptides, dtype=torch.float32),
  > })
  > np.testing.assert_array_almost_equal(test_predictions2.detach().numpy(), ref_predictions2.detach().numpy())
  > EOF

Test training model without CDR1-2.
  $ python -m tcr_antigen_prediction.models.apps.train_tcr_struct_map \
  > --log-level error \
  > --seed 123 \
  > --features-to-include cdr3_alpha cdr3_beta peptide mhc_pseudo \
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
  > with h5py.File(os.path.join(os.environ['TESTDIR'], 'data', 'contact_maps.h5')) as fh:
  >     cdr_peptide_contact_maps = (
  >         torch.tensor(fh['peptide']['cdr1_alpha'][:], dtype=torch.float32),
  >         torch.tensor(fh['peptide']['cdr2_alpha'][:], dtype=torch.float32),
  >         torch.tensor(fh['peptide']['cdr3_alpha'][:], dtype=torch.float32),
  >         torch.tensor(fh['peptide']['cdr1_beta'][:], dtype=torch.float32),
  >         torch.tensor(fh['peptide']['cdr2_beta'][:], dtype=torch.float32),
  >         torch.tensor(fh['peptide']['cdr3_beta'][:], dtype=torch.float32),
  >     )
  >     cdr_mhc_contact_maps = (
  >         torch.tensor(fh['mhc_pseudo']['cdr1_alpha'][:], dtype=torch.float32),
  >         torch.tensor(fh['mhc_pseudo']['cdr2_alpha'][:], dtype=torch.float32),
  >         torch.tensor(fh['mhc_pseudo']['cdr3_alpha'][:], dtype=torch.float32),
  >         torch.tensor(fh['mhc_pseudo']['cdr1_beta'][:], dtype=torch.float32),
  >         torch.tensor(fh['mhc_pseudo']['cdr2_beta'][:], dtype=torch.float32),
  >         torch.tensor(fh['mhc_pseudo']['cdr3_beta'][:], dtype=torch.float32),
  >     )
  > with h5py.File(os.path.join(os.environ['TESTDIR'], 'data', 'test.h5')) as fh:
  >     cdr3_alphas = fh['cdr3_alpha'][:]
  >     cdr3_betas = fh['cdr3_beta'][:]
  >     peptides = fh['peptide'][:]
  >     mhc_pseudos = fh['mhc_pseudo'][:]
  > ref_model1 = TCRStructMap(
  >     cdr_peptide_contact_maps,
  >     cdr_mhc_contact_maps,
  >     cdr1_alpha_length=0,
  >     cdr2_alpha_length=0,
  >     cdr1_beta_length=0,
  >     cdr2_beta_length=0,
  > )
  > ref_model1.load_state_dict(torch.load(
  >     os.path.join(os.environ['TESTDIR'], 'reference', 'no-cdr-1-and-2', 'model_1.pt'),
  >     map_location=torch.device('cpu'),
  >     weights_only=True,
  > ))
  > ref_model1.eval()
  > ref_predictions1 = ref_model1(**{
  >     'cdr3_alpha': torch.tensor(cdr3_alphas, dtype=torch.float32),
  >     'cdr3_beta': torch.tensor(cdr3_betas, dtype=torch.float32),
  >     'peptide': torch.tensor(peptides, dtype=torch.float32),
  >     'mhc_pseudo': torch.tensor(mhc_pseudos, dtype=torch.float32),
  > })
  > test_model1 = TCRStructMap(
  >     cdr_peptide_contact_maps,
  >     cdr_mhc_contact_maps,
  >     cdr1_alpha_length=0,
  >     cdr2_alpha_length=0,
  >     cdr1_beta_length=0,
  >     cdr2_beta_length=0,
  > )
  > test_model1.load_state_dict(torch.load(
  >     os.path.join('model', 'model_1.pt'),
  >     map_location=torch.device('cpu'),
  >     weights_only=True,
  > ))
  > test_model1.eval()
  > test_predictions1 = test_model1(**{
  >     'cdr3_alpha': torch.tensor(cdr3_alphas, dtype=torch.float32),
  >     'cdr3_beta': torch.tensor(cdr3_betas, dtype=torch.float32),
  >     'peptide': torch.tensor(peptides, dtype=torch.float32),
  >     'mhc_pseudo': torch.tensor(mhc_pseudos, dtype=torch.float32),
  > })
  > np.testing.assert_array_almost_equal(test_predictions1.detach().numpy(), ref_predictions1.detach().numpy())
  > ref_model2 = TCRStructMap(
  >     cdr_peptide_contact_maps,
  >     cdr_mhc_contact_maps,
  >     cdr1_alpha_length=0,
  >     cdr2_alpha_length=0,
  >     cdr1_beta_length=0,
  >     cdr2_beta_length=0,
  > )
  > ref_model2.load_state_dict(torch.load(
  >     os.path.join(os.environ['TESTDIR'], 'reference', 'no-cdr-1-and-2', 'model_2.pt'),
  >     map_location=torch.device('cpu'),
  >     weights_only=True,
  > ))
  > ref_model2.eval()
  > ref_predictions2 = ref_model2(**{
  >     'cdr3_alpha': torch.tensor(cdr3_alphas, dtype=torch.float32),
  >     'cdr3_beta': torch.tensor(cdr3_betas, dtype=torch.float32),
  >     'peptide': torch.tensor(peptides, dtype=torch.float32),
  >     'mhc_pseudo': torch.tensor(mhc_pseudos, dtype=torch.float32),
  > })
  > test_model2 = TCRStructMap(
  >     cdr_peptide_contact_maps,
  >     cdr_mhc_contact_maps,
  >     cdr1_alpha_length=0,
  >     cdr2_alpha_length=0,
  >     cdr1_beta_length=0,
  >     cdr2_beta_length=0,
  > )
  > test_model2.load_state_dict(torch.load(
  >     os.path.join('model', 'model_2.pt'),
  >     map_location=torch.device('cpu'),
  >     weights_only=True,
  > ))
  > test_model2.eval()
  > test_predictions2 = test_model2(**{
  >     'cdr3_alpha': torch.tensor(cdr3_alphas, dtype=torch.float32),
  >     'cdr3_beta': torch.tensor(cdr3_betas, dtype=torch.float32),
  >     'peptide': torch.tensor(peptides, dtype=torch.float32),
  >     'mhc_pseudo': torch.tensor(mhc_pseudos, dtype=torch.float32),
  > })
  > np.testing.assert_array_almost_equal(test_predictions2.detach().numpy(), ref_predictions2.detach().numpy())
  > EOF
