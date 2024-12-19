Test training model.
  $ python -m tcr_antigen_prediction.apps.train_tcr_contact_map_predictor \
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

TODO: compare output models to something
