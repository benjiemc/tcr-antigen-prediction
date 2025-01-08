Test training model.
  $ python -m tcr_antigen_prediction.models.apps.train_confidence_predictor \
  > --seed 123 \
  > -o test.onnx \
  > $TESTDIR/data/test.h5

TODO: compare output models to something
