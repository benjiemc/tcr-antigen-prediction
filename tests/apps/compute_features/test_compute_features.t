Compute pair features
  $ python -m tcr_antigen_prediction.apps.compute_features \
  > -o ./ \
  > --mode ppi_search \
  > $TESTDIR/data/1ao7_DECAB_CA.ply $TESTDIR/data/1ao7_DECAB_DE.ply
  */site-packages/numpy/core/_asarray.py:136: VisibleDeprecationWarning: Creating an ndarray from ragged nested sequences (which is a list-or-tuple of lists-or-tuples-or ndarrays with different lengths or shapes) is deprecated. If you meant to do this, you must specify 'dtype=object' when creating the ndarray (glob)
    return array(a, dtype, copy=False, order=order, subok=True)

  $ python $TESTDIR/compare_numpy_arrays.py $TESTDIR/reference ./