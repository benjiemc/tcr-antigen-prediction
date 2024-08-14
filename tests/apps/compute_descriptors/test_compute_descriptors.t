Test Example PPI
  $ python -m tcr_antigen_prediction.apps.compute_descriptors \
  > --model "$TESTDIR/../../../models/baseline_masif_ppi/model" \
  > -o ./ \
  > $TESTDIR/data  2> /dev/null

  $ python -c "import numpy as np; import os; arr1 = np.load('example-ppi/p1_desc_straight.npy'); arr2 = np.load(os.path.join(os.environ['TESTDIR'], 'reference/example-ppi/p1_desc_straight.npy')); np.testing.assert_array_almost_equal(arr1, arr2)"
  $ python -c "import numpy as np; import os; arr1 = np.load('example-ppi/p2_desc_straight.npy'); arr2 = np.load(os.path.join(os.environ['TESTDIR'], 'reference/example-ppi/p2_desc_straight.npy')); np.testing.assert_array_almost_equal(arr1, arr2)"
  $ python -c "import numpy as np; import os; arr1 = np.load('example-ppi/p1_desc_flipped.npy'); arr2 = np.load(os.path.join(os.environ['TESTDIR'], 'reference/example-ppi/p1_desc_flipped.npy')); np.testing.assert_array_almost_equal(arr1, arr2)"
  $ python -c "import numpy as np; import os; arr1 = np.load('example-ppi/p2_desc_flipped.npy'); arr2 = np.load(os.path.join(os.environ['TESTDIR'], 'reference/example-ppi/p2_desc_flipped.npy')); np.testing.assert_array_almost_equal(arr1, arr2)"