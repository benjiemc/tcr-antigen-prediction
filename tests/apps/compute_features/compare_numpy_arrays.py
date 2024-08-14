import glob
import os
import sys

import numpy as np

for file_ in glob.glob(os.path.join(sys.argv[1], '*.npy')):
    arr1 = np.load(file_, allow_pickle=True)
    arr2 = np.load(os.path.join(sys.argv[2], os.path.basename(file_)), allow_pickle=True)

    np.testing.assert_array_equal(arr1, arr2)