Test.
  $ python -m tcr_antigen_prediction.data.apps.process_interacting_residues \
  > --residue-column-names residue_name_tcr residue_name_pmhc \
  > -o standard.txt \
  > $TESTDIR/data/tcr_pmhc_contacts.csv

  $ python << EOF
  > import os
  > import numpy as np
  > test = np.loadtxt('standard.txt')
  > ref = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'standard.txt'))
  > np.testing.assert_array_almost_equal(test, ref)
  > EOF

Test querying antigen chain.
  $ python -m tcr_antigen_prediction.data.apps.process_interacting_residues \
  > --residue-column-names residue_name_tcr residue_name_pmhc \
  > --selection-query "chain_type_pmhc == 'antigen_chain'" \
  > -o query.txt \
  > $TESTDIR/data/tcr_pmhc_contacts.csv

  $ python << EOF
  > import os
  > import numpy as np
  > test = np.loadtxt('query.txt')
  > ref = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'query.txt'))
  > np.testing.assert_array_almost_equal(test, ref)
  > EOF
