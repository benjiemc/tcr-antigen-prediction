Test MHC1
  $ python -m tcr_antigen_prediction.data.apps.crop_tcr_pmhc \
  > -o tcr_mhc1_crop.pdb \
  > $TESTDIR/data/1ao7_DECAB.pdb \
  > --tcr-chains D E \
  > --mhc-chains A \
  > --antigen-chain C

  $ diff tcr_mhc1_crop.pdb $TESTDIR/reference/tcr_mhc1_crop.pdb

Test MHC1 without Het Atoms
  $ python -m tcr_antigen_prediction.data.apps.crop_tcr_pmhc \
  > -o tcr_mhc1_crop_no_het.pdb \
  > $TESTDIR/data/1ao7_DECAB.pdb \
  > --remove-het-atoms \
  > --tcr-chains D E \
  > --mhc-chains A \
  > --antigen-chain C

  $ diff tcr_mhc1_crop_no_het.pdb $TESTDIR/reference/tcr_mhc1_crop_no_het.pdb

Test MHC2 without Het Atoms
  $ python -m tcr_antigen_prediction.data.apps.crop_tcr_pmhc \
  > -o tcr_mhc2_crop_no_het.pdb \
  > $TESTDIR/data/2pxy_ABPCD.pdb \
  > --remove-het-atoms \
  > --tcr-chains A B \
  > --mhc-chains C D \
  > --antigen-chain P

  $ diff tcr_mhc2_crop_no_het.pdb $TESTDIR/reference/tcr_mhc2_crop_no_het.pdb
