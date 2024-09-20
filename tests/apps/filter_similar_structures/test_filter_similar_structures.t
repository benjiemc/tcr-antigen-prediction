  $ python -m tcr_antigen_prediction.apps.filter_similar_structures \
  > --structural-similarity-cutoff 2.0 \
  > --summary-csv "$TESTDIR/data/structures_summary.csv" \
  > "$TESTDIR/data/"
  name,Achain,Bchain,antigen_chain,mhc_chain1,mhc_type,CDR1alpha_sequence,CDR2alpha_sequence,CDR3alpha_sequence,CDR1beta_sequence,CDR2beta_sequence,CDR3beta_sequence,peptide_sequence,collated_cdrs,mhc_chain2
  7q9b_DECA,D,E,C,A,MH1,DRGSQS,IYSNGD,AVQKLV,MNHEY,SVGAGI,ASSYSFTEATYEQY,EAAGIGILTV,DRGSQS-IYSNGD-AVQKLV-MNHEY-SVGAGI-ASSYSFTEATYEQY,
  3qiw_CDEAB,C,D,E,A,MH2,TTMRA,LASGT,AAEPSSGQKLV,KGHPV,FQNQEV,ASSLNNANSDYT,ADLIAYLEQATKG,TTMRA-LASGT-AAEPSSGQKLV-KGHPV-FQNQEV-ASSLNNANSDYT,B
  3tf7_cCBA,c,C,B,A,MH1,YSATPY,YYSGDPVV,AVSAKGTGSKLS,NSHNY,SYGAGN,ASSDAPGQLY,QLSPFPFDL,YSATPY-YYSGDPVV-AVSAKGTGSKLS-NSHNY-SYGAGN-ASSDAPGQLY,