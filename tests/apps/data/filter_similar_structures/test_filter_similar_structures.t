  $ python -m tcr_antigen_prediction.data.apps.filter_similar_structures \
  > --structural-similarity-cutoff 2.0 \
  > --summary-csv "$TESTDIR/data/structures_summary.csv" \
  > "$TESTDIR/data/"
  name,Achain,Bchain,antigen_chain,mhc_chain1,mhc_chain2,mhc_type,cdr1_alpha,cdr2_alpha,cdr3_alpha,cdr1_beta,cdr2_beta,cdr3_beta,peptide
  7q9b_DECA.pdb,D,E,C,A,,MH1,DRGSQS,IYSNGD,AVQKLV,MNHEY,SVGAGI,ASSYSFTEATYEQY,EAAGIGILTV
  3qiw_CDEAB.pdb,C,D,E,A,B,MH2,TTMRA,LASGT,AAEPSSGQKLV,KGHPV,FQNQEV,ASSLNNANSDYT,ADLIAYLEQATKG
  3tf7_cCBA.pdb,c,C,B,A,,MH1,YSATPY,YYSGDPVV,AVSAKGTGSKLS,NSHNY,SYGAGN,ASSDAPGQLY,QLSPFPFDL