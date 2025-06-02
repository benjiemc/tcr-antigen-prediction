rule benchmark:
    input:
        "data/processed/TCRStructMap_predictions_on_sequences_pMHC_split.csv",
        "data/processed/TCRStructMap_sequence_only_predictions_on_sequences_pMHC_split.csv",
        "data/processed/TCRStructMap_cdrs_peptide_predictions_on_sequences_pMHC_split.csv",
        "data/processed/TCRStructMap_cdrs_mhc_pseudo_predictions_on_sequences_pMHC_split.csv",
        "data/processed/NetTCR_predictions_on_sequences_pMHC_split.csv",

rule run_tcr_struct_map_predictions:
    input:
        sequences="data/processed/sequences_{split_type}.h5",
        model="models/{model_name}"
    output: "data/processed/{model_name}_predictions_on_sequences_{split_type}.csv"
    resources:
        runtime="5m",
        mem="10GB",
        tasks=1
    script: "../scripts/run_tcr_struct_map_predictions.py"

rule run_tcr_struct_map_predictions_cdr_peptide:
    input:
        sequences="data/processed/sequences_pMHC_split.h5",
        contact_maps="data/processed/contact_maps.h5",
        model="models/TCRStructMap_cdrs_peptide"
    output: "data/processed/TCRStructMap_cdrs_peptide_predictions_on_sequences_pMHC_split.csv"
    resources:
        runtime="5m",
        mem="10GB",
        tasks=1
    script: "../scripts/run_tcr_struct_map_predictions_cdrs_peptide.py"

rule run_tcr_struct_map_predictions_cdr_mhc_pseudo:
    input:
        sequences="data/processed/sequences_pMHC_split.h5",
        contact_maps="data/processed/contact_maps.h5",
        model="models/TCRStructMap_cdrs_mhc_pseudo"
    output: "data/processed/TCRStructMap_cdrs_mhc_pseudo_predictions_on_sequences_pMHC_split.csv"
    resources:
        runtime="5m",
        mem="10GB",
        tasks=1
    script: "../scripts/run_tcr_struct_map_predictions_cdrs_mhc_pseudo.py"


rule run_nettcr_predictions:
    input:
        sequences="data/processed/sequences_right_pad_blosum.h5",
        model="models/NetTCR"
    output: "data/processed/NetTCR_predictions_on_sequences_pMHC_split.csv"
    resources:
        runtime="5m",
        mem="10GB",
        tasks=1
    script: "../scripts/run_nettcr_predictions.py"
