wildcard_constraints:
    model_type="(_sequence_only)?",
    split_type="(random)|(tcr)|(peptide)|(pMHC)|(levenshtein)"

rule benchmark:
    input:
        "data/processed/TCRStructMap_predictions_on_sequences_pMHC_split.csv",
        "data/processed/TCRStructMap_sequence_only_predictions_on_sequences_pMHC_split.csv",
        "data/processed/TCRStructMap_cdrs_peptide_predictions_on_sequences_pMHC_split.csv",
        "data/processed/TCRStructMap_cdrs_mhc_pseudo_predictions_on_sequences_pMHC_split.csv",
        "data/processed/NetTCR_predictions_on_sequences_pMHC_split.csv",

rule run_tcr_struct_map_predictions:
    input:
        sequences="data/processed/sequences_{split_type}_split.h5",
        model="models/TCRStructMap_{split_type}_split{model_type}"
    output: "data/processed/TCRStructMap{model_type}_predictions_on_sequences_{split_type}_split.csv"
    resources:
        runtime="5m",
        mem="10GB",
        gpu=1,
        tasks=1
    script: "../scripts/run_tcr_struct_map_predictions.py"

rule run_tcr_struct_map_predictions_cdr_peptide:
    input:
        sequences="data/processed/sequences_{split_type}_split.h5",
        contact_maps="data/processed/contact_maps.h5",
        model="models/TCRStructMap_{split_type}_split_cdrs_peptide"
    output: "data/processed/TCRStructMap_cdrs_peptide_predictions_on_sequences_{split_type}_split.csv"
    resources:
        runtime="5m",
        mem="10GB",
        gpu=1,
        tasks=1
    script: "../scripts/run_tcr_struct_map_predictions_cdrs_peptide.py"

rule run_tcr_struct_map_predictions_cdr_mhc_pseudo:
    input:
        sequences="data/processed/sequences_{split_type}_split.h5",
        contact_maps="data/processed/contact_maps.h5",
        model="models/TCRStructMap_{split_type}_split_cdrs_mhc_pseudo"
    output: "data/processed/TCRStructMap_cdrs_mhc_pseudo_predictions_on_sequences_{split_type}_split.csv"
    resources:
        runtime="5m",
        mem="10GB",
        gpu=1,
        tasks=1
    script: "../scripts/run_tcr_struct_map_predictions_cdrs_mhc_pseudo.py"


rule run_nettcr_predictions:
    input:
        sequences="data/processed/sequences_right_pad_blosum_{split_type}_split.h5",
        model="models/NetTCR_{split_type}_split"
    output: "data/processed/NetTCR_predictions_on_sequences_{split_type}_split.csv"
    resources:
        runtime="5m",
        mem="10GB",
        gpu=1,
        tasks=1
    script: "../scripts/run_nettcr_predictions.py"
