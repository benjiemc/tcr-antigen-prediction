rule models:
    input:
        "models/TCRen",
        "models/TCRStructMap_pMHC_split",
        "models/TCRStructMap_pMHC_split_sequence_only",
        "models/NetTCR_pMHC_split",
        "models/ConfidencePredictor_pMHC_split"

rule train_TCRen:
    input:
        summary="data/processed/structures_summary.csv",
        structures="data/processed/structures"
    output: directory("models/TCRen")
    resources:
        runtime="5m",
        mem="1GB",
        tasks=1
    shell:
        """
        mkdir -p {output}
        python -m tcr_antigen_prediction.models.apps.train_tcr_en \
            --log-level {config[log_level]} \
            -o {output}/TCRen_probabilities.csv \
            --summary-csv {input.summary} \
            {input.structures}/*.pdb
        """

rule train_TCRStructMap:
    input:
        data="data/processed/sequences_{split_type}_split.h5",
        contact_maps="data/processed/contact_maps.h5"
    output: directory("models/TCRStructMap_{split_type}_split")
    log: "data/logs/train_TCRStructMap_{split_type}_split.log"
    resources:
        runtime="5h",
        mem="5GB",
        gpu=1,
        tasks=1
    shell:
        """
        python -m tcr_antigen_prediction.models.apps.train_tcr_struct_map \
            --log-level {config[log_level]} \
            --log-file {log} \
            --seed {config[seed]} \
            -o {output} \
            --contact-maps {input.contact_maps} \
            {input.data}
        """

rule train_TCRStructMap_sequence_only:
    input: "data/processed/sequences_{split_type}_split.h5"
    output: directory("models/TCRStructMap_{split_type}_split_sequence_only")
    log: "data/logs/train_TCRStructMap_{split_type}_split_sequence_only.log"
    resources:
        runtime="5h",
        mem="5GB",
        gpu=1,
        tasks=1
    shell:
        """
        python -m tcr_antigen_prediction.models.apps.train_tcr_struct_map \
            --log-level {config[log_level]} \
            --log-file {log} \
            --seed {config[seed]} \
            -o {output} \
            {input}
        """

rule train_TCRStructMap_cdrs_peptide:
    input:
        data="data/processed/sequences_{split_type}_split.h5",
        contact_maps="data/processed/contact_maps.h5"
    output: directory("models/TCRStructMap_{split_type}_split_cdrs_peptide")
    log: "data/logs/train_TCRStructMap_{split_type}_split_cdrs_peptide.log"
    resources:
        runtime="5h",
        mem="5GB",
        gpu=1,
        tasks=1
    shell:
        """
        python -m tcr_antigen_prediction.models.apps.train_tcr_struct_map \
            --log-level {config[log_level]} \
            --log-file {log} \
            --seed {config[seed]} \
            -o {output} \
            --features-to-include cdr1_alpha cdr2_alpha cdr3_alpha cdr1_beta cdr2_beta cdr3_beta peptide \
            --contact-maps {input.contact_maps} \
            {input.data}
        """

rule train_TCRStructMap_cdrs_mhc_pseudo:
    input:
        data="data/processed/sequences_{split_type}_split.h5",
        contact_maps="data/processed/contact_maps.h5"
    output: directory("models/TCRStructMap_{split_type}_split_cdrs_mhc_pseudo")
    log: "data/logs/train_TCRStructMap_{split_type}_split_cdrs_mhc_pseudo.log"
    resources:
        runtime="5h",
        mem="5GB",
        gpu=1,
        tasks=1
    shell:
        """
        python -m tcr_antigen_prediction.models.apps.train_tcr_struct_map \
            --log-level {config[log_level]} \
            --log-file {log} \
            --seed {config[seed]} \
            -o {output} \
            --features-to-include cdr1_alpha cdr2_alpha cdr3_alpha cdr1_beta cdr2_beta cdr3_beta mhc_pseudo \
            --contact-maps {input.contact_maps} \
            {input.data}
        """

rule train_NetTCR:
    input: "data/processed/sequences_right_pad_blosum_{split_type}_split.h5"
    output: directory("models/NetTCR_{split_type}_split")
    log: "data/logs/train_NetTCR_{split_type}_split.log"
    resources:
        runtime="10h",
        mem="20GB",
        gpu=1,
        tasks=1
    shell:
        """
        python -m nettcr.apps.train_nettcr \
            --log-level {config[log_level]} \
            --log-file {log} \
            --seed {config[seed]} \
            -o {output} \
            {input}
        """

rule train_confidence_predictor:
    input: "data/processed/sequences_{split_type}_split.h5"
    output: directory("models/ConfidencePredictor_{split_type}_split")
    log: "data/logs/train_confidence_predictor_{split_type}_split.log"
    resources:
        runtime="5h",
        mem="150GB",
        tasks=1
    shell:
        """
        mkdir -p {output}
        python -m tcr_antigen_prediction.models.apps.train_confidence_predictor \
            --log-level {config[log_level]} \
            --log-file {log} \
            --seed {config[seed]} \
            -o {output}/model.onnx \
            {input}
        """
