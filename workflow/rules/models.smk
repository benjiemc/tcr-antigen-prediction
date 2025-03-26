rule models:
    input:
        "models/TCRen",
        "models/TCRStructMap",
        "models/TCRStructMap_sequence_only",
        "models/NetTCR",
        "models/ConfidencePredictor"

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
        data="data/processed/sequences.h5",
        contact_maps="data/processed/contact_maps.h5"
    output: directory("models/TCRStructMap")
    log: "data/logs/train_tcr_struct_map.log"
    resources:
        runtime="2h",
        mem="5GB",
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
    input: "data/processed/sequences.h5"
    output: directory("models/TCRStructMap_sequence_only")
    log: "data/logs/train_tcr_contact_map_predictor_sequence_only.log"
    resources:
        runtime="2h",
        mem="5GB",
        tasks=1
    shell:
        """
        python -m tcr_antigen_prediction.models.apps.train_tcr_contact_map_predictor \
            --log-level {config[log_level]} \
            --log-file {log} \
            --seed {config[seed]} \
            -o {output} \
            {input}
        """

rule train_NetTCR:
    input: "data/processed/sequences_right_pad_blosum.h5"
    output: directory("models/NetTCR")
    log: "data/logs/train_NetTCR.log"
    resources:
        runtime="10h",
        mem="20GB",
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
    input: "data/processed/sequences.h5"
    output: directory("models/ConfidencePredictor")
    log: "data/logs/train_confidence_predictor.log"
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
