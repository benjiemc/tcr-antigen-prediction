rule models:
    input: "models/TCRen", "models/TCRContactMapPredictor", "models/TCRContactMapPredictor_sequence_only"

rule train_TCRen:
    input: "data/processed/selected-stcrdab"
    output: directory("models/TCRen")
    resources:
        runtime="5m",
        mem="1GB",
        tasks=1
    shell:
        """
        @mkdir -p {output}
        @python -m tcr_antigen_prediction.models.apps.train_tcr_en \
            --log-level {config[log_level]} \
            -o {output}/TCRen_probabilities.csv \
            --summary-csv "{input}/stcrdab_split.csv" \
            $(cat "{input}/stcrdab_split.csv" | grep "train" | awk -F, -v dir="{input}" '{{ printf "%s/%s_%s%s%s%s%s.pdb ", dir, $1, $2, $3, $4, $5, $6 }}')
        """

rule train_TCRContactMapPredictor:
    input:
        data="data/processed/nettcr.h5",
        contact_maps="data/processed/contact_maps.h5"
    output: directory("models/TCRContactMapPredictor")
    log: "data/logs/train_tcr_contact_map_predictor.log"
    resources:
        runtime="1h",
        mem="5GB",
        tasks=1
    shell:
        """
        python -m tcr_antigen_prediction.models.apps.train_tcr_contact_map_predictor \
            --log-level {config[log_level]} \
            --log-file {log} \
            -o {output} \
            --contact-maps {input.contact_maps} \
            {input.data}
        """

rule train_TCRContactMapPredictor_sequence_only:
    input: "data/processed/nettcr.h5"
    output: directory("models/TCRContactMapPredictor_sequence_only")
    log: "data/logs/train_tcr_contact_map_predictor_sequence_only.log"
    resources:
        runtime="1h",
        mem="5GB",
        tasks=1
    shell:
        """
        python -m tcr_antigen_prediction.models.apps.train_tcr_contact_map_predictor \
            --log-level {config[log_level]} \
            --log-file {log} \
            -o {output} \
            {input}
        """
