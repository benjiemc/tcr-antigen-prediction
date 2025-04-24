rule visualisations:
    input:
        'report/figures/TCRStructMap_training.svg',
        'report/figures/TCRStructMap_sequence_only_training.svg',
        'report/figures/TCRStructMap_contact_probabilities_training.svg'
        'report/figures/TCRStructMap_cdrs_peptide_training.svg',
        'report/figures/TCRStructMap_cdrs_mhc_pseudo_training.svg',
        'report/figures/NetTCR_training.svg'

rule visualise_model_training:
    input: 'data/logs/train_{model_name}.log'
    output: 'report/figures/{model_name}_training.svg'
    resources:
        runtime="1m",
        mem="100MB",
        tasks=1
    shell:
        """
        python -m tcr_antigen_prediction.visualisations.apps.visualise_model_training \
            --title {wildcards.model_name} \
            -o {output} \
            {input}
        """
