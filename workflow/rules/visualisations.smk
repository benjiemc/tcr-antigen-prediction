rule visualisations:
    input:
        'report/figures/CDR1alpha_peptide_contacts.pdf',
        'report/figures/CDR2alpha_peptide_contacts.pdf',
        'report/figures/CDR3alpha_peptide_contacts.pdf',
        'report/figures/CDR1beta_peptide_contacts.pdf',
        'report/figures/CDR2beta_peptide_contacts.pdf',
        'report/figures/CDR3beta_peptide_contacts.pdf',
        'report/figures/CDR1alpha_mhc_contacts.pdf',
        'report/figures/CDR2alpha_mhc_contacts.pdf',
        'report/figures/CDR3alpha_mhc_contacts.pdf',
        'report/figures/CDR1beta_mhc_contacts.pdf',
        'report/figures/CDR2beta_mhc_contacts.pdf',
        'report/figures/CDR3beta_mhc_contacts.pdf',
        'report/figures/TCRStructMap_pMHC_split_training.svg',
        'report/figures/TCRStructMap_pMHC_split_sequence_only_training.svg',
        'report/figures/TCRStructMap_pMHC_split_cdrs_peptide_training.svg',
        'report/figures/TCRStructMap_pMHC_split_cdrs_mhc_pseudo_training.svg',
        'report/figures/NetTCR_pMHC_split_training.svg'

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
            --log-level {config[log_level]} \
            --title {wildcards.model_name} \
            -o {output} \
            {input}
        """

rule create_contact_map_plots:
    input:
        contacts="data/interim/tcr_pmhc_contacts.csv",
        mhc_pseudo_positions="data/interim/mhc_pseudo_seq_imgt_positions.json"
    resources:
        runtime="1m",
        mem="100MB",
        tasks=1
    output: expand("report/figures/{cdr_name}_{interaction}_contacts.pdf", cdr_name=[f'CDR{num}{chain}' for num in range(1, 4) for chain in ('alpha', 'beta')], interaction=['peptide', 'mhc'])
    shell:
        """
        python -m tcr_antigen_prediction.visualisations.apps.visualise_contact_maps \
            --log-level {config[log_level]} \
            --separate-plots \
            -o report/figures/peptide_contacts.pdf \
            --mhc-types MH1 \
            --interaction peptide \
            --pad \
            {input.contacts}

        python -m tcr_antigen_prediction.visualisations.apps.visualise_contact_maps \
            --log-level {config[log_level]} \
            --separate-plots \
            -o report/figures/mhc_contacts.pdf \
            --mhc-types MH1 \
            --interaction MHC \
            --mhc-pseudo-seq-positions {input.mhc_pseudo_positions} \
            --pad \
            {input.contacts}
        """
