rule run_visualise_structure_data_notebook:
    input: "data/processed/structures_summary.csv", "data/logs/select_stcrdab_tcr_pmhc_structures.log"
    resources:
        runtime="5m",
        mem="1GB",
        tasks=1
    shell:
        """
        papermill notebooks/visualise_structure_data.ipynb notebooks/visualise_structure_data.ipynb
        """

rule run_visualise_sequence_data_notebook:
    input: "data/interim/sequences.csv", "data/logs/collate_sequence_data.log"
    resources:
        runtime="5m",
        mem="1GB",
        tasks=1
    notebook: "../../notebooks/visualise_sequence_data.ipynb"

rule run_compairing_sequence_and_structure_data_notebook:
    input: "data/processed/structures_summary.csv", "data/interim/sequences.csv"
    resources:
        runtime="5m",
        mem="1GB",
        tasks=1
    shell:
        """
        papermill notebooks/compairing_sequence_and_structure_data.ipynb notebooks/compairing_sequence_and_structure_data.ipynb
        """

rule run_benchmark_model_performance_notebook:
    input:
        "data/interim/sequences_pMHC_split.csv",
        "data/interim/peptides.txt",
        "data/interim/peptide_distances.txt",
        "data/processed/TCRStructMap_predictions_on_sequences.csv",
        "data/processed/TCRStructMap_sequence_only_predictions_on_sequences.csv",
        "data/processed/NetTCR_predictions_on_sequences_pMHC_split.csv",
    resources:
        runtime="5m",
        mem="1GB",
        tasks=1
    shell:
        """
        papermill notebooks/benchmark_model_performance.ipynb notebooks/benchmark_model_performance.ipynb
        """

rule run_ablate_model_notebook:
    input:
        "data/interim/sequences_pMHC_split.csv",
        "data/processed/TCRStructMap_predictions_on_sequences.csv",
        "data/processed/TCRStructMap_cdrs_peptide_predictions_on_sequences.csv",
        "data/processed/TCRStructMap_cdrs_mhc_pseudo_predictions_on_sequences.csv",
        "data/interim/peptides.txt",
        "data/interim/peptide_distances.txt"
    resources:
        runtime="5m",
        mem="1GB",
        tasks=1
    shell:
        """
        papermill notebooks/ablate_model.ipynb notebooks/ablate_model.ipynb
        """

rule run_evaluate_using_structures_vs_sequence_notebook:
    input:
        "data/interim/sequences_random_split.csv",
        "data/interim/sequences_pMHC_split.csv",
        "data/processed/TCRStructMap_predictions_on_sequences_random_split.csv",
        "data/processed/TCRStructMap_sequence_only_predictions_on_sequences_random_split.csv",
        "data/processed/TCRStructMap_predictions_on_sequences_peptide_split.csv",
        "data/processed/TCRStructMap_sequence_only_predictions_on_sequences_peptide_split.csv",
        "data/processed/TCRStructMap_predictions_on_sequences_pMHC_split.csv",
        "data/processed/TCRStructMap_sequence_only_predictions_on_sequences_pMHC_split.csv",
    resources:
        runtime="5m",
        mem="1GB"
    shell:
        """
        papermill \
            notebooks/evaluate_using_structures_vs_sequence.ipynb \
            notebooks/evaluate_using_structures_vs_sequence.ipynb
        """

rule run_evaluate_model_on_data_splits_notebook:
    input:
        "data/interim/sequences_random_split.csv",
        "data/interim/sequences_peptide_split.csv",
        "data/interim/sequences_pMHC_split.csv",
        "data/processed/TCRStructMap_predictions_on_sequences_random_split.csv",
        "data/processed/TCRStructMap_cdrs_peptide_predictions_on_sequences_random_split.csv",
        "data/processed/TCRStructMap_cdrs_mhc_pseudo_predictions_on_sequences_random_split.csv",
        "data/processed/TCRStructMap_predictions_on_sequences_peptide_split.csv",
        "data/processed/TCRStructMap_cdrs_peptide_predictions_on_sequences_peptide_split.csv",
        "data/processed/TCRStructMap_cdrs_mhc_pseudo_predictions_on_sequences_peptide_split.csv",
        "data/processed/TCRStructMap_predictions_on_sequences_pMHC_split.csv",
        "data/processed/TCRStructMap_cdrs_peptide_predictions_on_sequences_pMHC_split.csv",
        "data/processed/TCRStructMap_cdrs_mhc_pseudo_predictions_on_sequences_pMHC_split.csv",
    resources:
        runtime="5m",
        mem="1GB"
    shell:
        """
        papermill notebooks/evaluate_model_on_data_splits.ipynb notebooks/evaluate_model_on_data_splits.ipynb
        """

rule run_evaluate_distance_based_confidence_predictor_performance_notebook:
    input: "data/logs/train_confidence_predictor.log"
    resources:
        runtime="1m",
        mem="1GB",
        tasks=1
    shell:
        """
        papermill notebooks/evaluate_distance_based_confidence_predictor_performance.ipynb notebooks/evaluate_distance_based_confidence_predictor_performance.ipynb
        """

rule run_evaluate_confidence_predictions_notebook:
    input:
        "models/ConfidencePredictor",
        "models/TCRStructMap",
        "data/processed/sequences_pMHC_split.h5",
        "data/interim/sequences_pMHC_split.csv",
        "data/interim/peptides.txt",
        "data/interim/peptide_distances.txt"
    resources:
        runtime="5m",
        mem="5GB",
        tasks=1
    shell:
        """
        papermill notebooks/evaluate_confidence_predictions.ipynb notebooks/evaluate_confidence_predictions.ipynb
        """

rule run_immrep_2025_notebook:
    input:
        "data/interim/test.csv",
        expand("data/external/hla_sequences/{mhc}.json", mhc=[
            'h2_d',
            'h2_k',
            'h2_l',
            'hla_a',
            'hla_b',
            'hla_c',
            'hla_e',
            'hla_f',
            'hla_g',
        ]),
        "data/interim/mhc_pseudo_seq_imgt_positions.json",
        "models/TCRStructMap"
    output: "data/processed/immrep_2025_submission.csv"
    resources:
        runtime="10m",
        mem="5GB",
        tasks=1
    shell:
        """
        papermill notebooks/immrep_2025.ipynb notebooks/immrep_2025.ipynb
        """
