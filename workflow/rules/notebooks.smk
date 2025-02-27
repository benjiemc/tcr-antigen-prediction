rule run_identify_contact_residues_notebook:
    input: "data/raw/stcrdab"
    output:
        "data/interim/tcr_mhc_class_I_contacts.csv",
        "data/interim/tcr_mhc_class_II_contacts.csv"
    resources:
        runtime="20m",
        mem="10GB",
        tasks=1
    notebook: "../../notebooks/Identify_contact_residues_on_MHC_molecules.ipynb"

rule run_visualise_structure_data_notebook:
    input: "data/interim/structures_summary.csv", "data/logs/select_stcrdab_tcr_pmhc_structures.log"
    resources:
        runtime="5m",
        mem="1GB",
        tasks=1
    notebook: "../../notebooks/visualise_structure_data.ipynb"

rule run_visualise_sequence_data_notebook:
    input: "data/interim/sequences.csv", "data/logs/collate_sequence_data.log"
    resources:
        runtime="5m",
        mem="1GB",
        tasks=1
    notebook: "../../notebooks/visualise_sequence_data.ipynb"

rule run_process_and_visualise_tcr_pmhc_contact_maps_notebook:
    input: "data/processed/tcr_pmhc_contacts.csv", "data/interim/mhc_pseudo_seq_imgt_positions.json"
    output: "data/processed/contact_maps.h5"
    resources:
        runtime="5m",
        mem="1GB",
        tasks=1
    notebook: "../../notebooks/process_and_visualise_tcr_pmhc_contact_maps.ipynb"

rule run_visualising_and_evaluating_tcr_struct_map_training_notebook:
    input: "data/logs/train_tcr_struct_map.log"
    resources:
        runtime="1m",
        mem="500MB",
        tasks=1
    notebook: "../../notebooks/visualising_and_evaluating_tcr_struct_map_training.ipynb"

rule run_visualising_and_evaluating_tcr_struct_map_sequence_only_training_notebook:
    input: "data/logs/train_tcr_struct_map_sequence_only.log"
    resources:
        runtime="1m",
        mem="500MB",
        tasks=1
    notebook: "../../notebooks/visualising_and_evaluating_tcr_struct_map_sequence_only_training.ipynb"

rule run_visualising_and_evaluating_nettcr_training_notebook:
    input: "data/logs/train_NetTCR.log"
    resources:
        runtime="1m",
        mem="500MB",
        tasks=1
    notebook: "../../notebooks/visualising_and_evaluating_nettcr_training.ipynb"

rule run_benchmark_model_performance_notebook:
    input:
        "data/interim/peptides.txt",
        "data/interim/peptide_distances.txt",
        "data/processed/sequences.h5",
        "models/TCRStructMap",
        "models/TCRStructMap_sequence_only"
    resources:
        runtime="5m",
        mem="1GB",
        tasks=1
    notebook: "../../notebooks/benchmark_model_performance.ipynb"

rule run_evaluate_distance_based_confidence_predictor_performance_notebook:
    input: "data/logs/train_confidence_predictor.log"
    resources:
        runtime="1m",
        mem="1GB",
        tasks=1
    notebook: "../../notebooks/evaluate_distance_based_confidence_predictor_performance.ipynb"

rule run_evaluate_confidence_predictions_notebook:
    input:
        "models/ConfidencePredictor",
        "models/TCRStructMap",
        "data/processed/sequences.h5",
        "data/interim/peptides.txt",
        "data/interim/peptide_distances.txt"
    resources:
        runtime="5m",
        mem="5GB",
        tasks=1
    notebook: "../../notebooks/evaluate_confidence_predictions.ipynb"

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
    notebook: "../../notebooks/immrep_2025.ipynb"
