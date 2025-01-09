rule run_identify_contact_residues_notebook:
    input: "data/raw/stcrdab"
    output:
        "data/interim/tcr_mhc_class_I_contacts.csv",
        "data/interim/tcr_mhc_class_II_contacts.csv"
    resources:
        runtime="20m",
        mem="10GB",
        tasks=1
    notebook: "notebooks/Identify_contact_residues_on_MHC_molecules.ipynb"

rule run_data_summary_notebook:
    input: "data/processed/selected-stcrdab", "data/logs/select_stcrdab_tcr_pmhc_structures.log"
    resources:
        runtime="5m",
        mem="1GB",
        tasks=1
    notebook: "notebooks/data_summary.ipynb"

rule run_process_and_visualise_tcr_pmhc_contact_maps_notebook:
    input: "data/processed/tcr_pmhc_contacts.csv", "data/interim/mhc_pseudo_seq_imgt_positions.json"
    output: "data/processed/contact_maps.h5"
    resources:
        runtime="5m",
        mem="1GB",
        tasks=1
    notebook: "notebooks/process_and_visualise_tcr_pmhc_contact_maps.ipynb"

rule run_visualising_and_evaluating_tcr_contact_map_predictor_training_notebook:
    input: "data/logs/train_tcr_contact_map_predictor.log"
    resources:
        runtime="1m",
        mem="500MB",
        tasks=1
    notebook: "notebooks/visualising_and_evaluating_tcr_contact_map_predictor_training.ipynb"
