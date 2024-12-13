default_config = {
    'log_level': 'info',
}

config = {**default_config, **config}

rule environment:
    shell:
        """
        mamba env create -f environment.yml
        mamba run -n tcr-antigen-prediction python -m pip install .
        git submodule init
        git submodule update --recursive
        mamba_prefix=$(mamba run -n tcr-antigen-prediction mamba info --json | jq '."env location"' | sed s/\\"//g)
        python_version=$(mamba run -n tcr-antigen-prediction python --version | cut -d " " -f2 | cut -d "." -f1-2)
        cp -r third_party/anarci $mamba_prefix/lib/python$python_version/site-packages
        """

rule data:
    input: "data/processed/selected-stcrdab", "data/processed/tcr_pmhc_contacts.csv"

rule download_stcrdab:
    output: directory("data/raw/stcrdab")
    resources:
        runtime="20m",
        mem="500MB",
        tasks=1
    shell: "python -m tcr_antigen_prediction.apps.download_stcrdab --log-level {config[log_level]} {output}"

rule select_stcrdab_structures:
    input:
        stcrdab_path="data/raw/stcrdab",
        tcr_mhc_class_I_contacts="data/interim/tcr_mhc_class_I_contacts.csv",
        tcr_mhc_class_II_contacts="data/interim/tcr_mhc_class_II_contacts.csv"
    output: directory("data/processed/selected-stcrdab")
    resources:
        runtime="1h",
        mem="500MB",
        tasks=1
    log: "data/logs/select_stcrdab_tcr_pmhc_structures.log"
    shell:
        """
        python -m tcr_antigen_prediction.apps.select_stcrdab_tcr_pmhc_structures \
            --log-level {config[log_level]} \
            --log-file {log} \
            --seed 123 \
            --tcr-types abTCR \
            --mhc-types MH1 MH2 \
            --antigen-types peptide \
            --mhc-class-I-tcr-contact-residues $(cut -d, -f2 {input.tcr_mhc_class_I_contacts} | sed 1d | sort | uniq | tr '\n' ' ') \
            --mhc-class-II-alpha-chain-tcr-contact-residues $(awk -F ',' '$2 == "mhc_chain1" {{ print $3 }}' {input.tcr_mhc_class_II_contacts} | sort | uniq | tr '\n' ' ') \
            --mhc-class-II-beta-chain-tcr-contact-residues $(awk -F ',' '$2 == "mhc_chain2" {{ print $3 }}' {input.tcr_mhc_class_II_contacts} | sort | uniq | tr '\n' ' ') \
            --crop-structures \
            --remove-het-atoms \
            --remove-structures-missing-residues \
            --fix-structures-missing-residues \
            --structural-similarity-cutoff 2.0 \
            -o {output} \
            {input.stcrdab_path} 2> {log}
        """

rule create_contact_maps:
    input: "data/processed/selected-stcrdab"
    output: "data/processed/tcr_pmhc_contacts.csv"
    resources:
        runtime="10m",
        mem="500MB",
        tasks=1
    shell:
        """
        python -m tcr_antigen_prediction.apps.create_contact_maps \
            --log-level {config[log_level]} \
            --tcr-norm relative_pos_centre \
            --mhc-norm imgt_number \
            --peptide-norm relative_pos_centre \
            -o {output} \
            --summary-csv {input}/stcrdab_split.csv \
            {input}/*.pdb \
        """

rule get_mhc_pseudo_sequence_imgt_numbers:
    input: "data/processed/tcr_pmhc_contacts.csv"
    output: "data/interim/mhc_pseudo_seq_imgt_positions.json"
    shell: "python -m tcr_antigen_prediction.apps.get_mhc_pseudo_sequence_imgt_numbers --log-level {config[log_level]} -o {output} {input}"

rule data_external:
    input: "data/processed/external_validation_data_selected"

rule renumber_external_structures:
    input: "data/external/ClassI_ternaries"
    output: directory("data/interim/external_validation_data_renumbered")
    resources:
        runtime="10m",
        mem="500MB",
        tasks=1
    shell:
        """
        mkdir -p "{output}"
        find "{input}" -name "*.pdb" | \
            xargs -I % bash -c 'python -m tcr_antigen_prediction.apps.renumber_tcr_pmhc_structure \
                --log-level {config[log_level]} \
                -o "{output}/$(basename "%")" "%"'
        """

rule identify_external_tcr_pmhc_interactions:
    input: "data/interim/external_validation_data_renumbered"
    output: directory("data/interim/external_validation_data_entities")
    resources:
        runtime="10m",
        mem="1GB",
        tasks=1
    shell:
        """
        mkdir -p {output}
        echo "name,Achain,Bchain,antigen_chain,mhc_chain1,mhc_chain2,mhc_type" > "{output}/structures_summary.csv"
        for file_name in "{input}"/*; do
            file_name_base=$(basename $file_name .pdb);
            python -m tcr_antigen_prediction.apps.identify_tcr_pmhc_interactions --log-level {config[log_level]} -o "/tmp/$file_name_base.csv" $file_name
            cat "/tmp/$file_name_base.csv" | sed 1d | cut -d, -f1-5 | tr ',' ' ' \
                | xargs -I % bash -c \
                'python -m tcr_antigen_prediction.apps.extract_chains_from_structure --chains % -o "${{3}}/${{1}}_$(echo "%" | tr -d " ").pdb" "$2"' \
                _ $file_name_base $file_name {output} || continue
            cat "/tmp/$file_name_base.csv" | sed 1d | xargs -I % bash -c 'echo ${{1}}_$(echo % | cut -d, -f1-5 | sed s/,//g),%' \
            _ $file_name_base \
            >> "{output}/structures_summary.csv"
        done
        """

rule crop_external_structures:
    input: "data/interim/external_validation_data_entities"
    output: directory("data/interim/external_validation_data_entities_crop")
    resources:
        runtime="5m",
        mem="100mb",
        tasks=1
    shell:
        """
        mkdir -p {output}
        head -n1 "{input}/structures_summary.csv" > "{output}/structures_summary.csv"
        num_lines=$(cat {input}/structures_summary.csv | wc -l)
        line=1
        while [ $line -lt $num_lines ]; do
            line=$(expr $line + 1)
            name=$(sed -n "${{line}}p" {input}/structures_summary.csv | cut -d, -f1)
            alpha_chain=$(sed -n "${{line}}p" {input}/structures_summary.csv | cut -d, -f2)
            beta_chain=$(sed -n "${{line}}p" {input}/structures_summary.csv | cut -d, -f3)
            antigen_chain=$(sed -n "${{line}}p" {input}/structures_summary.csv | cut -d, -f4)
            mhc_chain1=$(sed -n "${{line}}p" {input}/structures_summary.csv | cut -d, -f5)
            mhc_chain2=$(sed -n "${{line}}p" {input}/structures_summary.csv | cut -d, -f6)
            mhc_type=$(sed -n "${{line}}p" {input}/structures_summary.csv | cut -d, -f7)
            chains="${{alpha_chain}}${{beta_chain}}${{antigen_chain}}${{mhc_chain1}}${{mhc_chain2}}"
            file_name="${{name}}.pdb"
            echo "Working on $name..."
            if [ "$alpha_chain" = "" ] || [ "$beta_chain" = "" ]; then
                echo "Skipping entry without TCR"
                continue
            fi
            if [ "$mhc_type" = "MH1" ]; then
                python -m tcr_antigen_prediction.apps.crop_tcr_pmhc \
                    --log-level {config[log_level]} \
                    "{input}/$file_name" \
                    -o "{output}/$file_name" \
                    --tcr-chains $alpha_chain $beta_chain \
                    --mhc-chains $mhc_chain1 \
                    --antigen-chain $antigen_chain || continue
                echo "$(sed -n "${{line}}p" {input}/structures_summary.csv | cut -d, -f1-5),,$(sed -n "${{line}}p" {input}/structures_summary.csv | cut -d, -f7-)" >> "{output}/structures_summary.csv"
            elif [ "$mhc_type" = "MH2" ]; then
                python -m tcr_antigen_prediction.apps.crop_tcr_pmhc \
                    --log-level {config[log_level]} \
                    "{input}/$file_name" \
                    -o "{output}/$file_name" \
                    --tcr-chains $alpha_chain $beta_chain \
                    --mhc-chains $mhc_chain1 $mhc_chain2 \
                    --antigen-chain $antigen_chain || continue
                sed -n "${{line}}p" {input}/structures_summary.csv >> "{output}/structures_summary.csv"
            else
                echo "Skipping entry without MHC"
                continue
            fi
        done
        """

rule get_external_structures_sequences:
    input: "data/interim/external_validation_data_entities_crop"
    output: "data/interim/external_validation_data_annotated_sequences.csv"
    resources:
        runtime="5m",
        mem="100MB",
        tasks=1
    shell:
        """
        echo "$(head -n1 {input}/structures_summary.csv),CDR1alpha_sequence,CDR2alpha_sequence,CDR3alpha_sequence,CDR1beta_sequence,CDR2beta_sequence,CDR3beta_sequence,peptide_sequence" > {output}
        num_lines=$(cat {input}/structures_summary.csv | wc -l)
        line=1
        while [ $line -lt $num_lines ]; do
            line=$(expr $line + 1)
            name=$(sed -n "${{line}}p" {input}/structures_summary.csv | cut -d, -f1)
            alpha_chain=$(sed -n "${{line}}p" {input}/structures_summary.csv | cut -d, -f2)
            beta_chain=$(sed -n "${{line}}p" {input}/structures_summary.csv | cut -d, -f3)
            antigen_chain=$(sed -n "${{line}}p" {input}/structures_summary.csv | cut -d, -f4)
            mhc_chain1=$(sed -n "${{line}}p" {input}/structures_summary.csv | cut -d, -f5)
            mhc_chain2=$(sed -n "${{line}}p" {input}/structures_summary.csv | cut -d, -f6)
            chains="${{alpha_chain}}${{beta_chain}}${{antigen_chain}}${{mhc_chain1}}${{mhc_chain2}}"
            file_name="${{name}}.pdb"
            output_name="/tmp/$(basename $file_name .pdb).csv"
            python -m tcr_antigen_prediction.apps.annotate_tcr_pmhc_sequences \
                --log-level {config[log_level]} \
                --alpha-chain-id $alpha_chain \
                --beta-chain-id $beta_chain \
                --antigen-chain-id $antigen_chain \
                -o "$output_name" \
                "{input}/$file_name" || continue
            echo "$(sed -n "${{line}}p" {input}/structures_summary.csv),$(cat $output_name | sed 1d)" >> {output}
        done
        """

rule select_external_structures:
    input:
        data_dir="data/interim/external_validation_data_entities_crop",
        summary_file="data/interim/external_validation_data_annotated_sequences.csv"
    output: directory("data/processed/external_validation_data_selected")
    resources:
        runtime="10m",
        mem="1GB",
        tasks=1
    shell:
        """
        mkdir -p {output}
        python -m tcr_antigen_prediction.apps.filter_similar_structures \
            --log-level {config[log_level]} \
            --structural-similarity-cutoff 2.0 \
            --summary-csv {input.summary_file} \
            -o "{output}/structures_summary.csv" \
            {input.data_dir}
        cat "{output}/structures_summary.csv" | sed 1d | cut -d, -f1 | xargs -I % cp {input.data_dir}/%.pdb {output}/
        """

rule models:
    input: "data", "models/TCRen"

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
        @python -m tcr_antigen_prediction.apps.train_tcr_en \
            --log-level {config[log_level]} \
            -o {output}/TCRen_probabilities.csv \
            --summary-csv "{input}/stcrdab_split.csv" \
            $(cat "{input}/stcrdab_split.csv" | grep "train" | awk -F, -v dir="{input}" '{{ printf "%s/%s_%s%s%s%s%s.pdb ", dir, $1, $2, $3, $4, $5, $6 }}')
        """

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

rule run_visualise_tcr_pmhc_contact_maps_notebook:
    input: "data/processed/tcr_pmhc_contacts.csv"
    resources:
        runtime="5m",
        mem="1GB",
        tasks=1
    notebook: "notebooks/visualise_tcr_pmhc_contact_maps.ipynb"

rule lint:
    shell:
        """
        ruff check
        ruff format --check
        """

rule test:
    shell: "pytest tests"

rule docs:
    shell:
        """
        sphinx-apidoc -f -e -o docs/source src/tcr_antigen_prediction src/**/apps/*
        python docs/document_clis.py docs/source
        python docs/document_notebooks.py notebooks docs/source
        sphinx-build -b html ./docs ./docs/public
        """
