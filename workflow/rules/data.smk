rule data:
    input:
        "data/processed/selected-stcrdab",
        "data/processed/tcr_pmhc_contacts.csv",
        "data/processed/sequences.h5",
        "data/processed/contact_maps.h5"

rule download_stcrdab:
    output: directory("data/raw/stcrdab")
    resources:
        runtime="20m",
        mem="500MB",
        tasks=1
    shell: "python -m tcr_antigen_prediction.data.apps.download_stcrdab --log-level {config[log_level]} {output}"

rule select_stcrdab_structures:
    input: "data/raw/stcrdab"
    output:
        structures=directory("data/processed/selected-stcrdab"),
        summary="data/interim/structures_summary.csv"
    resources:
        runtime="1h",
        mem="500MB",
        tasks=1
    log: "data/logs/select_stcrdab_tcr_pmhc_structures.log"
    shell:
        """
        python -m tcr_antigen_prediction.data.apps.select_stcrdab_tcr_pmhc_structures \
            --log-level {config[log_level]} \
            --log-file {log} \
            --seed 123 \
            --tcr-types abTCR \
            --mhc-types MH1 MH2 \
            --antigen-types peptide \
            --crop-structures \
            --remove-het-atoms \
            --remove-structures-missing-residues \
            --fix-structures-missing-residues \
            --structural-similarity-cutoff 2.0 \
            --output-summary-csv {output.summary} \
            -o {output.structures} \
            {input}
        """

rule create_contact_maps:
    input:
        structures="data/processed/selected-stcrdab",
        summary="data/interim/structures_summary.csv"
    output: "data/processed/tcr_pmhc_contacts.csv"
    resources:
        runtime="10m",
        mem="500MB",
        tasks=1
    shell:
        """
        python -m tcr_antigen_prediction.data.apps.create_contact_maps \
            --log-level {config[log_level]} \
            --tcr-norm relative_pos_centre \
            --mhc-norm imgt_number \
            --peptide-norm relative_pos_centre \
            -o {output} \
            --summary-csv {input.summary} \
            {input.structures}/*.pdb \
        """

rule get_mhc_pseudo_sequence_imgt_numbers:
    input: "data/processed/tcr_pmhc_contacts.csv"
    output: "data/interim/mhc_pseudo_seq_imgt_positions.json"
    shell: "python -m tcr_antigen_prediction.data.apps.get_mhc_pseudo_sequence_imgt_numbers --log-level {config[log_level]} -o {output} {input}"

rule add_mhc_pseudo_sequences:
    input:
        structures="data/processed/selected-stcrdab",
        summary="data/interim/structures_summary.csv",
        mhc_pseudo_imgt="data/interim/mhc_pseudo_seq_imgt_positions.json"
    output: "data/processed/structures_summary.csv"
    resources:
        runtime="10m",
        mem="500MB",
        tasks=1
    shell:
        """
        python -m tcr_antigen_prediction.data.apps.annotate_mhc_pseudo_sequences \
            --log-level {config[log_level]} \
            -o {output} \
            --summary-csv {input.summary} \
            --mhc-pseudo-sequence-imgt-numbers {input.mhc_pseudo_imgt} \
            {input.structures}/*.pdb
        """

rule collate_sequence_data:
    input:
        iedb="data/raw/iedb.csv",
        vdjdb="data/raw/vdjdb.tsv",
        mcpas_tcr="data/raw/McPAS-TCR.csv",
        itrap="data/raw/itrap.csv",
        mhc_sequences=expand("data/external/hla_sequences/{mhc}.json", mhc=[
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
        mhc_pseudo_imgt="data/interim/mhc_pseudo_seq_imgt_positions.json"
    output: "data/interim/sequences.csv"
    log: "data/logs/collate_sequence_data.log"
    resources:
        runtime="20m",
        mem="500MB",
        tasks=1
    shell:
        """
        python -m tcr_antigen_prediction.data.apps.collate_sequence_data \
            --log-level {config[log_level]} \
            --log-file {log} \
            --mhc-sequences {input.mhc_sequences} \
            --mhc-pseudo-sequence-imgt-numbers {input.mhc_pseudo_imgt} \
            --iedb-path {input.iedb} \
            --vdjdb-path {input.vdjdb} \
            --itrap-path {input.itrap} \
            --mcpas-tcr-path {input.mcpas_tcr} \
            -o {output}
        """

rule calculate_peptide_distances:
    input: "data/interim/sequences.csv"
    output:
        distance_matrix="data/interim/peptide_distances.txt",
        peptides="data/interim/peptides.txt"
    resources:
        runtime="1m",
        mem="1GB",
        tasks=1
    shell:
        """
        cut -d"," -f7 {input} | sed 1d | sort | uniq > {output.peptides}
        python -m tcr_antigen_prediction.data.apps.compute_pw_distances \
            -o {output.distance_matrix} \
            $(cat {output.peptides} | tr '\n' ' ')
        """

rule process_sequence_data:
    input:
        sequences="data/interim/sequences.csv",
    output: "data/processed/sequences.h5"
    resources:
        runtime="20m",
        mem="2GB",
        tasks=1
    shell:
        """
        python -m tcr_antigen_prediction.data.apps.process_sequence_data \
            --log-level {config[log_level]} \
            --seed 123 \
            -o {output} \
            {input.sequences}
        """

rule process_sequence_data_for_nettcr:
    input: "data/interim/sequences.csv"
    output: "data/processed/sequences_right_pad_blosum.h5"
    resources:
        runtime="20m",
        mem="2GB",
        tasks=1
    shell:
        """
        python -m tcr_antigen_prediction.data.apps.process_sequence_data \
            --log-level {config[log_level]} \
            --seed 123 \
            --cdr1-alpha-length 7 \
            --cdr2-alpha-length 8 \
            --cdr3-alpha-length 22 \
            --cdr1-beta-length 6 \
            --cdr2-beta-length 7 \
            --cdr3-beta-length 23 \
            --peptide-length 12 \
            --pad-direction right \
            --encoding blosum50 \
            --normalisation-factor 5.0 \
            -o {output} \
            {input}
        """

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
            xargs -I % bash -c 'python -m tcr_antigen_prediction.data.apps.renumber_tcr_pmhc_structure \
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
            python -m tcr_antigen_prediction.data.apps.identify_tcr_pmhc_interactions --log-level {config[log_level]} -o "/tmp/$file_name_base.csv" $file_name
            cat "/tmp/$file_name_base.csv" | sed 1d | cut -d, -f1-5 | tr ',' ' ' \
                | xargs -I % bash -c \
                'python -m tcr_antigen_prediction.data.apps.extract_chains_from_structure --chains % -o "${{3}}/${{1}}_$(echo "%" | tr -d " ").pdb" "$2"' \
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
                python -m tcr_antigen_prediction.data.apps.crop_tcr_pmhc \
                    --log-level {config[log_level]} \
                    "{input}/$file_name" \
                    -o "{output}/$file_name" \
                    --tcr-chains $alpha_chain $beta_chain \
                    --mhc-chains $mhc_chain1 \
                    --antigen-chain $antigen_chain || continue
                echo "$(sed -n "${{line}}p" {input}/structures_summary.csv | cut -d, -f1-5),,$(sed -n "${{line}}p" {input}/structures_summary.csv | cut -d, -f7-)" >> "{output}/structures_summary.csv"
            elif [ "$mhc_type" = "MH2" ]; then
                python -m tcr_antigen_prediction.data.apps.crop_tcr_pmhc \
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
            python -m tcr_antigen_prediction.data.apps.annotate_tcr_pmhc_sequences \
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
        python -m tcr_antigen_prediction.data.apps.filter_similar_structures \
            --log-level {config[log_level]} \
            --structural-similarity-cutoff 2.0 \
            --summary-csv {input.summary_file} \
            -o "{output}/structures_summary.csv" \
            {input.data_dir}
        cat "{output}/structures_summary.csv" | sed 1d | cut -d, -f1 | xargs -I % cp {input.data_dir}/%.pdb {output}/
        """

rule download_immrep_2025_data:
    output: "data/external/immrep2025.zip"
    shell:
        """
        wget -O {output} "https://storage.googleapis.com/kaggle-competitions-data/kaggle-v2/90596/11270508/bundle/archive.zip?GoogleAccessId=web-data@kaggle-161607.iam.gserviceaccount.com&Expires=1741257712&Signature=KKcAzJCqhfISPZc2pOl6wUc02auUwRJ%2BUPRfo1o%2FSLPA1aK1PGLfJGvA7XgtoEmSS6gztZgjg4gwfJP%2BTgYiP%2B2Msyil0xcM2zRoEthK1kWEPuXXEolgWgSbwzukc%2FdOM8MBuV%2BdQl%2FoRADT45hZc4TaQpVvssczN1YChsPrrrlaQGJLxQeSxkgckmBahc%2Fi3E%2Ffg6yF8np2DqU5o07d9AucJCdtjxWzo8YXPS7d4Jy%2B1mQVQkYBWnWNlF%2F6hLfHd%2FaeosI9MOVRT4kASWJsa0iviMKHzn7OgOpGM88Op5t0Wxf1WyKWBtH8NqilhrLfOL9WJAAkYAAwIVUBayZeCA%3D%3D&response-content-disposition=attachment%3B+filename%3Dimmrep25.zip"
        """

rule uncompress_immrep_2025_data:
    input: "data/external/immrep2025.zip"
    output:
        "data/interim/iedb_positives.csv",
        "data/interim/sample_submission.csv",
        "data/interim/test.csv",
        "data/interim/vdjdb_positives.csv"
    shell: "unzip {input} -d data/interim/"
