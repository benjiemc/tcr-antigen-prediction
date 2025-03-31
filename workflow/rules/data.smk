rule data:
    input:
        "data/processed/structures",
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
        structures=directory("data/processed/structures"),
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
        structures="data/processed/structures",
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
        structures="data/processed/structures",
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

rule process_contact_maps:
    input:
        contacts="data/processed/tcr_pmhc_contacts.csv",
        mhc_pseudo_imgt="data/interim/mhc_pseudo_seq_imgt_positions.json"
    output: "data/processed/contact_maps.h5"
    resources:
        runtime="5m",
        mem="100MB",
        tasks=1
    shell:
        """
        python -m tcr_antigen_prediction.data.apps.process_contact_maps \
            --log-level {config[log_level]} \
            --mhc-types MH1 \
            --pad \
            -o {output} \
            --mhc-pseudo-sequence-imgt-numbers {input.mhc_pseudo_imgt} \
            {input.contacts}
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

rule calculate_levenshtein_distances:
    input: "data/interim/sequences.csv"
    output:
        distance_matrix="data/interim/{name}_distances.txt",
        names="data/interim/{name}s.txt"
    resources:
        runtime="1h",
        mem="1GB",
        tasks=1
    shell:
        """
        column_number=$(head -n1 {input} | tr ',' '\n' | nl | grep {wildcards.name} | cut -f1 | xargs)
        cut -d"," -f $column_number {input} | sed 1d | sort | uniq > {output.names}
        python -m tcr_antigen_prediction.data.apps.compute_pw_distances \
            -o {output.distance_matrix} \
            $(cat {output.names} | tr '\n' ' ')
        """

rule aggregate_levenshtein_distances:
    input:
        distance_matrices=expand("data/interim/{name}_distances.txt", name=[
            'cdr1_alpha',
            'cdr2_alpha',
            'cdr3_alpha',
            'cdr1_beta',
            'cdr2_beta',
            'cdr3_beta',
            'peptide',
            'mhc_pseudo',
        ]),
        names=expand("data/interim/{name}s.txt", name=[
            'cdr1_alpha',
            'cdr2_alpha',
            'cdr3_alpha',
            'cdr1_beta',
            'cdr2_beta',
            'cdr3_beta',
            'peptide',
            'mhc_pseudo',
        ])
    output: "data/interim/distance_matrices.h5"
    resources:
        runtime="5m",
        mem="1GB",
        tasks=1
    run:
        import os

        import h5py
        import numpy as np

        with h5py.File(output[0], 'w') as hdf5:
            for names, distance_matrix in zip(input.names, input.distance_matrices, strict=True):
                group = hdf5.create_group(os.path.basename(names).replace('s.txt', ''))

                with open(names, 'r') as fh:
                    group['names'] = [line.strip() for line in fh.readlines() if line]

                group['distance_matrix'] = np.loadtxt(distance_matrix, dtype=int)

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
            --seed {config[seed]} \
            -o {output} \
            {input.sequences}
        """

rule process_sequence_data_pmhc_split:
    input:
        sequences="data/interim/sequences.csv",
    output: "data/processed/sequences_pmhc_split.h5"
    resources:
        runtime="20m",
        mem="2GB",
        tasks=1
    shell:
        """
        python -m tcr_antigen_prediction.data.apps.process_sequence_data \
            --log-level {config[log_level]} \
            --seed {config[seed]} \
            --split-type pMHC \
            -o {output} \
            {input.sequences}
        """

rule process_sequence_data_tcr_split:
    input:
        sequences="data/interim/sequences.csv",
    output: "data/processed/sequences_tcr_split.h5"
    resources:
        runtime="20m",
        mem="2GB",
        tasks=1
    shell:
        """
        python -m tcr_antigen_prediction.data.apps.process_sequence_data \
            --log-level {config[log_level]} \
            --seed 123 \
            --split-type tcr \
            -o {output} \
            {input.sequences}
        """

rule process_sequence_data_tcr_levenshtein_split:
    input:
        sequences="data/interim/sequences.csv",
        distances="data/interim/distance_matrices.h5"
    output: "data/processed/sequences_tcr_levenshtein_split.h5"
    resources:
        runtime="20m",
        mem="275GB",
        tasks=1
    shell:
        """
        python -m tcr_antigen_prediction.data.apps.process_sequence_data \
            --log-level {config[log_level]} \
            --seed 123 \
            --split-type levenshtein \
            --split-distance 6 \
            --split-entities cdr1_alpha cdr2_alpha cdr3_alpha cdr1_beta cdr2_beta cdr3_beta \
            --distances {input.distances} \
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
            --seed {config[seed]} \
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
    input: "data/processed/external_structures_summary.csv",
           "data/processed/external_structures"

rule flatten_and_sanitize_external_structures_file_paths:
    input: "data/external/ClassI_ternaries"
    output: directory("data/interim/external_structures_sanitized")
    resources:
        runtime="1m",
        mem="500MB",
        tasks=1
    shell:
        """
        mkdir -p {output}
        find "{input}" -name "*.pdb" \
            | xargs -I % bash -c 'mv "%" {output}/$(basename "%" .pdb | tr '.' '_' | tr ',' '_' | tr '-' '_').pdb'
        """

rule renumber_external_structures:
    input: "data/interim/external_structures_sanitized"
    output: directory("data/interim/external_structures_renumbered")
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
    input: "data/interim/external_structures_renumbered"
    output: "data/interim/external_structures_summary.csv"
    resources:
        runtime="10m",
        mem="1GB",
        tasks=1
    shell:
        """
        echo "path,Achain,Bchain,antigen_chain,mhc_chain1,mhc_chain2,mhc_type" > {output}
        for file_name in "{input}"/*; do
            file_name_base=$(basename $file_name .pdb);
            python -m tcr_antigen_prediction.data.apps.identify_tcr_pmhc_interactions --log-level {config[log_level]} $file_name \
                | sed 1d \
                | tr '\t' ',' \
                | xargs -I % bash -c 'echo ${{1}}_$(echo % | cut -d, -f1-5 | sed s/,//g).pdb,%' _ $file_name_base \
                > /tmp/${{file_name_base}}.csv
            num_lines=$(cat /tmp/${{file_name_base}}.csv | wc -l)
            line=1
            while [ $line -le $num_lines ]; do
                alpha_chain=$(sed -n "${{line}}p" /tmp/${{file_name_base}}.csv | cut -d, -f2)
                beta_chain=$(sed -n "${{line}}p" /tmp/${{file_name_base}}.csv | cut -d, -f3)
                antigen_chain=$(sed -n "${{line}}p" /tmp/${{file_name_base}}.csv | cut -d, -f4)
                mhc_type=$(sed -n "${{line}}p" /tmp/${{file_name_base}}.csv | cut -d, -f7)
                if [ "$alpha_chain" != "" ] && [ "$beta_chain" != "" ] && [ "$antigen_chain" != "" ] && [ "$mhc_type" != "" ]; then
                   sed -n "${{line}}p" /tmp/${{file_name_base}}.csv >> {output}
                fi
                line=$(expr $line + 1)
            done
        done
        """

rule isolate_external_tcr_pmhc_complexes:
    input:
        summary="data/interim/external_structures_summary.csv",
        structures="data/interim/external_structures_renumbered"
    output: directory("data/interim/external_structures_entities")
    resources:
        runtime="10m",
        mem="1GB",
        tasks=1
    shell:
        """
        mkdir -p {output}
        cat "{input.summary}" | sed 1d | cut -d, -f1-6 \
            | xargs -I % bash -c \
                'python -m tcr_antigen_prediction.data.apps.extract_chains_from_structure \
                    --chains $(echo "%" | cut -d, -f2-6 | tr "," " ") \
                    -o {output}/$(echo "%" | cut -d, -f1) {input.structures}/$(echo "%" | cut -d, -f1 | rev | cut -d_ -f2- | rev ).pdb'
        """

rule crop_external_structures:
    input:
        summary="data/interim/external_structures_summary.csv",
        structures="data/interim/external_structures_entities"
    output:
        summary="data/interim/external_structures_summary_cropped.csv",
        structures=directory("data/interim/external_structures_crop")
    resources:
        runtime="5m",
        mem="100mb",
        tasks=1
    shell:
        """
        mkdir -p {output.structures}
        head -n1 {input.summary} > {output.summary}
        num_lines=$(cat {input.summary} | wc -l)
        line=1
        while [ $line -lt $num_lines ]; do
            line=$(expr $line + 1)
            file_name=$(sed -n "${{line}}p" {input.summary} | cut -d, -f1)
            alpha_chain=$(sed -n "${{line}}p" {input.summary} | cut -d, -f2)
            beta_chain=$(sed -n "${{line}}p" {input.summary} | cut -d, -f3)
            antigen_chain=$(sed -n "${{line}}p" {input.summary} | cut -d, -f4)
            mhc_chain1=$(sed -n "${{line}}p" {input.summary} | cut -d, -f5)
            mhc_chain2=$(sed -n "${{line}}p" {input.summary} | cut -d, -f6)
            mhc_type=$(sed -n "${{line}}p" {input.summary} | cut -d, -f7)
            chains="${{alpha_chain}}${{beta_chain}}${{antigen_chain}}${{mhc_chain1}}${{mhc_chain2}}"
            echo "Working on $file_name..."
            if [ "$mhc_type" = "MH1" ]; then
                python -m tcr_antigen_prediction.data.apps.crop_tcr_pmhc \
                    --log-level {config[log_level]} \
                    "{input.structures}/$file_name" \
                    -o "{output.structures}/$file_name" \
                    --tcr-chains $alpha_chain $beta_chain \
                    --mhc-chains $mhc_chain1 \
                    --antigen-chain $antigen_chain
                echo $(sed -n "${{line}}p" {input.summary} \
                    | cut -d, -f1-5),,$(sed -n "${{line}}p" {input.summary} \
                    | cut -d, -f7-) >> {output.summary}
            elif [ "$mhc_type" = "MH2" ]; then
                python -m tcr_antigen_prediction.data.apps.crop_tcr_pmhc \
                    --log-level {config[log_level]} \
                    "{input.structures}/$file_name" \
                    -o "{output.structures}/$file_name" \
                    --tcr-chains $alpha_chain $beta_chain \
                    --mhc-chains $mhc_chain1 $mhc_chain2 \
                    --antigen-chain $antigen_chain
                sed -n "${{line}}p" {input.summary} >> {output.summary}
            fi
        done
        """

rule get_external_structures_sequences:
    input:
        summary="data/interim/external_structures_summary_cropped.csv",
        structures="data/interim/external_structures_crop",
        mhc_pseudo_imgt="data/interim/mhc_pseudo_seq_imgt_positions.json"
    output: "data/interim/external_structures_summary_annotated_sequences.csv"
    resources:
        runtime="5m",
        mem="100MB",
        tasks=1
    shell:
        """
        echo "$(head -n1 {input.summary}),cdr1_alpha,cdr2_alpha,cdr3_alpha,cdr1_beta,cdr2_beta,cdr3_beta,peptide" > {output}
        num_lines=$(cat {input.summary} | wc -l)
        line=1
        while [ $line -lt $num_lines ]; do
            line=$(expr $line + 1)
            file_name=$(sed -n "${{line}}p" {input.summary} | cut -d, -f1)
            alpha_chain=$(sed -n "${{line}}p" {input.summary} | cut -d, -f2)
            beta_chain=$(sed -n "${{line}}p" {input.summary} | cut -d, -f3)
            antigen_chain=$(sed -n "${{line}}p" {input.summary} | cut -d, -f4)
            output_name="/tmp/$(basename $file_name .pdb).csv"
            python -m tcr_antigen_prediction.data.apps.annotate_tcr_pmhc_sequences \
                --log-level {config[log_level]} \
                --alpha-chain-id $alpha_chain \
                --beta-chain-id $beta_chain \
                --antigen-chain-id $antigen_chain \
                -o "$output_name" \
                {input.structures}/$file_name
            echo $(sed -n "${{line}}p" {input.summary}),$(cat $output_name | sed 1d) >> {output}
        done

        python -m tcr_antigen_prediction.data.apps.annotate_mhc_pseudo_sequences \
            --log-level {config[log_level]} \
            --summary-csv {output} \
            --mhc-pseudo-sequence-imgt-numbers {input.mhc_pseudo_imgt} \
            -o {output} \
            {input.structures}/*.pdb
        """

rule select_external_structures:
    input:
        data_dir="data/interim/external_structures_crop",
        summary_file="data/interim/external_structures_summary_annotated_sequences.csv"
    output:
        summary="data/processed/external_structures_summary.csv",
        structures=directory("data/processed/external_structures")
    resources:
        runtime="10m",
        mem="1GB",
        tasks=1
    shell:
        """
        mkdir -p {output.structures}
        python -m tcr_antigen_prediction.data.apps.filter_similar_structures \
            --log-level {config[log_level]} \
            --structural-similarity-cutoff 2.0 \
            --summary-csv {input.summary_file} \
            -o {output.summary} \
            {input.data_dir}
        cat {output.summary} | sed 1d | cut -d, -f1 | xargs -I % cp {input.data_dir}/% {output.structures}/
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
