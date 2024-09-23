.PHONY: all data models lint test docs


all: data models

data: \
	data/processed/selected-stcrdab_feats \
	data/processed/selected-stcrdab_ply

data/raw/stcrdab:
	@python -m tcr_antigen_prediction.apps.download_stcrdab $@
	@touch $@

data/interim/selected-stcrdab: data/raw/stcrdab data/external/masif_ppi_search_training_set.txt
	@python -m tcr_antigen_prediction.apps.select_stcrdab_tcr_pmhc_structures \
		--seed 123 \
		--tcr-types abTCR \
		--mhc-types MH1 MH2 \
		--antigen-types peptide \
		--remove-structures-missing-residues \
		--structural-similarity-cutoff 2.0 \
		--pdb-ids-to-exclude $$(cut -d _ -f 1 $(word 2,$^) | tr '[:upper:]' '[:lower:]' | sort | uniq | tr '\n' ' ') \
		-o $@ \
		$(word 1,$^)
	@touch $@

data/interim/selected-stcrdab_crop: data/interim/selected-stcrdab
	@mkdir -p $@
	@echo "Processing structures..."
	@head -n1 "$^/stcrdab_split.csv" > "$@/stcrdab_split.csv"
	@num_lines=$$(cat "$^/stcrdab_split.csv" | wc -l); \
	line=1; \
	while [ $$line -lt $$num_lines ]; do \
		line=$$(expr $$line + 1); \
		pdb_id=$$(sed -n "$${line}p" $^/stcrdab_split.csv | cut -d ',' -f1); \
		alpha_chain=$$(sed -n "$${line}p" $^/stcrdab_split.csv | cut -d ',' -f2); \
		beta_chain=$$(sed -n "$${line}p" $^/stcrdab_split.csv | cut -d ',' -f3); \
		antigen_chain=$$(sed -n "$${line}p" $^/stcrdab_split.csv | cut -d ',' -f4); \
		mhc_chain1=$$(sed -n "$${line}p" $^/stcrdab_split.csv | cut -d ',' -f5); \
		mhc_chain2=$$(sed -n "$${line}p" $^/stcrdab_split.csv | cut -d ',' -f6); \
		mhc_type=$$(sed -n "$${line}p" $^/stcrdab_split.csv | cut -d ',' -f7); \
		chains="$${alpha_chain}$${beta_chain}$${antigen_chain}$${mhc_chain1}$${mhc_chain2}"; \
		file_name="$${pdb_id}_$${chains}.pdb"; \
		echo "Working on $${file_name}..."; \
		if [ "$${mhc_type}" = "MH1" ]; then \
			python -m tcr_antigen_prediction.apps.crop_tcr_pmhc \
				"$^/$${file_name}" \
				-o "$@/$${file_name}" \
				--tcr-chains $$alpha_chain $$beta_chain \
				--mhc-chains $$mhc_chain1 \
				--antigen-chain $$antigen_chain || continue; \
		elif [ "$${mhc_type}" = "MH2" ]; then \
			python -m tcr_antigen_prediction.apps.crop_tcr_pmhc \
				"$^/$${file_name}" \
				-o "$@/$${file_name}" \
				--tcr-chains $$alpha_chain $$beta_chain \
				--mhc-chains $$mhc_chain1 $$mhc_chain2 \
				--antigen-chain $$antigen_chain || continue; \
		else \
			echo "INVALID MHC TYPE"; \
			continue; \
		fi; \
		sed -n "$${line}p" "$^/stcrdab_split.csv" >> "$@/stcrdab_split.csv"; \
		echo "Finished Structure"; \
	done
	@echo "All done."
	@touch $@

data/processed/selected-stcrdab_ply: data/interim/selected-stcrdab_crop
	@mkdir -p $@
	@head -n1 "$^/stcrdab_split.csv" > "$@/stcrdab_split.csv"
	@num_lines=$$(cat "$^/stcrdab_split.csv" | wc -l); \
	line=1; \
	while [ $$line -lt $$num_lines ]; do \
	    line=$$(expr $$line + 1); \
	    pdb_id=$$(sed -n "$${line}p" "$^/stcrdab_split.csv" | cut -d ',' -f1); \
	    alpha_chain=$$(sed -n "$${line}p" "$^/stcrdab_split.csv" | cut -d ',' -f2); \
	    beta_chain=$$(sed -n "$${line}p" "$^/stcrdab_split.csv" | cut -d ',' -f3); \
	    antigen_chain=$$(sed -n "$${line}p" "$^/stcrdab_split.csv" | cut -d ',' -f4); \
	    mhc_chain1=$$(sed -n "$${line}p" "$^/stcrdab_split.csv" | cut -d ',' -f5); \
	    mhc_chain2=$$(sed -n "$${line}p" "$^/stcrdab_split.csv" | cut -d ',' -f6); \
	    chains="$${alpha_chain}$${beta_chain}$${antigen_chain}$${mhc_chain1}$${mhc_chain2}"; \
	    echo "Working on $${pdb_id}_$${chains}..."; \
	    echo "Computing TCR (chains $$alpha_chain and $$beta_chain)"; \
	    python -m tcr_antigen_prediction.apps.prepare_structure "$^/$${pdb_id}_$${chains}.pdb" -o "$@" --chains $$alpha_chain $$beta_chain || continue; \
	    echo "Computing pMHC (chains $$antigen_chain, $$mhc_chain1, and $$mhc_chain2)"; \
	    python -m tcr_antigen_prediction.apps.prepare_structure "$^/$${pdb_id}_$${chains}.pdb" -o "$@" --chains $$antigen_chain $$mhc_chain1 $$mhc_chain2 || continue; \
	    sed -n "$${line}p" "$^/stcrdab_split.csv" >> "$@/stcrdab_split.csv"; \
	    echo "Finished Structure"; \
	done
	@echo "All done."
	@touch $@

data/processed/selected-stcrdab_feats: data/processed/selected-stcrdab_ply
	@mkdir -p $@
	@head -n1 "$^/stcrdab_split.csv" > "$@/stcrdab_split.csv"
	@num_lines=$$(sed -n '$$=' "$^/stcrdab_split.csv"); \
	line=1; \
	while [ $$line -lt $$num_lines ]; do \
	    line=$$(expr $$line + 1); \
	    pdb_id=$$(sed -n "$${line}p" "$^/stcrdab_split.csv" | cut -d ',' -f1); \
	    alpha_chain=$$(sed -n "$${line}p" "$^/stcrdab_split.csv" | cut -d ',' -f2); \
	    beta_chain=$$(sed -n "$${line}p" "$^/stcrdab_split.csv" | cut -d ',' -f3); \
	    antigen_chain=$$(sed -n "$${line}p" "$^/stcrdab_split.csv" | cut -d ',' -f4); \
	    mhc_chain1=$$(sed -n "$${line}p" "$^/stcrdab_split.csv" | cut -d ',' -f5); \
	    mhc_chain2=$$(sed -n "$${line}p" "$^/stcrdab_split.csv" | cut -d ',' -f6); \
	    chains="$${alpha_chain}$${beta_chain}$${antigen_chain}$${mhc_chain1}$${mhc_chain2}"; \
	    echo "Working on $$pdb_id chains $$chains..."; \
	    output_name="$@/$${pdb_id}_$${chains}"; \
		mkdir -p $$output_name; \
	    python -m tcr_antigen_prediction.apps.compute_features \
			--mode ppi_search \
			-o "$$output_name" \
			"$^/$${pdb_id}_$${chains}_$${alpha_chain}$${beta_chain}.ply" \
			"$^/$${pdb_id}_$${chains}_$${antigen_chain}$${mhc_chain1}$${mhc_chain2}.ply" || continue; \
	    sed -n "$${line}p" "$^/stcrdab_split.csv" >> "$@/stcrdab_split.csv"; \
	    echo "Finished pair"; \
	done
	@echo "All done."
	@touch $@

data/external/masif_ppi_search_training_set.txt:
	@wget -O $@ https://raw.githubusercontent.com/LPDI-EPFL/masif/master/data/masif_ppi_search/lists/training.txt

data/interim/external_validation_data_renumbered: data/external/ClassI_ternaries
	@mkdir -p "$@"
	@find "$^" -name "*.pdb" | xargs -I % bash -c 'python -m tcr_antigen_prediction.apps.renumber_tcr_pmhc_structure -o "$@/$$(basename "%")" "%"'

data/interim/external_validation_data_fix: data/interim/external_validation_data_renumbered
	@mkdir -p $@
	@for file_name in "$^"/*.pdb; do \
		file_name_base=$$(basename $$file_name .pdb); \
		python -m tcr_antigen_prediction.apps.fix_pdb -o "$@/$$(echo $$file_name_base | tr '.' '_').pdb" "$$file_name"; \
	done

data/interim/external_validation_data_entities: data/interim/external_validation_data_fix
	@mkdir -p $@
	@echo "name,Achain,Bchain,antigen_chain,mhc_chain1,mhc_chain2,mhc_type" > "$@/structures_summary.csv"
	@for file_name in "$^"/*; do \
		file_name_base=$$(basename $$file_name .pdb); \
		python -m tcr_antigen_prediction.apps.identify_tcr_pmhc_interactions -o "/tmp/$$file_name_base.csv" $$file_name; \
		cat "/tmp/$$file_name_base.csv" | sed 1d | cut -d, -f1-5 | tr ',' ' ' \
			| xargs -I % bash -c \
			'python -m tcr_antigen_prediction.apps.extract_chains_from_structure --chains % -o "$${3}/$${1}_$$(echo "%" | tr -d " ").pdb" "$$2"' \
			_ $$file_name_base $$file_name $@ || continue; \
		cat "/tmp/$$file_name_base.csv" | sed 1d | xargs -I % bash -c 'echo $${1}_$$(echo % | cut -d, -f1-5 | sed s/,//g),%' \
		_ $$file_name_base \
		>> "$@/structures_summary.csv"; \
	done

data/interim/external_validation_data_entities_crop: data/interim/external_validation_data_entities
	@mkdir -p $@
	@head -n1 "$^/structures_summary.csv" > "$@/structures_summary.csv"
	@num_lines=$$(cat $^/structures_summary.csv | wc -l); \
	line=1; \
	while [ $$line -lt $$num_lines ]; do \
		line=$$(expr $$line + 1); \
		name=$$(sed -n "$${line}p" $^/structures_summary.csv | cut -d, -f1); \
		alpha_chain=$$(sed -n "$${line}p" $^/structures_summary.csv | cut -d, -f2); \
		beta_chain=$$(sed -n "$${line}p" $^/structures_summary.csv | cut -d, -f3); \
		antigen_chain=$$(sed -n "$${line}p" $^/structures_summary.csv | cut -d, -f4); \
		mhc_chain1=$$(sed -n "$${line}p" $^/structures_summary.csv | cut -d, -f5); \
		mhc_chain2=$$(sed -n "$${line}p" $^/structures_summary.csv | cut -d, -f6); \
		mhc_type=$$(sed -n "$${line}p" $^/structures_summary.csv | cut -d, -f7); \
		chains="$${alpha_chain}$${beta_chain}$${antigen_chain}$${mhc_chain1}$${mhc_chain2}"; \
		file_name="$${name}.pdb"; \
		if [ "$$mhc_type" = "MH1" ]; then \
			python -m tcr_antigen_prediction.apps.crop_tcr_pmhc \
				"$^/$$file_name" \
				-o "$@/$$file_name" \
				--tcr-chains $$alpha_chain $$beta_chain \
				--mhc-chains $$mhc_chain1 \
				--antigen-chain $$antigen_chain || continue; \
				echo "$$(sed -n "$${line}p" $^/structures_summary.csv | cut -d, -f1-5),,$$(sed -n "$${line}p" $^/structures_summary.csv | cut -d, -f7-)" >> "$@/structures_summary.csv"; \
		elif [ "$$mhc_type" = "MH2" ]; then \
			python -m tcr_antigen_prediction.apps.crop_tcr_pmhc \
				"$^/$$file_name" \
				-o "$@/$$file_name" \
				--tcr-chains $$alpha_chain $$beta_chain \
				--mhc-chains $$mhc_chain1 $$mhc_chain2 \
				--antigen-chain $$antigen_chain || continue; \
				sed -n "$${line}p" $^/structures_summary.csv >> "$@/structures_summary.csv"; \
		else \
			continue; \
		fi; \
	done

data/interim/external_validation_data_annotated_sequences.csv: data/interim/external_validation_data_entities_crop
	@echo "$$(head -n1 $^/structures_summary.csv),CDR1alpha_sequence,CDR2alpha_sequence,CDR3alpha_sequence,CDR1beta_sequence,CDR2beta_sequence,CDR3beta_sequence,peptide_sequence" > $@
	@num_lines=$$(cat $^/structures_summary.csv | wc -l); \
	line=1; \
	while [ $$line -lt $$num_lines ]; do \
		line=$$(expr $$line + 1); \
		name=$$(sed -n "$${line}p" $^/structures_summary.csv | cut -d, -f1); \
		alpha_chain=$$(sed -n "$${line}p" $^/structures_summary.csv | cut -d, -f2); \
		beta_chain=$$(sed -n "$${line}p" $^/structures_summary.csv | cut -d, -f3); \
		antigen_chain=$$(sed -n "$${line}p" $^/structures_summary.csv | cut -d, -f4); \
		mhc_chain1=$$(sed -n "$${line}p" $^/structures_summary.csv | cut -d, -f5); \
		mhc_chain2=$$(sed -n "$${line}p" $^/structures_summary.csv | cut -d, -f6); \
		chains="$${alpha_chain}$${beta_chain}$${antigen_chain}$${mhc_chain1}$${mhc_chain2}"; \
		file_name="$${name}.pdb"; \
		output_name="/tmp/$$(basename $$file_name .pdb).csv"; \
		python -m tcr_antigen_prediction.apps.annotate_tcr_pmhc_sequences \
			--alpha-chain-id $$alpha_chain \
			--beta-chain-id $$beta_chain \
			--antigen-chain-id $$antigen_chain \
			-o "$$output_name" \
			"$^/$$file_name" || continue; \
		echo "$$(sed -n "$${line}p" $^/structures_summary.csv),$$(cat $$output_name | sed 1d)" >> $@; \
	done

data/interim/external_validation_data_selected: data/interim/external_validation_data_entities_crop data/interim/external_validation_data_annotated_sequences.csv
	@mkdir -p $@
	@python -m tcr_antigen_prediction.apps.filter_similar_structures -o "$@/structures_summary.csv" --structural-similarity-cutoff 2.0 --summary-csv $(word 2,$^) $(word 1,$^)
	@cat "$@/structures_summary.csv" | sed 1d | cut -d, -f1 | xargs -I % cp $(word 1,$^)/%.pdb $@/

data/processed/external_validation_data_selected_ply: data/interim/external_validation_data_selected
	@mkdir -p $@
	@head -n1 "$^/structures_summary.csv" > "$@/structures_summary.csv"
	@num_lines=$$(cat "$^/structures_summary.csv" | wc -l); \
	line=1; \
	while [ $$line -lt $$num_lines ]; do \
	    line=$$(expr $$line + 1); \
	    name=$$(sed -n "$${line}p" "$^/structures_summary.csv" | cut -d ',' -f1); \
	    alpha_chain=$$(sed -n "$${line}p" "$^/structures_summary.csv" | cut -d ',' -f2); \
	    beta_chain=$$(sed -n "$${line}p" "$^/structures_summary.csv" | cut -d ',' -f3); \
	    antigen_chain=$$(sed -n "$${line}p" "$^/structures_summary.csv" | cut -d ',' -f4); \
	    mhc_chain1=$$(sed -n "$${line}p" "$^/structures_summary.csv" | cut -d ',' -f5); \
	    mhc_chain2=$$(sed -n "$${line}p" "$^/structures_summary.csv" | cut -d ',' -f6); \
		mhc_type=$$(sed -n "$${line}p" "$^/structures_summary.csv" | cut -d ',' -f7); \
		echo "Working on $$name"; \
	    echo "Computing TCR (chains $$alpha_chain and $$beta_chain)"; \
	    python -m tcr_antigen_prediction.apps.prepare_structure "$^/$$name.pdb" -o "$@" --chains $$alpha_chain $$beta_chain || continue; \
		if [ "$$mhc_type" = "MH1" ]; then \
			echo "Computing pMHC-I (chains $$antigen_chain and $$mhc_chain1)"; \
			python -m tcr_antigen_prediction.apps.prepare_structure "$^/$$name.pdb" -o "$@" --chains $$antigen_chain $$mhc_chain1 || continue; \
		elif [ "$$mhc_type" = "MH2" ]; then \
			echo "Computing pMHC-II (chains $$antigen_chain, $$mhc_chain1, and $$mhc_chain2)"; \
			python -m tcr_antigen_prediction.apps.prepare_structure "$^/$$name.pdb" -o "$@" --chains $$antigen_chain $$mhc_chain1 $$mhc_chain2 || continue; \
		else \
			continue; \
		fi; \
	    sed -n "$${line}p" "$^/structures_summary.csv" >> "$@/structures_summary.csv"; \
	    echo "Finished Structure"; \
	done
	@echo "All done."

data/processed/external_validation_data_selected_feats: data/processed/external_validation_data_selected_ply
	@mkdir -p $@
	@head -n1 "$^/structures_summary.csv" > "$@/structures_summary.csv"
	@num_lines=$$(sed -n '$$=' "$^/structures_summary.csv"); \
	line=1; \
	while [ $$line -lt $$num_lines ]; do \
	    line=$$(expr $$line + 1); \
	    name=$$(sed -n "$${line}p" "$^/structures_summary.csv" | cut -d ',' -f1); \
	    alpha_chain=$$(sed -n "$${line}p" "$^/structures_summary.csv" | cut -d ',' -f2); \
	    beta_chain=$$(sed -n "$${line}p" "$^/structures_summary.csv" | cut -d ',' -f3); \
	    antigen_chain=$$(sed -n "$${line}p" "$^/structures_summary.csv" | cut -d ',' -f4); \
	    mhc_chain1=$$(sed -n "$${line}p" "$^/structures_summary.csv" | cut -d ',' -f5); \
	    mhc_chain2=$$(sed -n "$${line}p" "$^/structures_summary.csv" | cut -d ',' -f6); \
		mhc_type=$$(sed -n "$${line}p" "$^/structures_summary.csv" | cut -d ',' -f7); \
	    echo "Working on $$name..."; \
	    output_name="$@/$$name"; \
		mkdir -p $$output_name; \
		if [ "$$mhc_type" = "MH1" ]; then \
			python -m tcr_antigen_prediction.apps.compute_features \
				--mode ppi_search \
				-o "$$output_name" \
				"$^/$${name}_$${alpha_chain}$${beta_chain}.ply" \
				"$^/$${name}_$${antigen_chain}$${mhc_chain1}.ply" || continue; \
		elif [ "$$mhc_type" = "MH2" ]; then \
			python -m tcr_antigen_prediction.apps.compute_features \
				--mode ppi_search \
				-o "$$output_name" \
				"$^/$${name}_$${alpha_chain}$${beta_chain}.ply" \
				"$^/$${name}_$${antigen_chain}$${mhc_chain1}$${mhc_chain2}.ply" || continue; \
		else \
			continue; \
		fi; \
	    sed -n "$${line}p" "$^/structures_summary.csv" >> "$@/structures_summary.csv"; \
	    echo "Finished pair"; \
	done
	@echo "All done."

models: data models/baseline_masif_ppi models/finetune_masif_ppi

models/finetune_masif_ppi: data/processed/selected-stcrdab_feats data/processed/selected-stcrdab_ply  models/baseline_masif_ppi
	@python -m tcr_antigen_prediction.apps.train_masif \
		-o $@ \
		--validation-data $$(cat "$(word 1,$^)/stcrdab_split.csv" | grep "validation" | awk -F, -v dir="$(word 1,$^)" '{ printf "%s/%s_%s%s%s%s%s ", dir, $$1, $$2, $$3, $$4, $$5, $$6 }') \
		--ply-dir "$(word 2,$^)" \
		--model "$(word 3,$^)/model" \
		--seed 123 \
		--binder-name p2 \
		--positive-name p1 \
		--num_iterations 10_000 \
		--num_iter_eval 250 \
		--batch_size 32 \
		--validation_batch_size 1000 \
		--feat-mask 1.0 1.0 1.0 1.0 1.0 1.0 \
		--contact-distance 5.0 \
		--sc-max-cutoff 1.0 \
		--sc-min-cutoff 0.5 \
		--pos-surf-accept-probability 1.0 \
		$$(cat "$(word 1,$^)/stcrdab_split.csv" | grep "training" | awk -F, -v dir="$(word 1,$^)" '{ printf "%s/%s_%s%s%s%s%s ", dir, $$1, $$2, $$3, $$4, $$5, $$6 }')
	@rm "$@/checkpoint"
	@touch $@

lint:
	@FL_STATUS=0; PY_STATUS=0; \
	flake8 src || FL_STATUS=$$?; \
	pylint src || PY_STATUS=$$?; \
	exit $$(($$FL_STATUS | $$PY_STATUS))

test:
	@pytest tests/
	@prysk tests/apps

docs:
	@sphinx-apidoc -f -e -o docs/source src/tcr_antigen_prediction src/**/apps/*
	@python docs/document_clis.py docs/source
	@sphinx-build -b html ./docs ./docs/public