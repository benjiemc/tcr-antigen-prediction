.PHONY: all data lint test docs


all: data

data: \
	data/processed/selected-stcrdab_feats \
	data/processed/selected-stcrdab_ply

data/raw/stcrdab:
	@python -m tcr_antigen_prediction.apps.download_stcrdab $@

data/interim/selected-stcrdab: data/raw/stcrdab data/external/masif_ppi_search_training_set.txt
	@python -m tcr_antigen_prediction.apps.select_tcr_pmhc_structures \
		--seed 123 \
		--remove-structures-missing-residues \
		--structural-similarity-cutoff 2.0 \
		--pdb-ids-to-exclude $$(cut -d _ -f 1 $(word 2,$^) | tr '[:upper:]' '[:lower:]' | sort | uniq | tr '\n' ' ') \
		-o $@ \
		$(word 1,$^)

data/interim/selected-stcrdab_crop: data/interim/selected-stcrdab
	@mkdir -p $@
	@echo "Processing structures..."
	@head -n1 "$^/stcrdab_split.csv" > "$@/stcrdab_split.csv"

	num_lines=$$(cat "$^/stcrdab_split.csv" | wc -l); \
	line=1; \
	while [ $$line -lt $$num_lines ]; do \
		let line++; \

		pdb_id=$$(sed -n "$${line}p" $^/stcrdab_split.csv | cut -d ',' -f1); \
		alpha_chain=$$(sed -n "$${line}p" $^/stcrdab_split.csv | cut -d ',' -f2); \
		beta_chain=$$(sed -n "$${line}p" $^/stcrdab_split.csv | cut -d ',' -f3); \
		antigen_chain=$$(sed -n "$${line}p" $^/stcrdab_split.csv | cut -d ',' -f4); \
		mhc_chain1=$$(sed -n "$${line}p" $^/stcrdab_split.csv | cut -d ',' -f5); \
		mhc_chain2=$$(sed -n "$${line}p" $^/stcrdab_split.csv | cut -d ',' -f6); \
		mhc_type=$$(sed -n "$${line}p" $^/stcrdab_split.csv | cut -d ',' -f7); \

		chains="$${alpha_chain}$${beta_chain}$${antigen_chain}$${mhc_chain1}$${mhc_chain2}"; \
		file_name="$${pdb_id}_$${chains}.pdb"; \

		@echo "Working on $${file_name}..."; \

		if [ "$${mhc_type}" = "MH1" ]; then \
			@python -m tcr_antigen_prediction.apps.crop_tcr_pmhc \
				"$^/$${file_name}" \
				-o "$@/$${file_name}" \
				--tcr-chains $$alpha_chain $$beta_chain \
				--mhc-chains $$mhc_chain1 \
				--antigen-chain $$antigen_chain; \

		elif [ "$${mhc_type}" = "MH2" ]; then \
			@python -m tcr_antigen_prediction.apps.crop_tcr_pmhc \
				"$^/$${file_name}" \
				-o "$@/$${file_name}" \
				--tcr-chains $$alpha_chain $$beta_chain \
				--mhc-chains $$mhc_chain1 $$mhc_chain2 \
				--antigen-chain $$antigen_chain; \

		else \
			@echo "INVALID MHC TYPE"; \
			exit 1; \

		fi; \
		@sed -n "$${line}p" "$^/stcrdab_split.csv" >> "$@/stcrdab_split.csv"; \
		@echo "Finished Structure"; \

	done

data/processed/selected-stcrdab_ply: data/interim/selected-stcrdab_crop
	@mkdir -p $@
	head -n1 "$^/stcrdab_split.csv" > "$@/stcrdab_split.csv"

	num_lines=$$(cat "$^/stcrdab_split.csv" | wc -l); \
	line=1; \
	while [ $$line -lt $$num_lines ]; do \
	    let line++; \
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

data/processed/selected-stcrdab_feats: data/processed/selected-stcrdab_ply
	mkdir -p $@
	head -n1 "$^/stcrdab_split.csv" > "$@/stcrdab_split.csv"

	num_lines=$$(sed -n '$$=' "$^/stcrdab_split.csv"); \
	line=1; \
	while [ $$line -lt $$num_lines ]; do \
	    let line++; \
	    pdb_id=$$(sed -n "$${line}p" "$^/stcrdab_split.csv" | cut -d ',' -f1); \
	    alpha_chain=$$(sed -n "$${line}p" "$^/stcrdab_split.csv" | cut -d ',' -f2); \
	    beta_chain=$$(sed -n "$${line}p" "$^/stcrdab_split.csv" | cut -d ',' -f3); \
	    antigen_chain=$$(sed -n "$${line}p" "$^/stcrdab_split.csv" | cut -d ',' -f4); \
	    mhc_chain1=$$(sed -n "$${line}p" "$^/stcrdab_split.csv" | cut -d ',' -f5); \
	    mhc_chain2=$$(sed -n "$${line}p" "$^/stcrdab_split.csv" | cut -d ',' -f6); \
	    chains="$${alpha_chain}$${beta_chain}$${antigen_chain}$${mhc_chain1}$${mhc_chain2}"; \

	    echo "Working on $$pdb_id chains $$chains..."; \
	    output_name="$@/$${pdb_id}_$${chains}"; \

	    if [ ! -d $$output_name ]; then \
	        mkdir $$output_name; \
	    fi; \

	    python -m tcr_antigen_prediction.apps.compute_features --mode ppi_search -o "$$output_name" "$^/$${pdb_id}_$${chains}_$${alpha_chain}$${beta_chain}.ply" "$^/$${pdb_id}_$${chains}_$${antigen_chain}$${mhc_chain1}$${mhc_chain2}.ply" || continue; \
	    sed -n "$${line}p" "$^/stcrdab_split.csv" >> "$@/stcrdab_split.csv"; \
	    echo "Finished pair"; \
	done

	@echo "All done."

data/external/masif_ppi_search_training_set.txt:
	@wget -O $@ https://raw.githubusercontent.com/LPDI-EPFL/masif/master/data/masif_ppi_search/lists/training.txt

lint:
	@FL_STATUS=0; PY_STATUS=0; \
	flake8 src || FL_STATUS=$$?; \
	pylint src || PY_STATUS=$$?; \
	exit $$(($$FL_STATUS | $$PY_STATUS))

test:
	@pytest tests/
	@prysk tests/apps

docs:
	@sphinx-apidoc -f -e -o docs/source src/tcr_antigen_prediction
	@sphinx-build -b html ./docs public