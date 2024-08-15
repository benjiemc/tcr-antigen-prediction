.PHONY: all lint test docs


all:
	@echo TODO

data/raw/stcrdab:
	@python -m tcr_antigen_prediction.apps.download_stcrdab $@

data/interim/selected-stcrdab: data/raw/stcrdab
	@python -m tcr_antigen_prediction.apps.select_tcr_pmhc_structures --seed 123 -o $@ $^

lint:
	@flake8 src
	@pylint src

test:
	@pytest tests/
	@prysk tests/apps

docs:
	@sphinx-apidoc -f -e -o docs/source src/tcr_antigen_prediction
	@sphinx-build -b html ./docs public