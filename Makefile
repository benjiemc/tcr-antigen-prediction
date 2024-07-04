.PHONY: all lint

all:
	@echo TODO

data/raw/stcrdab:
	@python -m tcr_antigen_prediction.apps.download_stcrdab $@

lint:
	@flake8 src/
	@pylint src