.PHONY: all lint test docs


all:
	@echo TODO

data/raw/stcrdab:
	@python -m tcr_antigen_prediction.apps.download_stcrdab $@

lint:
	@flake8 src
	@pylint src

test:
	@prysk tests/apps

docs:
	@sphinx-apidoc -f -e -o docs/source src/tcr_antigen_prediction
	@sphinx-build -b html ./docs public