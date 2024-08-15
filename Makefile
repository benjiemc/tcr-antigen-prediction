.PHONY: all lint test docs


all:
	@echo TODO

data/raw/stcrdab:
	@python -m tcr_antigen_prediction.apps.download_stcrdab $@

data/interim/selected-stcrdab: data/raw/stcrdab
	@python -m tcr_antigen_prediction.apps.select_tcr_pmhc_structures --seed 123 -o $@ $^

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