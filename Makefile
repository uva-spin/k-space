PYTHON ?= python
.PHONY: check test policy list
check:
	$(PYTHON) tools/check_repository_policy.py
	$(PYTHON) run_all_checks.py --clean
test:
	$(PYTHON) -m pytest -q
policy:
	$(PYTHON) tools/check_repository_policy.py
list:
	$(PYTHON) run_all_checks.py --list

# Preserve the earlier repository command aliases.
.PHONY: checks linebreak-check clean
checks: check
linebreak-check:
	$(PYTHON) tools/check_line_breaks.py .
clean:
	$(PYTHON) -c "from pathlib import Path; import shutil; [shutil.rmtree(p) for p in [Path('outputs/latest'), Path('outputs/latest_test'), Path('__pycache__'), Path('validators/__pycache__')] if p.exists() and not p.is_symlink()]"

.PHONY: docs examples
docs:
	$(PYTHON) tools/check_documentation.py
examples:
	$(PYTHON) examples/canonical_evolution.py
	$(PYTHON) examples/momentum_grid.py --reference
