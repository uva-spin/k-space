.PHONY: checks clean linebreak-check

checks:
	python run_all_checks.py --clean

linebreak-check:
	python tools/check_line_breaks.py .

clean:
	rm -rf outputs/latest outputs/latest_test __pycache__ validators/__pycache__
