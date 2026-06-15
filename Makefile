.PHONY: checks clean

checks:
	python run_all_checks.py --clean

clean:
	rm -rf outputs/latest outputs/latest_test __pycache__ validators/__pycache__
