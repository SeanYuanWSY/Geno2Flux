.PHONY: setup test

setup:
	bash run.sh setup

test:
	uv run --frozen --extra test python -m pytest -q
