.PHONY: install install-dev lint format test test-unit test-fast test-integration test-coverage clean help

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-18s\033[0m %s\n", $$1, $$2}'

install:  ## Install package with core dependencies
	pip install -e .

install-dev:  ## Install with dev, ml, and openmm extras
	pip install -e ".[dev,ml,openmm]"

install-minimal:  ## Install with core dependencies only (no ML/OpenMM)
	pip install -e ".[minimal]"

lint:  ## Run ruff linter and format check
	python -m ruff check CondenSimAdapter/ tests/
	python -m ruff format --check CondenSimAdapter/ tests/

format:  ## Auto-format code with ruff
	python -m ruff format CondenSimAdapter/ tests/
	python -m ruff check --fix CondenSimAdapter/ tests/

test:  ## Run all tests
	pytest -v

test-unit:  ## Run unit tests only
	pytest tests/unit/ -v

test-fast:  ## Run fast unit tests (no GPU/GROMACS/slow)
	pytest tests/unit/ -v -m "not slow and not gpu and not gmx" --no-header -q

test-integration:  ## Run integration tests only
	pytest tests/integration/ -v

test-coverage:  ## Run tests with coverage report
	pytest --cov=CondenSimAdapter --cov-report=term-missing --cov-report=html -v

clean:  ## Remove build artifacts and caches
	rm -rf build/ dist/ *.egg-info .pytest_cache .ruff_cache htmlcov .coverage
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
