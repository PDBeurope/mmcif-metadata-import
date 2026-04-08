# Changelog

All notable changes to this project will be documented in this file.

## [Released]

## [0.1.0] - 2026-04-08

### Added
- PyPI publishing workflow at `.github/workflows/pypi-publish.yml`, triggered by version tags (`v*`).
- Python packaging metadata in `pyproject.toml` so the project can be built and published as a package.
- Console entry point `mmcif-metadata-import` mapped to `import_metadata:main`.
- Distribution of specification CSV files from `specs/*.csv` in built artifacts.
- Spec file path resolution fallback in `import_metadata.py` to support both repository and installed-package usage.
