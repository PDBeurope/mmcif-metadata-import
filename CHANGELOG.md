# Changelog

All notable changes to this project will be documented in this file.

## [Released]

## [0.1.0] - 2026-04-08

### Changed
- README and `docs/user-tutorial.html`: document PyPI install and the `mmcif-metadata-import` CLI; note `python import_metadata.py` for runs from a cloned repo.

### Added
- `docs/index.html`: redirect from the GitHub Pages site root to `user-tutorial.html`.
- PyPI publishing workflow at `.github/workflows/pypi-publish.yml`, triggered by version tags (`v*`).
- Python packaging metadata in `pyproject.toml` so the project can be built and published as a package.
- Console entry point `mmcif-metadata-import` mapped to `import_metadata:main`.
- Distribution of specification CSV files from `specs/*.csv` in built artifacts.
- Spec file path resolution fallback in `import_metadata.py` to support both repository and installed-package usage.
