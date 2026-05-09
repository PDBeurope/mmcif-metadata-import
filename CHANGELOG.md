# Changelog

All notable changes to this project will be documented in this file.

## [Released]

## [0.2.0] - 2026-05-09

### Added
- Macromolecule merge safeguards (`polymer_safeguards.py`): when merging with `MACROMOLECULES.csv`, validate reference vs target polymer chains (`_entity` / `_struct_asym` when both qualify, otherwise forced `_atom_site` fallback on both sides) before copying macromolecule categories.
- CLI flag `--no-macromolecule-safeguards` to skip validation.
- Exit code **2** when the merge completes but macromolecule categories were omitted due to failed safeguards (exit **1** remains for hard failures).
- Import log section **MACROMOLECULE SAFEGUARDS** with structured failure details when checks run.
- Web app warning when macromolecules are skipped; `alert-warning` styling in the template.
- Unit tests in `tests/test_polymer_safeguards.py`.
- Design notes in `dev/macromolecules-import-safeguards-plan.md`.

### Changed
- **`import_metadata()`** now returns **`ImportMetadataOutcome`** (`ok`, `exit_status`, optional `safeguard_result`) instead of a bare boolean. Callers should use `outcome.ok` (and optionally `outcome.exit_status`).
- `metadata_import.ipynb` updated for the new return type.

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
