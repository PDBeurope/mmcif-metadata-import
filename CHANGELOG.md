# Changelog

All notable changes to this project will be documented in this file.

## [Released]

## [0.3.0] - 2026-05-10

### Added
- **Merge overwrite mode**: `--overwrite-existing` (CLI) and `import_metadata(..., overwrite_existing=True)` (Python). When merging, conflicting pairs and loops in the merge target’s first data block are removed and replaced with imported metadata; additional `data_` blocks are preserved. Requires `--merge_to_file` on the CLI (error if used alone).
- **`MergeMetadataResult`** returned by **`merge_metadata_to_file()`**: `success`, `skipped_categories`, `skipped_items`, `overwritten_categories`, `overwritten_items`.
- Import log sections **CATEGORIES OVERWRITTEN IN MERGE TARGET** and **ITEMS OVERWRITTEN IN MERGE TARGET** when overwrite mode is used; summary counts for overwritten categories/items.
- Web app: checkbox **Overwrite existing metadata in merge target** (used when a merge file is uploaded).

### Changed
- Default merge behavior is unchanged: without `--overwrite-existing`, categories/items already present in the target are still skipped and listed under “not imported” in the log.

## [0.2.0] - 2026-05-09

### Added
- [`docs/macromolecule-safeguards.md`](docs/macromolecule-safeguards.md): user reference for macromolecule merge checks and log rule codes (`ALIGN-*`).
- Macromolecule merge safeguards (`polymer_safeguards.py`): when merging with `MACROMOLECULES.csv`, validate reference vs target polymer chains (`_entity` / `_struct_asym` when both qualify, otherwise forced `_atom_site` fallback on both sides) before copying macromolecule categories.
- CLI flag `--no-macromolecule-safeguards` to skip validation.
- Exit code **2** when the merge completes but macromolecule categories were omitted due to failed safeguards (exit **1** remains for hard failures).
- Import log section **MACROMOLECULE SAFEGUARDS** with structured failure details when checks run.
- Web app warning when macromolecules are skipped; `alert-warning` styling in the template.
- Unit tests in `tests/test_polymer_safeguards.py`.

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
