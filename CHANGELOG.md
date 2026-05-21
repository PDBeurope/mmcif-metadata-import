# Changelog

All notable changes to this project will be documented in this file.

## [Released]

## [Unreleased]

## [0.5.5] - 2026-05-21

### Changed
- EM sample / entity assembly categories moved from `specs/EM.csv` to method-dependent `specs/MACROMOLECULES_EM.csv`, appended when `--macromolecules` is selected and the profile deposition is electron microscopy (same profile rules as conditional author specs).
- Added EM macromolecule categories: `em_entity_assembly_molwt`, `em_entity_assembly_recombinant`, `em_virus_entity`, `em_virus_natural_host`, `em_virus_shell` (with existing `em_entity_assembly` and `em_entity_assembly_naturalsource`).
- `MACROMOLECULE_CATEGORIES` in `polymer_safeguards.py` extended to include EM sample categories so selective macromolecule overwrite and merge remapping treat them with the rest of the macromolecule bundle.

### Added
- `profile_indicates_em_method()`, `resolve_macromolecules_em_supplement_spec_path()`, and `_resolve_macromolecules_specs_for_import()` in `import_metadata.py`.

## [0.5.4] - 2026-05-19

### Fixed
- **Frame-format** macromolecule categories (e.g. archive-style **`_entity_src_gen`**, **`_entity_poly`**) now participate in polymer **`entity_id`** remapping when content-aligned chain pairing maps a reference entity to one or more target entities (e.g. reference entity **`1`** → target **`A`** and **`B`**). Single-target remaps stay as frames; multi-target remaps emit a **`loop_`** with duplicated rows so **`entity_id`** matches target **`_entity.id`**.
- Loop construction for remapped macromolecule tables uses gemmi **`init_mmcif_loop`** so values with commas (e.g. **`pdbx_strand_id` `A,B`**) or semicolon multiline sequences round-trip correctly.

### Added
- Unit tests for frame **`_entity_src_gen`** and **`_entity_poly`** remapping in **`tests/test_polymer_safeguards.py`**.

## [0.5.3] - 2026-05-19

### Fixed
- Polymer **`loop_` `_entity`** rows are preserved and remapped to target **`entity.id`** values with reference **`pdbx_description`**, so content-aligned merges no longer drop polymer names when upload entity ids differ from the reference (e.g. letter ids vs numeric).
- **`needs_polymer_metadata_id_remap()`** triggers entity/asym remapping when chains align by content but entity ids still need alignment, even if **`label_asym_id`** sets already match.
- Macromolecule merge uses **two-phase** replacement (non-macromolecule category splice, then gemmi macromolecule overwrite) to avoid duplicating categories such as **`_software`**.
- Reference **`loop_`** blocks are skipped when **any** tag from the loop already exists in the target, avoiding partial-loop splices.

### Added
- Unit tests for polymer **`_entity`** / **`_entity_poly`** remapping in **`tests/test_polymer_safeguards.py`**.

## [0.5.2] - 2026-05-15

### Changed
- After successful content-aligned chain remapping, macromolecule categories in the merge target are **always replaced** (selective overwrite of `MACROMOLECULES.csv` categories only), so bare target `_entity` / `_entity_poly` rows do not block a clean merge when `--overwrite-existing` is not set.

## [0.5.1] - 2026-05-15

### Fixed
- Content-aligned macromolecule merges now **remap `entity_id` references** to the target model and **reconcile polymer `_struct_asym`** with `_atom_site`, so merged mmCIF files do not mix reference chain names (e.g. `A`/`B`) with target coordinates (`Axp`/`Bxp`).

## [0.5.0] - 2026-05-15

### Changed
- Macromolecule merge safeguards: when reference and target polymer **`label_asym_id`** sets differ, align chains by **coordinate-derived sequence** (unique 1:1 content match) instead of failing immediately. Failure rule **`ALIGN-1-CONTENT-MISMATCH`** replaces **`ALIGN-1-ASYMM-SET`** for unalignable cases. Successful content alignment is logged as **`content_aligned`** and **`chain_pairing`** in safeguard JSON.

### Added
- **`pair_polymer_chains_by_content()`** and **`polymer_profiles_by_asym()`** in `polymer_safeguards.py`; unit tests for renamed-chain and content-mismatch cases.

### Documentation
- [`docs/macromolecule-safeguards.md`](docs/macromolecule-safeguards.md) and README updated for content alignment and new rule code.

## [0.4.0] - 2026-05-10

### Added
- **Conditional author specifications**: `specs/AUTHORS_EM_MAP_ONLY.csv`, `AUTHORS_EM_WITH_ATOM_SITE.csv`, and `AUTHORS_DEFAULT.csv`. With **`--authors`**, the importer chooses one of these using a **profile** mmCIF (merge target when **`--merge_to_file`** is set, otherwise the input file): electron microscopy without an `_atom_site` loop uses the map-only set; EM with `_atom_site` uses the full EM set; other methods use audit + contact only. If `_exptl.method` (and EM `database_2` where required) cannot be determined, the tool falls back to **`specs/AUTHORS.csv`** (union of all three author categories).
- **`resolve_authors_spec_path()`**, **`block_has_atom_site()`**, and **`EM_METHOD_CODES`** in `import_metadata.py`.
- **Web app** and **`metadata_import.ipynb`**: authors checkbox uses the same resolution rules as the CLI.
- **Developer demos**: `dev/temp_test/author_demo/run_author_demos.py` plus minimal mmCIF fixtures under `dev/temp_test/author_demo/` (outputs go to `author_demo/output/`).

### Changed
- CLI help text for **`--authors`** now refers to the dynamic `AUTHORS*.csv` selection.
- **README**: documents author profiles, fallback, Python helpers, and the author demo script.

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
