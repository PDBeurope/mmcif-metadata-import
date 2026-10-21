# Macromolecule merge safeguards — user reference

This document describes the **automatic checks** that run when you merge metadata using **`--macromolecules`** together with **`--merge_to_file`**. Checks are **per polymer entity**, not all-or-nothing for the macromolecule spec. Matching entities are imported; entities that fail alignment are left as they are on the target. Other selected categories are still merged. The CLI exits with **code 2** only when **no** polymer entity is eligible (see [README](../README.md)). A partial import exits **0** and sets **`partial`: true** in the safeguard JSON.

The authoritative implementation is [`polymer_safeguards.py`](../polymer_safeguards.py) in the repository root.

---

## Terms

| Term | Meaning |
|------|---------|
| **Reference** | The mmCIF file given as the positional **`input_file`** (metadata is read from here). |
| **Target** | The file passed to **`--merge_to_file`** (metadata is merged into the first `data_` block of this file). |
| **Safeguard failure JSON** | When checks fail, details appear in the import **log** under **MACROMOLECULE SAFEGUARDS** and in `ImportMetadataOutcome.safeguard_result` (programmatic use). Each failure object includes a **`rule`** string (the codes below). |

---

## When do safeguards run?

- **`--macromolecules`** is selected **and**
- **`--merge_to_file`** is set **and**
- **`--no-macromolecule-safeguards`** is **not** used.

They apply only to the macromolecule category set from `MACROMOLECULES.csv`; they do not replace normal method validation (`--xray` / `--em` / …).

### EM map-only exceptions

| Reference | Merge target | Macromolecules behaviour |
|-----------|--------------|---------------------------|
| **EM map-only** | **EM map-only** | **Blind copy** — safeguards skipped; categories replaced on target (`mode`: `blind_copy_em_map_only`, exit **0**) |
| **EM map+model** / **model-only** (any ref with coordinates) | **EM map-only** | **Blind copy** on target (same as above) |
| **EM map-only** | **Any target with `_atom_site`** | **Blocked** — macromolecule categories not merged; exit **1** (`mode`: `blocked_map_only_reference_to_model_target`) |
| **EM map+model** / **model-only** | **EM map+model** / **model-only** | **Full safeguards** (unchanged) |

Map-only detection uses **`detect_method_from_input`** on the merge target / reference (`EM_MAP_ONLY` = WWPDB + EMDB, no PDB), not “missing `_atom_site`” alone.

---


## What is being compared?

1. **Which polymer chains exist**  
   The tool builds the set of **`_atom_site.label_asym_id`** values that count as **polymer chains** on reference and target, using the same **mode** on both sides (see next section).

2. **Chain alignment (names or content)**  
   - If the polymer **`label_asym_id`** sets are **identical** on both sides, each reference chain is compared to the target chain with the **same** ID.  
   - If the sets **differ**, each shared `label_asym_id` is still compared to the chain with the **same** ID. Only chains present on one side are paired by a **unique** residue count and coordinate-derived sequence. A same-ID sequence mismatch is not reassigned to a different chain. Unmatched extra chains are skipped. Duplicate profiles that are not a unique 1:1 match are **not** paired (`ALIGN-1-AMBIGUOUS-PROFILE`). When an eligible pair uses different chain names, the log includes **`content_aligned`: true** and a **`chain_pairing`** map for the **eligible** pairs only (reference ID → target ID).  
   - If **no** chain can be paired, the run fails with **`ALIGN-1-CONTENT-MISMATCH`** and macromolecule categories are not merged.

3. **ID remapping on merge**  
   When chain names differ, or entity ids differ for an eligible pair, the importer **rewrites** macromolecule metadata before merge so `entity_id` values match the target. A **full** content-aligned match still **replaces** macromolecule categories on the target (even without `--overwrite-existing`). A **partial** match replaces only rows for eligible entity ids and leaves every other row in those categories in place. Other requested categories still use the default skip-if-present merge behaviour.
   - **`entity_id`** fields point at the target’s polymer **`entity.id`** values (e.g. reference entity `1` shared by chains `A`/`B` may expand to target entities `A` and `B`). This applies to **`loop_`** and **frame** macromolecule categories (e.g. archive-style **`_entity_src_gen`**); when one reference row maps to several target polymers, the importer may emit a **`loop_`** with one row per target entity and duplicated field values.
   - Reference **`_entity`** polymer rows that would conflict with the target are omitted (the target’s polymer entity rows are kept).
   - After a full content-aligned category replacement, **`_struct_asym`** polymer rows are reconciled so **`_struct_asym.id`** matches **`_atom_site.label_asym_id`** and **`_struct_asym.entity_id`** matches **`_atom_site.label_entity_id`** for each polymer chain. Partial import does not rewrite `_struct_asym`, so rows for skipped chains stay as they were.

4. **Per aligned chain pair**  
   For each reference polymer chain and its aligned target chain:
   - **Residue count**: number of distinct sequence positions from **`_atom_site.label_seq_id`**, or **`auth_seq_id`** when `label_seq_id` is missing / `.` / `?`.
   - **Sequence from coordinates**: a string derived from **`_atom_site.label_comp_id`** per position, in wwPDB-style one-letter form (standard amino acids as one letter; non-standard residues as **`(COMP_ID)`** in parentheses).
   - **Entity mode only**: that atom-site sequence is also compared to **`_entity_poly.pdbx_seq_one_letter_code`** for the entity linked to the chain on that side, and the reference and target entity-poly strings are compared when both are present.

An entity is eligible only when **every** polymer chain of that entity, on both the reference and the target after pairing, passes ALIGN-2 and (in `entity` mode) ALIGN-3. Metadata for that entity is imported. A mismatched chain does not block entities that fully match.

Eligible rows replace the target’s rows for those entity ids. Ligand, water, and other non-polymer `_entity` rows are kept, as are polymer rows for entities that were not eligible. EM assembly and virus categories are imported only when every polymer entity they reference is eligible; otherwise that category is omitted (`omitted_categories` in the JSON) and the target’s copy of it is left in place.

If no polymer entity is eligible, none of the macromolecule categories are merged (exit **2**).

---

## Identification mode (`mode` in the JSON)

The log / JSON includes **`mode`**:

| Value | Meaning |
|-------|---------|
| **`entity`** | Both files have usable **`_entity`** and **`_struct_asym`** data. Polymer chains are those whose **`_struct_asym.entity_id`** points to an **`_entity.type`** in the polymer list (polypeptide(L/D), polydeoxyribonucleotide, polyribonucleotide, polysaccharide). **`branched`** entity types are **not** included. |
| **`forced_fallback`** | At least one side does not qualify for the entity path, so **both** sides use **`_atom_site`** heuristics only: chains with **≥ 2** distinct sequence positions and **at least one** `ATOM` row (so chains that are only `HETATM` are not inferred as polymer here). Sequence positions use **`label_seq_id`** when set; otherwise **`auth_seq_id`** (PDB-style mmCIF where `label_seq_id` is `.`). |

---

## Failure rule codes

These are the values of the **`rule`** field in each object inside **`failures`** in the JSON.

### `ALIGN-0-NO-POLYMER-CHAINS`

**Meaning:** No polymer chains could be identified on the reference (and typically neither side) for safeguard comparison, so macromolecule categories are **not** merged. This replaces the old vacuous pass when both polymer sets were empty.

**Typical cause:** Missing or unusable coordinates for polymer detection (e.g. session-generated reference without `_atom_site`), or `_atom_site` present but neither `label_seq_id` nor `auth_seq_id` usable. **Exception:** EM map-only merge targets use blind copy upstream and do not hit this rule.

**Fields (among others):** `mode`, `reference_asym_ids`, `target_asym_ids`, `note`.

---

### `ALIGN-1-CONTENT-MISMATCH`

**Meaning:** Polymer **`label_asym_id`** sets differ **and** no chain could be paired by a unique coordinate-derived sequence. Extra chains alone do not produce this rule when at least one unique pair exists.

**Typical cause:** No shared polymer profile, or every shared profile is ambiguous (see `ALIGN-1-AMBIGUOUS-PROFILE`).

**Fields (among others):** `reference_asym_ids`, `target_asym_ids`, `reference_chain_count`, `target_chain_count`, `mode`.

**Note:** Different chain **names** alone (e.g. `A` vs `Axp`) do **not** fail if a unique content pair exists. Unmatched extra chains are skipped rather than failing the whole macromolecule spec. See **Successful check** for `content_aligned` / `chain_pairing`, and **Partial import** for `partial`.

---

### `ALIGN-1-AMBIGUOUS-PROFILE`

**Meaning:** This chain’s coordinate-derived sequence is not unique on both sides, so it cannot be paired. The chain (and its entity) is skipped. Other uniquely paired entities may still be imported.

**Fields (among others):** `label_asym_id`, `note`.

---

### `ALIGN-ENTITY-INCOMPLETE`

**Meaning:** This chain’s coordinates matched, but another polymer chain of the same entity did not. Shared entity metadata is not imported for that entity.

**Fields (among others):** `label_asym_id`, `entity_id`, `note`.

---

### `ALIGN-2-LENGTH-OR-ATOM-SEQ`

**Meaning:** For a given **`label_asym_id`** that exists on both sides, either the **number of distinct sequence positions** (`label_seq_id`, else `auth_seq_id`) differs, or the **coordinate-derived sequence string** differs.

**Typical cause:** Different construct, missing residues, different numbering coverage, or real sequence mismatch.

**Fields (among others):** `label_asym_id`, `reference_count`, `target_count`, `reference_atom_site_seq`, `target_atom_site_seq`.

---

### `ALIGN-3-REF-ATOM-VS-ENTITY-POLY`

**Meaning:** In **`entity`** mode, on the **reference** file, the sequence built from **`_atom_site`** for this chain does not match **`_entity_poly.pdbx_seq_one_letter_code`** for the linked entity.

**Typical cause:** Internal inconsistency in the reference mmCIF (coordinates vs deposited sequence).

**Fields (among others):** `label_asym_id`, `entity_id`, `atom_site_seq`, `entity_poly_seq`.

---

### `ALIGN-3-TGT-ATOM-VS-ENTITY-POLY`

**Meaning:** Same as above, but for the **target** file.

**Typical cause:** Internal inconsistency in the target mmCIF.

**Fields (among others):** `label_asym_id`, `entity_id`, `atom_site_seq`, `entity_poly_seq`.

---

### `ALIGN-3-ENTITY-POLY-MISMATCH`

**Meaning:** In **`entity`** mode, for a shared **`label_asym_id`**, the **`_entity_poly.pdbx_seq_one_letter_code`** strings for the reference entity and the target entity differ.

**Typical cause:** Reference and target describe genuinely different polymers while reusing the same chain ID.

**Fields (among others):** `label_asym_id`, `reference_entity_poly`, `target_entity_poly`.

---

## Successful check (`ok: true`)

If **`ok`** is **true**, at least one polymer entity was eligible and its macromolecule metadata was merged. The log may still list **`checked_asym_ids`** and **`mode`** for traceability.

When alignment used content matching, the JSON also includes **`content_aligned`: true** and **`chain_pairing`** (reference `label_asym_id` → target `label_asym_id`) for the eligible pairs.

### Partial import (`partial: true`)

`partial` is true when some polymer entities were imported and some were not. Exit status is still **0**. The JSON then also includes:

| Field | Meaning |
|-------|---------|
| **`imported_asym_ids`** | Reference polymer chains whose entity was imported |
| **`skipped_asym_ids`** | Reference chains that were not imported, plus target chains that were not the image of an imported pair |
| **`imported_entity_ids`** | Merge-target `_entity.id` values that received reference metadata |
| **`skipped_entity_ids`** | Merge-target polymer `_entity.id` values left unchanged |
| **`omitted_categories`** | Requested macromolecule categories dropped because they referenced a skipped polymer (for example an EM assembly `entity_id_list`) |

---

## Disabling safeguards

**`--no-macromolecule-safeguards`** skips all of the above. Use only when you accept the risk of inconsistent or dictionary-invalid macromolecule metadata after merge.

---

## Where to look in output

- **CLI / `--log`:** section **MACROMOLECULE SAFEGUARDS** in the `.log` file (includes JSON-style details).
- **Python:** `import_metadata(...)` returns **`ImportMetadataOutcome`**; inspect **`safeguard_result`** (dataclass-like object with **`ok`**, **`mode`**, **`failures`**, **`checked_asym_ids`**, and when applicable **`partial`**, **`content_aligned`**, **`chain_pairing`**, **`imported_asym_ids`**, **`skipped_asym_ids`**, **`imported_entity_ids`**, **`skipped_entity_ids`**) and **`to_json_fragment()`** if needed.

---

*This reference is meant to stay aligned with `polymer_safeguards.py`. If behaviour changes, update this file in the same change.*
