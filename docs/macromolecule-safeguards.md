# Macromolecule merge safeguards — user reference

This document describes the **automatic checks** that run when you merge metadata using **`--macromolecules`** together with **`--merge_to_file`**. If checks fail, categories from `specs/MACROMOLECULES.csv` are **not** copied; other selected categories are still merged. The CLI exits with **code 2** in that case (see [README](../README.md)).

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
   - If the sets **differ**, the tool tries **content alignment**: same number of chains, and a **unique 1:1 match** where each pair has the same residue count and coordinate-derived sequence (order of chain names does not matter). When this succeeds, the log includes **`content_aligned`: true** and a **`chain_pairing`** map (reference ID → target ID).  
   - If content alignment cannot be established, the run fails with **`ALIGN-1-CONTENT-MISMATCH`** (macromolecule categories are not merged).

3. **ID remapping on merge (content alignment only)**  
   When chain names differ but content alignment succeeds, the importer **rewrites** macromolecule metadata before merge and **replaces** any existing macromolecule categories in the merge target (even without `--overwrite-existing`). Other requested categories still use the default skip-if-present merge behaviour.
   - **`entity_id`** fields point at the target’s polymer **`entity.id`** values (e.g. reference entity `1` shared by chains `A`/`B` may expand to target entities `A` and `B`). This applies to **`loop_`** and **frame** macromolecule categories (e.g. archive-style **`_entity_src_gen`**); when one reference row maps to several target polymers, the importer may emit a **`loop_`** with one row per target entity and duplicated field values.
   - Reference **`_entity`** polymer rows that would conflict with the target are omitted (the target’s polymer entity rows are kept).
   - After merge, **`_struct_asym`** polymer rows are reconciled so **`_struct_asym.id`** matches **`_atom_site.label_asym_id`** and **`_struct_asym.entity_id`** matches **`_atom_site.label_entity_id`** for each polymer chain.

4. **Per aligned chain pair**  
   For each reference polymer chain and its aligned target chain:
   - **Residue count**: number of distinct **`_atom_site.label_seq_id`** values (polymer positions).
   - **Sequence from coordinates**: a string derived from **`_atom_site.label_comp_id`** per position, in wwPDB-style one-letter form (standard amino acids as one letter; non-standard residues as **`(COMP_ID)`** in parentheses).
   - **Entity mode only**: that atom-site sequence is also compared to **`_entity_poly.pdbx_seq_one_letter_code`** for the entity linked to the chain on that side, and the reference and target entity-poly strings are compared when both are present.

If every required check passes, macromolecule metadata is merged. If any check fails, **none** of the macromolecule categories from the safeguard set are merged for that run (all-or-nothing for those categories).

---

## Identification mode (`mode` in the JSON)

The log / JSON includes **`mode`**:

| Value | Meaning |
|-------|---------|
| **`entity`** | Both files have usable **`_entity`** and **`_struct_asym`** data. Polymer chains are those whose **`_struct_asym.entity_id`** points to an **`_entity.type`** in the polymer list (polypeptide(L/D), polydeoxyribonucleotide, polyribonucleotide, polysaccharide). **`branched`** entity types are **not** included. |
| **`forced_fallback`** | At least one side does not qualify for the entity path, so **both** sides use **`_atom_site`** heuristics only: chains with **≥ 2** distinct `label_seq_id` values and **at least one** `ATOM` row (so chains that are only `HETATM` are not inferred as polymer here). |

---

## Failure rule codes

These are the values of the **`rule`** field in each object inside **`failures`** in the JSON.

### `ALIGN-1-CONTENT-MISMATCH`

**Meaning:** Polymer **`label_asym_id`** sets differ **and** chains cannot be aligned by content: different chain counts, no matching multiset of coordinate-derived sequences, or ambiguous pairing (e.g. two chains with identical sequence on one side).

**Typical cause:** Genuinely different structures (extra/missing chain, different sequence), or duplicate polymer profiles that prevent a unique match.

**Fields (among others):** `reference_asym_ids`, `target_asym_ids`, `reference_chain_count`, `target_chain_count`, `mode`.

**Note:** Different chain **names** alone (e.g. `A` vs `Axp`) do **not** fail if content alignment succeeds; see **Successful check** for `content_aligned` / `chain_pairing` in the log.

---

### `ALIGN-2-LENGTH-OR-ATOM-SEQ`

**Meaning:** For a given **`label_asym_id`** that exists on both sides, either the **number of distinct `label_seq_id`** values differs, or the **coordinate-derived sequence string** differs.

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

If **`failures`** is empty and **`ok`** is **true**, macromolecule categories were allowed by safeguards. The log may still list **`checked_asym_ids`** and **`mode`** for traceability.

When alignment used content matching, the JSON also includes **`content_aligned`: true** and **`chain_pairing`** (reference `label_asym_id` → target `label_asym_id`).

---

## Disabling safeguards

**`--no-macromolecule-safeguards`** skips all of the above. Use only when you accept the risk of inconsistent or dictionary-invalid macromolecule metadata after merge.

---

## Where to look in output

- **CLI / `--log`:** section **MACROMOLECULE SAFEGUARDS** in the `.log` file (includes JSON-style details).
- **Python:** `import_metadata(...)` returns **`ImportMetadataOutcome`**; inspect **`safeguard_result`** (dataclass-like object with **`ok`**, **`mode`**, **`failures`**, **`checked_asym_ids`**, and when applicable **`content_aligned`**, **`chain_pairing`**) and **`to_json_fragment()`** if needed.

---

*This reference is meant to stay aligned with `polymer_safeguards.py`. If behaviour changes, update this file in the same change.*
