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

---

## What is being compared?

1. **Which polymer chains exist**  
   The tool builds the set of **`_atom_site.label_asym_id`** values that count as **polymer chains** on reference and target, using the same **mode** on both sides (see next section).

2. **Per shared chain**  
   For each `label_asym_id` that appears on **both** sides:
   - **Residue count**: number of distinct **`_atom_site.label_seq_id`** values (polymer positions).
   - **Sequence from coordinates**: a string derived from **`_atom_site.label_comp_id`** per position, in wwPDB-style one-letter form (standard amino acids as one letter; non-standard residues as **`(COMP_ID)`** in parentheses).
   - **Entity mode only**: that atom-site sequence is also compared to **`_entity_poly.pdbx_seq_one_letter_code`** for the entity linked to the chain, and the reference and target entity-poly strings are compared when both are present.

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

### `ALIGN-1-ASYMM-SET`

**Meaning:** The set of polymer **`label_asym_id`** values on the reference does not match the set on the target (same IDs required; order does not matter).

**Typical cause:** Different chain naming between files (e.g. `A` vs `Axp`), or one file has extra/missing polymer chains.

**Fields (among others):** `reference_asym_ids`, `target_asym_ids`, `mode`.

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

---

## Disabling safeguards

**`--no-macromolecule-safeguards`** skips all of the above. Use only when you accept the risk of inconsistent or dictionary-invalid macromolecule metadata after merge.

---

## Where to look in output

- **CLI / `--log`:** section **MACROMOLECULE SAFEGUARDS** in the `.log` file (includes JSON-style details).
- **Python:** `import_metadata(...)` returns **`ImportMetadataOutcome`**; inspect **`safeguard_result`** (dataclass-like object with **`ok`**, **`mode`**, **`failures`**, **`checked_asym_ids`**) and **`to_json_fragment()`** if needed.

---

*This reference is meant to stay aligned with `polymer_safeguards.py`. If behaviour changes, update this file in the same change.*
