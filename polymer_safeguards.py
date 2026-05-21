"""
Macromolecule merge safeguards (reference vs target mmCIF).

See docs/macromolecule-safeguards.md for user-facing documentation of checks and rule codes.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

import gemmi

# Must match specs/MACROMOLECULES.csv and specs/MACROMOLECULES_EM.csv (category column, should_import Y).
MACROMOLECULE_CATEGORIES: Set[str] = frozenset(
    {
        "_entity",
        "_entity_name_com",
        "_entity_poly",
        "_entity_poly_seq",
        "_entity_src_nat",
        "_entity_src_gen",
        "_pdbx_entity_src_syn",
        "_struct_ref",
        "_struct_ref_seq",
        "_struct_ref_seq_dif",
        "_em_entity_assembly",
        "_em_entity_assembly_naturalsource",
        "_em_entity_assembly_molwt",
        "_em_entity_assembly_recombinant",
        "_em_virus_entity",
        "_em_virus_natural_host",
        "_em_virus_shell",
    }
)

POLYMER_ENTITY_TYPES: Set[str] = frozenset(
    {
        "polypeptide(L)",
        "polypeptide(D)",
        "polydeoxyribonucleotide",
        "polyribonucleotide",
        "polysaccharide",
    }
)

# Standard protein residues → one-letter (same convention as typical pdbx one-letter strings).
AA_3TO1: Dict[str, str] = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "ASX": "B",
    "GLX": "Z",
    "UNK": "X",
}

# Standard DNA/RNA one-letter (entity_poly style for simple polymers).
NA_3TO1: Dict[str, str] = {
    "DA": "A",
    "DC": "C",
    "DG": "G",
    "DT": "T",
    "DU": "U",
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "U",
}


def _cif_string_raw(s: str) -> str:
    if not s:
        return ""
    t = s.strip()
    if (t.startswith("'") and t.endswith("'")) or (t.startswith('"') and t.endswith('"')):
        return t[1:-1]
    return t


def _is_missing_cif_value(s: str) -> bool:
    t = (s or "").strip()
    return t in ("?", ".", "")


def _get_loop(block: gemmi.cif.Block, tag: str) -> Optional[gemmi.cif.Loop]:
    """Return full Loop for a column tag (gemmi find_loop may return Column or Loop)."""
    col = block.find_loop(tag)
    if col is None:
        return None
    if isinstance(col, gemmi.cif.Loop):
        return col
    return col.get_loop()


def _loop_as_table(loop: gemmi.cif.Loop) -> Tuple[List[str], List[List[str]]]:
    """Return (tags, rows) with raw string cell values."""
    tags = list(loop.tags)
    w = loop.width()
    vals = list(loop.values)
    rows: List[List[str]] = []
    for i in range(0, len(vals), w):
        rows.append([str(vals[i + j]) for j in range(w)])
    return tags, rows


def _find_column(tags: List[str], *candidates: str) -> Optional[int]:
    for c in candidates:
        if c in tags:
            return tags.index(c)
    return None


def _entity_type_by_id(block: gemmi.cif.Block) -> Dict[str, str]:
    loop = _get_loop(block, "_entity.id")
    if not loop:
        return {}
    tags, rows = _loop_as_table(loop)
    ic = _find_column(tags, "_entity.id")
    it = _find_column(tags, "_entity.type")
    if ic is None or it is None:
        return {}
    out: Dict[str, str] = {}
    for row in rows:
        eid = _cif_string_raw(row[ic])
        et = _cif_string_raw(row[it])
        if eid and not _is_missing_cif_value(et):
            out[eid] = et
    return out


def _struct_asym_map(block: gemmi.cif.Block) -> Dict[str, str]:
    """label_asym_id -> entity_id (empty if loop missing)."""
    loop = _get_loop(block, "_struct_asym.id")
    if not loop:
        return {}
    tags, rows = _loop_as_table(loop)
    ia = _find_column(tags, "_struct_asym.id")
    ie = _find_column(tags, "_struct_asym.entity_id")
    if ia is None or ie is None:
        return {}
    out: Dict[str, str] = {}
    for row in rows:
        asym = _cif_string_raw(row[ia])
        eid = _cif_string_raw(row[ie])
        if asym and eid:
            out[asym] = eid
    return out


def primary_polymer_path_usable(block: gemmi.cif.Block) -> bool:
    """R-POLY-1 availability: _entity + _struct_asym with resolvable polymer asym ids."""
    et = _entity_type_by_id(block)
    sm = _struct_asym_map(block)
    if not et or not sm:
        return False
    for asym, eid in sm.items():
        if eid not in et:
            return False
        if et[eid] in POLYMER_ENTITY_TYPES:
            return True
    return False


def polymer_asym_ids_entity_path(block: gemmi.cif.Block) -> Set[str]:
    """Non-branched polymer asym ids from _struct_asym + _entity (R-POLY-1)."""
    et = _entity_type_by_id(block)
    sm = _struct_asym_map(block)
    out: Set[str] = set()
    for asym, eid in sm.items():
        t = et.get(eid)
        if t in POLYMER_ENTITY_TYPES:
            out.add(asym)
    return out


def _atom_site_columns(block: gemmi.cif.Block) -> Optional[Tuple[List[str], List[List[str]]]]:
    loop = _get_loop(block, "_atom_site.label_asym_id")
    if not loop:
        return None
    return _loop_as_table(loop)


def polymer_asym_ids_fallback(block: gemmi.cif.Block) -> Set[str]:
    """R-POLY-2: _atom_site heuristic (length >= 2 distinct label_seq_id, >= one ATOM)."""
    tbl = _atom_site_columns(block)
    if not tbl:
        return set()
    tags, rows = tbl
    ia = _find_column(tags, "_atom_site.label_asym_id")
    iseq = _find_column(tags, "_atom_site.label_seq_id")
    ig = _find_column(tags, "_atom_site.group_PDB")
    if ia is None or iseq is None or ig is None:
        return set()

    by_asym: Dict[str, Set[str]] = {}
    all_het: Dict[str, bool] = {}

    for row in rows:
        asym = _cif_string_raw(row[ia])
        if not asym:
            continue
        seq = _cif_string_raw(row[iseq])
        g = _cif_string_raw(row[ig]).upper()
        by_asym.setdefault(asym, set())
        all_het.setdefault(asym, True)
        if g == "ATOM":
            all_het[asym] = False
        if not _is_missing_cif_value(seq) and seq not in (".",):
            by_asym[asym].add(seq)

    out: Set[str] = set()
    for asym, seqs in by_asym.items():
        if len(seqs) >= 2 and not all_het.get(asym, True):
            out.add(asym)
    return out


def identification_mode(ref: gemmi.cif.Block, tgt: gemmi.cif.Block) -> str:
    """Return 'entity' or 'forced_fallback' (R-POLY-3)."""
    pr = primary_polymer_path_usable(ref)
    pt = primary_polymer_path_usable(tgt)
    if pr and pt:
        return "entity"
    return "forced_fallback"


def polymer_asym_ids_for_mode(block: gemmi.cif.Block, mode: str) -> Set[str]:
    if mode == "entity":
        return polymer_asym_ids_entity_path(block)
    return polymer_asym_ids_fallback(block)


def _entity_poly_one_letter(block: gemmi.cif.Block, entity_id: str) -> Optional[str]:
    loop = _get_loop(block, "_entity_poly.entity_id")
    if not loop:
        return None
    tags, rows = _loop_as_table(loop)
    ie = _find_column(tags, "_entity_poly.entity_id")
    isq = _find_column(tags, "_entity_poly.pdbx_seq_one_letter_code")
    if ie is None or isq is None:
        return None
    for row in rows:
        if _cif_string_raw(row[ie]) == entity_id:
            raw = _cif_string_raw(row[isq])
            return raw if raw and not _is_missing_cif_value(raw) else None
    return None


def _comp_id_to_seq_segment(comp_id: str, entity_type: str) -> str:
    c = comp_id.strip().upper()
    if not c:
        return ""
    if entity_type in ("polydeoxyribonucleotide", "polyribonucleotide"):
        if c in NA_3TO1:
            return NA_3TO1[c]
        return f"({c})"
    if c in AA_3TO1:
        return AA_3TO1[c]
    return f"({c})"


def _entity_type_for_asym(block: gemmi.cif.Block, asym: str) -> str:
    sm = _struct_asym_map(block)
    et = _entity_type_by_id(block)
    eid = sm.get(asym, "")
    return et.get(eid, "polypeptide(L)")


def atom_site_sequence_and_count(
    block: gemmi.cif.Block, asym_id: str
) -> Tuple[int, str]:
    """
    Distinct label_seq_id count and canonical pdbx-style one-letter string from _atom_site.
    One residue per label_seq_id; first non-missing label_comp_id wins (microheterogeneity).
    """
    tbl = _atom_site_columns(block)
    if not tbl:
        return 0, ""
    tags, rows = tbl
    ia = _find_column(tags, "_atom_site.label_asym_id")
    iseq = _find_column(tags, "_atom_site.label_seq_id")
    ic = _find_column(tags, "_atom_site.label_comp_id")
    if ia is None or iseq is None or ic is None:
        return 0, ""

    seq_to_comp: Dict[str, str] = {}
    for row in rows:
        if _cif_string_raw(row[ia]) != asym_id:
            continue
        sid = _cif_string_raw(row[iseq])
        if _is_missing_cif_value(sid) or sid in (".",):
            continue
        comp = _cif_string_raw(row[ic]).upper()
        if not comp or _is_missing_cif_value(comp):
            continue
        if sid not in seq_to_comp:
            seq_to_comp[sid] = comp

    def _seq_sort_key(x: str) -> Tuple[int, str]:
        try:
            return int(x), ""
        except ValueError:
            return -1, x

    ordered = sorted(seq_to_comp.keys(), key=_seq_sort_key)
    etype = _entity_type_for_asym(block, asym_id)
    parts = [_comp_id_to_seq_segment(seq_to_comp[s], etype) for s in ordered]
    return len(ordered), "".join(parts)


PolymerProfile = Tuple[int, str]


def polymer_profiles_by_asym(
    block: gemmi.cif.Block, mode: str
) -> Dict[str, PolymerProfile]:
    """label_asym_id -> (residue_count, atom_site_sequence) for each polymer chain."""
    asym_ids = polymer_asym_ids_for_mode(block, mode)
    return {a: atom_site_sequence_and_count(block, a) for a in asym_ids}


def pair_polymer_chains_by_content(
    reference_profiles: Dict[str, PolymerProfile],
    target_profiles: Dict[str, PolymerProfile],
) -> Optional[Dict[str, str]]:
    """
    Map reference label_asym_id -> target label_asym_id when chains match by content.

    Requires equal chain counts, matching multisets of (count, sequence), and a unique
    1:1 pairing. Returns None if alignment cannot be established.
    """
    if len(reference_profiles) != len(target_profiles):
        return None
    if not reference_profiles and not target_profiles:
        return {}

    ref_multiset = sorted(reference_profiles.values())
    tgt_multiset = sorted(target_profiles.values())
    if ref_multiset != tgt_multiset:
        return None

    pairing: Dict[str, str] = {}
    used_target: Set[str] = set()
    for ref_asym, ref_profile in reference_profiles.items():
        candidates = [
            tgt_asym
            for tgt_asym, tgt_profile in target_profiles.items()
            if tgt_asym not in used_target and tgt_profile == ref_profile
        ]
        if len(candidates) != 1:
            return None
        pairing[ref_asym] = candidates[0]
        used_target.add(candidates[0])
    return pairing


def entity_id_for_asym(block: gemmi.cif.Block, asym_id: str) -> str:
    """Resolve entity_id for label_asym_id from _struct_asym or _atom_site."""
    sm = _struct_asym_map(block)
    if asym_id in sm:
        return sm[asym_id]
    tbl = _atom_site_columns(block)
    if not tbl:
        return ""
    tags, rows = tbl
    ia = _find_column(tags, "_atom_site.label_asym_id")
    ie = _find_column(tags, "_atom_site.label_entity_id")
    ig = _find_column(tags, "_atom_site.group_PDB")
    if ia is None or ie is None:
        return ""
    for row in rows:
        if _cif_string_raw(row[ia]) != asym_id:
            continue
        if ig is not None and _cif_string_raw(row[ig]).upper() == "ATOM":
            eid = _cif_string_raw(row[ie])
            if eid and not _is_missing_cif_value(eid):
                return eid
    for row in rows:
        if _cif_string_raw(row[ia]) == asym_id:
            eid = _cif_string_raw(row[ie])
            if eid and not _is_missing_cif_value(eid):
                return eid
    return ""


def build_polymer_entity_remapping(
    chain_pairing: Dict[str, str],
    reference_block: gemmi.cif.Block,
    target_block: gemmi.cif.Block,
) -> Tuple[Dict[str, List[str]], Set[str]]:
    """
    Map reference polymer entity_id -> ordered target entity_id list (one per paired chain).

    Also returns the set of reference polymer entity ids involved in the pairing.
    """
    ref_entity_to_targets: Dict[str, List[str]] = {}
    ref_polymer_entities: Set[str] = set()
    for ref_asym, tgt_asym in chain_pairing.items():
        re = entity_id_for_asym(reference_block, ref_asym)
        te = entity_id_for_asym(target_block, tgt_asym)
        if not re or not te:
            continue
        ref_polymer_entities.add(re)
        bucket = ref_entity_to_targets.setdefault(re, [])
        if te not in bucket:
            bucket.append(te)
    for re in ref_entity_to_targets:
        ref_entity_to_targets[re] = sorted(ref_entity_to_targets[re])
    return ref_entity_to_targets, ref_polymer_entities


def needs_content_id_remap(chain_pairing: Dict[str, str]) -> bool:
    return any(ref_a != tgt_a for ref_a, tgt_a in chain_pairing.items())


def needs_polymer_metadata_id_remap(
    chain_pairing: Dict[str, str],
    reference_block: gemmi.cif.Block,
    target_block: gemmi.cif.Block,
) -> bool:
    """True when polymer entity/asym ids in metadata must be rewritten for the merge target."""
    if needs_content_id_remap(chain_pairing):
        return True
    ref_entity_to_targets, _ = build_polymer_entity_remapping(
        chain_pairing, reference_block, target_block
    )
    for ref_eid, target_eids in ref_entity_to_targets.items():
        if len(target_eids) != 1 or target_eids[0] != ref_eid:
            return True
    return False


def _is_entity_id_tag(tag: str) -> bool:
    return tag.endswith(".entity_id") or tag == "_entity.id"


def _is_label_asym_tag(tag: str) -> bool:
    if _is_entity_id_tag(tag):
        return False
    return "label_asym_id" in tag or tag.endswith(".asym_id")


def _loop_item_from_table(tags: List[str], rows: List[List[str]]):
    """Return a gemmi.cif.Item wrapping a loop built from tags and rows."""
    if not tags:
        raise ValueError("empty loop tags")
    w = len(tags)
    for row in rows:
        if len(row) != w:
            raise ValueError(f"row width {len(row)} != {w}")
    category = get_category_from_tag(tags[0])
    short_tags = [
        tag.split(".", 1)[1] if tag.startswith("_") and "." in tag else tag
        for tag in tags
    ]
    tmp = gemmi.cif.Block("_loop_build")
    loop = tmp.init_mmcif_loop(category, short_tags)
    for row in rows:
        loop.add_row(row)
    for item in tmp:
        if item.loop is not None:
            return item
    raise ValueError("failed to build loop from table")


def _remap_cell_asym(value: str, chain_pairing: Dict[str, str]) -> str:
    raw = _cif_string_raw(value)
    if raw in chain_pairing:
        return chain_pairing[raw]
    return value


def _expand_rows_for_entity_id(
    row: List[str],
    tags: List[str],
    ref_entity_to_targets: Dict[str, List[str]],
) -> List[List[str]]:
    entity_cols = [i for i, t in enumerate(tags) if _is_entity_id_tag(t)]
    if not entity_cols:
        return [row]
    ref_eid = _cif_string_raw(row[entity_cols[0]])
    targets = ref_entity_to_targets.get(ref_eid)
    if not targets:
        return [row]
    if len(targets) == 1:
        new_row = list(row)
        for ci in entity_cols:
            new_row[ci] = targets[0]
        return [new_row]
    out: List[List[str]] = []
    for te in targets:
        new_row = list(row)
        for ci in entity_cols:
            new_row[ci] = te
        out.append(new_row)
    return out


def _frame_pairs_have_entity_id(pairs: List[Tuple[str, str]]) -> bool:
    return any(_is_entity_id_tag(tag) for tag, _ in pairs)


def _emit_remapped_frame_pairs(
    out: gemmi.cif.Block,
    pairs: List[Tuple[str, str]],
    ref_entity_to_targets: Dict[str, List[str]],
) -> None:
    """Write frame pairs to ``out``, remapping entity_id and expanding to a loop when needed."""
    tags = [tag for tag, _ in pairs]
    row = [val for _, val in pairs]
    new_rows = _expand_rows_for_entity_id(row, tags, ref_entity_to_targets)
    if len(new_rows) == 1:
        for tag, val in zip(tags, new_rows[0]):
            out.set_pair(tag, val)
        return
    out.add_item(_loop_item_from_table(tags, new_rows))


def _flush_frame_pair_buffer(
    out: gemmi.cif.Block,
    pairs: List[Tuple[str, str]],
    category: Optional[str],
    ref_entity_to_targets: Dict[str, List[str]],
) -> None:
    if not pairs:
        return
    if (
        category in MACROMOLECULE_CATEGORIES
        and _frame_pairs_have_entity_id(pairs)
    ):
        _emit_remapped_frame_pairs(out, pairs, ref_entity_to_targets)
        return
    for tag, val in pairs:
        out.set_pair(tag, val)


def _remap_loop_rows(
    tags: List[str],
    rows: List[List[str]],
    category: str,
    chain_pairing: Dict[str, str],
    ref_entity_to_targets: Dict[str, List[str]],
    ref_polymer_entities: Set[str],
) -> List[List[str]]:
    if category == "_entity":
        ie = _find_column(tags, "_entity.id")
        out_rows: List[List[str]] = []
        for row in rows:
            new_row = list(row)
            if ie is not None:
                eid = _cif_string_raw(row[ie])
                if eid in ref_polymer_entities:
                    out_rows.extend(
                        _expand_rows_for_entity_id(new_row, tags, ref_entity_to_targets)
                    )
                    continue
            out_rows.append(new_row)
        return out_rows

    out_rows: List[List[str]] = []
    asym_cols = [i for i, t in enumerate(tags) if _is_label_asym_tag(t)]
    for row in rows:
        new_row = list(row)
        for ci in asym_cols:
            new_row[ci] = _remap_cell_asym(new_row[ci], chain_pairing)
        expanded = _expand_rows_for_entity_id(new_row, tags, ref_entity_to_targets)
        out_rows.extend(expanded)
    return out_rows


def remap_macromolecule_metadata_for_target(
    metadata_block: gemmi.cif.Block,
    reference_block: gemmi.cif.Block,
    target_block: gemmi.cif.Block,
    chain_pairing: Dict[str, str],
) -> gemmi.cif.Block:
    """
    Rewrite macromolecule metadata so entity and asym ids match the merge target coordinates.

    Used after content-based chain pairing (reference asym names differ from target).
    """
    if not needs_polymer_metadata_id_remap(chain_pairing, reference_block, target_block):
        return metadata_block

    ref_entity_to_targets, ref_polymer_entities = build_polymer_entity_remapping(
        chain_pairing, reference_block, target_block
    )

    out = gemmi.cif.Block(metadata_block.name)
    frame_pairs: List[Tuple[str, str]] = []
    frame_category: Optional[str] = None
    for item in metadata_block:
        if item.pair is not None:
            cat = get_category_from_tag(item.pair[0])
            if frame_category is not None and cat != frame_category:
                _flush_frame_pair_buffer(
                    out, frame_pairs, frame_category, ref_entity_to_targets
                )
                frame_pairs = []
                frame_category = None
            if frame_category is None:
                frame_category = cat
            frame_pairs.append((item.pair[0], item.pair[1]))
            continue
        _flush_frame_pair_buffer(
            out, frame_pairs, frame_category, ref_entity_to_targets
        )
        frame_pairs = []
        frame_category = None
        if item.loop is not None and item.loop.tags:
            cat = get_category_from_tag(item.loop.tags[0])
            if cat not in MACROMOLECULE_CATEGORIES:
                out.add_item(item)
                continue
            tags, rows = _loop_as_table(item.loop)
            new_rows = _remap_loop_rows(
                tags,
                rows,
                cat,
                chain_pairing,
                ref_entity_to_targets,
                ref_polymer_entities,
            )
            if new_rows:
                out.add_item(_loop_item_from_table(tags, new_rows))
        else:
            out.add_item(item)
    _flush_frame_pair_buffer(
        out, frame_pairs, frame_category, ref_entity_to_targets
    )
    return out


def reconcile_polymer_struct_asym_in_block(block: gemmi.cif.Block) -> Tuple[gemmi.cif.Block, bool]:
    """
    Align _struct_asym polymer rows with _atom_site label_asym_id / label_entity_id.

    Returns (block, changed). Replaces or adds rows for polymer chains inferred from coordinates.
    """
    desired: List[Tuple[str, str]] = []
    for asym in sorted(polymer_asym_ids_fallback(block)):
        eid = entity_id_for_asym(block, asym)
        if eid:
            desired.append((asym, eid))
    if not desired:
        return block, False

    polymer_asym_set = {a for a, _ in desired}
    loop = _get_loop(block, "_struct_asym.id")
    if loop is None:
        new_item = _loop_item_from_table(
            ["_struct_asym.id", "_struct_asym.entity_id"],
            [[asym, eid] for asym, eid in desired],
        )
        return _block_with_replaced_loop_item(block, "_struct_asym", new_item), True

    tags, rows = _loop_as_table(loop)
    ia = _find_column(tags, "_struct_asym.id")
    ie = _find_column(tags, "_struct_asym.entity_id")
    if ia is None or ie is None:
        new_item = _loop_item_from_table(
            ["_struct_asym.id", "_struct_asym.entity_id"],
            [[asym, eid] for asym, eid in desired],
        )
        return _block_with_replaced_loop_item(block, "_struct_asym", new_item), True

    kept: List[List[str]] = []
    for row in rows:
        asym = _cif_string_raw(row[ia])
        if asym in polymer_asym_set:
            continue
        kept.append(list(row))

    for asym, eid in desired:
        new_row = ["?"] * len(tags)
        new_row[ia] = asym
        new_row[ie] = eid
        kept.append(new_row)

    return _block_with_replaced_loop_item(
        block, "_struct_asym", _loop_item_from_table(tags, kept)
    ), True


def _block_with_replaced_loop_item(block: gemmi.cif.Block, category: str, new_item) -> gemmi.cif.Block:
    """Return a copy of ``block`` with the loop for ``category`` replaced by ``new_item``."""
    cat_tag = category if category.startswith("_") else f"_{category}"
    out = gemmi.cif.Block(block.name)
    replaced = False
    for item in block:
        if item.loop is not None and item.loop.tags:
            if get_category_from_tag(item.loop.tags[0]) == cat_tag:
                if not replaced:
                    out.add_item(new_item)
                    replaced = True
                continue
        if item.pair is not None:
            out.set_pair(item.pair[0], item.pair[1])
        else:
            out.add_item(item)
    if not replaced:
        out.add_item(new_item)
    return out


def reconcile_polymer_struct_asym_in_mmcif_file(mmcif_path: str) -> bool:
    """Load mmCIF, reconcile polymer _struct_asym in the first data block, write back if changed."""
    doc = gemmi.cif.read(mmcif_path)
    if not doc:
        return False
    new_block, changed = reconcile_polymer_struct_asym_in_block(doc[0])
    if not changed:
        return False
    out_doc = gemmi.cif.Document()
    out_doc.add_copied_block(new_block)
    for bi in range(1, len(doc)):
        out_doc.add_copied_block(doc[bi])
    out_doc.write_file(mmcif_path)
    return True


@dataclass
class SafeguardResult:
    ok: bool
    mode: str
    failures: List[Dict[str, Any]] = field(default_factory=list)
    checked_asym_ids: List[str] = field(default_factory=list)
    chain_pairing: Dict[str, str] = field(default_factory=dict)
    content_aligned: bool = False

    def to_json_fragment(self) -> str:
        payload: Dict[str, Any] = {
            "ok": self.ok,
            "mode": self.mode,
            "failures": self.failures,
            "checked_asym_ids": self.checked_asym_ids,
        }
        if self.content_aligned:
            payload["content_aligned"] = True
            payload["chain_pairing"] = self.chain_pairing
        return json.dumps(payload, indent=2)


def validate_macromolecule_merge(reference_block: gemmi.cif.Block, target_block: gemmi.cif.Block) -> SafeguardResult:
    """
    R-ALIGN / R-POLY checks between reference and target first data blocks.
    """
    mode = identification_mode(reference_block, target_block)
    ref_asym = polymer_asym_ids_for_mode(reference_block, mode)
    tgt_asym = polymer_asym_ids_for_mode(target_block, mode)

    failures: List[Dict[str, Any]] = []
    content_aligned = False
    chain_pairing: Dict[str, str] = {}

    if ref_asym == tgt_asym:
        chain_pairing = {a: a for a in ref_asym}
    else:
        ref_profiles = polymer_profiles_by_asym(reference_block, mode)
        tgt_profiles = polymer_profiles_by_asym(target_block, mode)
        paired = pair_polymer_chains_by_content(ref_profiles, tgt_profiles)
        if paired is None:
            failures.append(
                {
                    "rule": "ALIGN-1-CONTENT-MISMATCH",
                    "mode": mode,
                    "reference_asym_ids": sorted(ref_asym),
                    "target_asym_ids": sorted(tgt_asym),
                    "reference_chain_count": len(ref_asym),
                    "target_chain_count": len(tgt_asym),
                }
            )
            return SafeguardResult(
                False,
                mode,
                failures,
                sorted(ref_asym | tgt_asym),
                {},
                False,
            )
        chain_pairing = paired
        content_aligned = True

    asym_ids = sorted(ref_asym)
    if not asym_ids:
        return SafeguardResult(True, mode, [], [], chain_pairing, content_aligned)

    sm_r = _struct_asym_map(reference_block)
    sm_t = _struct_asym_map(target_block)

    for ref_asym_id in asym_ids:
        tgt_asym_id = chain_pairing[ref_asym_id]
        cr, sr = atom_site_sequence_and_count(reference_block, ref_asym_id)
        ct, st = atom_site_sequence_and_count(target_block, tgt_asym_id)
        if cr != ct or sr != st:
            entry: Dict[str, Any] = {
                "rule": "ALIGN-2-LENGTH-OR-ATOM-SEQ",
                "label_asym_id": ref_asym_id,
                "reference_count": cr,
                "target_count": ct,
                "reference_atom_site_seq": sr,
                "target_atom_site_seq": st,
            }
            if content_aligned:
                entry["target_label_asym_id"] = tgt_asym_id
            failures.append(entry)
            continue

        if mode == "entity":
            er = sm_r.get(ref_asym_id)
            et = sm_t.get(tgt_asym_id)
            if er and et:
                pr = _entity_poly_one_letter(reference_block, er)
                pt = _entity_poly_one_letter(target_block, et)
                if pr is not None and sr != pr:
                    failures.append(
                        {
                            "rule": "ALIGN-3-REF-ATOM-VS-ENTITY-POLY",
                            "label_asym_id": ref_asym_id,
                            "target_label_asym_id": tgt_asym_id,
                            "entity_id": er,
                            "atom_site_seq": sr,
                            "entity_poly_seq": pr,
                        }
                    )
                if pt is not None and st != pt:
                    failures.append(
                        {
                            "rule": "ALIGN-3-TGT-ATOM-VS-ENTITY-POLY",
                            "label_asym_id": ref_asym_id,
                            "target_label_asym_id": tgt_asym_id,
                            "entity_id": et,
                            "atom_site_seq": st,
                            "entity_poly_seq": pt,
                        }
                    )
                if pr is not None and pt is not None and pr != pt:
                    failures.append(
                        {
                            "rule": "ALIGN-3-ENTITY-POLY-MISMATCH",
                            "label_asym_id": ref_asym_id,
                            "target_label_asym_id": tgt_asym_id,
                            "reference_entity_poly": pr,
                            "target_entity_poly": pt,
                        }
                    )

    return SafeguardResult(
        len(failures) == 0, mode, failures, asym_ids, chain_pairing, content_aligned
    )


def strip_macromolecule_categories(block: gemmi.cif.Block) -> gemmi.cif.Block:
    """Return a new block without categories in MACROMOLECULE_CATEGORIES."""
    out = gemmi.cif.Block(block.name)
    for item in block:
        if item.pair is not None:
            cat = get_category_from_tag(item.pair[0])
            if cat in MACROMOLECULE_CATEGORIES:
                continue
            out.set_pair(item.pair[0], item.pair[1])
        elif item.loop is not None:
            loop = item.loop
            if not loop.tags:
                continue
            cat = get_category_from_tag(loop.tags[0])
            if cat in MACROMOLECULE_CATEGORIES:
                continue
            out.add_item(item)
    return out


def get_category_from_tag(item_name: str) -> str:
    if "." in item_name:
        return item_name.split(".", 1)[0]
    return item_name


def macromolecule_spec_in_use(spec_files: Iterable[str]) -> bool:
    for p in spec_files:
        if "MACROMOLECULES" in Path(p).name.upper():
            return True
    return False


def block_has_items(block: gemmi.cif.Block) -> bool:
    for _ in block:
        return True
    return False
