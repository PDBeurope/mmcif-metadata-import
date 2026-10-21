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


def _atom_site_seq_id_from_row(
    row: List[str],
    label_seq_idx: Optional[int],
    auth_seq_idx: Optional[int],
) -> str:
    """
    Prefer ``_atom_site.label_seq_id`` when present; otherwise ``auth_seq_id``.

    PDB-style / early-upload mmCIF often has ``label_seq_id`` as ``.`` while residue
    numbers live only in ``auth_seq_id``.
    """
    if label_seq_idx is not None:
        seq = _cif_string_raw(row[label_seq_idx])
        if not _is_missing_cif_value(seq):
            return seq
    if auth_seq_idx is not None:
        seq = _cif_string_raw(row[auth_seq_idx])
        if not _is_missing_cif_value(seq):
            return seq
    return ""


def polymer_asym_ids_fallback(block: gemmi.cif.Block) -> Set[str]:
    """R-POLY-2: _atom_site heuristic (length >= 2 distinct seq ids, >= one ATOM).

    Sequence positions use ``label_seq_id`` when set; otherwise ``auth_seq_id``.
    """
    tbl = _atom_site_columns(block)
    if not tbl:
        return set()
    tags, rows = tbl
    ia = _find_column(tags, "_atom_site.label_asym_id")
    iseq = _find_column(tags, "_atom_site.label_seq_id")
    iauth = _find_column(tags, "_atom_site.auth_seq_id")
    ig = _find_column(tags, "_atom_site.group_PDB")
    if ia is None or ig is None:
        return set()
    if iseq is None and iauth is None:
        return set()

    by_asym: Dict[str, Set[str]] = {}
    all_het: Dict[str, bool] = {}

    for row in rows:
        asym = _cif_string_raw(row[ia])
        if not asym:
            continue
        seq = _atom_site_seq_id_from_row(row, iseq, iauth)
        g = _cif_string_raw(row[ig]).upper()
        by_asym.setdefault(asym, set())
        all_het.setdefault(asym, True)
        if g == "ATOM":
            all_het[asym] = False
        if seq:
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
    Distinct residue-position count and canonical pdbx-style one-letter string from _atom_site.

    Positions use ``label_seq_id`` when set; otherwise ``auth_seq_id`` (PDB-style mmCIF).
    One residue per seq id; first non-missing label_comp_id wins (microheterogeneity).
    """
    tbl = _atom_site_columns(block)
    if not tbl:
        return 0, ""
    tags, rows = tbl
    ia = _find_column(tags, "_atom_site.label_asym_id")
    iseq = _find_column(tags, "_atom_site.label_seq_id")
    iauth = _find_column(tags, "_atom_site.auth_seq_id")
    ic = _find_column(tags, "_atom_site.label_comp_id")
    if ia is None or ic is None:
        return 0, ""
    if iseq is None and iauth is None:
        return 0, ""

    seq_to_comp: Dict[str, str] = {}
    for row in rows:
        if _cif_string_raw(row[ia]) != asym_id:
            continue
        sid = _atom_site_seq_id_from_row(row, iseq, iauth)
        if not sid:
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


def pair_polymer_chains_partial(
    reference_profiles: Dict[str, PolymerProfile],
    target_profiles: Dict[str, PolymerProfile],
) -> Tuple[Dict[str, str], List[str], List[str]]:
    """
    Pair chains that have a unique matching polymer profile.

    Returns ``(pairing, ambiguous_reference_asyms, ambiguous_target_asyms)``.
    A profile with exactly one reference chain and one target chain is paired.
    A profile present on both sides but not uniquely 1:1 is ambiguous (those chains
    are not paired). Chains whose profile exists on only one side are left unmatched
    and are not listed as ambiguous.
    """
    ref_by: Dict[PolymerProfile, List[str]] = {}
    tgt_by: Dict[PolymerProfile, List[str]] = {}
    for asym, profile in reference_profiles.items():
        ref_by.setdefault(profile, []).append(asym)
    for asym, profile in target_profiles.items():
        tgt_by.setdefault(profile, []).append(asym)

    pairing: Dict[str, str] = {}
    ambiguous_ref: List[str] = []
    ambiguous_tgt: List[str] = []
    for profile in set(ref_by) | set(tgt_by):
        refs = ref_by.get(profile, [])
        tgts = tgt_by.get(profile, [])
        if len(refs) == 1 and len(tgts) == 1:
            pairing[refs[0]] = tgts[0]
        elif refs and tgts:
            ambiguous_ref.extend(refs)
            ambiguous_tgt.extend(tgts)
    return pairing, ambiguous_ref, ambiguous_tgt


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


def _is_entity_id_list_tag(tag: str) -> bool:
    return tag.endswith("entity_id_list")


def _is_entity_assembly_fk_tag(tag: str) -> bool:
    return tag.endswith(".entity_assembly_id")


def _split_id_list(value: str) -> List[str]:
    raw = _cif_string_raw(value)
    if _is_missing_cif_value(raw):
        return []
    out: List[str] = []
    for chunk in raw.replace(";", ",").split(","):
        part = chunk.strip()
        if part and not _is_missing_cif_value(part):
            out.append(part)
    return out


def _remap_entity_id_list_value(
    value: str, ref_entity_to_targets: Dict[str, List[str]]
) -> str:
    ids = _split_id_list(value)
    if not ids:
        return value
    rewritten: List[str] = []
    for eid in ids:
        targets = ref_entity_to_targets.get(eid)
        if targets:
            for te in targets:
                if te not in rewritten:
                    rewritten.append(te)
        elif eid not in rewritten:
            rewritten.append(eid)
    return ", ".join(rewritten)


def _rewrite_entity_id_lists(
    row: List[str], tags: List[str], ref_entity_to_targets: Dict[str, List[str]]
) -> None:
    for i, tag in enumerate(tags):
        if _is_entity_id_list_tag(tag):
            row[i] = _remap_entity_id_list_value(row[i], ref_entity_to_targets)


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
    for expanded in new_rows:
        _rewrite_entity_id_lists(expanded, tags, ref_entity_to_targets)
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
        for expanded_row in expanded:
            _rewrite_entity_id_lists(expanded_row, tags, ref_entity_to_targets)
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


def _polymer_entity_groups(
    block: gemmi.cif.Block, asym_ids: Iterable[str]
) -> Dict[str, Set[str]]:
    """Polymer asym ids grouped by entity_id. Chains with no entity id are their own group."""
    groups: Dict[str, Set[str]] = {}
    for asym in asym_ids:
        eid = entity_id_for_asym(block, asym)
        key = eid if eid else f"__asym__:{asym}"
        groups.setdefault(key, set()).add(asym)
    return groups


def _add_table_row(
    tables: Dict[str, Tuple[List[str], List[List[str]]]],
    category: str,
    tags: List[str],
    row: List[str],
) -> None:
    if category not in tables:
        tables[category] = (list(tags), [list(row)])
        return
    old_tags, rows = tables[category]
    if list(old_tags) == list(tags):
        rows.append(list(row))
        return
    union = list(old_tags)
    for tag in tags:
        if tag not in union:
            union.append(tag)

    def _align(src_tags: List[str], src_row: List[str]) -> List[str]:
        idx = {tag: i for i, tag in enumerate(src_tags)}
        return [src_row[idx[tag]] if tag in idx else "?" for tag in union]

    tables[category] = (
        union,
        [_align(old_tags, existing) for existing in rows] + [_align(list(tags), row)],
    )


def _macro_tables(block: gemmi.cif.Block) -> Dict[str, Tuple[List[str], List[List[str]]]]:
    """Macromolecule categories as (tags, rows). Frame pairs become a single row."""
    tables: Dict[str, Tuple[List[str], List[List[str]]]] = {}
    frame_cat: Optional[str] = None
    frame_pairs: List[Tuple[str, str]] = []

    def _flush() -> None:
        nonlocal frame_cat, frame_pairs
        if frame_pairs and frame_cat in MACROMOLECULE_CATEGORIES:
            tags = [tag for tag, _ in frame_pairs]
            row = [val for _, val in frame_pairs]
            _add_table_row(tables, frame_cat, tags, row)
        frame_pairs = []
        frame_cat = None

    for item in block:
        if item.pair is not None:
            cat = get_category_from_tag(item.pair[0])
            if cat in MACROMOLECULE_CATEGORIES:
                if frame_cat is not None and cat != frame_cat:
                    _flush()
                frame_cat = cat
                frame_pairs.append((item.pair[0], item.pair[1]))
            else:
                _flush()
            continue
        _flush()
        if item.loop is not None and item.loop.tags:
            cat = get_category_from_tag(item.loop.tags[0])
            if cat in MACROMOLECULE_CATEGORIES:
                tags, rows = _loop_as_table(item.loop)
                for row in rows:
                    _add_table_row(tables, cat, tags, row)
    _flush()
    return tables


def _row_entity_ids(tags: List[str], row: List[str]) -> Set[str]:
    ids: Set[str] = set()
    for i, tag in enumerate(tags):
        if _is_entity_id_tag(tag):
            value = _cif_string_raw(row[i])
            if value and not _is_missing_cif_value(value):
                ids.add(value)
        elif _is_entity_id_list_tag(tag):
            ids.update(_split_id_list(row[i]))
    return ids


def _row_mentions_entities(tags: List[str], row: List[str], entity_ids: Set[str]) -> bool:
    return bool(_row_entity_ids(tags, row) & entity_ids)


def _fresh_numeric_id(used: Set[str]) -> str:
    n = 1
    while str(n) in used:
        n += 1
    return str(n)


def _rekey_import_ids(
    import_rows: List[List[str]],
    tags: List[str],
    id_tag: str,
    reserved: Set[str],
) -> Dict[str, str]:
    """Rename import ids that collide with ``reserved``. Returns old -> new for changed ids."""
    col = _find_column(tags, id_tag)
    if col is None:
        return {}
    mapping: Dict[str, str] = {}
    used = set(reserved)
    for row in import_rows:
        current = _cif_string_raw(row[col])
        if not current or _is_missing_cif_value(current):
            continue
        if current in mapping:
            row[col] = mapping[current]
            continue
        if current in used:
            new_id = _fresh_numeric_id(used)
            mapping[current] = new_id
            row[col] = new_id
            used.add(new_id)
        else:
            used.add(current)
    return mapping


def _apply_id_map(
    rows: List[List[str]], tags: List[str], id_tag: str, mapping: Dict[str, str]
) -> None:
    if not mapping:
        return
    col = _find_column(tags, id_tag)
    if col is None:
        return
    for row in rows:
        current = _cif_string_raw(row[col])
        if current in mapping:
            row[col] = mapping[current]


def _align_rows_to_tags(
    dst_tags: List[str], src_tags: List[str], src_rows: List[List[str]]
) -> List[List[str]]:
    idx = {tag: i for i, tag in enumerate(src_tags)}
    aligned: List[List[str]] = []
    for row in src_rows:
        aligned.append([row[idx[tag]] if tag in idx else "?" for tag in dst_tags])
    return aligned


def _union_tags(left: List[str], right: List[str]) -> List[str]:
    out = list(left)
    for tag in right:
        if tag not in out:
            out.append(tag)
    return out


def categories_in_block(block: gemmi.cif.Block) -> Set[str]:
    cats: Set[str] = set()
    for item in block:
        if item.pair is not None:
            cats.add(get_category_from_tag(item.pair[0]))
        elif item.loop is not None and item.loop.tags:
            cats.add(get_category_from_tag(item.loop.tags[0]))
    return cats


def _entity_row_decision(
    tags: List[str],
    row: List[str],
    eligible_entities: Set[str],
    eligible_asyms: Set[str],
    polymer_entities: Set[str],
) -> Optional[bool]:
    """True = keep, False = drop, None = row has no entity or asym key."""
    saw = False
    for i, tag in enumerate(tags):
        if _is_entity_id_list_tag(tag):
            saw = True
            for eid in _split_id_list(row[i]):
                if eid in polymer_entities and eid not in eligible_entities:
                    return False
        elif _is_entity_id_tag(tag):
            saw = True
            eid = _cif_string_raw(row[i])
            if _is_missing_cif_value(eid) or eid not in eligible_entities:
                return False
        elif _is_label_asym_tag(tag):
            saw = True
            asym = _cif_string_raw(row[i])
            if _is_missing_cif_value(asym) or asym not in eligible_asyms:
                return False
    if not saw:
        return None
    return True


def restrict_macromolecule_metadata_to_entities(
    block: gemmi.cif.Block,
    eligible_ref_entity_ids: Iterable[str],
    eligible_ref_asym_ids: Iterable[str],
    ref_polymer_entity_ids: Iterable[str],
) -> gemmi.cif.Block:
    """
    Keep macromolecule rows for eligible reference polymer entities only.

    Rows that cannot be tied to an entity or chain are omitted (fail-closed), including
    EM assembly rows whose ``entity_id_list`` names a polymer entity that is not eligible.
    Child assembly rows are kept only when their parent assembly row was kept.
    ``_struct_ref_seq`` / ``_struct_ref_seq_dif`` follow the kept ``_struct_ref`` rows.
    Non-macromolecule categories are copied unchanged.
    """
    eligible = {eid for eid in eligible_ref_entity_ids if eid and not str(eid).startswith("__asym__:")}
    asyms = set(eligible_ref_asym_ids)
    polymers = {eid for eid in ref_polymer_entity_ids if eid and not str(eid).startswith("__asym__:")}
    tables = _macro_tables(block)
    kept: Dict[str, Tuple[List[str], List[List[str]]]] = {}

    for cat, (tags, rows) in tables.items():
        if cat in ("_struct_ref_seq", "_struct_ref_seq_dif"):
            continue
        filtered: List[List[str]] = []
        keyed = False
        for row in rows:
            decision = _entity_row_decision(tags, row, eligible, asyms, polymers)
            if decision is None:
                continue
            keyed = True
            if decision:
                filtered.append(row)
        if keyed:
            kept[cat] = (tags, filtered)

    kept_ref_ids: Set[str] = set()
    if "_struct_ref" in kept:
        tags, rows = kept["_struct_ref"]
        id_col = _find_column(tags, "_struct_ref.id")
        if id_col is not None:
            kept_ref_ids = {
                _cif_string_raw(row[id_col])
                for row in rows
                if not _is_missing_cif_value(_cif_string_raw(row[id_col]))
            }
    if "_struct_ref_seq" in tables and kept_ref_ids:
        tags, rows = tables["_struct_ref_seq"]
        ref_col = _find_column(tags, "_struct_ref_seq.ref_id")
        if ref_col is not None:
            seq_rows = [row for row in rows if _cif_string_raw(row[ref_col]) in kept_ref_ids]
            kept["_struct_ref_seq"] = (tags, seq_rows)
            align_col = _find_column(tags, "_struct_ref_seq.align_id")
            align_ids = {
                _cif_string_raw(row[align_col])
                for row in seq_rows
                if align_col is not None and not _is_missing_cif_value(_cif_string_raw(row[align_col]))
            }
            if "_struct_ref_seq_dif" in tables and align_ids:
                dtags, drows = tables["_struct_ref_seq_dif"]
                dif_col = _find_column(dtags, "_struct_ref_seq_dif.align_id")
                if dif_col is not None:
                    kept["_struct_ref_seq_dif"] = (
                        dtags,
                        [row for row in drows if _cif_string_raw(row[dif_col]) in align_ids],
                    )

    assembly_ids: Set[str] = set()
    if "_em_entity_assembly" in kept:
        tags, rows = kept["_em_entity_assembly"]
        id_col = _find_column(tags, "_em_entity_assembly.id")
        if id_col is not None:
            assembly_ids = {
                _cif_string_raw(row[id_col])
                for row in rows
                if not _is_missing_cif_value(_cif_string_raw(row[id_col]))
            }
    for cat, (tags, rows) in tables.items():
        if cat in kept or cat.startswith("_struct_ref"):
            continue
        fk_cols = [i for i, tag in enumerate(tags) if _is_entity_assembly_fk_tag(tag)]
        if not fk_cols or not assembly_ids:
            continue
        filtered = []
        for row in rows:
            fks = [_cif_string_raw(row[i]) for i in fk_cols]
            if fks and all(fk in assembly_ids for fk in fks):
                filtered.append(row)
        if filtered:
            kept[cat] = (tags, filtered)

    out = gemmi.cif.Block(block.name)
    emitted: Set[str] = set()

    def _emit(cat: str) -> None:
        if cat in emitted:
            return
        emitted.add(cat)
        packed = kept.get(cat)
        if not packed:
            return
        tags, rows = packed
        if rows:
            out.add_item(_loop_item_from_table(tags, rows))

    for item in block:
        if item.pair is not None:
            cat = get_category_from_tag(item.pair[0])
            if cat in MACROMOLECULE_CATEGORIES:
                _emit(cat)
                continue
            out.set_pair(item.pair[0], item.pair[1])
            continue
        if item.loop is not None and item.loop.tags:
            cat = get_category_from_tag(item.loop.tags[0])
            if cat in MACROMOLECULE_CATEGORIES:
                _emit(cat)
                continue
            out.add_item(item)
            continue
        out.add_item(item)
    return out


def merge_partial_macromolecule_rows(
    target_block: gemmi.cif.Block,
    metadata_block: gemmi.cif.Block,
    eligible_target_entity_ids: Iterable[str],
) -> gemmi.cif.Block:
    """
    Replace target macromolecule rows for ``eligible_target_entity_ids`` with metadata rows.

    Other polymer rows, ligand/water rows, and categories that are not in the metadata
    block stay on the target. ``_struct_ref`` ids and EM assembly ids that would collide
    with a kept target row are renumbered, and child rows follow.
    """
    eligible = {eid for eid in eligible_target_entity_ids if eid}
    target_tables = _macro_tables(target_block)
    import_tables = _macro_tables(metadata_block)
    if not import_tables:
        return target_block

    def _keep_target_row(
        cat: str,
        tags: List[str],
        row: List[str],
        drop_ref_ids: Set[str],
        drop_align_ids: Set[str],
        drop_asm_ids: Set[str],
    ) -> bool:
        if _row_mentions_entities(tags, row, eligible):
            return False
        if cat == "_struct_ref_seq" and drop_ref_ids:
            col = _find_column(tags, "_struct_ref_seq.ref_id")
            if col is not None and _cif_string_raw(row[col]) in drop_ref_ids:
                return False
        if cat == "_struct_ref_seq_dif" and drop_align_ids:
            col = _find_column(tags, "_struct_ref_seq_dif.align_id")
            if col is not None and _cif_string_raw(row[col]) in drop_align_ids:
                return False
        if drop_asm_ids:
            for i, tag in enumerate(tags):
                if _is_entity_assembly_fk_tag(tag) and _cif_string_raw(row[i]) in drop_asm_ids:
                    return False
        return True

    drop_ref_ids: Set[str] = set()
    if "_struct_ref" in import_tables and "_struct_ref" in target_tables:
        tags, rows = target_tables["_struct_ref"]
        ent_col = _find_column(tags, "_struct_ref.entity_id")
        id_col = _find_column(tags, "_struct_ref.id")
        if ent_col is not None and id_col is not None:
            for row in rows:
                if _cif_string_raw(row[ent_col]) in eligible:
                    drop_ref_ids.add(_cif_string_raw(row[id_col]))

    drop_align_ids: Set[str] = set()
    if drop_ref_ids and "_struct_ref_seq" in target_tables:
        tags, rows = target_tables["_struct_ref_seq"]
        ref_col = _find_column(tags, "_struct_ref_seq.ref_id")
        align_col = _find_column(tags, "_struct_ref_seq.align_id")
        if ref_col is not None:
            for row in rows:
                if _cif_string_raw(row[ref_col]) in drop_ref_ids and align_col is not None:
                    drop_align_ids.add(_cif_string_raw(row[align_col]))

    drop_asm_ids: Set[str] = set()
    if "_em_entity_assembly" in import_tables and "_em_entity_assembly" in target_tables:
        tags, rows = target_tables["_em_entity_assembly"]
        id_col = _find_column(tags, "_em_entity_assembly.id")
        if id_col is not None:
            for row in rows:
                if _row_mentions_entities(tags, row, eligible):
                    drop_asm_ids.add(_cif_string_raw(row[id_col]))

    # Copy import rows so rekeying does not mutate the metadata block's collected tables
    # if the caller reuses them. _macro_tables already copied cells.
    import_rows_by_cat: Dict[str, Tuple[List[str], List[List[str]]]] = {
        cat: (tags, [list(row) for row in rows]) for cat, (tags, rows) in import_tables.items()
    }
    if not eligible:
        # No target entity id to replace. Adding rows into an existing category would duplicate it.
        import_rows_by_cat = {
            cat: packed
            for cat, packed in import_rows_by_cat.items()
            if cat not in target_tables and packed[1]
        }

    if "_struct_ref" in import_rows_by_cat:
        tags, rows = import_rows_by_cat["_struct_ref"]
        reserved: Set[str] = set()
        if "_struct_ref" in target_tables:
            tt, tr = target_tables["_struct_ref"]
            id_col = _find_column(tt, "_struct_ref.id")
            if id_col is not None:
                reserved = {
                    _cif_string_raw(row[id_col])
                    for row in tr
                    if _keep_target_row("_struct_ref", tt, row, drop_ref_ids, drop_align_ids, drop_asm_ids)
                }
        mapping = _rekey_import_ids(rows, tags, "_struct_ref.id", reserved)
        if "_struct_ref_seq" in import_rows_by_cat:
            st, sr = import_rows_by_cat["_struct_ref_seq"]
            _apply_id_map(sr, st, "_struct_ref_seq.ref_id", mapping)

    if "_struct_ref_seq" in import_rows_by_cat:
        tags, rows = import_rows_by_cat["_struct_ref_seq"]
        reserved = set()
        if "_struct_ref_seq" in target_tables:
            tt, tr = target_tables["_struct_ref_seq"]
            align_col = _find_column(tt, "_struct_ref_seq.align_id")
            if align_col is not None:
                reserved = {
                    _cif_string_raw(row[align_col])
                    for row in tr
                    if _keep_target_row("_struct_ref_seq", tt, row, drop_ref_ids, drop_align_ids, drop_asm_ids)
                }
        mapping = _rekey_import_ids(rows, tags, "_struct_ref_seq.align_id", reserved)
        if "_struct_ref_seq_dif" in import_rows_by_cat:
            dt, dr = import_rows_by_cat["_struct_ref_seq_dif"]
            _apply_id_map(dr, dt, "_struct_ref_seq_dif.align_id", mapping)

    if "_em_entity_assembly" in import_rows_by_cat:
        tags, rows = import_rows_by_cat["_em_entity_assembly"]
        reserved = set()
        if "_em_entity_assembly" in target_tables:
            tt, tr = target_tables["_em_entity_assembly"]
            id_col = _find_column(tt, "_em_entity_assembly.id")
            if id_col is not None:
                reserved = {
                    _cif_string_raw(row[id_col])
                    for row in tr
                    if _keep_target_row("_em_entity_assembly", tt, row, drop_ref_ids, drop_align_ids, drop_asm_ids)
                }
        mapping = _rekey_import_ids(rows, tags, "_em_entity_assembly.id", reserved)
        for cat, (ctags, crows) in import_rows_by_cat.items():
            if cat == "_em_entity_assembly":
                continue
            for tag in ctags:
                if _is_entity_assembly_fk_tag(tag):
                    _apply_id_map(crows, ctags, tag, mapping)

    merged: Dict[str, Tuple[List[str], List[List[str]]]] = {}
    cats_to_merge = set(import_rows_by_cat)
    if drop_ref_ids and "_struct_ref_seq" in target_tables:
        cats_to_merge.add("_struct_ref_seq")
    if drop_align_ids and "_struct_ref_seq_dif" in target_tables:
        cats_to_merge.add("_struct_ref_seq_dif")
    if drop_asm_ids:
        for cat, (tags, _rows) in target_tables.items():
            if any(_is_entity_assembly_fk_tag(tag) for tag in tags):
                cats_to_merge.add(cat)

    for cat in cats_to_merge:
        import_packed = import_rows_by_cat.get(cat)
        target_packed = target_tables.get(cat)
        kept_rows: List[List[str]] = []
        target_tags: List[str] = []
        if target_packed:
            target_tags, target_rows = target_packed
            kept_rows = [
                list(row)
                for row in target_rows
                if _keep_target_row(cat, target_tags, row, drop_ref_ids, drop_align_ids, drop_asm_ids)
            ]
        if import_packed:
            import_tags, import_rows = import_packed
        else:
            import_tags, import_rows = [], []
        if not import_rows and len(kept_rows) == (len(target_packed[1]) if target_packed else 0):
            continue
        tags = _union_tags(target_tags, import_tags)
        rows = _align_rows_to_tags(tags, target_tags, kept_rows) if target_tags else []
        if import_rows:
            rows.extend(_align_rows_to_tags(tags, import_tags, import_rows))
        merged[cat] = (tags, rows)

    if not merged:
        return target_block

    out = gemmi.cif.Block(target_block.name)
    emitted: Set[str] = set()

    def _emit(cat: str) -> None:
        if cat in emitted:
            return
        emitted.add(cat)
        packed = merged.get(cat)
        if not packed:
            return
        tags, rows = packed
        if rows:
            out.add_item(_loop_item_from_table(tags, rows))

    for item in target_block:
        if item.pair is not None:
            cat = get_category_from_tag(item.pair[0])
            if cat in merged:
                _emit(cat)
                continue
            out.set_pair(item.pair[0], item.pair[1])
            continue
        if item.loop is not None and item.loop.tags:
            cat = get_category_from_tag(item.loop.tags[0])
            if cat in merged:
                _emit(cat)
                continue
            out.add_item(item)
            continue
        out.add_item(item)
    for cat in import_rows_by_cat:
        if cat not in emitted and cat in merged:
            _emit(cat)
    return out


@dataclass
class SafeguardResult:
    ok: bool
    mode: str
    failures: List[Dict[str, Any]] = field(default_factory=list)
    checked_asym_ids: List[str] = field(default_factory=list)
    chain_pairing: Dict[str, str] = field(default_factory=dict)
    content_aligned: bool = False
    partial: bool = False
    imported_asym_ids: List[str] = field(default_factory=list)
    skipped_asym_ids: List[str] = field(default_factory=list)
    imported_entity_ids: List[str] = field(default_factory=list)
    skipped_entity_ids: List[str] = field(default_factory=list)
    # Reference polymer entity ids eligible for import (not written to the log JSON).
    eligible_ref_entity_ids: List[str] = field(default_factory=list)
    ref_polymer_entity_ids: List[str] = field(default_factory=list)
    omitted_categories: List[str] = field(default_factory=list)

    def to_json_fragment(self) -> str:
        payload: Dict[str, Any] = {
            "ok": self.ok,
            "mode": self.mode,
            "failures": self.failures,
            "checked_asym_ids": self.checked_asym_ids,
        }
        if self.partial:
            payload["partial"] = True
        if self.content_aligned:
            payload["content_aligned"] = True
            payload["chain_pairing"] = self.chain_pairing
        if (
            self.partial
            or self.imported_asym_ids
            or self.skipped_asym_ids
            or self.imported_entity_ids
            or self.skipped_entity_ids
        ):
            payload["imported_asym_ids"] = self.imported_asym_ids
            payload["skipped_asym_ids"] = self.skipped_asym_ids
            payload["imported_entity_ids"] = self.imported_entity_ids
            payload["skipped_entity_ids"] = self.skipped_entity_ids
        if self.omitted_categories:
            payload["omitted_categories"] = self.omitted_categories
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

    if not ref_asym and not tgt_asym:
        chain_pairing = {}
    elif ref_asym == tgt_asym:
        chain_pairing = {a: a for a in ref_asym}
    else:
        # Shared chain IDs are always compared to themselves. Content pairing is only
        # for chains that exist on one side, so a same-ID mismatch cannot be reassigned
        # to a different chain that happens to have the same sequence.
        shared = ref_asym & tgt_asym
        chain_pairing = {a: a for a in shared}
        ref_profiles = polymer_profiles_by_asym(reference_block, mode)
        tgt_profiles = polymer_profiles_by_asym(target_block, mode)
        extra_ref = {a: ref_profiles[a] for a in ref_asym - shared}
        extra_tgt = {a: tgt_profiles[a] for a in tgt_asym - shared}
        extra_pairs, ambiguous_ref, _ambiguous_tgt = pair_polymer_chains_partial(
            extra_ref, extra_tgt
        )
        for ambiguous_asym in sorted(ambiguous_ref):
            failures.append(
                {
                    "rule": "ALIGN-1-AMBIGUOUS-PROFILE",
                    "label_asym_id": ambiguous_asym,
                    "note": "Polymer profile is not unique; chain was not paired.",
                }
            )
        chain_pairing.update(extra_pairs)
        if not chain_pairing:
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
                False,
                [],
                sorted(ref_asym | tgt_asym),
            )
        content_aligned = any(ref_id != tgt_id for ref_id, tgt_id in chain_pairing.items())

    asym_ids = sorted(ref_asym)
    if not asym_ids:
        # Fail closed: empty polymer sets must not vacuous-pass and allow wholesale
        # macromolecule overwrite (e.g. blank label_seq_id before auth_seq_id fallback,
        # or session CIF without usable coordinates). EM map-only blind copy is handled
        # upstream in import_metadata and never reaches this function.
        failures.append(
            {
                "rule": "ALIGN-0-NO-POLYMER-CHAINS",
                "mode": mode,
                "reference_asym_ids": sorted(ref_asym),
                "target_asym_ids": sorted(tgt_asym),
                "note": (
                    "No polymer chains detected for safeguard comparison; "
                    "macromolecule import blocked."
                ),
            }
        )
        return SafeguardResult(False, mode, failures, [], {}, False)

    sm_r = _struct_asym_map(reference_block)
    sm_t = _struct_asym_map(target_block)
    passed_pairs: Dict[str, str] = {}

    for ref_asym_id in sorted(chain_pairing):
        tgt_asym_id = chain_pairing[ref_asym_id]
        chain_failed = False
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
            ref_ent = entity_id_for_asym(reference_block, ref_asym_id)
            tgt_ent = entity_id_for_asym(target_block, tgt_asym_id)
            if ref_ent:
                entry["reference_entity_id"] = ref_ent
            if tgt_ent:
                entry["target_entity_id"] = tgt_ent
            failures.append(entry)
            chain_failed = True

        if not chain_failed and mode == "entity":
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
                    chain_failed = True
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
                    chain_failed = True
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
                    chain_failed = True
        if not chain_failed:
            passed_pairs[ref_asym_id] = tgt_asym_id

    ref_groups = _polymer_entity_groups(reference_block, ref_asym)
    tgt_groups = _polymer_entity_groups(target_block, tgt_asym)
    ref_asym_entity = {asym: eid for eid, chains in ref_groups.items() for asym in chains}
    tgt_asym_entity = {asym: eid for eid, chains in tgt_groups.items() for asym in chains}
    ref_polymer_entity_ids = sorted(
        eid for eid in ref_groups if not str(eid).startswith("__asym__:")
    )

    eligible_pairs: Dict[str, str] = {}
    eligible_ref_entities: List[str] = []
    for ref_eid, ref_chains in ref_groups.items():
        if any(chain not in passed_pairs for chain in ref_chains):
            continue
        target_eids = {tgt_asym_entity[passed_pairs[chain]] for chain in ref_chains}
        entity_ok = True
        for tgt_eid in target_eids:
            for tgt_chain in tgt_groups.get(tgt_eid, set()):
                sources = [ref_chain for ref_chain, paired_tgt in passed_pairs.items() if paired_tgt == tgt_chain]
                if len(sources) != 1 or ref_asym_entity.get(sources[0]) != ref_eid:
                    entity_ok = False
                    break
            if not entity_ok:
                break
        if not entity_ok:
            continue
        eligible_ref_entities.append(ref_eid)
        for chain in ref_chains:
            eligible_pairs[chain] = passed_pairs[chain]

    for ref_chain, tgt_chain in sorted(passed_pairs.items()):
        if ref_chain in eligible_pairs:
            continue
        failures.append(
            {
                "rule": "ALIGN-ENTITY-INCOMPLETE",
                "label_asym_id": ref_chain,
                "target_label_asym_id": tgt_chain,
                "entity_id": ref_asym_entity.get(ref_chain, ""),
                "note": (
                    "Polymer matched, but another chain of the same entity did not; "
                    "entity metadata was not imported."
                ),
            }
        )

    def _public_entity_id(entity_id: str) -> bool:
        return bool(entity_id) and not str(entity_id).startswith("__asym__:")

    imported_entity_ids = sorted(
        {
            tgt_asym_entity[tgt_chain]
            for tgt_chain in eligible_pairs.values()
            if _public_entity_id(tgt_asym_entity.get(tgt_chain, ""))
        }
    )
    skipped_entity_ids = sorted(
        {
            tgt_asym_entity[asym]
            for asym in tgt_asym
            if _public_entity_id(tgt_asym_entity.get(asym, ""))
            and tgt_asym_entity[asym] not in imported_entity_ids
        }
    )
    imported_asym_ids = sorted(eligible_pairs)
    skipped_asym_ids = sorted((ref_asym - set(eligible_pairs)) | (tgt_asym - set(eligible_pairs.values())))
    ok = bool(eligible_pairs)
    partial = ok and bool(skipped_asym_ids or failures)
    checked = sorted(ref_asym | tgt_asym)
    # Full success keeps the previous log shape (no imported/skipped lists).
    report_subset = partial or not ok
    return SafeguardResult(
        ok,
        mode,
        failures,
        checked,
        eligible_pairs if ok else {},
        content_aligned and ok,
        partial,
        imported_asym_ids if report_subset and ok else [],
        skipped_asym_ids if report_subset else [],
        imported_entity_ids if report_subset and ok else [],
        skipped_entity_ids if report_subset else [],
        sorted(eid for eid in eligible_ref_entities if _public_entity_id(eid)) if ok else [],
        ref_polymer_entity_ids,
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
