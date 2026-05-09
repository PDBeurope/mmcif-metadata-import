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

# Must match specs/MACROMOLECULES.csv (category column, should_import Y).
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


@dataclass
class SafeguardResult:
    ok: bool
    mode: str
    failures: List[Dict[str, Any]] = field(default_factory=list)
    checked_asym_ids: List[str] = field(default_factory=list)

    def to_json_fragment(self) -> str:
        return json.dumps(
            {"ok": self.ok, "mode": self.mode, "failures": self.failures, "checked_asym_ids": self.checked_asym_ids},
            indent=2,
        )


def validate_macromolecule_merge(reference_block: gemmi.cif.Block, target_block: gemmi.cif.Block) -> SafeguardResult:
    """
    R-ALIGN / R-POLY checks between reference and target first data blocks.
    """
    mode = identification_mode(reference_block, target_block)
    ref_asym = polymer_asym_ids_for_mode(reference_block, mode)
    tgt_asym = polymer_asym_ids_for_mode(target_block, mode)

    failures: List[Dict[str, Any]] = []

    if ref_asym != tgt_asym:
        failures.append(
            {
                "rule": "ALIGN-1-ASYMM-SET",
                "mode": mode,
                "reference_asym_ids": sorted(ref_asym),
                "target_asym_ids": sorted(tgt_asym),
            }
        )
        return SafeguardResult(False, mode, failures, sorted(ref_asym | tgt_asym))

    asym_ids = sorted(ref_asym)
    if not asym_ids:
        return SafeguardResult(True, mode, [], [])

    sm_r = _struct_asym_map(reference_block)
    sm_t = _struct_asym_map(target_block)

    for asym in asym_ids:
        cr, sr = atom_site_sequence_and_count(reference_block, asym)
        ct, st = atom_site_sequence_and_count(target_block, asym)
        if cr != ct or sr != st:
            failures.append(
                {
                    "rule": "ALIGN-2-LENGTH-OR-ATOM-SEQ",
                    "label_asym_id": asym,
                    "reference_count": cr,
                    "target_count": ct,
                    "reference_atom_site_seq": sr,
                    "target_atom_site_seq": st,
                }
            )
            continue

        if mode == "entity":
            er = sm_r.get(asym)
            et = sm_t.get(asym)
            if er and et:
                pr = _entity_poly_one_letter(reference_block, er)
                pt = _entity_poly_one_letter(target_block, et)
                if pr is not None and sr != pr:
                    failures.append(
                        {
                            "rule": "ALIGN-3-REF-ATOM-VS-ENTITY-POLY",
                            "label_asym_id": asym,
                            "entity_id": er,
                            "atom_site_seq": sr,
                            "entity_poly_seq": pr,
                        }
                    )
                if pt is not None and st != pt:
                    failures.append(
                        {
                            "rule": "ALIGN-3-TGT-ATOM-VS-ENTITY-POLY",
                            "label_asym_id": asym,
                            "entity_id": et,
                            "atom_site_seq": st,
                            "entity_poly_seq": pt,
                        }
                    )
                if pr is not None and pt is not None and pr != pt:
                    failures.append(
                        {
                            "rule": "ALIGN-3-ENTITY-POLY-MISMATCH",
                            "label_asym_id": asym,
                            "reference_entity_poly": pr,
                            "target_entity_poly": pt,
                        }
                    )

    return SafeguardResult(len(failures) == 0, mode, failures, asym_ids)


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
