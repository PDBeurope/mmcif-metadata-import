"""Unit tests for polymer_safeguards."""

import unittest
from typing import List

import gemmi

from polymer_safeguards import (
    atom_site_sequence_and_count,
    build_polymer_entity_remapping,
    pair_polymer_chains_by_content,
    reconcile_polymer_struct_asym_in_block,
    remap_macromolecule_metadata_for_target,
    strip_macromolecule_categories,
    validate_macromolecule_merge,
    _struct_asym_map,
)


def _read(s: str) -> gemmi.cif.Block:
    doc = gemmi.cif.read_string(s)
    return doc[0]


class PolymerSafeguardsTest(unittest.TestCase):
    def _minimal_poly_cif(self, seq_one_letter: str, atom_rows: str) -> str:
        """seq_one_letter for entity_poly; atom_rows lines like: ATOM A 1 ALA"""
        return f"""data_test
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
#
loop_
_entity.id
_entity.type
1 polypeptide(L)
#
loop_
_entity_poly.entity_id
_entity_poly.pdbx_seq_one_letter_code
1 {seq_one_letter}
#
loop_
_atom_site.group_PDB
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.label_comp_id
{atom_rows}
"""

    def test_validate_identical_ok(self):
        cif = self._minimal_poly_cif(
            "AG",
            "ATOM A 1 ALA\nATOM A 2 GLY\n",
        )
        ref = _read(cif)
        tgt = _read(cif)
        r = validate_macromolecule_merge(ref, tgt)
        self.assertTrue(r.ok)
        self.assertEqual(r.mode, "entity")

    def test_validate_sequence_mismatch(self):
        ref = _read(self._minimal_poly_cif("AG", "ATOM A 1 ALA\nATOM A 2 GLY\n"))
        tgt = _read(self._minimal_poly_cif("AA", "ATOM A 1 ALA\nATOM A 2 ALA\n"))
        r = validate_macromolecule_merge(ref, tgt)
        self.assertFalse(r.ok)
        self.assertTrue(any(f.get("rule") == "ALIGN-2-LENGTH-OR-ATOM-SEQ" for f in r.failures))

    def test_mse_encoding_length(self):
        cif = self._minimal_poly_cif(
            "A(MSE)G",
            "ATOM A 1 ALA\nATOM A 2 MSE\nHETATM A 2 MSE\nATOM A 3 GLY\n",
        )
        b = _read(cif)
        n, s = atom_site_sequence_and_count(b, "A")
        self.assertEqual(n, 3)
        self.assertEqual(s, "A(MSE)G")

    def test_strip_macromolecule_categories(self):
        cif = """data_test
_entity.id 1
_entity.type polypeptide(L)
_software.name foo
"""
        b = _read(cif)
        out = strip_macromolecule_categories(b)
        cats = set()
        for it in out:
            if it.pair:
                cats.add(it.pair[0].split(".")[0])
            elif it.loop and it.loop.tags:
                cats.add(it.loop.tags[0].split(".")[0])
        self.assertNotIn("_entity", cats)
        self.assertIn("_software", cats)

    def test_forced_fallback_both(self):
        """No _struct_asym -> forced_fallback; same atom_site on both."""
        cif = """data_test
loop_
_atom_site.group_PDB
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.label_comp_id
ATOM A 1 ALA
ATOM A 2 GLY
ATOM B 1 SER
ATOM B 2 THR
"""
        ref = _read(cif)
        tgt = _read(cif)
        r = validate_macromolecule_merge(ref, tgt)
        self.assertTrue(r.ok)
        self.assertEqual(r.mode, "forced_fallback")
        self.assertFalse(r.content_aligned)

    def test_content_aligned_renamed_chains(self):
        """Same polymer content, different label_asym_id (refmac-style)."""
        ref_cif = """data_test
loop_
_atom_site.group_PDB
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.label_comp_id
ATOM A 1 ALA
ATOM A 2 GLY
ATOM B 1 SER
ATOM B 2 THR
"""
        tgt_cif = """data_test
loop_
_atom_site.group_PDB
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.label_comp_id
ATOM Axp 1 ALA
ATOM Axp 2 GLY
ATOM Bxp 1 SER
ATOM Bxp 2 THR
"""
        r = validate_macromolecule_merge(_read(ref_cif), _read(tgt_cif))
        self.assertTrue(r.ok)
        self.assertTrue(r.content_aligned)
        self.assertEqual(r.chain_pairing, {"A": "Axp", "B": "Bxp"})

    def test_content_mismatch_extra_chain(self):
        ref = _read(
            """data_test
loop_
_atom_site.group_PDB
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.label_comp_id
ATOM A 1 ALA
ATOM A 2 GLY
ATOM B 1 SER
ATOM B 2 THR
"""
        )
        tgt = _read(
            """data_test
loop_
_atom_site.group_PDB
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.label_comp_id
ATOM A 1 ALA
ATOM A 2 GLY
"""
        )
        r = validate_macromolecule_merge(ref, tgt)
        self.assertFalse(r.ok)
        self.assertTrue(any(f.get("rule") == "ALIGN-1-CONTENT-MISMATCH" for f in r.failures))

    def test_content_mismatch_different_sequence(self):
        ref = _read(
            self._minimal_poly_cif("AG", "ATOM A 1 ALA\nATOM A 2 GLY\n")
        )
        tgt = _read(
            """data_test
loop_
_atom_site.group_PDB
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.label_comp_id
ATOM X 1 ALA
ATOM X 2 ALA
"""
        )
        r = validate_macromolecule_merge(ref, tgt)
        self.assertFalse(r.ok)
        self.assertTrue(any(f.get("rule") == "ALIGN-1-CONTENT-MISMATCH" for f in r.failures))

    def test_pair_polymer_chains_by_content(self):
        ref = {"A": (2, "AG"), "B": (2, "ST")}
        tgt = {"Axp": (2, "AG"), "Bxp": (2, "ST")}
        self.assertEqual(pair_polymer_chains_by_content(ref, tgt), {"A": "Axp", "B": "Bxp"})
        self.assertIsNone(pair_polymer_chains_by_content(ref, {"Axp": (2, "AA")}))

    def test_remap_entity_poly_to_target_entity_ids(self):
        ref_cif = """data_test
loop_
_entity.id
_entity.type
1 polymer
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
B 1
loop_
_entity_poly.entity_id
_entity_poly.pdbx_seq_one_letter_code
1 AG
loop_
_atom_site.group_PDB
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.label_comp_id
_atom_site.label_entity_id
ATOM A 1 ALA A
ATOM A 2 GLY A
ATOM B 1 SER B
ATOM B 2 THR B
"""
        tgt_cif = """data_test
loop_
_entity.id
_entity.type
A polymer
B polymer
loop_
_atom_site.group_PDB
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.label_comp_id
_atom_site.label_entity_id
ATOM Axp 1 ALA A
ATOM Axp 2 GLY A
ATOM Bxp 1 SER B
ATOM Bxp 2 THR B
"""
        ref = _read(ref_cif)
        tgt = _read(tgt_cif)
        meta = _read("""data_test
loop_
_entity_poly.entity_id
_entity_poly.pdbx_seq_one_letter_code
1 AG
""")
        pairing = {"A": "Axp", "B": "Bxp"}
        out = remap_macromolecule_metadata_for_target(meta, ref, tgt, pairing)
        self.assertEqual(_entity_poly_ids_from_block(out), ["A", "B"])

    def test_reconcile_struct_asym_from_atom_site(self):
        cif = """data_test
loop_
_entity.id
_entity.type
A polymer
B polymer
loop_
_atom_site.group_PDB
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.label_comp_id
_atom_site.label_entity_id
ATOM Axp 1 ALA A
ATOM Axp 2 GLY A
ATOM Bxp 1 SER B
ATOM Bxp 2 THR B
"""
        b = _read(cif)
        new_b, changed = reconcile_polymer_struct_asym_in_block(b)
        self.assertTrue(changed)
        sm = _struct_asym_map(new_b)
        self.assertEqual(sm.get("Axp"), "A")
        self.assertEqual(sm.get("Bxp"), "B")


def _entity_poly_ids_from_block(block) -> List[str]:
    from polymer_safeguards import _get_loop, _loop_as_table, _find_column, _cif_string_raw

    loop = _get_loop(block, "_entity_poly.entity_id")
    if not loop:
        return []
    tags, rows = _loop_as_table(loop)
    ie = _find_column(tags, "_entity_poly.entity_id")
    if ie is None:
        return []
    return sorted({_cif_string_raw(row[ie]) for row in rows})


if __name__ == "__main__":
    unittest.main()
