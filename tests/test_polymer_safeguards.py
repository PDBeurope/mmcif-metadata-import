"""Unit tests for polymer_safeguards."""

import unittest

import gemmi

from polymer_safeguards import (
    atom_site_sequence_and_count,
    pair_polymer_chains_by_content,
    strip_macromolecule_categories,
    validate_macromolecule_merge,
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


if __name__ == "__main__":
    unittest.main()
