"""Unit tests for polymer_safeguards."""

import unittest
from typing import List

import gemmi

from polymer_safeguards import (
    atom_site_sequence_and_count,
    build_polymer_entity_remapping,
    merge_partial_macromolecule_rows,
    pair_polymer_chains_by_content,
    reconcile_polymer_struct_asym_in_block,
    remap_macromolecule_metadata_for_target,
    restrict_macromolecule_metadata_to_entities,
    strip_macromolecule_categories,
    validate_macromolecule_merge,
    _get_loop,
    _loop_as_table,
    _find_column,
    _cif_string_raw,
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

    def test_auth_seq_id_fallback_when_label_seq_id_missing(self):
        """PDB-style mmCIF: label_seq_id is '.' but auth_seq_id has residue numbers."""
        cif = """data_test
loop_
_atom_site.group_PDB
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.label_comp_id
ATOM A . 1 ALA
ATOM A . 2 GLY
ATOM B . 10 SER
ATOM B . 11 THR
"""
        ref = _read(cif)
        tgt = _read(cif)
        r = validate_macromolecule_merge(ref, tgt)
        self.assertTrue(r.ok)
        self.assertEqual(r.mode, "forced_fallback")
        self.assertEqual(r.checked_asym_ids, ["A", "B"])
        n, s = atom_site_sequence_and_count(ref, "A")
        self.assertEqual(n, 2)
        self.assertEqual(s, "AG")

    def test_no_vacuous_pass_when_no_polymer_chains(self):
        """Blank label_seq_id and no auth_seq_id must not vacuous-pass."""
        cif = """data_test
loop_
_atom_site.group_PDB
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.label_comp_id
ATOM A . ALA
ATOM A . GLY
ATOM B . SER
ATOM B . THR
"""
        r = validate_macromolecule_merge(_read(cif), _read(cif))
        self.assertFalse(r.ok)
        self.assertTrue(
            any(f.get("rule") == "ALIGN-0-NO-POLYMER-CHAINS" for f in r.failures)
        )

    def test_auth_seq_id_detects_mismatch_vs_session_style_empty_ref(self):
        """Target with auth_seq_id polymers vs empty ref -> ALIGN-1, not vacuous pass."""
        ref = _read(
            """data_test
_software.name foo
"""
        )
        tgt = _read(
            """data_test
loop_
_atom_site.group_PDB
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.label_comp_id
ATOM Axp . 1 ALA
ATOM Axp . 2 GLY
"""
        )
        r = validate_macromolecule_merge(ref, tgt)
        self.assertFalse(r.ok)
        self.assertTrue(
            any(f.get("rule") == "ALIGN-1-CONTENT-MISMATCH" for f in r.failures)
        )

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
        self.assertTrue(r.ok)
        self.assertTrue(r.partial)
        self.assertEqual(r.imported_asym_ids, ["A"])
        self.assertEqual(r.skipped_asym_ids, ["B"])
        self.assertFalse(any(f.get("rule") == "ALIGN-1-CONTENT-MISMATCH" for f in r.failures))

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

    def test_remap_entity_loop_keeps_polymer_description(self):
        ref_cif = """data_test
loop_
_entity.id
_entity.type
_entity.pdbx_description
1 polymer EAL-protein
2 non-polymer GLYCEROL
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
B 1
loop_
_atom_site.group_PDB
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.label_comp_id
_atom_site.label_entity_id
ATOM A 1 ALA A
ATOM B 1 SER B
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
ATOM A 1 ALA A
ATOM B 1 SER B
"""
        ref = _read(ref_cif)
        tgt = _read(tgt_cif)
        meta = _read("""data_test
loop_
_entity.id
_entity.type
_entity.pdbx_description
1 polymer EAL-protein
2 non-polymer GLYCEROL
""")
        out = remap_macromolecule_metadata_for_target(meta, ref, tgt, {"A": "A", "B": "B"})
        from polymer_safeguards import _get_loop, _loop_as_table, _find_column, _cif_string_raw

        loop = _get_loop(out, "_entity.id")
        self.assertIsNotNone(loop)
        tags, rows = _loop_as_table(loop)
        ie = _find_column(tags, "_entity.id")
        idesc = _find_column(tags, "_entity.pdbx_description")
        by_id = {_cif_string_raw(row[ie]): _cif_string_raw(row[idesc]) for row in rows}
        self.assertEqual(by_id.get("A"), "EAL-protein")
        self.assertEqual(by_id.get("B"), "EAL-protein")
        self.assertEqual(by_id.get("2"), "GLYCEROL")

    def test_remap_frame_entity_src_gen_one_to_two_target_entities(self):
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
_atom_site.group_PDB
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.label_comp_id
_atom_site.label_entity_id
ATOM A 1 ALA A
ATOM B 1 SER B
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
ATOM Bxp 1 SER B
"""
        ref = _read(ref_cif)
        tgt = _read(tgt_cif)
        meta = _read("""data_test
_entity_src_gen.entity_id 1
_entity_src_gen.pdbx_gene_src_gene PA0861
_entity_src_gen.pdbx_gene_src_scientific_name "Pseudomonas aeruginosa PAO1"
_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id 208964
_entity_src_gen.pdbx_host_org_scientific_name "Escherichia coli BL21(DE3)"
_entity_src_gen.pdbx_host_org_ncbi_taxonomy_id 469008
_entity_src_gen.plasmid_name pET28a
""")
        pairing = {"A": "Axp", "B": "Bxp"}
        out = remap_macromolecule_metadata_for_target(meta, ref, tgt, pairing)
        ids = _entity_src_gen_ids_from_block(out)
        self.assertEqual(ids, ["A", "B"])
        from polymer_safeguards import _get_loop, _loop_as_table, _find_column, _cif_string_raw

        loop = _get_loop(out, "_entity_src_gen.entity_id")
        self.assertIsNotNone(loop)
        tags, rows = _loop_as_table(loop)
        ig = _find_column(tags, "_entity_src_gen.pdbx_gene_src_gene")
        self.assertEqual(_cif_string_raw(rows[0][ig]), "PA0861")
        self.assertEqual(_cif_string_raw(rows[1][ig]), "PA0861")

    def test_remap_frame_entity_poly_to_target_entity_ids(self):
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
_entity_poly.entity_id 1
_entity_poly.type polypeptide(L)
_entity_poly.pdbx_strand_id A,B
_entity_poly.pdbx_seq_one_letter_code AG
""")
        pairing = {"A": "Axp", "B": "Bxp"}
        out = remap_macromolecule_metadata_for_target(meta, ref, tgt, pairing)
        self.assertEqual(_entity_poly_ids_from_block(out), ["A", "B"])

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

    def _chains_cif(self, chains, descriptions=None, assembly_list=None) -> str:
        """chains: (asym, entity_id, comp_ids). comp_ids is a list of 3-letter codes."""
        asym_rows = "\n".join(f"{asym} {eid}" for asym, eid, _comps in chains)
        entity_ids = []
        for _asym, eid, _comps in chains:
            if eid not in entity_ids:
                entity_ids.append(eid)
        desc = descriptions or {}
        ent_rows = "\n".join(
            f"{eid} polypeptide(L) '{desc.get(eid, eid)}'" for eid in entity_ids
        )
        atom_rows = []
        for asym, eid, comps in chains:
            for i, comp in enumerate(comps, start=1):
                atom_rows.append(f"ATOM {asym} {i} {comp} {eid}")
        assembly = ""
        if assembly_list is not None:
            joined = ", ".join(assembly_list)
            assembly = f"""
_em_entity_assembly.id 1
_em_entity_assembly.entity_id_list '{joined}'
_em_entity_assembly_naturalsource.entity_assembly_id 1
_em_entity_assembly_naturalsource.organism 'Escherichia coli'
"""
        return f"""data_test
loop_
_struct_asym.id
_struct_asym.entity_id
{asym_rows}
loop_
_entity.id
_entity.type
_entity.pdbx_description
{ent_rows}
loop_
_atom_site.group_PDB
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.label_comp_id
_atom_site.label_entity_id
{chr(10).join(atom_rows)}
{assembly}
"""

    def test_partial_import_skips_only_failed_entity(self):
        ref = _read(
            self._chains_cif(
                [("A", "1", ["ALA", "GLY"]), ("H", "2", ["GLY", "TRP"])],
                {"1": "ProteinA-ref", "2": "PeptideH-ref"},
            )
        )
        tgt = _read(
            self._chains_cif(
                [("A", "1", ["ALA", "GLY"]), ("H", "2", ["GLY", "TYR"])],
                {"1": "ProteinA-old", "2": "PeptideH-old"},
            )
        )
        result = validate_macromolecule_merge(ref, tgt)
        self.assertTrue(result.ok)
        self.assertTrue(result.partial)
        self.assertEqual(result.imported_asym_ids, ["A"])
        self.assertEqual(result.skipped_asym_ids, ["H"])
        self.assertEqual(result.imported_entity_ids, ["1"])
        self.assertEqual(result.skipped_entity_ids, ["2"])
        self.assertTrue(any(f.get("rule") == "ALIGN-2-LENGTH-OR-ATOM-SEQ" for f in result.failures))
        self.assertEqual(result.chain_pairing, {"A": "A"})

    def test_multichain_entity_skipped_when_one_chain_mismatches(self):
        ref = _read(
            self._chains_cif(
                [
                    ("A", "1", ["ALA", "GLY"]),
                    ("B", "1", ["SER", "THR"]),
                    ("C", "2", ["ALA", "ALA"]),
                ]
            )
        )
        tgt = _read(
            self._chains_cif(
                [
                    ("A", "1", ["ALA", "GLY"]),
                    ("B", "1", ["SER", "MET"]),
                    ("C", "2", ["ALA", "ALA"]),
                ]
            )
        )
        result = validate_macromolecule_merge(ref, tgt)
        self.assertTrue(result.ok)
        self.assertTrue(result.partial)
        self.assertEqual(result.imported_asym_ids, ["C"])
        self.assertEqual(result.skipped_asym_ids, ["A", "B"])
        self.assertEqual(result.imported_entity_ids, ["2"])
        self.assertEqual(result.skipped_entity_ids, ["1"])
        self.assertTrue(any(f.get("rule") == "ALIGN-ENTITY-INCOMPLETE" for f in result.failures))

    def test_no_eligible_entity_is_not_partial(self):
        ref = _read(self._chains_cif([("A", "1", ["ALA", "GLY"]), ("B", "2", ["SER", "THR"])]))
        tgt = _read(self._chains_cif([("A", "1", ["ALA", "TRP"]), ("B", "2", ["SER", "MET"])]))
        result = validate_macromolecule_merge(ref, tgt)
        self.assertFalse(result.ok)
        self.assertFalse(result.partial)
        self.assertEqual(result.imported_asym_ids, [])
        self.assertEqual(result.chain_pairing, {})

    def test_shared_chain_id_mismatch_is_not_reassigned(self):
        """A same-ID sequence mismatch stays on that ID even if another chain has the sequence."""
        ref = _read(
            self._chains_cif(
                [("A", "1", ["ALA", "GLY"]), ("H", "2", ["GLY", "TRP"])],
                {"1": "keep-A", "2": "peptide-H"},
            )
        )
        tgt = _read(
            self._chains_cif(
                [
                    ("A", "1", ["ALA", "GLY"]),
                    ("H", "2", ["GLY", "TYR"]),
                    ("Z", "3", ["GLY", "TRP"]),
                ],
                {"1": "old-A", "2": "old-H", "3": "extra"},
            )
        )
        result = validate_macromolecule_merge(ref, tgt)
        self.assertTrue(result.ok)
        self.assertTrue(result.partial)
        self.assertEqual(result.chain_pairing, {"A": "A"})
        self.assertIn("H", result.skipped_asym_ids)
        self.assertIn("Z", result.skipped_asym_ids)
        self.assertTrue(any(f.get("label_asym_id") == "H" and f.get("rule") == "ALIGN-2-LENGTH-OR-ATOM-SEQ" for f in result.failures))
        self.assertNotIn("Z", result.chain_pairing.values())

    def test_partial_content_pairing_imports_unique_profile_only(self):
        ref = _read(
            self._chains_cif(
                [("A", "1", ["ALA", "GLY"]), ("B", "2", ["SER", "THR"])],
                {"1": "from-ref", "2": "skip-ref"},
            )
        )
        tgt = _read(
            self._chains_cif(
                [("Axp", "9", ["ALA", "GLY"]), ("Bxp", "8", ["TRP", "TRP"])],
                {"9": "old-A", "8": "old-B"},
            )
        )
        result = validate_macromolecule_merge(ref, tgt)
        self.assertTrue(result.ok)
        self.assertTrue(result.partial)
        self.assertTrue(result.content_aligned)
        self.assertEqual(result.chain_pairing, {"A": "Axp"})
        self.assertEqual(result.imported_entity_ids, ["9"])
        self.assertIn("B", result.skipped_asym_ids)
        self.assertIn("Bxp", result.skipped_asym_ids)

    def test_ambiguous_profiles_do_not_block_a_unique_chain(self):
        ref = _read(
            self._chains_cif(
                [
                    ("A", "1", ["ALA", "GLY"]),
                    ("B", "2", ["ALA", "GLY"]),
                    ("C", "3", ["SER", "THR"]),
                ]
            )
        )
        tgt = _read(
            self._chains_cif(
                [
                    ("X", "1", ["ALA", "GLY"]),
                    ("Y", "2", ["ALA", "GLY"]),
                    ("Z", "3", ["SER", "THR"]),
                ]
            )
        )
        result = validate_macromolecule_merge(ref, tgt)
        self.assertTrue(result.ok)
        self.assertTrue(result.partial)
        self.assertEqual(result.chain_pairing, {"C": "Z"})
        self.assertTrue(any(f.get("rule") == "ALIGN-1-AMBIGUOUS-PROFILE" for f in result.failures))

    def test_em_assembly_omitted_when_it_names_a_skipped_entity(self):
        ref = _read(
            self._chains_cif(
                [("A", "1", ["ALA", "GLY"]), ("H", "2", ["GLY", "TRP"])],
                assembly_list=["1", "2"],
            )
        )
        meta = restrict_macromolecule_metadata_to_entities(
            ref, ["1"], ["A"], ["1", "2"]
        )
        cats = set()
        for item in meta:
            if item.pair is not None:
                cats.add(item.pair[0].split(".", 1)[0])
            elif item.loop is not None and item.loop.tags:
                cats.add(item.loop.tags[0].split(".", 1)[0])
        self.assertIn("_entity", cats)
        self.assertNotIn("_em_entity_assembly", cats)
        self.assertNotIn("_em_entity_assembly_naturalsource", cats)
        descriptions = _descriptions_by_entity(meta)
        self.assertIn("1", descriptions)
        self.assertNotIn("2", descriptions)

    def test_row_merge_keeps_ligand_and_failed_entity(self):
        target = _read(
            """data_test
loop_
_entity.id
_entity.type
_entity.pdbx_description
1 polymer old-A
2 polymer old-H
3 non-polymer GLYCEROL
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
H 2
L 3
"""
        )
        incoming = _read(
            """data_test
loop_
_entity.id
_entity.type
_entity.pdbx_description
1 polymer new-A
loop_
_entity_src_gen.entity_id
_entity_src_gen.pdbx_gene_src_scientific_name
1 'Homo sapiens'
"""
        )
        merged = merge_partial_macromolecule_rows(target, incoming, ["1"])
        descriptions = _descriptions_by_entity(merged)
        self.assertEqual(descriptions["1"], "new-A")
        self.assertEqual(descriptions["2"], "old-H")
        self.assertEqual(descriptions["3"], "GLYCEROL")
        self.assertEqual(_struct_asym_map(merged)["H"], "2")
        src = _entity_src_gen_ids_from_block(merged)
        self.assertEqual(src, ["1"])


def _descriptions_by_entity(block) -> dict:
    loop = _get_loop(block, "_entity.id")
    self_tags, rows = _loop_as_table(loop)
    ie = _find_column(self_tags, "_entity.id")
    idesc = _find_column(self_tags, "_entity.pdbx_description")
    return {
        _cif_string_raw(row[ie]): _cif_string_raw(row[idesc]) for row in rows
    }


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


def _entity_src_gen_ids_from_block(block) -> List[str]:
    from polymer_safeguards import _get_loop, _loop_as_table, _find_column, _cif_string_raw

    loop = _get_loop(block, "_entity_src_gen.entity_id")
    if loop:
        tags, rows = _loop_as_table(loop)
        ie = _find_column(tags, "_entity_src_gen.entity_id")
        if ie is None:
            return []
        return [_cif_string_raw(row[ie]) for row in rows]
    for item in block:
        if item.pair is not None and item.pair[0] == "_entity_src_gen.entity_id":
            return [_cif_string_raw(item.pair[1])]
    return []


if __name__ == "__main__":
    unittest.main()
