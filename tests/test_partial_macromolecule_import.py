"""Merge-mode tests for partial macromolecule import."""

import tempfile
import unittest
from pathlib import Path

import gemmi

from import_metadata import import_metadata


SPEC = Path(__file__).resolve().parents[1] / "specs" / "MACROMOLECULES.csv"
EM_SPEC = SPEC.parent / "MACROMOLECULES_EM.csv"


def _write(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def _descriptions(path: Path) -> dict:
    block = gemmi.cif.read(str(path))[0]
    loop = block.find_loop("_entity.id")
    lp = loop if isinstance(loop, gemmi.cif.Loop) else loop.get_loop()
    tags = list(lp.tags)
    width = lp.width()
    values = list(lp.values)
    ie = tags.index("_entity.id")
    idesc = tags.index("_entity.pdbx_description")

    def raw(value: str) -> str:
        text = str(value).strip()
        if len(text) >= 2 and text[0] == text[-1] and text[0] in "'\"":
            return text[1:-1]
        return text

    out = {}
    for i in range(0, len(values), width):
        row = [str(values[i + j]) for j in range(width)]
        out[raw(row[ie])] = raw(row[idesc])
    return out


class PartialMacromoleculeImportTest(unittest.TestCase):
    def test_matching_entity_imported_when_another_polymer_differs(self):
        reference = """data_ref
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
H 2
L 3
loop_
_entity.id
_entity.type
_entity.pdbx_description
1 polymer ProteinA-from-ref
2 polymer PeptideH-from-ref
3 non-polymer LIG-from-ref
loop_
_entity_src_gen.entity_id
_entity_src_gen.pdbx_gene_src_scientific_name
1 'Homo sapiens'
2 'Homo sapiens'
_em_entity_assembly.id 1
_em_entity_assembly.entity_id_list '1, 2'
_em_entity_assembly_naturalsource.entity_assembly_id 1
_em_entity_assembly_naturalsource.organism 'Escherichia coli'
loop_
_atom_site.group_PDB
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.label_comp_id
_atom_site.label_entity_id
ATOM A 1 ALA 1
ATOM A 2 GLY 1
ATOM H 1 GLY 2
ATOM H 2 TRP 2
HETATM L 1 LIG 3
"""
        target = """data_tgt
loop_
_struct_asym.id
_struct_asym.entity_id
A 10
H 20
L 30
loop_
_entity.id
_entity.type
_entity.pdbx_description
10 polymer ProteinA-old
20 polymer PeptideH-old
30 non-polymer GLYCEROL
_em_entity_assembly.id 7
_em_entity_assembly.entity_id_list 20
_em_entity_assembly.name KEEP-ME
loop_
_atom_site.group_PDB
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.label_comp_id
_atom_site.label_entity_id
ATOM A 1 ALA 10
ATOM A 2 GLY 10
ATOM H 1 GLY 20
ATOM H 2 TYR 20
HETATM L 1 LIG 30
"""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ref_path = root / "ref.cif"
            tgt_path = root / "tgt.cif"
            out_path = root / "merged.cif"
            log_path = root / "merged.log"
            _write(ref_path, reference)
            _write(tgt_path, target)
            outcome = import_metadata(
                str(ref_path),
                [str(SPEC), str(EM_SPEC)],
                str(root / "unused.cif"),
                merge_to_file=str(tgt_path),
                merge_output_file=str(out_path),
                log_file=str(log_path),
            )
            self.assertTrue(outcome.ok)
            self.assertEqual(outcome.exit_status, 0)
            self.assertIsNotNone(outcome.safeguard_result)
            self.assertTrue(outcome.safeguard_result.partial)
            self.assertEqual(outcome.safeguard_result.imported_entity_ids, ["10"])
            self.assertEqual(outcome.safeguard_result.skipped_entity_ids, ["20"])
            self.assertIn("_em_entity_assembly", outcome.safeguard_result.omitted_categories)
            descriptions = _descriptions(out_path)
            self.assertEqual(descriptions["10"], "ProteinA-from-ref")
            self.assertEqual(descriptions["20"], "PeptideH-old")
            self.assertEqual(descriptions["30"], "GLYCEROL")
            self.assertNotIn("1", descriptions)
            self.assertNotIn("2", descriptions)
            block = gemmi.cif.read(str(out_path))[0]
            src = block.find_values("_entity_src_gen.entity_id")
            src_ids = [str(v).strip().strip("'\"") for v in src]
            self.assertEqual(src_ids, ["10"])
            assembly_names = [
                str(v).strip().strip("'\"")
                for v in block.find_values("_em_entity_assembly.name")
            ]
            self.assertEqual(assembly_names, ["KEEP-ME"])
            log_text = log_path.read_text(encoding="utf-8")
            self.assertIn("left unchanged", log_text)
            self.assertIn('"partial": true', log_text)

    def test_no_matching_polymer_exits_2(self):
        reference = """data_ref
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
loop_
_entity.id
_entity.type
_entity.pdbx_description
1 polymer from-ref
loop_
_atom_site.group_PDB
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.label_comp_id
ATOM A 1 ALA
ATOM A 2 GLY
"""
        target = """data_tgt
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
loop_
_entity.id
_entity.type
_entity.pdbx_description
1 polymer keep-target
loop_
_atom_site.group_PDB
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.label_comp_id
ATOM A 1 ALA
ATOM A 2 TRP
"""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ref_path = root / "ref.cif"
            tgt_path = root / "tgt.cif"
            out_path = root / "merged.cif"
            _write(ref_path, reference)
            _write(tgt_path, target)
            outcome = import_metadata(
                str(ref_path),
                [str(SPEC)],
                str(root / "unused.cif"),
                merge_to_file=str(tgt_path),
                merge_output_file=str(out_path),
            )
            self.assertTrue(outcome.ok)
            self.assertEqual(outcome.exit_status, 2)
            self.assertFalse(outcome.safeguard_result.partial)
            self.assertEqual(_descriptions(out_path)["1"], "keep-target")


if __name__ == "__main__":
    unittest.main()
