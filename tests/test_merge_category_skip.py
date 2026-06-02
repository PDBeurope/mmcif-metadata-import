"""Tests for default-merge category-level skip (DAOTHER-10789)."""

import tempfile
import unittest
from pathlib import Path

import gemmi

import import_metadata as im


class MergeCategorySkipTest(unittest.TestCase):
    def test_skips_incoming_pairs_when_category_exists_as_frame(self):
        target = """data_target
_em_ctf_correction.details ?
_em_ctf_correction.id 1
#
loop_
_atom_site.group_PDB
_atom_site.id
ATOM 1
"""
        incoming = """data_in
_em_ctf_correction.amplitude_correction 0.5
_em_ctf_correction.id 2
"""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            target_file = tmp_path / "target.cif"
            target_file.write_text(target, encoding="utf-8")
            meta_block = gemmi.cif.read_string(incoming)[0]
            out_file = tmp_path / "merged.cif"
            result = im.merge_metadata_to_file(
                meta_block, str(target_file), str(out_file)
            )
            self.assertTrue(result.success)
            merged = out_file.read_text(encoding="utf-8")
            self.assertNotIn("amplitude_correction", merged)
            self.assertIn("_em_ctf_correction.details", merged)

    def test_skips_incoming_loop_when_category_exists_as_pairs(self):
        target = """data_target
_em_ctf_correction.details ?
_em_ctf_correction.id 1
#
"""
        incoming = """data_in
loop_
_em_ctf_correction.id
_em_ctf_correction.amplitude_correction
2 0.5
"""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            target_file = tmp_path / "target.cif"
            target_file.write_text(target, encoding="utf-8")
            meta_block = gemmi.cif.read_string(incoming)[0]
            out_file = tmp_path / "merged.cif"
            result = im.merge_metadata_to_file(
                meta_block, str(target_file), str(out_file)
            )
            self.assertTrue(result.success)
            merged = out_file.read_text(encoding="utf-8")
            self.assertNotIn("amplitude_correction", merged)
            self.assertEqual(merged.count("loop_"), 0)


if __name__ == "__main__":
    unittest.main()
