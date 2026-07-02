"""Tests for EM reconstruction-method gating in merge mode (DAOTHER-10724)."""

import tempfile
import unittest
from pathlib import Path
from typing import Optional

import import_metadata as im


def _em_reference_cif(reconstruction_method: str) -> str:
    return f"""data_ref
_exptl.method 'ELECTRON MICROSCOPY'
_em_experiment.reconstruction_method '{reconstruction_method}'
#
loop_
_database_2.database_id
WWPDB
EMDB
PDB
"""


def _em_target_cif(reconstruction_method: Optional[str]) -> str:
    recon_line = ""
    if reconstruction_method is not None:
        recon_line = f"_em_experiment.reconstruction_method '{reconstruction_method}'\n"
    return f"""data_tgt
_exptl.method 'ELECTRON MICROSCOPY'
{recon_line}#
loop_
_database_2.database_id
WWPDB
EMDB
PDB
"""


class EmSubmethodGateTest(unittest.TestCase):
    def test_match_includes_em_csv(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ref = tmp_path / "ref.cif"
            tgt = tmp_path / "tgt.cif"
            ref.write_text(_em_reference_cif("TOMOGRAPHY"), encoding="utf-8")
            tgt.write_text(_em_target_cif("TOMOGRAPHY"), encoding="utf-8")
            skipped = []
            out = im.gate_em_csv_for_merge_submethod(
                ["specs/EM.csv"],
                input_file=str(ref),
                merge_to_file=str(tgt),
                skipped_specs=skipped,
            )
            self.assertEqual(out, ["specs/EM.csv"])
            self.assertEqual(skipped, [])

    def test_mismatch_skips_em_csv(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ref = tmp_path / "ref.cif"
            tgt = tmp_path / "tgt.cif"
            ref.write_text(_em_reference_cif("SINGLE PARTICLE"), encoding="utf-8")
            tgt.write_text(_em_target_cif("TOMOGRAPHY"), encoding="utf-8")
            skipped = []
            out = im.gate_em_csv_for_merge_submethod(
                ["specs/EM.csv"],
                input_file=str(ref),
                merge_to_file=str(tgt),
                skipped_specs=skipped,
            )
            self.assertEqual(out, [])
            self.assertEqual(len(skipped), 1)
            self.assertIn("do not match", skipped[0][1])

    def test_target_missing_uses_hint_match(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ref = tmp_path / "ref.cif"
            tgt = tmp_path / "tgt.cif"
            ref.write_text(_em_reference_cif("SINGLE PARTICLE"), encoding="utf-8")
            tgt.write_text(_em_target_cif(None), encoding="utf-8")
            skipped = []
            out = im.gate_em_csv_for_merge_submethod(
                ["specs/EM.csv"],
                input_file=str(ref),
                merge_to_file=str(tgt),
                merge_target_em_reconstruction_method="SINGLE PARTICLE",
                skipped_specs=skipped,
            )
            self.assertEqual(out, ["specs/EM.csv"])
            self.assertEqual(skipped, [])

    def test_hint_mismatch_skips(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ref = tmp_path / "ref.cif"
            tgt = tmp_path / "tgt.cif"
            ref.write_text(_em_reference_cif("TOMOGRAPHY"), encoding="utf-8")
            tgt.write_text(_em_target_cif(None), encoding="utf-8")
            skipped = []
            out = im.gate_em_csv_for_merge_submethod(
                ["specs/EM.csv"],
                input_file=str(ref),
                merge_to_file=str(tgt),
                merge_target_em_reconstruction_method="SINGLE PARTICLE",
                skipped_specs=skipped,
            )
            self.assertEqual(out, [])
            self.assertIn("do not match", skipped[0][1])

    def test_unresolved_both_sides_skips(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ref = tmp_path / "ref.cif"
            tgt = tmp_path / "tgt.cif"
            ref.write_text(_em_target_cif(None), encoding="utf-8")
            tgt.write_text(_em_target_cif(None), encoding="utf-8")
            skipped = []
            out = im.gate_em_csv_for_merge_submethod(
                ["specs/EM.csv"],
                input_file=str(ref),
                merge_to_file=str(tgt),
                skipped_specs=skipped,
            )
            self.assertEqual(out, [])
            self.assertIn("could not be resolved", skipped[0][1])

    def test_bypass_flag_keeps_em_csv(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ref = tmp_path / "ref.cif"
            tgt = tmp_path / "tgt.cif"
            ref.write_text(_em_target_cif(None), encoding="utf-8")
            tgt.write_text(_em_target_cif(None), encoding="utf-8")
            skipped = []
            out = im.gate_em_csv_for_merge_submethod(
                ["specs/EM.csv"],
                input_file=str(ref),
                merge_to_file=str(tgt),
                ignore_em_submethod_mismatch=True,
                skipped_specs=skipped,
            )
            self.assertEqual(out, ["specs/EM.csv"])
            self.assertEqual(skipped, [])

    def test_no_merge_to_file_no_gate(self):
        skipped = []
        out = im.gate_em_csv_for_merge_submethod(
            ["specs/EM.csv"],
            input_file="ref.cif",
            merge_to_file=None,
            skipped_specs=skipped,
        )
        self.assertEqual(out, ["specs/EM.csv"])


if __name__ == "__main__":
    unittest.main()
