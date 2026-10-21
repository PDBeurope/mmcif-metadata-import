"""Tests for specification CSV path resolution (repository and pip install layouts)."""

import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import import_metadata as im


class SpecPathResolutionTest(unittest.TestCase):
    def test_get_spec_dir_finds_repo_specs(self):
        spec_dir = im.get_spec_dir()
        self.assertIsNotNone(spec_dir)
        self.assertTrue((spec_dir / "EM.csv").exists())

    def test_resolve_spec_path_from_repo(self):
        resolved = im.resolve_spec_path("specs/EM.csv")
        self.assertTrue(resolved.exists())
        self.assertEqual(resolved.name, "EM.csv")

    def test_resolve_spec_path_prefers_cwd_override(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            custom = tmp_path / "specs"
            custom.mkdir()
            marker = custom / "EM.csv"
            marker.write_text("custom\n", encoding="utf-8")
            previous = os.getcwd()
            try:
                os.chdir(tmp_path)
                # Resolve before restoring cwd. A relative path would otherwise
                # point at the repository specs/ once the test returns there.
                resolved = im.resolve_spec_path("specs/EM.csv").resolve()
            finally:
                os.chdir(previous)
            self.assertEqual(resolved, marker.resolve())

    def test_resolve_spec_path_from_pip_install_prefix(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            site_packages = tmp_path / "site-packages"
            site_packages.mkdir()
            prefix_specs = tmp_path / "specs"
            prefix_specs.mkdir()
            em_csv = prefix_specs / "EM.csv"
            em_csv.write_text("category,item,should_import,type\n", encoding="utf-8")
            empty_cwd = tmp_path / "work"
            empty_cwd.mkdir()

            previous = os.getcwd()
            try:
                os.chdir(empty_cwd)
                with patch.object(im, "PROJECT_ROOT", site_packages), patch.object(
                    sys, "prefix", str(tmp_path)
                ):
                    self.assertEqual(im.get_spec_dir(), prefix_specs)
                    resolved = im.resolve_spec_path("specs/EM.csv")
            finally:
                os.chdir(previous)
            self.assertEqual(resolved, em_csv)
            self.assertTrue(resolved.exists())


if __name__ == "__main__":
    unittest.main()
