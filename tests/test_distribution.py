"""Privacy and integrity checks, using temporary dummy files only."""
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("verify_package", ROOT / "verify_package.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class DistributionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        (self.root / "README.md").write_text("Code package\n")
        raw = (self.root / "README.md").read_bytes()
        (self.root / "RELEASE_MANIFEST.json").write_text(json.dumps({"files": [
            {"path": "README.md", "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}]}))

    def tearDown(self):
        self.temp.cleanup()

    def test_unchanged_code(self):
        mod.verify(self.root)

    def test_extra_observations_rejected(self):
        (self.root / "Data").mkdir()
        (self.root / "Data" / "full_RTG_data.csv").write_text("dummy\n")
        with self.assertRaises(ValueError):
            mod.verify(self.root)

    def test_modified_source_rejected(self):
        (self.root / "README.md").write_text("Fake package\n")
        with self.assertRaises(ValueError):
            mod.verify(self.root)

    def test_link_to_external_file_rejected(self):
        external = tempfile.TemporaryDirectory()
        self.addCleanup(external.cleanup)
        target = Path(external.name) / "external.txt"
        target.write_text("Code package\n")
        (self.root / "README.md").unlink()
        (self.root / "README.md").symlink_to(target)
        with self.assertRaises(ValueError):
            mod.verify(self.root)

    def test_workspace_inside_release_rejected_before_data_copy(self):
        destination = ROOT / "SHOULD_NOT_BE_CREATED"
        result = subprocess.run([sys.executable, str(ROOT / "run_analysis.py"),
            "--data-dir", str(self.root), "--work-dir", str(destination)], capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(b"private work directory must be outside the release directory", result.stderr)
        self.assertFalse(destination.exists())

    def test_small_bootstrap_cannot_render(self):
        destination = self.root / "private-run"
        result = subprocess.run([sys.executable, str(ROOT / "run_analysis.py"),
            "--data-dir", str(self.root), "--work-dir", str(destination),
            "--bootstrap", "80"], capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(b"Small bootstrap checks cannot render", result.stderr)
        self.assertFalse(destination.exists())

    def test_normal_clone_metadata_allowed(self):
        (self.root / ".git").mkdir()
        (self.root / ".git" / "HEAD").write_text("ref: refs/heads/main\n")
        mod.verify(self.root)

    def test_linked_git_metadata_rejected(self):
        target = self.root / "outside-git"
        target.mkdir()
        (self.root / ".git").symlink_to(target, target_is_directory=True)
        with self.assertRaises(ValueError):
            mod.verify(self.root)

    def test_private_workspace_does_not_copy_git_history(self):
        clone = self.root / "clone"
        shutil.copytree(ROOT, clone, ignore=shutil.ignore_patterns(".git", "__pycache__"))
        (clone / ".git").mkdir()
        (clone / ".git" / "HEAD").write_text("ref: refs/heads/main\n")
        inputs = self.root / "inputs"
        inputs.mkdir()
        for name in ("full_RTG_data.csv", "demographics.csv"):
            (inputs / name).write_text("invented\n")
        archive = self.root / "archive"
        (archive / "results" / "HMM").mkdir(parents=True)
        (archive / "modsData").mkdir()
        work = self.root / "private-work"
        result = subprocess.run([sys.executable, str(clone / "run_analysis.py"),
            "--data-dir", str(inputs), "--work-dir", str(work), "--mode", "archived",
            "--private-archive", str(archive), "--no-render"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((work / "Data" / "full_RTG_data.csv").is_file())
        self.assertFalse((work / ".git").exists())


if __name__ == "__main__":
    unittest.main()
