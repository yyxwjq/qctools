import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
VASPJOB = PROJECT_ROOT / "qctools" / "vaspjob"


class VaspJobTests(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.workdir = Path(self.tmpdir.name)
        self.script = self.workdir / "vaspjob"
        shutil.copy2(VASPJOB, self.script)
        self.script.chmod(0o755)

    def tearDown(self):
        self.tmpdir.cleanup()

    def run_vaspjob(self, *args):
        return subprocess.run(
            [str(self.script), *args],
            cwd=self.workdir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

    def write_task(self, name, incar, oszicar=None, outcar=None):
        task = self.workdir / name
        task.mkdir(parents=True)
        (task / "INCAR").write_text(incar, encoding="utf-8")
        (task / "POSCAR").write_text("POSCAR\n", encoding="utf-8")
        if oszicar is not None:
            (task / "OSZICAR").write_text(oszicar, encoding="utf-8")
        if outcar is not None:
            (task / "OUTCAR").write_text(outcar, encoding="utf-8")
        return task

    def test_extracts_incar_parameters_and_classifies_ion_opt(self):
        self.write_task(
            "opt",
            "NSW = 10\nISIF = 2\nISPIN = 2\n",
            oszicar=" 1 F= -.123456 E0= -.123 d E =0\n",
            outcar=" free  energy   TOTEN  =       -123.456789 eV\n",
        )

        result = self.run_vaspjob("-v")

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("./opt", result.stdout)
        self.assertIn("ion_opt", result.stdout)
        self.assertIn("1/10", result.stdout)
        self.assertIn("↑↓", result.stdout)
        self.assertNotIn("integer expression expected", result.stderr)

    def test_identifies_neb_parent_directory(self):
        neb = self.workdir / "neb"
        for image in ["00", "01", "02"]:
            image_dir = neb / image
            image_dir.mkdir(parents=True)
            (image_dir / "POSCAR").write_text("POSCAR\n", encoding="utf-8")
        (neb / "INCAR").write_text("IMAGES = 3\nNSW = 20\n", encoding="utf-8")
        (neb / "01" / "OSZICAR").write_text(
            " 1 F= -.100000 E0= -.1 d E =0\n", encoding="utf-8"
        )

        result = self.run_vaspjob("-v")

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("./neb", result.stdout)
        self.assertIn("neb", result.stdout)
        self.assertIn("1/20", result.stdout)

    def test_nbands_alone_does_not_force_cohp_classification(self):
        self.write_task(
            "nbands",
            "NSW = 0\nNBANDS = 128\n",
            oszicar=" 1 F= -.500000 E0= -.5 d E =0\n",
        )

        result = self.run_vaspjob("-v")

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("./nbands", result.stdout)
        self.assertIn("single_point", result.stdout)
        self.assertNotIn("cohp", result.stdout)

    def test_wavecar_list_and_removal_handle_spaces(self):
        wavecar_dir = self.workdir / "path with spaces"
        wavecar_dir.mkdir()
        wavecar = wavecar_dir / "WAVECAR"
        wavecar.write_text("wavecar\n", encoding="utf-8")

        list_result = self.run_vaspjob("-m")
        self.assertEqual(list_result.returncode, 0, list_result.stderr)
        self.assertIn("./path with spaces/WAVECAR", (self.workdir / "wavecar_list.txt").read_text())

        remove_result = self.run_vaspjob("-r")
        self.assertEqual(remove_result.returncode, 0, remove_result.stderr)
        self.assertFalse(wavecar.exists())


if __name__ == "__main__":
    unittest.main()
