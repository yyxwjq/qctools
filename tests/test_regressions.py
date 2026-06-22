import importlib
import os
import tempfile
import unittest

import numpy as np
from ase import Atoms
from ase.io import read, write


class LoggingRegressionTests(unittest.TestCase):
    def test_importing_analysis_module_does_not_delete_existing_log(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            old_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                with open("qctools.log", "w", encoding="utf-8") as handle:
                    handle.write("keep me")

                import qctools.rdf

                importlib.reload(qctools.rdf)

                with open("qctools.log", encoding="utf-8") as handle:
                    self.assertEqual(handle.read(), "keep me")
            finally:
                os.chdir(old_cwd)


class ADFRegressionTests(unittest.TestCase):
    def test_adf_respects_pair_cutoff(self):
        from qctools.adf import get_angular_distribution

        image = Atoms(
            "HOH",
            positions=[(-5.0, 0.0, 0.0), (0.0, 0.0, 0.0), (5.0, 0.0, 0.0)],
            cell=[20.0, 20.0, 20.0],
            pbc=False,
        )

        result = get_angular_distribution([image], ["H", "O"], rcut=3.0, bin_size=10.0)

        self.assertEqual(result[("H", "O", "H")].sum(), 0)

    def test_adf_places_180_degree_angle_in_last_bin(self):
        from qctools.adf import get_angular_distribution

        image = Atoms(
            "HOH",
            positions=[(-1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
            cell=[10.0, 10.0, 10.0],
            pbc=False,
        )

        result = get_angular_distribution([image], ["H", "O"], rcut=2.0, bin_size=10.0)

        self.assertEqual(result[("H", "O", "H")][-1], 1)


class RDFRegressionTests(unittest.TestCase):
    def test_rdf_normalization_is_frame_averaged(self):
        from qctools.rdf import get_radial_distribution, normalize_rdf

        image = Atoms(
            "HO",
            positions=[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
            cell=[10.0, 10.0, 10.0],
            pbc=False,
        )

        one_frame_counts, one_frame_meta = get_radial_distribution(
            [image], ["H", "O"], cutoff=2.0, bin_size=0.5, collect_metadata=True
        )
        two_frame_counts, two_frame_meta = get_radial_distribution(
            [image, image], ["H", "O"], cutoff=2.0, bin_size=0.5, collect_metadata=True
        )

        one_frame_rdf = normalize_rdf(
            one_frame_counts, one_frame_meta, cutoff=2.0, bin_size=0.5
        )[("H", "O")]
        two_frame_rdf = normalize_rdf(
            two_frame_counts, two_frame_meta, cutoff=2.0, bin_size=0.5
        )[("H", "O")]

        np.testing.assert_allclose(one_frame_rdf[:, 1], two_frame_rdf[:, 1])


class CoordinationRegressionTests(unittest.TestCase):
    def test_group_coordnum_uses_switching_limit_at_r0(self):
        from qctools.coord import _group_coordnum_serial

        image = Atoms(
            "HH",
            positions=[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
            cell=[10.0, 10.0, 10.0],
            pbc=False,
        )

        self.assertAlmostEqual(
            _group_coordnum_serial(image, [0], [1], r0=1.0, en=6, ed=12),
            0.5,
        )


class EditAtomsRegressionTests(unittest.TestCase):
    def test_remove_preserves_original_order_of_remaining_atoms(self):
        from qctools.edit_atoms import remove

        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = os.path.join(tmpdir, "mixed.vasp")
            atoms = Atoms(
                "HOHC",
                positions=[
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 1.0),
                    (0.0, 1.0, 0.0),
                    (1.0, 0.0, 0.0),
                ],
                cell=[5.0, 5.0, 5.0],
                pbc=True,
            )
            write(input_file, atoms, format="vasp", vasp5=True, direct=True)

            remove(input_file, ["H"])

            output = read(os.path.join(tmpdir, "mixed_remove.vasp"))
            self.assertEqual(output.get_chemical_symbols(), ["O", "C"])


class MLErrorAnnotationRegressionTests(unittest.TestCase):
    def test_energy_error_structures_use_image_ids_not_energy_rows(self):
        from qctools.ml.error_img import err_structure_finding

        images = [
            Atoms("H", positions=[(0.0, 0.0, 0.0)], cell=[5.0, 5.0, 5.0], pbc=False),
            Atoms("He", positions=[(0.0, 0.0, 0.0)], cell=[5.0, 5.0, 5.0], pbc=False),
            Atoms("Li", positions=[(0.0, 0.0, 0.0)], cell=[5.0, 5.0, 5.0], pbc=False),
        ]
        for index, atoms in enumerate(images):
            atoms.info["source_index"] = index

        with tempfile.TemporaryDirectory() as tmpdir:
            old_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                np.savetxt(
                    "energy.data",
                    np.array(
                        [
                            [2, 0.0, 10.0],
                            [0, 0.0, 0.0],
                            [1, 0.0, 0.0],
                        ]
                    ),
                    fmt="%s",
                )

                err_structure_finding(
                    error_bar=0.5,
                    images=images,
                    fontsize=10,
                    replace_atom="Au",
                    Cutimg=True,
                    comment=False,
                    show_marginals=False,
                )

                err_images = read("Err-energy.xyz", ":")
                leave_images = read("leave-E-img.xyz", ":")
                self.assertEqual([atoms.info["source_index"] for atoms in err_images], [2])
                self.assertEqual(
                    [atoms.info["source_index"] for atoms in leave_images], [0, 1]
                )
            finally:
                os.chdir(old_cwd)

    def test_energy_error_annotations_accept_current_four_column_file(self):
        from qctools.ml.error_img import plot_scatter_with_marginals

        with tempfile.TemporaryDirectory() as tmpdir:
            old_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                np.savetxt(
                    "energy_error.txt",
                    np.array([[0.0, -1.0, -1.5, 0.5]]),
                    fmt="%s",
                )

                plot_scatter_with_marginals(
                    x=np.array([[-1.0], [-2.0]]),
                    y=np.array([[-1.5], [-2.1]]),
                    fontsize=10,
                    plot_type="energy",
                    Err_comment=True,
                    show_marginals=False,
                )

                self.assertTrue(os.path.exists("energy_error_analysis.png"))
            finally:
                os.chdir(old_cwd)

    def test_force_error_annotations_accept_component_column(self):
        from qctools.ml.error_img import plot_scatter_with_marginals

        with tempfile.TemporaryDirectory() as tmpdir:
            old_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                with open("force_error.txt", "w", encoding="utf-8") as handle:
                    handle.write("0 1 0.2 0.8 0.6 x\n")

                plot_scatter_with_marginals(
                    x=np.array([[0.2, 0.0, 0.0], [0.1, 0.2, 0.3]]),
                    y=np.array([[0.8, 0.0, 0.0], [0.2, 0.1, 0.4]]),
                    fontsize=10,
                    plot_type="force",
                    Err_comment=True,
                    show_marginals=False,
                )

                self.assertTrue(os.path.exists("force_error_analysis.png"))
            finally:
                os.chdir(old_cwd)


if __name__ == "__main__":
    unittest.main()
