import json
import os
import unittest

from ase.build import bulk

from pyiron_rdm import concept_dict


class TestConceptDict(unittest.TestCase):
    def test_get_unit_cell_parameters(self):
        Al = bulk("Al", cubic=True)
        self.assertDictEqual(
            concept_dict.get_unit_cell_parameters(Al),
            {
                "a": 4.05,
                "alpha": 90.0,
                "beta": 90.0,
                "gamma": 90.0,
                "volume": 66.4301,
                "space_group": "Fm-3m",
                "space_group_number": 225,
                "bravais_lattice": "fcc",
            },
        )
        Fe = bulk("Fe")
        self.assertDictEqual(
            concept_dict.get_unit_cell_parameters(Fe),
            {
                "a": 2.87,
                "alpha": 90.0,
                "beta": 90.0,
                "gamma": 90.0,
                "volume": 23.6399,
                "space_group": "Im-3m",
                "space_group_number": 229,
                "bravais_lattice": "bcc",
            },
        )
        Mg = bulk("Mg")
        self.assertEqual(
            concept_dict.get_unit_cell_parameters(Mg)["bravais_lattice"],
            "hcp",
        )

    def test_process_structure_crystal(self):
        Fe = bulk("Fe", cubic=True)
        result = concept_dict.process_structure_crystal(
            "path", "name", Fe, "structure_name", "structure_path"
        )
        self.assertDictEqual(
            result,
            {
                "@context": {
                    "path": "http://purls.helmholtz-metadaten.de/cmso/hasPath",
                    "unit_cell": "http://purls.helmholtz-metadaten.de/cmso/UnitCell",
                    "atoms": "http://purls.helmholtz-metadaten.de/cmso/Atom",
                    "molecules": "http://purls.helmholtz-metadaten.de/cmso/Molecule",
                    "bravais_lattice": "http://purls.helmholtz-metadaten.de/cmso/hasBravaisLattice",
                    "chemical_species": "http://purls.helmholtz-metadaten.de/cmso/ChemicalSpecies",
                    "simulation_cell": "http://www.w3.org/2000/01/rdf-schema#label",
                    "label": "http://www.w3.org/2000/01/rdf-schema#label",
                    "unit": "http://purls.helmholtz-metadaten.de/cmso/hasUnit",
                    "value": "http://purls.helmholtz-metadaten.de/asmo/hasValue",
                    "vector": "http://purls.helmholtz-metadaten.de/cmso/Vector",
                    "job_details": "http://id-from-pmdco-pending",
                    "workflow_manager": "http://demo.fiz-karlsruhe.de/matwerk/E457491",
                    "lattice_parameter_a": "http://purls.helmholtz-metadaten.de/cmso/hasLatticeParameter",
                    "lattice_angle_alpha": "http://purls.helmholtz-metadaten.de/cmso/hasAngle",
                    "lattice_angle_beta": "http://purls.helmholtz-metadaten.de/cmso/hasAngle",
                    "lattice_angle_gamma": "http://purls.helmholtz-metadaten.de/cmso/hasAngle",
                    "lattice_volume": "http://purls.helmholtz-metadaten.de/asmo/Volume",
                    "space_group": "http://purls.helmholtz-metadaten.de/cmso/hasSpaceGroup",
                    "simulation_cell_lengths": "http://purls.helmholtz-metadaten.de/cmso/hasLength",
                    "simulation_cell_vectors": "http://purls.helmholtz-metadaten.de/cmso/hasVector",
                    "simulation_cell_volume": "http://purls.helmholtz-metadaten.de/cmso/hasVolume",
                    "simulation_cell_angle": "http://purls.helmholtz-metadaten.de/cmso/hasAngle",
                },
                "atoms": [
                    {"value": 2, "label": "Fe"},
                    {"value": 2, "label": "total_number_atoms"},
                ],
                "simulation_cell": [
                    {
                        "value": "[2.87, 2.87, 2.87]",
                        "unit": "ANGSTROM",
                        "label": "simulation_cell_lengths",
                    },
                    {
                        "value": "[array([2.87, 0.  , 0.  ]), array([0.  , 2.87, 0.  ]), array([0.  , 0.  , 2.87])]",
                        "unit": "ANGSTROM",
                        "label": "simulation_cell_vectors",
                    },
                    {
                        "value": "[90.0, 90.0, 90.0]",
                        "unit": "DEGREES",
                        "label": "simulation_cell_angles",
                    },
                    {
                        "value": 23.6399,
                        "unit": "ANGSTROM3",
                        "label": "simulation_cell_volume",
                    },
                ],
                "workflow_manager": {"label": ""},
                "job_details": [
                    {"label": "structure_name", "value": "structure_name"},
                    {"label": "project_name", "value": "name"},
                ],
                "path": "structure_path",
            },
        )

    def test_flatten_cdict(self):
        d = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(d, "..", "static", "lammps_concept_dict.json")) as f:
            concept_dicts = json.load(f)
        flattened_dicts = {
            key: concept_dict.flatten_cdict(value) for key, value in concept_dicts.items()
        }
        with open(os.path.join(d, "..", "static", "lammps_flattened_dict.json")) as f:
            ref_flattened_dicts = json.load(f)
        self.assertDictEqual(flattened_dicts, ref_flattened_dicts)


if __name__ == "__main__":
    unittest.main()
