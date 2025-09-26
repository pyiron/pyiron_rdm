import unittest

from ase.build import bulk

from pyiron_rdm.concept_dict import get_unit_cell_parameters


class TestVersion(unittest.TestCase):
    def test_get_unit_cell_parameters(self):
        Al = bulk("Al", cubic=True)
        self.assertEqual(
            get_unit_cell_parameters(Al),
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
        self.assertEqual(
            get_unit_cell_parameters(Fe),
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
            get_unit_cell_parameters(Mg)["bravais_lattice"],
            "hcp",
        )


if __name__ == "__main__":
    unittest.main()
