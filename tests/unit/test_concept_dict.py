import unittest

from ase.build import bulk

from pyiron_rdm.concept_dict import get_unit_cell_parameters


class TestVersion(unittest.TestCase):
    def test_get_unit_cell_parameters(self):
        structure = bulk("Al", cubic=True)
        self.assertEqual(
            get_unit_cell_parameters(structure),
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


if __name__ == "__main__":
    unittest.main()
