import json
import unittest

from pyiron_rdm import ob_cfg_sfb1394
from pyiron_rdm.concept_dict import flatten_cdict


class TestObCFG(unittest.TestCase):
    def test_map_cdict_to_ob(self):
        with open("../static/lammps_concept_dict.json") as f:
            concept_dicts = json.load(f)
        for concept_dict in concept_dicts.values():
            cdict = flatten_cdict(concept_dict)
            props_dict = ob_cfg_sfb1394.map_cdict_to_ob(
                user_name="test",
                cdict=cdict,
                concept_dict=concept_dict,
            )


if __name__ == "__main__":
    unittest.main()
