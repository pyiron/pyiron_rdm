import json
import os
import unittest

from pyiron_rdm import ob_cfg_sfb1394
from pyiron_rdm.concept_dict import flatten_cdict


class TestObCFG(unittest.TestCase):
    def test_map_cdict_to_ob(self):
        d = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(d, "..", "static", "lammps_concept_dict.json")) as f:
            concept_dicts = json.load(f)
        props_dict = {}
        for key, concept_dict in concept_dicts.items():
            cdict = flatten_cdict(concept_dict)
            p_dict = ob_cfg_sfb1394.map_cdict_to_ob(
                user_name="test",
                cdict=cdict,
                concept_dict=concept_dict,
            )
            p_dict.pop("pyiron_conceptual_dictionary")
            p_dict.pop("description_multiline")
            if "date" in p_dict:
                p_dict.pop("date")
            props_dict[key] = p_dict
        with open(os.path.join(d, "..", "static", "lammps_props_dict_sfb.json")) as f:
            props_dict_ref = json.load(f)
        self.assertDictEqual(props_dict, props_dict_ref)


if __name__ == "__main__":
    unittest.main()
