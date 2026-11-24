import unittest

from pyiron_rdm import ob_OT_bam, ob_OT_sfb1394


class TestObOT(unittest.TestCase):
    def test_validate_options(self):
        masterdata_dict = {"allowed_defects": {"surface"}}
        with self.assertRaises(TypeError):
            ob_OT_bam.validate_options(**{"my_key": 123})
        with self.assertRaises(TypeError):
            ob_OT_sfb1394.validate_options(
                **{"my_key": 123}, masterdata_dict=masterdata_dict
            )
        self.assertIsNone(
            ob_OT_bam.validate_options(**{"materials": "abc", "comments": "hello"})
        )
        self.assertIsNone(
            ob_OT_sfb1394.validate_options(
                **{"materials": "abc", "comments": "hello"},
                masterdata_dict=masterdata_dict,
            )
        )
        with self.assertRaises(TypeError):
            ob_OT_sfb1394.validate_options(
                **{"defects": 123}, masterdata_dict=masterdata_dict
            )
        with self.assertRaises(ValueError):
            ob_OT_sfb1394.validate_options(
                **{"defects": ["hello"]}, masterdata_dict=masterdata_dict
            )
        self.assertIsNone(
            ob_OT_sfb1394.validate_options(
                **{"defects": ["surface"]}, masterdata_dict=masterdata_dict
            )
        )


if __name__ == "__main__":
    unittest.main()
