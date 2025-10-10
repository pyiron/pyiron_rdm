import unittest

from pyiron_rdm import ob_OT_sfb1394, ob_OT_bam


class TestObOT(unittest.TestCase):
    def test_validate_options(self):
        with self.assertRaises(TypeError):
            ob_OT_bam.validate_options(**{"my_key": 123})
        with self.assertRaises(TypeError):
            ob_OT_sfb1394.validate_options(**{"my_key": 123})
        self.assertIsNone(
            ob_OT_bam.validate_options(**{"materials": "abc", "comments": "hello"})
        )
        self.assertIsNone(
            ob_OT_sfb1394.validate_options(**{"materials": "abc", "comments": "hello"})
        )
        with self.assertRaises(TypeError):
            ob_OT_sfb1394.validate_options(**{"defects": 123})
        with self.assertRaises(ValueError):
            ob_OT_sfb1394.validate_options(**{"defects": ["hello"]})
        self.assertIsNone(ob_OT_sfb1394.validate_options(**{"defects": ["surface"]}))


if __name__ == "__main__":
    unittest.main()
