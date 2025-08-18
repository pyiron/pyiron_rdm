import unittest

import pyiron_rdm


class TestVersion(unittest.TestCase):
    def test_version(self):
        version = pyiron_rdm.__version__
        print(version)
        self.assertTrue(version.startswith("0"))
