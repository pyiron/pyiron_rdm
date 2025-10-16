import os
import unittest


class TestConf(unittest.TestCase):
    def test_conf_present(self):
        print(os.getcwd())
        self.assertTrue(os.isfile('some.conf'))
