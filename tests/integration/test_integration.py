import os
import unittest


class TestConf(unittest.TestCase):
    def test_env_present(self):
        self.assertEqual(os.env[pyironautouser]), 'some')

    def test_conf_present(self):
        print(os.getcwd())
        self.assertTrue(os.path.isfile('some.conf'))
