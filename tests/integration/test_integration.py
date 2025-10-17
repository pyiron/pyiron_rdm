import json
import os
import unittest

from pybis_aixtended import OpenbisAixTended


class TestConf(unittest.TestCase):
    def test_conf_present(self):
        print(os.getcwd())
        self.assertTrue(os.path.isfile("some.conf"))

    def test_login_sfb1394_instance(self):
        ob = OpenbisAixTended.OpenbisWithS3(
            "https://openbis.imm.rwth-aachen.de/openbis/webapp/eln-lims/",
            s3_config_path="some.conf",
        )
        with open("pyironautouser") as f:
            ob.login("pyironautouser", f.readline())

        self.assertTrue(ob.is_session_active())
