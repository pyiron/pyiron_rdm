
import os
import unittest

from datetime import datetime
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


class TestOpenBISinteractions(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.created_permIDs = []
        ob = OpenbisAixTended.OpenbisWithS3(
            "https://openbis.imm.rwth-aachen.de/openbis/webapp/eln-lims/",
            s3_config_path="some.conf",
        )
        with open("pyironautouser") as f:
            ob.login("pyironautouser", f.readline())
        cls.pybis_w_s3 = ob
        date = datetime.now()
        date_code = f'{str(date.date()).replace("-","")}{date.hour}{date.minute}{date.second}'
        cls.pybis_pr = ob.pybis_w_s3.new_project(
        space       = 'PYIRONAUTOUSER',
        code        = 'CiProject' + date_code,
        description = 'Temporary CI project from pyiron_rdm'
        )
        cls.pybis_pr.save()
        super().setUpClass()

    @classmethod
    def tearDownClass(cls):
        # ToDo remove all created_permIDs
        cls.pybis_w_s3.logout()
        super().tearDownClass()

    def test_nothing(self):
        self.assertTrue(True)
