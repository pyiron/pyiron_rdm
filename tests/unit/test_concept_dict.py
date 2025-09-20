import unittest
import os

import numpy as np

from pyiron_rdm.concept_dict import process_lammps_job
from pyiron_atomistics import Project


class TestLammps(unittest.TestCase):
    def setUp(self):
        self._project = Project("test_lammps")

    def tearDown(self):
        self._project.remove_jobs(recursive=True, silently=True)
        self._project.remove(enable=True)

    def test_calc_minimize_pressure_zero(self):
        pr = Project("test_project")
        job = pr.create.job.Lammps("lmp")
        job.structure = pr.create.structure.bulk("Al", cubic=True)
        job.calc_minimize(pressure=0.0)
        job.run()
        job_dict = process_lammps_job(job)
        print(job_dict)
        self.assertEqual(job_dict["workflow_manager"]["label"], "pyiron")
        self.assertTrue(job_dict["software"][0]["label"].startswith("LAMMPS"))
        for item in job_dict["job_details"]:
            if item["label"] == "job_name":
                self.assertEqual(item["value"], job.job_name)
            elif item["label"] == "project_name":
                self.assertEqual(item["value"], os.path.basename(os.path.dirname(job.project.path)))
            elif item["label"] == "job_type":
                self.assertEqual(item["value"], "Lammps")
            elif item["label"] == "job_status":
                self.assertEqual(item["value"], "finished")
            elif item["label"] == "sim_coretime_hours":
                self.assertTrue(item["value"] < 0.2)
            elif item["label"] == "number_cores":
                self.assertEqual(item["value"], job.server.cores)
            else:
                self.assertTrue(item["label"] in ["job_starttime", "job_stoptime", "host"])
        self.assertEqual(job_dict["path"], os.path.abspath(os.path.join(job.working_directory, "../..", job.job_name)))
        context_dict = {
            'sample': 'http://purls.helmholtz-metadaten.de/cmso/AtomicScaleSample',
            'path': 'http://purls.helmholtz-metadaten.de/cmso/hasPath',
            'dof': 'http://purls.helmholtz-metadaten.de/asmo/hasRelaxationDOF',
            'inputs': 'http://purls.helmholtz-metadaten.de/asmo/hasInputParameter',
            'label': 'http://www.w3.org/2000/01/rdf-schema#label',
            'unit': 'http://purls.helmholtz-metadaten.de/asmo/hasUnit',
            'value': 'http://purls.helmholtz-metadaten.de/asmo/hasValue',
            'outputs': 'http://purls.helmholtz-metadaten.de/cmso/hasCalculatedProperty',
            'workflow_manager': 'http://demo.fiz-karlsruhe.de/matwerk/E457491',
            'molecular_dynamics': 'http://purls.helmholtz-metadaten.de/asmo/MolecularDynamics',
            'molecular_statics': 'http://purls.helmholtz-metadaten.de/asmo/MolecularStatics',
            'ensemble': 'http://purls.helmholtz-metadaten.de/asmo/hasStatisticalEnsemble',
            'job_details': 'http://id-from-pmdco-pending',
            'periodic_boundary_condition': 'http://purls.helmholtz-metadaten.de/asmo/PeriodicBoundaryCondition',
            'initial_temperature': 'http://purls.helmholtz-metadaten.de/asmo/Temperature',
            'initial_pressure': 'http://purls.helmholtz-metadaten.de/asmo/Pressure',
            'target_temperature': 'http://purls.helmholtz-metadaten.de/asmo/Temperature',
            'target_pressure': 'http://purls.helmholtz-metadaten.de/asmo/Pressure',
            'ionic_energy_tolerance': 'http://purls.helmholtz-metadaten.de/asmo/InputParameter',
            'force_tolerace': 'http://purls.helmholtz-metadaten.de/asmo/InputParameter',
            'maximum_iterations': 'http://purls.helmholtz-metadaten.de/asmo/InputParameter',
            'potential': 'http://purls.helmholtz-metadaten.de/asmo/EmbeddedAtomModel',
            'average_temperature': 'http://purls.helmholtz-metadaten.de/asmo/Temperature',
            'average_pressure': 'http://purls.helmholtz-metadaten.de/asmo/Pressure',
            'average_total_energy': 'http://purls.helmholtz-metadaten.de/asmo/TotalEnergy',
            'final_total_energy': 'http://purls.helmholtz-metadaten.de/asmo/TotalEnergy',
            'final_potential_energy': 'http://purls.helmholtz-metadaten.de/asmo/PotentialEnergy',
            'average_total_volume': 'http://purls.helmholtz-metadaten.de/asmo/Volume',
            'final_total_volume': 'http://purls.helmholtz-metadaten.de/asmo/Volume',
            'final_maximum_force': 'http://purls.helmholtz-metadaten.de/asmo/Force',
            'number_ionic_steps': 'http://purls.helmholtz-metadaten.de/asmo/NumberOfIonicSteps',
            'time_step': 'http://purls.helmholtz-metadaten.de/asmo/TimeStep',
            'simulation_time': 'http://purls.helmholtz-metadaten.de/asmo/Time',
            'LAMMPS': 'http://demo.fiz-karlsruhe.de/matwerk/E447986'
        }
        for k,v in context_dict.items():
            self.assertEqual(job_dict["@context"][k], v)
        for item in job_dict["outputs"]:
            if item["label"] == 'average_total_energy':
                self.assertEqual(item["unit"], "EV")
                self.assertEqual(item["value"], np.round(job["output/generic/energy_tot"][-1], 2))
            elif item["label"] == 'average_total_volume':
                self.assertEqual(item["unit"], "ANGSTROM3")
                self.assertEqual(item["value"], np.round(job["output/generic/volume"][-1], 4))
            elif item["label"] == 'final_total_volume':
                self.assertEqual(item["unit"], "ANGSTROM3")
                self.assertEqual(item["value"], np.round(job["output/generic/volume"][-1], 4))
            elif item["label"] == 'final_total_energy':
                self.assertEqual(item["unit"], "EV")
                self.assertEqual(item["value"], np.round(job["output/generic/energy_tot"][-1], 2))
            elif item["label"] == 'final_potential_energy':
                self.assertEqual(item["unit"], "EV")
                self.assertEqual(item["value"], np.round(job["output/generic/energy_pot"][-1], 2))
            elif item["label"] == 'number_ionic_steps':
                self.assertEqual(item["value"], 1)
            elif item["label"] == 'final_maximum_force':
                self.assertEqual(item["unit"], "EV-PER-ANGSTROM")
                self.assertEqual(item["value"], 1.4e-15)
            else:
                raise ValueError()
        for item in job_dict["dof"]:
            self.assertTrue(item in ['http://purls.helmholtz-metadaten.de/asmo/AtomicPositionRelaxation', 'http://purls.helmholtz-metadaten.de/asmo/CellVolumeRelaxation'])
        self.assertEqual(job_dict["molecular_statics"]['minimization_algorithm'], "cg")
        self.assertTrue(job_dict["molecular_statics"]['periodicity_in_x'])
        self.assertTrue(job_dict["molecular_statics"]['periodicity_in_y'])
        self.assertTrue(job_dict["molecular_statics"]['periodicity_in_z'])
        for item in job_dict["molecular_statics"]["inputs"]:
            if item["label"] == 'initial_temperature':
                self.assertEqual(item["unit"], "K")
                self.assertIsNone(item["value"])
            elif item["label"] == 'target_temperature':
                self.assertEqual(item["unit"], "K")
                self.assertIsNone(item["value"])
            elif item["label"] == 'initial_pressure':
                self.assertEqual(item["unit"], "GigaPA")
                self.assertIsNone(item["value"])
            elif item["label"] == 'target_pressure':
                self.assertEqual(item["unit"], "GigaPA")
                self.assertEqual(item["value"], 0.0)
            elif item["label"] == 'ionic_energy_tolerance':
                self.assertEqual(item["unit"], "EV")
                self.assertEqual(item["value"], 0.0)
            elif item["label"] == 'force_tolerance':
                self.assertEqual(item["unit"], "EV-PER-ANGSTROM")
                self.assertEqual(item["value"], 0.0001)
            elif item["label"] == 'maximum_iterations':
                self.assertEqual(item["value"], 100000)
            elif item["label"] == 'timestep':
                self.assertIsNone(item["value"])
                self.assertEqual(item["unit"], "PICOSECOND")
            elif item["label"] == 'simulation_time':
                self.assertIsNone(item["value"])
                self.assertEqual(item["unit"], "PICOSECOND")
            else:
                raise ValueError()
        self.assertIsNone(job_dict["molecular_statics"]['ensemble'])
        self.assertEqual(job_dict["molecular_statics"]["potential"]["label"], '1995--Angelo-J-E--Ni-Al-H--LAMMPS--ipr1')
        self.assertEqual(job_dict["molecular_statics"]["potential"]["url"], 'https://doi.org/10.1088%2F0965-0393%2F3%2F3%2F001')
