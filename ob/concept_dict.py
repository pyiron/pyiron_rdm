"""
LAMMPS specific functions for parsing

Use this a reference for specific implementations

TODO:
- add vasp module
- Instead of job.to_dict() use pyiron_dataclass - convert from hdf5 to dataclass - convert to conc. dict.



"""
import os
import numpy as np
import ast
import json
from typing import Optional

def process_general_job(job):
    method_dict = {}
    add_simulation_software(job, method_dict)
    get_simulation_folder(job, method_dict)
    file_name = job.path + '_concept_dict.json'
    with open(file_name, 'w') as f:
        json.dump(method_dict, f, indent=2)
    return method_dict

def process_lammps_job(job):
    method_dict = {}
    add_lammps_contexts(method_dict)
    #get_structures(job, method_dict)
    identify_lammps_method(job, method_dict)
    extract_lammps_calculated_quantities(job, method_dict)
    add_simulation_software(job, method_dict)
    get_simulation_folder(job, method_dict)
    file_name = job.path + '_concept_dict.json'
    with open(file_name, 'w') as f:
        json.dump(method_dict, f, indent=2)
    return method_dict

def process_structure_crystal(pr, structure, structure_name, structure_path, structure_parameters: dict = None):
    sample_dict = {}
    add_structure_contexts(sample_dict)
    get_chemical_species(structure, sample_dict)
    identify_structure_parameters(structure_parameters, sample_dict)
    get_simulation_cell(structure, sample_dict)
    add_structure_software(pr, structure_name, sample_dict)
    get_structure_folder(structure_path, sample_dict)
    json_file_name = structure_path + structure_name +'_concept_dict.json'
    with open(json_file_name, 'w') as f:
        json.dump(sample_dict, f, indent=2)
    return sample_dict

def process_murnaghan_job(job):
    method_dict = {}
    add_murnaghan_contexts(method_dict)
    identify_murnaghan_method(job, method_dict)
    extract_murnaghan_calculated_quantities(job, method_dict)
    add_simulation_software(job, method_dict)
    get_simulation_folder(job, method_dict)
    file_name = job.path + '_concept_dict.json'
    with open(file_name, 'w') as f:
        json.dump(method_dict, f, indent=2)
    return method_dict

def process_vasp_job(job):
    method_dict = {}
    add_vasp_contexts(method_dict)
    identify_vasp_method(job, method_dict)
    extract_vasp_calculated_quantities(job, method_dict)
    add_simulation_software(job, method_dict)
    get_simulation_folder(job, method_dict)
    file_name = job.path + '_concept_dict.json'
    with open(file_name, 'w') as f:
        json.dump(method_dict, f, indent=2)
    return method_dict

def add_lammps_contexts(method_dict):
    method_dict['@context'] = {}
    method_dict['@context']['sample'] = 'http://purls.helmholtz-metadaten.de/cmso/AtomicScaleSample'
    method_dict['@context']['path'] = 'http://purls.helmholtz-metadaten.de/cmso/hasPath'
    method_dict['@context']['dof'] = 'http://purls.helmholtz-metadaten.de/asmo/hasRelaxationDOF'
    method_dict['@context']['inputs'] = 'http://purls.helmholtz-metadaten.de/asmo/hasInputParameter'
    method_dict['@context']['label'] = 'http://www.w3.org/2000/01/rdf-schema#label'
    method_dict['@context']['unit'] = 'http://purls.helmholtz-metadaten.de/asmo/hasUnit'
    method_dict['@context']['value'] = 'http://purls.helmholtz-metadaten.de/asmo/hasValue'
    method_dict['@context']['outputs'] = 'http://purls.helmholtz-metadaten.de/cmso/hasCalculatedProperty'
    method_dict['@context']['workflow_manager'] = 'http://demo.fiz-karlsruhe.de/matwerk/E457491'
    #method_dict['@context']['software'] = ''
    method_dict['@context']['molecular_dynamics'] = 'http://purls.helmholtz-metadaten.de/asmo/MolecularDynamics'
    method_dict['@context']['molecular_statics'] = 'http://purls.helmholtz-metadaten.de/asmo/MolecularStatics'
    method_dict['@context']['ensemble'] = 'http://purls.helmholtz-metadaten.de/asmo/hasStatisticalEnsemble'
    method_dict['@context']['job_details'] = 'id-from-pmdco-pending'
    method_dict['@context']['periodic_boundary_condition'] = 'http://purls.helmholtz-metadaten.de/asmo/PeriodicBoundaryCondition'
    method_dict['@context']['temperature'] = 'http://purls.helmholtz-metadaten.de/asmo/Temperature'
    method_dict['@context']['pressure'] = 'http://purls.helmholtz-metadaten.de/asmo/Pressure'
    method_dict['@context']['ionic_energy_tolerance'] = 'http://purls.helmholtz-metadaten.de/asmo/InputParameter'
    method_dict['@context']['force_tolerace'] = 'http://purls.helmholtz-metadaten.de/asmo/InputParameter'
    method_dict['@context']['maximum_iterations'] = 'http://purls.helmholtz-metadaten.de/asmo/InputParameter'
    method_dict['@context']['potential'] = "http://purls.helmholtz-metadaten.de/asmo/InteratomicPotential"
    method_dict['@context']['average_total_energy'] = 'http://purls.helmholtz-metadaten.de/asmo/TotalEnergy'
    method_dict['@context']['final_total_energy'] = 'http://purls.helmholtz-metadaten.de/asmo/TotalEnergy'
    method_dict['@context']['final_potential_energy'] = 'http://purls.helmholtz-metadaten.de/asmo/PotentialEnergy'
    method_dict['@context']['average_total_volume'] = 'http://purls.helmholtz-metadaten.de/asmo/Volume'
    method_dict['@context']['final_total_volume'] = 'http://purls.helmholtz-metadaten.de/asmo/Volume'
    method_dict['@context']['final_maximum_force'] = 'http://purls.helmholtz-metadaten.de/asmo/Force'
    method_dict['@context']['number_ionic_steps'] = 'http://purls.helmholtz-metadaten.de/asmo/NumberOfIonicSteps'
    method_dict['@context']['LAMMPS'] = 'http://demo.fiz-karlsruhe.de/matwerk/E447986'

def get_structures(job, method_dict):
    initial_pyiron_structure = job.structure
    final_pyiron_structure = job.get_structure(frame=-1)
    
    method_dict['sample'] =  {'initial':initial_pyiron_structure, 
                            'final': final_pyiron_structure}

def identify_lammps_method(job, method_dict):
    job_dict = job.input.to_dict()
    input_dict = {
        job_dict["control_inp/data_dict"]["Parameter"][x]: job_dict[
            "control_inp/data_dict"
        ]["Value"][x]
        for x in range(len(job_dict["control_inp/data_dict"]["Parameter"]))
    }
    dof = []
    temp = None
    press = None
    md_method = None
    ensemble = None

    if "min_style" in input_dict.keys():
        dof.append("http://purls.helmholtz-metadaten.de/asmo/AtomicPositionRelaxation")
        if job.input.control['fix___ensemble'] != 'all nve': dof.append("http://purls.helmholtz-metadaten.de/asmo/CellVolumeRelaxation")
        md_method = "molecular_statics"
        e_tol = float(input_dict['minimize'].split()[0])
        f_tol = float(input_dict['minimize'].split()[1])
        maxiter = int(input_dict['minimize'].split()[2])

    elif "nve" in input_dict["fix___ensemble"]:
        if int(input_dict["run"]) == 0:
            md_method = "molecular_statics"
            ensemble = "http://purls.helmholtz-metadaten.de/asmo/MicrocanonicalEnsemble"

        elif int(input_dict["run"]) > 0:
            dof.append("http://purls.helmholtz-metadaten.de/asmo/AtomicPositionRelaxation")
            md_method = "molecular_dynamics"
            ensemble = "http://purls.helmholtz-metadaten.de/asmo/MicrocanonicalEnsemble"

    elif "nvt" in input_dict["fix___ensemble"]:
        raw = input_dict["fix___ensemble"].split()
        temp = float(raw[3])
        dof.append("http://purls.helmholtz-metadaten.de/asmo/AtomicPositionRelaxation")
        md_method = "molecular_dynamics"
        ensemble = "http://purls.helmholtz-metadaten.de/asmo/CanonicalEnsemble"

    elif "npt" in input_dict["fix___ensemble"]:
        dof.append("http://purls.helmholtz-metadaten.de/asmo/AtomicPositionRelaxation")
        dof.append("http://purls.helmholtz-metadaten.de/asmo/CellVolumeRelaxation")
        if "aniso" in input_dict["fix___ensemble"]:
            dof.append("http://purls.helmholtz-metadaten.de/asmo/CellShapeRelaxation")
        md_method = "molecular_dynamics"
        raw = input_dict["fix___ensemble"].split()
        temp = float(raw[3])
        press = float(raw[7])
        ensemble = "http://purls.helmholtz-metadaten.de/asmo/IsothermalIsobaricEnsemble"

    method_dict[md_method] = {}

    if md_method == "molecular_statics":
        method_dict[md_method]['minimization_algorithm'] = input_dict['min_style']

    input_dict = {
        job_dict["control_inp/data_dict"]["Parameter"][x]: job_dict[
            "control_inp/data_dict"
        ]["Value"][x]
        for x in range(len(job_dict["control_inp/data_dict"]["Parameter"]))
                }
    pb = []
    pb.append(input_dict['boundary'])
    if (pb[0][0] == "p"):
        method_dict[md_method]['periodicity_in_x'] = True
    else:
        method_dict[md_method]['periodicity_in_x'] = False
    if (pb[0][2] == "p"):
        method_dict[md_method]['periodicity_in_y'] = True
    else:
        method_dict[md_method]['periodicity_in_y'] = False
    if (pb[0][4] == "p"):
        method_dict[md_method]['periodicity_in_z'] = True
    else:
        method_dict[md_method]['periodicity_in_z'] = False

    method_dict[md_method]['inputs'] = []

    temperature = {}
    temperature["value"] = temp
    temperature["unit"] = "K"
    temperature["label"] = "temperature"

    method_dict[md_method]['inputs'].append(temperature)

    pressure = {}
    pressure["value"] = press
    pressure["unit"] = "GigaPA"
    pressure["label"] = "pressure"

    method_dict[md_method]['inputs'].append(pressure)

    energy_tol = {}
    energy_tol["value"] = e_tol
    energy_tol["unit"] = "EV"
    energy_tol["label"] = "ionic_energy_tolerance"

    method_dict[md_method]['inputs'].append(energy_tol)

    force_tol = {}
    force_tol["value"] = f_tol
    force_tol["unit"] = "EV-PER-ANGSTROM"
    force_tol["label"] = "force_tolerance"

    method_dict[md_method]['inputs'].append(force_tol)

    maximum_iterations = {}
    maximum_iterations["value"] = maxiter
    maximum_iterations["label"] = "maximum_iterations"

    method_dict[md_method]['inputs'].append(maximum_iterations)

    method_dict[md_method]["ensemble"] = ensemble
    
    method_dict["dof"] = dof


    # now process potential
    inpdict = job.input.to_dict()
    ps = inpdict["potential_inp/data_dict"]["Value"][0]
    name = inpdict["potential_inp/potential/Name"]
    potstr = job.input.to_dict()["potential_inp/potential/Citations"]
    potdict = ast.literal_eval(potstr[1:-1])
    url = None
    if "url" in potdict[list(potdict.keys())[0]].keys():
        url = potdict[list(potdict.keys())[0]]["url"]
    if 'meam' in ps:
        method_dict['@context']['potential'] = "http://purls.helmholtz-metadaten.de/asmo/ModifiedEmbeddedAtomModel"
    elif 'eam' in ps:
        method_dict['@context']['potential'] = "http://purls.helmholtz-metadaten.de/asmo/EmbeddedAtomModel"
    elif 'lj' in ps:
        method_dict['@context']['potential'] = "http://purls.helmholtz-metadaten.de/asmo/LennardJonesPotential"
    elif 'ace' in ps:
        method_dict['@context']['potential'] = "http://purls.helmholtz-metadaten.de/asmo/MachineLearningPotential"

    method_dict[md_method]["potential"] = {}
    method_dict[md_method]["potential"]["label"] = name
    if url is not None:
        method_dict[md_method]["potential"]["url"] = url

def extract_lammps_calculated_quantities(job, method_dict):
    """
    Extracts calculated quantities from a job.

    Parameters
    ----------
    job : pyiron.Job
        The job object containing the calculated quantities.

    Returns
    -------
    list
        A list of dictionaries, each containing the label, value, unit, and associate_to_sample of a calculated quantity.

    """
    aen = np.mean(job.output.energy_tot)
    fen = job.output.energy_tot[-1]
    fpe = job.output.energy_tot[-1] # TODO: are these always supposed to be the same?
    avol = np.mean(job.output.volume)
    fvol = job.output.volume[-1]
    fmax = job.output.force_max[-1]
    nionic = len(job.output.steps)-1
    outputs = []
    outputs.append(
        {
            "label": "average_total_energy",
            "value": np.round(aen, decimals=4),
            "unit": "EV",
        }
    )
    outputs.append(
        {
            "label": "final_total_energy",
            "value": np.round(fen, decimals=4),
            "unit": "EV",
        }
    )
    outputs.append(
        {
            "label": "final_potential_energy",
            "value": np.round(fpe, decimals=4),
            "unit": "EV",
        }
    )
    outputs.append(
        {
            "label": "average_total_volume",
            "value": np.round(avol, decimals=4),
            "unit": "ANGSTROM3",
        }
    )
    outputs.append(
        {
            "label": "final_total_volume",
            "value": np.round(fvol, decimals=4),
            "unit": "ANGSTROM3",
        }
    )
    outputs.append(
        {
            "label": "final_maximum_force",
            "value": np.round(fmax, decimals=16),
            "unit": "EV-PER-ANGSTROM",
        }
    )
    outputs.append(
        {
            "label": "number_ionic_steps",
            "value": nionic,
        }
    )
    method_dict['outputs'] =  outputs

def add_simulation_software(job, method_dict):
    method_dict["workflow_manager"] = {}
    import subprocess
    import platform
    try:
        if "Windows" in platform.system():
            output1 = subprocess.check_output(['findstr', 'pyiron_atomistics', job.path.replace('/', '\\') + '_environment.yml'])
        else:
            output1 = subprocess.check_output(['grep', 'pyiron_atomistics', job.path + '_environment.yml'])
        s1 = str((output1.decode('utf-8')))
    except:
        s1 = ''
    try:
        if "Windows" in platform.system():
            output2 = subprocess.check_output(['findstr', 'pyiron_workflow', job.path.replace('/', '\\') + '_environment.yml'])
        else:
            output2 = subprocess.check_output(['grep', 'pyiron_workflow', job.path + '_environment.yml'])
        s2 = str((output2.decode('utf-8')))
    except:
        s2 = ''
    try:
        if "Windows" in platform.system():
            output3 = subprocess.check_output(['findstr', 'pyironflow', job.path.replace('/', '\\') + '_environment.yml'])
        else:
            output3 = subprocess.check_output(['grep', 'pyironflow', job.path + '_environment.yml'])
        s3 = str((output3.decode('utf-8')))
    except:
        s3 = ''
    try:
        if "Windows" in platform.system():
            output4 = subprocess.check_output(['findstr', 'executorlib', job.path.replace('/', '\\') + '_environment.yml'])
        else:
            output4 = subprocess.check_output(['grep', 'executorlib', job.path + '_environment.yml'])
        s4 = str((output4.decode('utf-8')))
    except:
        s4 = ''

    #hdf_ver = job.to_dict()['HDF_VERSION']
    try:
        st1 = 'p' + s1.split('=')[0].split('p')[1] + "=" + s1.split('=')[1] + ', '
    except:
        st1 = ''
    try:
        st2 = 'p' + s2.split('=')[0].split('p')[1] + "=" + s2.split('=')[1] + ', '
    except:
        st2 = ''
    try:
        st3 = 'p' + s3.split('=')[0].split('p')[1] + "=" + s3.split('=')[1] + ', '
    except:
        st3 = ''
    try:
        st4 = s4.split('=')[0].split('- ')[1] + "=" + s4.split('=')[1]
    except:
        st4 = ''
    st = st1 + st2 + st3 + st4
    
    #+ ', pyiron_HDF_version=' + hdf_ver
    method_dict["workflow_manager"]["label"] = st

    pyiron_job_details = []
    pyiron_job_details.append(
        {
            "label": "job_name",
            "value": job.name,
        }
    )
    pyiron_job_details.append(
        {
            "label": "project_name",
            "value": job.project.name,
        }
    )
    pyiron_job_details.append(
        {
            "label": "job_type",
            "value": job.database_entry.hamilton,
        }
    )
    pyiron_job_details.append(
        {
            "label": "job_status",
            "value": str(job.status),
        }
    )
    pyiron_job_details.append(
        {
            "label": "job_starttime",
            "value": str(job.database_entry.timestart.strftime("%Y-%m-%d %H:%M:%S")),
        }
    )
    pyiron_job_details.append(
        {
            "label": "job_stoptime",
            "value": str(job.database_entry.timestop.strftime("%Y-%m-%d %H:%M:%S")),
        }
    )
    pyiron_job_details.append(
        {
            "label": "sim_coretime_hours",
            "value": np.round(job.database_entry.totalcputime/3600, 6),
        }
    )
    pyiron_job_details.append(
        {
            "label": "number_cores",
            "value": job.to_dict()['server']['cores'],
        }
    )
    pyiron_job_details.append(
        {
            "label": "host",
            "value": job.to_dict()['server']['host'],
        }
    )
 
    try:
        software = {
            "label": job.to_dict()['executable']['name'].upper() + ' ' + job.to_dict()['executable']['version'],
        }
        method_dict["software"] = [software]
    except KeyError:
        pass

    method_dict["job_details"] = pyiron_job_details

def get_simulation_folder(job, method_dict):
    method_dict['path'] = job.path

def add_murnaghan_contexts(method_dict):
    method_dict['@context'] = {}
    method_dict['@context']['sample'] = 'http://purls.helmholtz-metadaten.de/cmso/AtomicScaleSample'
    method_dict['@context']['path'] = 'http://purls.helmholtz-metadaten.de/cmso/hasPath'
    method_dict['@context']['equation_of_state_fit'] = 'http://purls.helmholtz-metadaten.de/asmo/EquationOfStateFit'
    method_dict['@context']['inputs'] = 'http://purls.helmholtz-metadaten.de/asmo/hasInputParameter'
    method_dict['@context']['label'] = 'http://www.w3.org/2000/01/rdf-schema#label'
    method_dict['@context']['unit'] = 'http://purls.helmholtz-metadaten.de/asmo/hasUnit'
    method_dict['@context']['value'] = 'http://purls.helmholtz-metadaten.de/asmo/hasValue'
    method_dict['@context']['outputs'] = 'http://purls.helmholtz-metadaten.de/cmso/hasCalculatedProperty'
    method_dict['@context']['workflow_manager'] = 'http://demo.fiz-karlsruhe.de/matwerk/E457491'
    #method_dict['@context']['software'] = ''
    method_dict['@context']['job_details'] = 'id-from-pmdco-pending'
    method_dict['@context']['strain_axes'] = 'http://purls.helmholtz-metadaten.de/asmo/InputParameter'
    method_dict['@context']['number_of_data_points'] = 'http://purls.helmholtz-metadaten.de/asmo/InputParameter'
    method_dict['@context']['volume_range'] = 'http://purls.helmholtz-metadaten.de/asmo/VolumeRange'
    method_dict['@context']['equilibrium_bulk_modulus'] = 'http://purls.helmholtz-metadaten.de/asmo/BulkModulus'
    method_dict['@context']['equilibrium_total_energy'] = 'http://purls.helmholtz-metadaten.de/asmo/TotalEnergy'
    method_dict['@context']['equilibrium_volume'] = 'http://purls.helmholtz-metadaten.de/asmo/Volume'
    
def identify_murnaghan_method(job, method_dict):

    if job.input['fit_type'] == 'birchmurnaghan':
        method_dict['equation_of_state_fit'] = "http://purls.helmholtz-metadaten.de/asmo/BirchMurnaghan"
    elif job.input['fit_type'] == 'murnaghan':
        method_dict['equation_of_state_fit'] = "http://purls.helmholtz-metadaten.de/asmo/Murnaghan"
    elif job.input['fit_type'] == 'vinet':
        method_dict['equation_of_state_fit'] = "http://purls.helmholtz-metadaten.de/asmo/Vinet"
    else:
        method_dict['equation_of_state_fit'] = "http://purls.helmholtz-metadaten.de/asmo/PolynomialFit"

    inputs = []
    if job.input['fit_type'] == 'polynomial':
        inputs.append(
            {
                "label": "fit_order",
                "value": job.input['fit_order'],
            }
        )
        method_dict['@context']['fit_order'] = 'http://purls.helmholtz-metadaten.de/asmo/InputParameter'

    inputs.append(
        {
            "label": "strain_axes",
            "value": ','.join(job.input['axes']),
        }
    )
    inputs.append(
        {
            "label": "number_of_data_points",
            "value": job.input['num_points'],
        }
    )
    inputs.append(
        {
            "label": "volume_range",
            "value": job.input['vol_range'],
        }
    )

    method_dict["inputs"] = inputs
    
    
def extract_murnaghan_calculated_quantities(job, method_dict):
    outputs = []
    outputs.append(
        {
            "label": "equilibrium_bulk_modulus",
            "value": np.round(job.output_to_pandas()['equilibrium_bulk_modulus'][0],4),
            "unit": "GigaPA",
        }
    )
    outputs.append(
        {
            "label": "equilibrium_total_energy",
            "value": np.round(job.output_to_pandas()['equilibrium_energy'][0],4),
            "unit": "EV",
        }
    )
    outputs.append(
        {
            "label": "equilibrium_volume",
            "value": np.round(job.output_to_pandas()['equilibrium_volume'][0],4),
            "unit": "ANGSTROM3",
        }
    )

    method_dict["outputs"] = outputs

def get_unit_cell_parameters(structure):
    if structure.get_symmetry().spacegroup['InternationalTableSymbol'] == "Im-3m":
        if structure.get_number_of_atoms() == 1:
            structure_parameters = {
                'a': np.round(structure.cell[1][0]*2, 4),
                'alpha': 90.0,
                'beta': 90.0,
                'gamma': 90.0,
                'volume': np.round((structure.cell[1][0]*2)**3, 4),
                'space_group': structure.get_symmetry().spacegroup['InternationalTableSymbol'],
                'space_group_number': structure.get_symmetry().spacegroup['Number'],
                'bravais_lattice': 'bcc'
            }
        else:
            structure_parameters = {
                'a': np.round(structure.cell[1][1], 4),
                'alpha': 90.0,
                'beta': 90.0,
                'gamma': 90.0,
                'volume': np.round(structure.get_volume(),4),
                'space_group': structure.get_symmetry().spacegroup['InternationalTableSymbol'],
                'space_group_number': structure.get_symmetry().spacegroup['Number'],
                'bravais_lattice': 'bcc'
            }
    elif structure.get_symmetry().spacegroup['InternationalTableSymbol'] == "Fm-3m":
        if structure.get_number_of_atoms() == 1:
            structure_parameters = {
                'a': np.round(structure.cell[1][0]*2, 4),
                'alpha': 90.0,
                'beta': 90.0,
                'gamma': 90.0,
                'volume': np.round((structure.cell[1][0]*2)**3, 4),
                'space_group': structure.get_symmetry().spacegroup['InternationalTableSymbol'],
                'space_group_number': structure.get_symmetry().spacegroup['Number'],
                'bravais_lattice': 'fcc'
            }
        else:
            structure_parameters = {
                'a': np.round(structure.cell[1][1], 4),
                'alpha': 90.0,
                'beta': 90.0,
                'gamma': 90.0,
                'volume': np.round(structure.get_volume(),4),
                'space_group': structure.get_symmetry().spacegroup['InternationalTableSymbol'],
                'space_group_number': structure.get_symmetry().spacegroup['Number'],
                'bravais_lattice': 'fcc'
            }
    elif structure.get_symmetry().spacegroup['InternationalTableSymbol'] == "P6_3/mmc":
        if structure.get_number_of_atoms() == 2:
            structure_parameters = {
                'a': np.round(structure.cell[0][0], 4),
                'c': np.round(structure.cell[2][2], 4),
                'alpha': 90.0,
                'beta': 90.0,
                'gamma': 120.0,
                'volume': np.round(structure.get_volume()*3,4),
                'space_group': structure.get_symmetry().spacegroup['InternationalTableSymbol'],
                'space_group_number': structure.get_symmetry().spacegroup['Number'],
                'bravais_lattice': 'hcp'
            }
        else:
            structure_parameters = {
                'a': np.round(structure.cell[0][0], 4),
                'c': np.round(structure.cell[2][2], 4),
                'alpha': 90.0,
                'beta': 90.0,
                'gamma': 120.0,
                'volume': np.round(structure.get_volume(),4),
                'space_group': structure.get_symmetry().spacegroup['InternationalTableSymbol'],
                'space_group_number': structure.get_symmetry().spacegroup['Number'],
                'bravais_lattice': 'hcp'
            }

    return structure_parameters

def add_structure_contexts(sample_dict):
    sample_dict['@context'] = {}
    sample_dict['@context']['path'] = 'http://purls.helmholtz-metadaten.de/cmso/hasPath'
    sample_dict['@context']['unit_cell'] = 'http://purls.helmholtz-metadaten.de/cmso/UnitCell'
    sample_dict['@context']['atoms'] = 'http://purls.helmholtz-metadaten.de/cmso/Atom'
    sample_dict['@context']['molecules'] = 'http://purls.helmholtz-metadaten.de/cmso/Molecule'
    sample_dict['@context']['bravais_lattice'] = 'http://purls.helmholtz-metadaten.de/cmso/hasBravaisLattice'
    sample_dict['@context']['chemical_species'] = 'http://purls.helmholtz-metadaten.de/cmso/ChemicalSpecies'
    sample_dict['@context']['simulation_cell'] = 'http://www.w3.org/2000/01/rdf-schema#label'
    sample_dict['@context']['label'] = 'http://www.w3.org/2000/01/rdf-schema#label'
    sample_dict['@context']['unit'] = 'http://purls.helmholtz-metadaten.de/cmso/hasUnit'
    sample_dict['@context']['value'] = 'http://purls.helmholtz-metadaten.de/asmo/hasValue'
    sample_dict['@context']['vector'] = 'http://purls.helmholtz-metadaten.de/cmso/Vector'
    sample_dict['@context']['job_details'] = 'id-from-pmdco-pending'
    sample_dict['@context']['workflow_manager'] = 'http://demo.fiz-karlsruhe.de/matwerk/E457491'
    #sample_dict['@context']['software'] = ''
    sample_dict['@context']['lattice_parameter_a'] = 'http://purls.helmholtz-metadaten.de/cmso/hasLatticeParameter'
    sample_dict['@context']['lattice_angle_alpha'] = "http://purls.helmholtz-metadaten.de/cmso/hasAngle"
    sample_dict['@context']['lattice_angle_beta'] = "http://purls.helmholtz-metadaten.de/cmso/hasAngle"
    sample_dict['@context']['lattice_angle_gamma'] = "http://purls.helmholtz-metadaten.de/cmso/hasAngle"
    sample_dict['@context']['lattice_volume'] = 'http://purls.helmholtz-metadaten.de/asmo/Volume'
    sample_dict['@context']['space_group'] = 'http://purls.helmholtz-metadaten.de/cmso/hasSpaceGroup'
    sample_dict['@context']['bravais_lattice'] = 'http://purls.helmholtz-metadaten.de/cmso/hasBravaisLattice'
    sample_dict['@context']['simulation_cell_lengths'] = 'http://purls.helmholtz-metadaten.de/cmso/hasLength'
    sample_dict['@context']['simulation_cell_vectors'] = 'http://purls.helmholtz-metadaten.de/cmso/hasVector'
    sample_dict["@context"]['simulation_cell_volume'] = "http://purls.helmholtz-metadaten.de/cmso/hasVolume"  
    sample_dict["@context"]['simulation_cell_angle'] = "http://purls.helmholtz-metadaten.de/cmso/hasAngle"

def identify_structure_parameters(structure_parameters, sample_dict):
    if not structure_parameters:
        return
    else:
        
        unit_cell_details = []

        cell_parameter_a_dict = {}
        cell_parameter_a_dict["value"] = structure_parameters['a']
        cell_parameter_a_dict["unit"] = "ANGSTROM"
        cell_parameter_a_dict["label"] = 'lattice_parameter_a'
        unit_cell_details.append(cell_parameter_a_dict)
        if 'b' in structure_parameters.keys():
            cell_parameter_b_dict = {}
            cell_parameter_b_dict["value"] = structure_parameters['b']
            cell_parameter_b_dict["unit"] = "ANGSTROM"
            cell_parameter_b_dict["label"] = 'lattice_parameter_b'
            sample_dict['@context']['lattice_parameter_b'] = 'http://purls.helmholtz-metadaten.de/cmso/hasLatticeParameter'

            unit_cell_details.append(cell_parameter_b_dict)
        if 'c' in structure_parameters.keys():
            cell_parameter_c_dict = {}
            cell_parameter_c_dict["value"] = structure_parameters['c']
            cell_parameter_c_dict["unit"] = "ANGSTROM"
            cell_parameter_c_dict["label"] = 'lattice_parameter_c'
            sample_dict['@context']['lattice_parameter_c'] = 'http://purls.helmholtz-metadaten.de/cmso/hasLatticeParameter'
            unit_cell_details.append(cell_parameter_c_dict)
        if 'c_over_a' in structure_parameters.keys():
            cell_parameter_c_over_a_dict = {}
            cell_parameter_c_over_a_dict["value"] = structure_parameters['c_over_a']
            cell_parameter_c_over_a_dict["unit"] = "ANGSTROM"
            cell_parameter_c_over_a_dict["label"] = 'lattice_parameter_c_over_a'
            sample_dict['@context']['lattice_parameter_c_over_a'] = 'http://purls.helmholtz-metadaten.de/cmso/hasLatticeParameter'
            unit_cell_details.append(cell_parameter_c_over_a_dict)
            
        cell_parameter_alpha_dict = {}
        cell_parameter_alpha_dict["value"] = structure_parameters['alpha']
        cell_parameter_alpha_dict["unit"] = "DEGREE"
        cell_parameter_alpha_dict["label"] = 'lattice_angle_alpha'
        unit_cell_details.append(cell_parameter_alpha_dict)
        cell_parameter_beta_dict = {}
        cell_parameter_beta_dict["value"] = structure_parameters['beta']
        cell_parameter_beta_dict["unit"] = "DEGREE"
        cell_parameter_beta_dict["label"] = 'lattice_angle_beta'
        unit_cell_details.append(cell_parameter_beta_dict)
        cell_parameter_gamma_dict = {}
        cell_parameter_gamma_dict["value"] = structure_parameters['gamma']
        cell_parameter_gamma_dict["unit"] = "DEGREE"
        cell_parameter_gamma_dict["label"] = 'lattice_angle_gamma'
        unit_cell_details.append(cell_parameter_gamma_dict)

        cell_parameter_vol_dict = {}
        cell_parameter_vol_dict["value"] = structure_parameters['volume']
        cell_parameter_vol_dict["unit"] = "ANGSTROM3"
        cell_parameter_vol_dict["label"] = 'lattice_volume'
        unit_cell_details.append(cell_parameter_vol_dict)

        cell_parameter_spg_dict = {}
        cell_parameter_spg_dict["value"] = structure_parameters['space_group']
        cell_parameter_spg_dict["label"] = 'space_group'
        unit_cell_details.append(cell_parameter_spg_dict)

        cell_parameter_bsl_dict = {}
        cell_parameter_bsl_dict["value"] = structure_parameters['bravais_lattice']
        cell_parameter_bsl_dict["label"] = 'bravais_lattice'
        unit_cell_details.append(cell_parameter_bsl_dict)

        sample_dict['unit_cell'] = unit_cell_details
            

def get_chemical_species(structure, sample_dict):
    structure = structure
    natoms = structure.get_number_of_atoms()
    species_dict = dict(structure.get_number_species_atoms())
    atoms_list = []
    for k in species_dict.keys():
        element = {}
        element["value"] = species_dict[k]
        element["label"] = k
        atoms_list.append(element)

    atoms_list.append({'value': natoms, 'label': 'total_number_atoms'})
        
    sample_dict['atoms'] = atoms_list

def get_simulation_cell(structure, sample_dict):
    structure = structure
    cell_lengths = str([np.round(structure.cell.cellpar()[0],4),np.round(structure.cell.cellpar()[1],4),np.round(structure.cell.cellpar()[2],4)])
    cell_vectors = str([np.round(structure.cell[0],4),np.round(structure.cell[1],4),np.round(structure.cell[2],4)])
    cell_angles = str([np.round(structure.cell.cellpar()[3],4),np.round(structure.cell.cellpar()[4],4),np.round(structure.cell.cellpar()[5],4)])
    cell_volume = structure.get_volume()
    
    simulation_cell_details = []

    cell_lengths_dict = {}
    cell_lengths_dict["value"] = cell_lengths
    cell_lengths_dict["unit"] = "ANGSTROM"
    cell_lengths_dict["label"] = 'simulation_cell_lengths'
    simulation_cell_details.append(cell_lengths_dict)
    
    cell_vector_dict = {}
    cell_vector_dict["value"] = cell_vectors
    cell_vector_dict["unit"] = "ANGSTROM"
    cell_vector_dict["label"] = 'simulation_cell_vectors'
    simulation_cell_details.append(cell_vector_dict)

    cell_angles_dict = {}
    cell_angles_dict["value"] = cell_angles
    cell_angles_dict["unit"] = "DEGREES"
    cell_angles_dict["label"] = 'simulation_cell_angles'
    simulation_cell_details.append(cell_angles_dict)

    cell_volume_dict = {}
    cell_volume_dict["value"] = np.round(cell_volume, decimals=4)
    cell_volume_dict["unit"] = "ANGSTROM3"
    cell_volume_dict["label"] = 'simulation_cell_volume'
    
    simulation_cell_details.append(cell_volume_dict)

    sample_dict['simulation_cell'] = simulation_cell_details

def add_structure_software(pr, structure_name, sample_dict):
    sample_dict["workflow_manager"] = {}
    import subprocess
    import platform
    try:
        if "Windows" in platform.system():
            output1 = subprocess.check_output(['findstr', 'pyiron_atomistics', pr.path + '\\' + pr.name + '_environment.yml'])
        else:
            output1 = subprocess.check_output(['grep', 'pyiron_atomistics', pr.path + pr.name + '_environment.yml'])
        s1 = str((output1.decode('utf-8')))
    except:
        s1 = ''
    try:
        if "Windows" in platform.system():
            output2 = subprocess.check_output(['findstr', 'pyiron_workflow', pr.path + '\\' + pr.name + '_environment.yml'])
        else:
            output2 = subprocess.check_output(['grep', 'pyiron_workflow', pr.path + pr.name + '_environment.yml'])
        s2 = str((output2.decode('utf-8')))
    except:
        s2 = ''
    try:
        if "Windows" in platform.system():
            output3 = subprocess.check_output(['findstr', 'pyironflow', pr.path + '\\' + pr.name + '_environment.yml'])
        else:
            output3 = subprocess.check_output(['grep', 'pyironflow', pr.path + pr.name + '_environment.yml'])
        s3 = str((output3.decode('utf-8')))
    except:
        s3 = ''
    try:
        if "Windows" in platform.system():
            output4 = subprocess.check_output(['findstr', 'executorlib', pr.path + '\\' + pr.name + '_environment.yml'])
        else:
            output4 = subprocess.check_output(['grep', 'executorlib', pr.path + pr.name + '_environment.yml'])
        s4 = str((output4.decode('utf-8')))
    except:
        s4 = ''

    #hdf_ver = job.to_dict()['HDF_VERSION']
    try:
        st1 = 'p' + s1.split('=')[0].split('p')[1] + "=" + s1.split('=')[1] + ', '
    except:
        st1 = ''
    try:
        st2 = 'p' + s2.split('=')[0].split('p')[1] + "=" + s2.split('=')[1] + ', '
    except:
        st2 = ''
    try:
        st3 = 'p' + s3.split('=')[0].split('p')[1] + "=" + s3.split('=')[1] + ', '
    except:
        st3 = ''
    try:
        st4 = 'e' + s3.split('=')[0].split('e')[1] + "=" + s4.split('=')[1] + ', '
    except:
        st4 = ''
    st = st1 + st2 + st3 + st4
    
    #+ ', pyiron_HDF_version=' + hdf_ver
    sample_dict["workflow_manager"]["label"] = st
    
    pyiron_job_details = []
    pyiron_job_details.append(
        {
            "label": "structure_name",
            "value": structure_name,
        }
    )
    pyiron_job_details.append(
        {
            "label": "project_name",
            "value": pr.name,
        }
    )
    sample_dict["job_details"] = pyiron_job_details

def get_structure_folder(structure_path, sample_dict):
    import subprocess
    import platform
    if "Windows" in platform.system():
        s_clean = os.getcwd().replace('\\', '/')
    else:
        output = subprocess.check_output(['pwd'])
        s = str((output.decode('utf-8')))
        s_clean = s.replace("\n", "")
    sample_dict['path'] = s_clean + '/' + structure_path

def add_vasp_contexts(method_dict):
    # TODO: check, expand on
    method_dict['@context'] = {}
    method_dict['@context']['sample'] = 'http://purls.helmholtz-metadaten.de/cmso/AtomicScaleSample'
    method_dict['@context']['path'] = 'http://purls.helmholtz-metadaten.de/cmso/hasPath'
    method_dict['@context']['dof'] = 'http://purls.helmholtz-metadaten.de/asmo/hasRelaxationDOF'
    method_dict['@context']['inputs'] = 'http://purls.helmholtz-metadaten.de/asmo/hasInputParameter'
    method_dict['@context']['label'] = 'http://www.w3.org/2000/01/rdf-schema#label'
    method_dict['@context']['unit'] = 'http://purls.helmholtz-metadaten.de/asmo/hasUnit'
    method_dict['@context']['value'] = 'http://purls.helmholtz-metadaten.de/asmo/hasValue'
    method_dict['@context']['outputs'] = 'http://purls.helmholtz-metadaten.de/cmso/hasCalculatedProperty'
    method_dict['@context']['workflow_manager'] = 'http://demo.fiz-karlsruhe.de/matwerk/E457491'
    #method_dict['@context']['software'] = '"@id": "https://www.vasp.at/","label": "VASP"'
    method_dict['@context']['dft'] = "http://purls.helmholtz-metadaten.de/asmo/DensityFunctionalTheory"
    method_dict['@context']['xc_functional'] = 'https://w3id.org/mdo/calculation/hasXCFunctional'

def identify_vasp_method(job, method_dict):
    # copy-pasted from pyiron-conceptual-dictionary
    indf = job.input.incar.to_dict()
    params = indf['data_dict']['Parameter']
    vals = indf['data_dict']['Value']
    mlist = []
    for p,v in zip(params, vals):
        mlist.append(p + '=' + v)
    mstring = ';'.join(mlist)
    raw = mstring.split(';')
    mdict = {}
    for r in raw:
        rsplit = r.split('=')
        if len(rsplit) == 2:
            mdict[rsplit[0].replace(' ','')] = rsplit[1].replace(' ','')
    dof = []
    if 'ISIF' in mdict.keys():
        if mdict['ISIF'] in ['0', '1', '2']:
            dof.append('AtomicPositionRelaxation')
        elif mdict['ISIF'] == '3':
            dof.append('AtomicPositionRelaxation')
            dof.append('CellShapeRelaxation')
            dof.append('CellVolumeRelaxation')
        elif mdict['ISIF'] == '4':
            dof.append('AtomicPositionRelaxation')
            dof.append('CellShapeRelaxation')
        elif mdict['ISIF'] == '5':
            dof.append('CellShapeRelaxation')
        elif mdict['ISIF'] == '6':
            dof.append('CellShapeRelaxation')
            dof.append('CellVolumeRelaxation')
        elif mdict['ISIF'] == '7':
            dof.append('CellVolumeRelaxation')
        elif mdict['ISIF'] == '8':
            dof.append('AtomicPositionRelaxation')
            dof.append('CellVolumeRelaxation')
    if 'NSW' in mdict.keys():
        if mdict['NSW'] == '0':
            dof = []

    method_dict['dft'] = {}
    method_dict['dof'] = dof

    method_dict['dft']['inputs'] = []

    encut_dict = {}
    encut_dict['value'] = mdict['ENCUT']
    encut_dict['label'] = 'energy_cutoff'
    encut_dict['unit'] = 'EV'
    method_dict['dft']['inputs'].append(encut_dict)

    indf = job.input.to_dict()['kpoints/data_dict']
    params = indf['Parameter']
    vals = indf['Value']   

    kpoint_type = vals[2]
    kpoint_grid = vals[3]

    kpoint_dict = {}
    kpoint_dict['label'] = f'kpoint_{kpoint_type}'
    kpoint_dict['value'] = kpoint_grid
    method_dict['dft']['inputs'].append(kpoint_dict)

    indf = job.input.to_dict()['potcar/data_dict']
    xc = indf['Value'][0]
    method_dict['xc_functional'] = xc

def extract_vasp_calculated_quantities(job, method_dict):
    # mostly copy-pasted from pyiron_concept_dictionary // lammps
    # TODO how many decimals? lammps does 4, not 5
    """
    Extracts calculated quantities from a job.

    Parameters
    ----------
    job : pyiron.Job
        The job object containing the calculated quantities.

    Returns
    -------
    list
        A list of dictionaries, each containing the label, value, unit, and associate_to_sample of a calculated quantity.

    """
    outputs = []
    outputs.append(
        {
            "label": "average_total_energy",
            "value": np.round(job.output.energy_tot[-1], decimals=5),
            "unit": "EV",
        }
    )
    outputs.append(
        {
            "label": "final_total_volume",
            "value": np.round(job.output.volume[-1], decimals=5),
            "unit": "ANGSTROM3",
        }
    )

def export_env(path):
    """Exports to path+_environment.yml"""
    import os
    import platform
    if "Windows" in platform.system():
        os.system('conda env export | findstr -v "^prefix: " > ' + path + '_environment.yml')
    else:
        os.system('conda env export | grep -v "^prefix: " > ' + path + '_environment.yml')