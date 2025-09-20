"""
LAMMPS specific functions for parsing

Use this a reference for specific implementations

TODO:
- add vasp module
- Instead of job.to_dict() use pyiron_dataclass - convert from hdf5 to dataclass - convert to conc. dict.



"""

import ast
import json
import os
import warnings
from typing import Optional

import numpy as np
from ase import units


def process_general_job(job):
    method_dict = {}
    _add_simulation_software(job, method_dict)
    _get_simulation_folder(job, method_dict)
    file_name = job.path + "_concept_dict.json"
    with open(file_name, "w") as f:
        json.dump(method_dict, f, indent=2)
    return method_dict


def process_lammps_job(job):
    method_dict = {}
    _add_lammps_contexts(method_dict)
    # get_structures(job, method_dict)
    _identify_lammps_method(job, method_dict)
    _extract_lammps_calculated_quantities(job, method_dict)
    _add_simulation_software(job, method_dict)
    _get_simulation_folder(job, method_dict)
    file_name = job.path + "_concept_dict.json"
    with open(file_name, "w") as f:
        json.dump(method_dict, f, indent=2)
    return method_dict


def process_structure_crystal(
    pr,
    structure,
    structure_name,
    structure_path,
    structure_parameters: dict = None,
    options=None,
):
    sample_dict = {}
    _add_structure_contexts(sample_dict)
    _get_chemical_species(structure, sample_dict)
    _identify_structure_parameters(structure_parameters, sample_dict)
    _get_simulation_cell(structure, sample_dict)
    _add_structure_software(pr, structure_name, sample_dict)
    sample_dict["path"] = structure_path
    if options.get("defects"):
        sample_dict["defects"] = options["defects"]
    if options.get("comments"):
        sample_dict["comments"] = options["comments"]
    json_file_name = structure_path + structure_name + "_concept_dict.json"
    with open(json_file_name, "w") as f:
        json.dump(sample_dict, f, indent=2)
    return sample_dict


def process_murnaghan_job(job):
    method_dict = {}
    _add_murnaghan_contexts(method_dict)
    _identify_murnaghan_method(job, method_dict)
    _extract_murnaghan_calculated_quantities(job, method_dict)
    _add_simulation_software(job, method_dict)
    _get_simulation_folder(job, method_dict)
    file_name = job.path + "_concept_dict.json"
    with open(file_name, "w") as f:
        json.dump(method_dict, f, indent=2)
    return method_dict


def process_vasp_job(job):
    method_dict = {}
    _add_vasp_contexts(method_dict)
    _identify_vasp_method(job, method_dict)
    _extract_vasp_calculated_quantities(job, method_dict)
    _add_simulation_software(job, method_dict)
    _get_simulation_folder(job, method_dict)
    file_name = job.path + "_concept_dict.json"
    with open(file_name, "w") as f:
        json.dump(method_dict, f, indent=2)
    return method_dict


def _add_lammps_contexts(method_dict):
    method_dict["@context"] = {}
    method_dict["@context"][
        "sample"
    ] = "http://purls.helmholtz-metadaten.de/cmso/AtomicScaleSample"
    method_dict["@context"]["path"] = "http://purls.helmholtz-metadaten.de/cmso/hasPath"
    method_dict["@context"][
        "dof"
    ] = "http://purls.helmholtz-metadaten.de/asmo/hasRelaxationDOF"
    method_dict["@context"][
        "inputs"
    ] = "http://purls.helmholtz-metadaten.de/asmo/hasInputParameter"
    method_dict["@context"]["label"] = "http://www.w3.org/2000/01/rdf-schema#label"
    method_dict["@context"]["unit"] = "http://purls.helmholtz-metadaten.de/asmo/hasUnit"
    method_dict["@context"][
        "value"
    ] = "http://purls.helmholtz-metadaten.de/asmo/hasValue"
    method_dict["@context"][
        "outputs"
    ] = "http://purls.helmholtz-metadaten.de/cmso/hasCalculatedProperty"
    method_dict["@context"][
        "workflow_manager"
    ] = "http://demo.fiz-karlsruhe.de/matwerk/E457491"
    # method_dict['@context']['software'] = ''
    method_dict["@context"][
        "molecular_dynamics"
    ] = "http://purls.helmholtz-metadaten.de/asmo/MolecularDynamics"
    method_dict["@context"][
        "molecular_statics"
    ] = "http://purls.helmholtz-metadaten.de/asmo/MolecularStatics"
    method_dict["@context"][
        "ensemble"
    ] = "http://purls.helmholtz-metadaten.de/asmo/hasStatisticalEnsemble"
    method_dict["@context"]["job_details"] = "http://id-from-pmdco-pending"
    method_dict["@context"][
        "periodic_boundary_condition"
    ] = "http://purls.helmholtz-metadaten.de/asmo/PeriodicBoundaryCondition"
    method_dict["@context"][
        "initial_temperature"
    ] = "http://purls.helmholtz-metadaten.de/asmo/Temperature"
    method_dict["@context"][
        "initial_pressure"
    ] = "http://purls.helmholtz-metadaten.de/asmo/Pressure"
    method_dict["@context"][
        "target_temperature"
    ] = "http://purls.helmholtz-metadaten.de/asmo/Temperature"
    method_dict["@context"][
        "target_pressure"
    ] = "http://purls.helmholtz-metadaten.de/asmo/Pressure"
    method_dict["@context"][
        "ionic_energy_tolerance"
    ] = "http://purls.helmholtz-metadaten.de/asmo/InputParameter"
    method_dict["@context"][
        "force_tolerace"
    ] = "http://purls.helmholtz-metadaten.de/asmo/InputParameter"
    method_dict["@context"][
        "maximum_iterations"
    ] = "http://purls.helmholtz-metadaten.de/asmo/InputParameter"
    method_dict["@context"][
        "potential"
    ] = "http://purls.helmholtz-metadaten.de/asmo/InteratomicPotential"
    method_dict["@context"][
        "average_temperature"
    ] = "http://purls.helmholtz-metadaten.de/asmo/Temperature"
    method_dict["@context"][
        "average_pressure"
    ] = "http://purls.helmholtz-metadaten.de/asmo/Pressure"
    method_dict["@context"][
        "average_total_energy"
    ] = "http://purls.helmholtz-metadaten.de/asmo/TotalEnergy"
    method_dict["@context"][
        "final_total_energy"
    ] = "http://purls.helmholtz-metadaten.de/asmo/TotalEnergy"
    method_dict["@context"][
        "final_potential_energy"
    ] = "http://purls.helmholtz-metadaten.de/asmo/PotentialEnergy"
    method_dict["@context"][
        "average_total_volume"
    ] = "http://purls.helmholtz-metadaten.de/asmo/Volume"
    method_dict["@context"][
        "final_total_volume"
    ] = "http://purls.helmholtz-metadaten.de/asmo/Volume"
    method_dict["@context"][
        "final_maximum_force"
    ] = "http://purls.helmholtz-metadaten.de/asmo/Force"
    method_dict["@context"][
        "number_ionic_steps"
    ] = "http://purls.helmholtz-metadaten.de/asmo/NumberOfIonicSteps"
    method_dict["@context"][
        "time_step"
    ] = "http://purls.helmholtz-metadaten.de/asmo/TimeStep"
    method_dict["@context"][
        "simulation_time"
    ] = "http://purls.helmholtz-metadaten.de/asmo/Time"
    method_dict["@context"]["LAMMPS"] = "http://demo.fiz-karlsruhe.de/matwerk/E447986"


def get_structures(job, method_dict):
    initial_pyiron_structure = job.structure
    final_pyiron_structure = job.get_structure(frame=-1)

    method_dict["sample"] = {
        "initial": initial_pyiron_structure,
        "final": final_pyiron_structure,
    }


def _identify_lammps_method(job, method_dict):
    job_dict = job.input.to_dict()
    input_dict = {
        job_dict["control_inp/data_dict"]["Parameter"][x]: job_dict[
            "control_inp/data_dict"
        ]["Value"][x]
        for x in range(len(job_dict["control_inp/data_dict"]["Parameter"]))
    }
    dof = []
    maxiter = None
    e_tol = None
    f_tol = None
    temp = None
    press = None
    md_method = None
    ensemble = None
    temp_target = None
    press_target = None
    timestep = None
    simulation_time = None

    if "min_style" in input_dict.keys():
        dof.append("http://purls.helmholtz-metadaten.de/asmo/AtomicPositionRelaxation")
        if job.input.control["fix___ensemble"] != "all nve":
            dof.append("http://purls.helmholtz-metadaten.de/asmo/CellVolumeRelaxation")
        raw = input_dict["fix___ensemble"].split()
        md_method = "molecular_statics"
        e_tol = float(input_dict["minimize"].split()[0])
        f_tol = float(input_dict["minimize"].split()[1])
        maxiter = int(input_dict["minimize"].split()[2])
        if job.input.control["fix___ensemble"] != "all nve":
            press_target = float(raw[3]) * 0.0001

    elif "nve" in input_dict["fix___ensemble"]:
        if int(input_dict["run"]) == 0:
            md_method = "molecular_statics"
            ensemble = "http://purls.helmholtz-metadaten.de/asmo/MicrocanonicalEnsemble"

        elif int(input_dict["run"]) > 0:
            dof.append(
                "http://purls.helmholtz-metadaten.de/asmo/AtomicPositionRelaxation"
            )
            md_method = "molecular_dynamics"
            ensemble = "http://purls.helmholtz-metadaten.de/asmo/MicrocanonicalEnsemble"
            timestep = float(input_dict["timestep"].split()[0])
            simulation_time = float(input_dict["run"].split()[0]) * float(
                input_dict["timestep"].split()[0]
            )

    elif "nvt" in input_dict["fix___ensemble"]:
        raw = input_dict["fix___ensemble"].split()
        temp = float(raw[3])
        temp_target = float(raw[4])
        dof.append("http://purls.helmholtz-metadaten.de/asmo/AtomicPositionRelaxation")
        md_method = "molecular_dynamics"
        ensemble = "http://purls.helmholtz-metadaten.de/asmo/CanonicalEnsemble"
        timestep = float(input_dict["timestep"].split()[0])
        simulation_time = float(input_dict["run"].split()[0]) * float(
            input_dict["timestep"].split()[0]
        )

    elif "npt" in input_dict["fix___ensemble"]:
        dof.append("http://purls.helmholtz-metadaten.de/asmo/AtomicPositionRelaxation")
        dof.append("http://purls.helmholtz-metadaten.de/asmo/CellVolumeRelaxation")
        if "aniso" in input_dict["fix___ensemble"]:
            dof.append("http://purls.helmholtz-metadaten.de/asmo/CellShapeRelaxation")
        md_method = "molecular_dynamics"
        raw = input_dict["fix___ensemble"].split()
        temp = float(raw[3])
        temp_target = float(raw[4])
        press = float(raw[7]) * 0.0001
        press_target = float(raw[8]) * 0.0001
        ensemble = "http://purls.helmholtz-metadaten.de/asmo/IsothermalIsobaricEnsemble"
        timestep = float(input_dict["timestep"].split()[0])
        simulation_time = float(input_dict["run"].split()[0]) * float(
            input_dict["timestep"].split()[0]
        )

    method_dict[md_method] = {}

    if md_method == "molecular_statics":
        method_dict[md_method]["minimization_algorithm"] = input_dict["min_style"]

    input_dict = {
        job_dict["control_inp/data_dict"]["Parameter"][x]: job_dict[
            "control_inp/data_dict"
        ]["Value"][x]
        for x in range(len(job_dict["control_inp/data_dict"]["Parameter"]))
    }
    pb = []
    pb.append(input_dict["boundary"])
    if pb[0][0] == "p":
        method_dict[md_method]["periodicity_in_x"] = True
    else:
        method_dict[md_method]["periodicity_in_x"] = False
    if pb[0][2] == "p":
        method_dict[md_method]["periodicity_in_y"] = True
    else:
        method_dict[md_method]["periodicity_in_y"] = False
    if pb[0][4] == "p":
        method_dict[md_method]["periodicity_in_z"] = True
    else:
        method_dict[md_method]["periodicity_in_z"] = False

    method_dict[md_method]["inputs"] = []

    temperature = {}
    temperature["value"] = temp
    temperature["unit"] = "K"
    temperature["label"] = "initial_temperature"

    method_dict[md_method]["inputs"].append(temperature)

    target_temperature = {}
    target_temperature["value"] = temp_target
    target_temperature["unit"] = "K"
    target_temperature["label"] = "target_temperature"

    method_dict[md_method]["inputs"].append(target_temperature)

    pressure = {}
    pressure["value"] = press
    pressure["unit"] = "GigaPA"
    pressure["label"] = "initial_pressure"

    method_dict[md_method]["inputs"].append(pressure)

    target_pressure = {}
    target_pressure["value"] = press_target
    target_pressure["unit"] = "GigaPA"
    target_pressure["label"] = "target_pressure"

    method_dict[md_method]["inputs"].append(target_pressure)

    energy_tol = {}
    energy_tol["value"] = e_tol
    energy_tol["unit"] = "EV"
    energy_tol["label"] = "ionic_energy_tolerance"

    method_dict[md_method]["inputs"].append(energy_tol)

    force_tol = {}
    force_tol["value"] = f_tol
    force_tol["unit"] = "EV-PER-ANGSTROM"
    force_tol["label"] = "force_tolerance"

    method_dict[md_method]["inputs"].append(force_tol)

    maximum_iterations = {}
    maximum_iterations["value"] = maxiter
    maximum_iterations["label"] = "maximum_iterations"

    method_dict[md_method]["inputs"].append(maximum_iterations)

    timestep_dict = {}
    timestep_dict["value"] = timestep
    timestep_dict["unit"] = "PICOSECOND"
    timestep_dict["label"] = "timestep"

    method_dict[md_method]["inputs"].append(timestep_dict)

    sim_time_dict = {}
    sim_time_dict["value"] = simulation_time
    sim_time_dict["unit"] = "PICOSECOND"
    sim_time_dict["label"] = "simulation_time"

    method_dict[md_method]["inputs"].append(sim_time_dict)

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
    if "meam" in ps:
        method_dict["@context"][
            "potential"
        ] = "http://purls.helmholtz-metadaten.de/asmo/ModifiedEmbeddedAtomModel"
    elif "eam" in ps:
        method_dict["@context"][
            "potential"
        ] = "http://purls.helmholtz-metadaten.de/asmo/EmbeddedAtomModel"
    elif "lj" in ps:
        method_dict["@context"][
            "potential"
        ] = "http://purls.helmholtz-metadaten.de/asmo/LennardJonesPotential"
    elif "ace" in ps:
        method_dict["@context"][
            "potential"
        ] = "http://purls.helmholtz-metadaten.de/asmo/MachineLearningPotential"

    method_dict[md_method]["potential"] = {}
    method_dict[md_method]["potential"]["label"] = name
    if url is not None:
        method_dict[md_method]["potential"]["url"] = url


def _extract_lammps_calculated_quantities(job, method_dict):
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
    fpe = job.output.energy_pot[-1]
    avol = np.mean(job.output.volume)
    fvol = job.output.volume[-1]
    fmax = job.output.force_max[-1]
    nionic = len(job.output.steps) - 1
    atemp = np.mean(job.output.temperature)
    apress = np.mean(
        np.array(
            [
                (1 / 3 * (tensor[0, 0] + tensor[1, 1] + tensor[2, 2]))
                for tensor in job.output.pressures
            ]
        )
    )
    apot = np.mean(job.output.energy_pot)
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
    if "molecular_statics" in method_dict.keys():
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
    else:
        outputs.append(
            {
                "label": "average_potential_energy",
                "value": np.round(apot, decimals=4),
                "unit": "EV",
            }
        )
        outputs.append(
            {
                "label": "average_temperature",
                "value": np.round(atemp, decimals=4),
                "unit": "K",
            }
        )
        outputs.append(
            {
                "label": "average_pressure",
                "value": np.round(apress, decimals=4),
                "unit": "GigaPA",
            }
        )
    method_dict["outputs"] = outputs


def _add_simulation_software(job, method_dict):
    method_dict["workflow_manager"] = {}
    import platform
    import subprocess

    try:
        if "Windows" in platform.system():
            output1 = subprocess.check_output(
                [
                    "findstr",
                    "pyiron_atomistics",
                    job.path.replace("/", "\\") + "_environment.yml",
                ]
            )
        else:
            output1 = subprocess.check_output(
                ["grep", "pyiron_atomistics", job.path + "_environment.yml"]
            )
        s1 = str((output1.decode("utf-8")))
        st1 = "p" + s1.split("=")[0].split("p")[1] + "=" + s1.split("=")[1] + ", "
    except:
        st1 = ""
    try:
        if "Windows" in platform.system():
            output2 = subprocess.check_output(
                [
                    "findstr",
                    "pyiron_workflow",
                    job.path.replace("/", "\\") + "_environment.yml",
                ]
            )
        else:
            output2 = subprocess.check_output(
                ["grep", "pyiron_workflow", job.path + "_environment.yml"]
            )
        s2 = str((output2.decode("utf-8")))
        st2 = "p" + s2.split("=")[0].split("p")[1] + "=" + s2.split("=")[1] + ", "
    except:
        st2 = ""
    try:
        if "Windows" in platform.system():
            output3 = subprocess.check_output(
                [
                    "findstr",
                    "pyironflow",
                    job.path.replace("/", "\\") + "_environment.yml",
                ]
            )
        else:
            output3 = subprocess.check_output(
                ["grep", "pyironflow", job.path + "_environment.yml"]
            )
        s3 = str((output3.decode("utf-8")))
        st3 = "p" + s3.split("=")[0].split("p")[1] + "=" + s3.split("=")[1] + ", "
    except:
        st3 = ""
    try:
        if "Windows" in platform.system():
            output4 = subprocess.check_output(
                [
                    "findstr",
                    "executorlib",
                    job.path.replace("/", "\\") + "_environment.yml",
                ]
            )
        else:
            output4 = subprocess.check_output(
                ["grep", "executorlib", job.path + "_environment.yml"]
            )
        s4 = str((output4.decode("utf-8")))
        st4 = s4.split("=")[0].split("- ")[1] + "=" + s4.split("=")[1]
    except:
        st4 = ""

    # hdf_ver = job.to_dict()['HDF_VERSION']
    st = st1 + st2 + st3 + st4

    # + ', pyiron_HDF_version=' + hdf_ver
    if st:
        method_dict["workflow_manager"]["label"] = st
    else:
        method_dict["workflow_manager"]["label"] = "pyiron"

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
    try:
        pyiron_job_details.append(
            {
                "label": "job_stoptime",
                "value": str(job.database_entry.timestop.strftime("%Y-%m-%d %H:%M:%S")),
            }
        )
    except (
        AttributeError
    ):  # if job.database_entry.timestop = None - TODO better error handling
        pass
    try:
        pyiron_job_details.append(
            {
                "label": "sim_coretime_hours",
                "value": np.round(job.database_entry.totalcputime / 3600, 6),
            }
        )
    except (
        TypeError
    ):  # if job.database_entry.totalcputime = None - TODO better error handling
        pass

    server = job.to_dict()["server"]
    pyiron_job_details.append(
        {
            "label": "number_cores",
            "value": server["cores"],
        }
    )
    pyiron_job_details.append(
        {
            "label": "host",
            "value": server["host"],
        }
    )
    if server[
        "queue"
    ]:  # TODO: testif we could skip conditionals and not get 'None' on openBIS
        pyiron_job_details.append({"label": "queue", "value": server["queue"]})
    if server["qid"]:
        pyiron_job_details.append({"label": "queue id", "value": str(server["qid"])})

    try:
        software = {
            "label": job.to_dict()["executable"]["name"].upper()
            + " "
            + job.to_dict()["executable"]["version"],
        }
        method_dict["software"] = [software]
    except (
        TypeError
    ):  # if version None - will throw an error when uploading but we have some info
        software = {
            "label": job.to_dict()["executable"]["name"].upper(),
        }
        # software['label'] += ' 5.4.4' #TODO remove !!!! workaround only!!!!
        method_dict["software"] = [software]
    except KeyError:
        pass

    method_dict["job_details"] = pyiron_job_details


def _get_simulation_folder(job, method_dict):
    method_dict["path"] = job.path


def _add_murnaghan_contexts(method_dict):
    method_dict["@context"] = {}
    method_dict["@context"][
        "sample"
    ] = "http://purls.helmholtz-metadaten.de/cmso/AtomicScaleSample"
    method_dict["@context"]["path"] = "http://purls.helmholtz-metadaten.de/cmso/hasPath"
    method_dict["@context"][
        "equation_of_state_fit"
    ] = "http://purls.helmholtz-metadaten.de/asmo/EquationOfStateFit"
    method_dict["@context"][
        "inputs"
    ] = "http://purls.helmholtz-metadaten.de/asmo/hasInputParameter"
    method_dict["@context"]["label"] = "http://www.w3.org/2000/01/rdf-schema#label"
    method_dict["@context"]["unit"] = "http://purls.helmholtz-metadaten.de/asmo/hasUnit"
    method_dict["@context"][
        "value"
    ] = "http://purls.helmholtz-metadaten.de/asmo/hasValue"
    method_dict["@context"][
        "outputs"
    ] = "http://purls.helmholtz-metadaten.de/cmso/hasCalculatedProperty"
    method_dict["@context"][
        "workflow_manager"
    ] = "http://demo.fiz-karlsruhe.de/matwerk/E457491"
    # method_dict['@context']['software'] = ''
    method_dict["@context"]["job_details"] = "http://id-from-pmdco-pending"
    method_dict["@context"][
        "strain_axes"
    ] = "http://purls.helmholtz-metadaten.de/asmo/InputParameter"
    method_dict["@context"][
        "number_of_data_points"
    ] = "http://purls.helmholtz-metadaten.de/asmo/InputParameter"
    method_dict["@context"][
        "volume_range"
    ] = "http://purls.helmholtz-metadaten.de/asmo/VolumeRange"
    method_dict["@context"][
        "equilibrium_bulk_modulus"
    ] = "http://purls.helmholtz-metadaten.de/asmo/BulkModulus"
    method_dict["@context"][
        "equilibrium_total_energy"
    ] = "http://purls.helmholtz-metadaten.de/asmo/TotalEnergy"
    method_dict["@context"][
        "equilibrium_volume"
    ] = "http://purls.helmholtz-metadaten.de/asmo/Volume"


def _identify_murnaghan_method(job, method_dict):

    if job.input["fit_type"] == "birchmurnaghan":
        method_dict["equation_of_state_fit"] = (
            "http://purls.helmholtz-metadaten.de/asmo/BirchMurnaghan"
        )
    elif job.input["fit_type"] == "murnaghan":
        method_dict["equation_of_state_fit"] = (
            "http://purls.helmholtz-metadaten.de/asmo/Murnaghan"
        )
    elif job.input["fit_type"] == "vinet":
        method_dict["equation_of_state_fit"] = (
            "http://purls.helmholtz-metadaten.de/asmo/Vinet"
        )
    else:
        method_dict["equation_of_state_fit"] = (
            "http://purls.helmholtz-metadaten.de/asmo/PolynomialFit"
        )

    inputs = []
    if job.input["fit_type"] == "polynomial":
        inputs.append(
            {
                "label": "fit_order",
                "value": job.input["fit_order"],
            }
        )
        method_dict["@context"][
            "fit_order"
        ] = "http://purls.helmholtz-metadaten.de/asmo/InputParameter"

    inputs.append(
        {
            "label": "strain_axes",
            "value": ",".join(job.input["axes"]),
        }
    )
    inputs.append(
        {
            "label": "number_of_data_points",
            "value": job.input["num_points"],
        }
    )
    inputs.append(
        {
            "label": "volume_range",
            "value": job.input["vol_range"],
        }
    )

    method_dict["inputs"] = inputs


def _extract_murnaghan_calculated_quantities(job, method_dict):
    outputs = []
    outputs.append(
        {
            "label": "equilibrium_bulk_modulus",
            "value": np.round(job.output_to_pandas()["equilibrium_bulk_modulus"][0], 4),
            "unit": "GigaPA",
        }
    )
    outputs.append(
        {
            "label": "equilibrium_total_energy",
            "value": np.round(job.output_to_pandas()["equilibrium_energy"][0], 4),
            "unit": "EV",
        }
    )
    outputs.append(
        {
            "label": "equilibrium_volume",
            "value": np.round(job.output_to_pandas()["equilibrium_volume"][0], 4),
            "unit": "ANGSTROM3",
        }
    )

    method_dict["outputs"] = outputs


def get_unit_cell_parameters(structure):
    if structure.get_symmetry().spacegroup["InternationalTableSymbol"] == "Im-3m":
        if structure.get_number_of_atoms() == 1:
            structure_parameters = {
                "a": np.round(structure.cell[1][0] * 2, 4),
                "alpha": 90.0,
                "beta": 90.0,
                "gamma": 90.0,
                "volume": np.round((structure.cell[1][0] * 2) ** 3, 4),
                "space_group": structure.get_symmetry().spacegroup[
                    "InternationalTableSymbol"
                ],
                "space_group_number": structure.get_symmetry().spacegroup["Number"],
                "bravais_lattice": "bcc",
            }
        else:
            structure_parameters = {
                "a": np.round(structure.cell[1][1], 4),
                "alpha": 90.0,
                "beta": 90.0,
                "gamma": 90.0,
                "volume": np.round(structure.get_volume(), 4),
                "space_group": structure.get_symmetry().spacegroup[
                    "InternationalTableSymbol"
                ],
                "space_group_number": structure.get_symmetry().spacegroup["Number"],
                "bravais_lattice": "bcc",
            }
    elif structure.get_symmetry().spacegroup["InternationalTableSymbol"] == "Fm-3m":
        if structure.get_number_of_atoms() == 1:
            structure_parameters = {
                "a": np.round(structure.cell[1][0] * 2, 4),
                "alpha": 90.0,
                "beta": 90.0,
                "gamma": 90.0,
                "volume": np.round((structure.cell[1][0] * 2) ** 3, 4),
                "space_group": structure.get_symmetry().spacegroup[
                    "InternationalTableSymbol"
                ],
                "space_group_number": structure.get_symmetry().spacegroup["Number"],
                "bravais_lattice": "fcc",
            }
        else:
            structure_parameters = {
                "a": np.round(structure.cell[1][1], 4),
                "alpha": 90.0,
                "beta": 90.0,
                "gamma": 90.0,
                "volume": np.round(structure.get_volume(), 4),
                "space_group": structure.get_symmetry().spacegroup[
                    "InternationalTableSymbol"
                ],
                "space_group_number": structure.get_symmetry().spacegroup["Number"],
                "bravais_lattice": "fcc",
            }
    elif structure.get_symmetry().spacegroup["InternationalTableSymbol"] == "P6_3/mmc":
        if structure.get_number_of_atoms() == 2:
            structure_parameters = {
                "a": np.round(structure.cell[0][0], 4),
                "c": np.round(structure.cell[2][2], 4),
                "alpha": 90.0,
                "beta": 90.0,
                "gamma": 120.0,
                "volume": np.round(structure.get_volume() * 3, 4),
                "space_group": structure.get_symmetry().spacegroup[
                    "InternationalTableSymbol"
                ],
                "space_group_number": structure.get_symmetry().spacegroup["Number"],
                "bravais_lattice": "hcp",
            }
        else:
            structure_parameters = {
                "a": np.round(structure.cell[0][0], 4),
                "c": np.round(structure.cell[2][2], 4),
                "alpha": 90.0,
                "beta": 90.0,
                "gamma": 120.0,
                "volume": np.round(structure.get_volume(), 4),
                "space_group": structure.get_symmetry().spacegroup[
                    "InternationalTableSymbol"
                ],
                "space_group_number": structure.get_symmetry().spacegroup["Number"],
                "bravais_lattice": "hcp",
            }

    return structure_parameters


def _add_structure_contexts(sample_dict):
    sample_dict["@context"] = {}
    sample_dict["@context"]["path"] = "http://purls.helmholtz-metadaten.de/cmso/hasPath"
    sample_dict["@context"][
        "unit_cell"
    ] = "http://purls.helmholtz-metadaten.de/cmso/UnitCell"
    sample_dict["@context"]["atoms"] = "http://purls.helmholtz-metadaten.de/cmso/Atom"
    sample_dict["@context"][
        "molecules"
    ] = "http://purls.helmholtz-metadaten.de/cmso/Molecule"
    sample_dict["@context"][
        "bravais_lattice"
    ] = "http://purls.helmholtz-metadaten.de/cmso/hasBravaisLattice"
    sample_dict["@context"][
        "chemical_species"
    ] = "http://purls.helmholtz-metadaten.de/cmso/ChemicalSpecies"
    sample_dict["@context"][
        "simulation_cell"
    ] = "http://www.w3.org/2000/01/rdf-schema#label"
    sample_dict["@context"]["label"] = "http://www.w3.org/2000/01/rdf-schema#label"
    sample_dict["@context"]["unit"] = "http://purls.helmholtz-metadaten.de/cmso/hasUnit"
    sample_dict["@context"][
        "value"
    ] = "http://purls.helmholtz-metadaten.de/asmo/hasValue"
    sample_dict["@context"][
        "vector"
    ] = "http://purls.helmholtz-metadaten.de/cmso/Vector"
    sample_dict["@context"]["job_details"] = "http://id-from-pmdco-pending"
    sample_dict["@context"][
        "workflow_manager"
    ] = "http://demo.fiz-karlsruhe.de/matwerk/E457491"
    # sample_dict['@context']['software'] = ''
    sample_dict["@context"][
        "lattice_parameter_a"
    ] = "http://purls.helmholtz-metadaten.de/cmso/hasLatticeParameter"
    sample_dict["@context"][
        "lattice_angle_alpha"
    ] = "http://purls.helmholtz-metadaten.de/cmso/hasAngle"
    sample_dict["@context"][
        "lattice_angle_beta"
    ] = "http://purls.helmholtz-metadaten.de/cmso/hasAngle"
    sample_dict["@context"][
        "lattice_angle_gamma"
    ] = "http://purls.helmholtz-metadaten.de/cmso/hasAngle"
    sample_dict["@context"][
        "lattice_volume"
    ] = "http://purls.helmholtz-metadaten.de/asmo/Volume"
    sample_dict["@context"][
        "space_group"
    ] = "http://purls.helmholtz-metadaten.de/cmso/hasSpaceGroup"
    sample_dict["@context"][
        "bravais_lattice"
    ] = "http://purls.helmholtz-metadaten.de/cmso/hasBravaisLattice"
    sample_dict["@context"][
        "simulation_cell_lengths"
    ] = "http://purls.helmholtz-metadaten.de/cmso/hasLength"
    sample_dict["@context"][
        "simulation_cell_vectors"
    ] = "http://purls.helmholtz-metadaten.de/cmso/hasVector"
    sample_dict["@context"][
        "simulation_cell_volume"
    ] = "http://purls.helmholtz-metadaten.de/cmso/hasVolume"
    sample_dict["@context"][
        "simulation_cell_angle"
    ] = "http://purls.helmholtz-metadaten.de/cmso/hasAngle"


def _identify_structure_parameters(structure_parameters, sample_dict):
    if not structure_parameters:
        return
    else:

        unit_cell_details = []

        cell_parameter_a_dict = {}
        cell_parameter_a_dict["value"] = structure_parameters["a"]
        cell_parameter_a_dict["unit"] = "ANGSTROM"
        cell_parameter_a_dict["label"] = "lattice_parameter_a"
        unit_cell_details.append(cell_parameter_a_dict)
        if "b" in structure_parameters.keys():
            cell_parameter_b_dict = {}
            cell_parameter_b_dict["value"] = structure_parameters["b"]
            cell_parameter_b_dict["unit"] = "ANGSTROM"
            cell_parameter_b_dict["label"] = "lattice_parameter_b"
            sample_dict["@context"][
                "lattice_parameter_b"
            ] = "http://purls.helmholtz-metadaten.de/cmso/hasLatticeParameter"

            unit_cell_details.append(cell_parameter_b_dict)
        if "c" in structure_parameters.keys():
            cell_parameter_c_dict = {}
            cell_parameter_c_dict["value"] = structure_parameters["c"]
            cell_parameter_c_dict["unit"] = "ANGSTROM"
            cell_parameter_c_dict["label"] = "lattice_parameter_c"
            sample_dict["@context"][
                "lattice_parameter_c"
            ] = "http://purls.helmholtz-metadaten.de/cmso/hasLatticeParameter"
            unit_cell_details.append(cell_parameter_c_dict)
        if "c_over_a" in structure_parameters.keys():
            cell_parameter_c_over_a_dict = {}
            cell_parameter_c_over_a_dict["value"] = structure_parameters["c_over_a"]
            cell_parameter_c_over_a_dict["unit"] = "ANGSTROM"
            cell_parameter_c_over_a_dict["label"] = "lattice_parameter_c_over_a"
            sample_dict["@context"][
                "lattice_parameter_c_over_a"
            ] = "http://purls.helmholtz-metadaten.de/cmso/hasLatticeParameter"
            unit_cell_details.append(cell_parameter_c_over_a_dict)

        cell_parameter_alpha_dict = {}
        cell_parameter_alpha_dict["value"] = structure_parameters["alpha"]
        cell_parameter_alpha_dict["unit"] = "DEGREE"
        cell_parameter_alpha_dict["label"] = "lattice_angle_alpha"
        unit_cell_details.append(cell_parameter_alpha_dict)
        cell_parameter_beta_dict = {}
        cell_parameter_beta_dict["value"] = structure_parameters["beta"]
        cell_parameter_beta_dict["unit"] = "DEGREE"
        cell_parameter_beta_dict["label"] = "lattice_angle_beta"
        unit_cell_details.append(cell_parameter_beta_dict)
        cell_parameter_gamma_dict = {}
        cell_parameter_gamma_dict["value"] = structure_parameters["gamma"]
        cell_parameter_gamma_dict["unit"] = "DEGREE"
        cell_parameter_gamma_dict["label"] = "lattice_angle_gamma"
        unit_cell_details.append(cell_parameter_gamma_dict)

        cell_parameter_vol_dict = {}
        cell_parameter_vol_dict["value"] = structure_parameters["volume"]
        cell_parameter_vol_dict["unit"] = "ANGSTROM3"
        cell_parameter_vol_dict["label"] = "lattice_volume"
        unit_cell_details.append(cell_parameter_vol_dict)

        cell_parameter_spg_dict = {}
        cell_parameter_spg_dict["value"] = structure_parameters["space_group"]
        cell_parameter_spg_dict["label"] = "space_group"
        unit_cell_details.append(cell_parameter_spg_dict)

        cell_parameter_bsl_dict = {}
        cell_parameter_bsl_dict["value"] = structure_parameters["bravais_lattice"]
        cell_parameter_bsl_dict["label"] = "bravais_lattice"
        unit_cell_details.append(cell_parameter_bsl_dict)

        sample_dict["unit_cell"] = unit_cell_details


def _get_chemical_species(structure, sample_dict):
    structure = structure
    natoms = structure.get_number_of_atoms()
    species_dict = dict(structure.get_number_species_atoms())
    atoms_list = []
    for k in species_dict.keys():
        element = {}
        element["value"] = species_dict[k]
        element["label"] = k
        atoms_list.append(element)

    atoms_list.append({"value": natoms, "label": "total_number_atoms"})

    sample_dict["atoms"] = atoms_list


def _get_simulation_cell(structure, sample_dict):
    structure = structure
    cell_lengths = str(
        [
            np.round(structure.cell.cellpar()[0], 4),
            np.round(structure.cell.cellpar()[1], 4),
            np.round(structure.cell.cellpar()[2], 4),
        ]
    )
    cell_vectors = str(
        [
            np.round(structure.cell[0], 4),
            np.round(structure.cell[1], 4),
            np.round(structure.cell[2], 4),
        ]
    )
    cell_angles = str(
        [
            np.round(structure.cell.cellpar()[3], 4),
            np.round(structure.cell.cellpar()[4], 4),
            np.round(structure.cell.cellpar()[5], 4),
        ]
    )
    cell_volume = structure.get_volume()

    simulation_cell_details = []

    cell_lengths_dict = {}
    cell_lengths_dict["value"] = cell_lengths
    cell_lengths_dict["unit"] = "ANGSTROM"
    cell_lengths_dict["label"] = "simulation_cell_lengths"
    simulation_cell_details.append(cell_lengths_dict)

    cell_vector_dict = {}
    cell_vector_dict["value"] = cell_vectors
    cell_vector_dict["unit"] = "ANGSTROM"
    cell_vector_dict["label"] = "simulation_cell_vectors"
    simulation_cell_details.append(cell_vector_dict)

    cell_angles_dict = {}
    cell_angles_dict["value"] = cell_angles
    cell_angles_dict["unit"] = "DEGREES"
    cell_angles_dict["label"] = "simulation_cell_angles"
    simulation_cell_details.append(cell_angles_dict)

    cell_volume_dict = {}
    cell_volume_dict["value"] = np.round(cell_volume, decimals=4)
    cell_volume_dict["unit"] = "ANGSTROM3"
    cell_volume_dict["label"] = "simulation_cell_volume"

    simulation_cell_details.append(cell_volume_dict)

    sample_dict["simulation_cell"] = simulation_cell_details


def _add_structure_software(pr, structure_name, sample_dict):
    sample_dict["workflow_manager"] = {}
    import platform
    import subprocess

    try:
        if "Windows" in platform.system():
            output1 = subprocess.check_output(
                [
                    "findstr",
                    "pyiron_atomistics",
                    pr.path + "\\" + pr.name + "_environment.yml",
                ]
            )
        else:
            output1 = subprocess.check_output(
                ["grep", "pyiron_atomistics", pr.path + pr.name + "_environment.yml"]
            )
        s1 = str((output1.decode("utf-8")))
    except:
        s1 = ""
    try:
        if "Windows" in platform.system():
            output2 = subprocess.check_output(
                [
                    "findstr",
                    "pyiron_workflow",
                    pr.path + "\\" + pr.name + "_environment.yml",
                ]
            )
        else:
            output2 = subprocess.check_output(
                ["grep", "pyiron_workflow", pr.path + pr.name + "_environment.yml"]
            )
        s2 = str((output2.decode("utf-8")))
    except:
        s2 = ""
    try:
        if "Windows" in platform.system():
            output3 = subprocess.check_output(
                ["findstr", "pyironflow", pr.path + "\\" + pr.name + "_environment.yml"]
            )
        else:
            output3 = subprocess.check_output(
                ["grep", "pyironflow", pr.path + pr.name + "_environment.yml"]
            )
        s3 = str((output3.decode("utf-8")))
    except:
        s3 = ""
    try:
        if "Windows" in platform.system():
            output4 = subprocess.check_output(
                [
                    "findstr",
                    "executorlib",
                    pr.path + "\\" + pr.name + "_environment.yml",
                ]
            )
        else:
            output4 = subprocess.check_output(
                ["grep", "executorlib", pr.path + pr.name + "_environment.yml"]
            )
        s4 = str((output4.decode("utf-8")))
    except:
        s4 = ""

    # hdf_ver = job.to_dict()['HDF_VERSION']
    try:
        st1 = "p" + s1.split("=")[0].split("p")[1] + "=" + s1.split("=")[1] + ", "
    except:
        st1 = ""
    try:
        st2 = "p" + s2.split("=")[0].split("p")[1] + "=" + s2.split("=")[1] + ", "
    except:
        st2 = ""
    try:
        st3 = "p" + s3.split("=")[0].split("p")[1] + "=" + s3.split("=")[1] + ", "
    except:
        st3 = ""
    try:
        st4 = "e" + s3.split("=")[0].split("e")[1] + "=" + s4.split("=")[1] + ", "
    except:
        st4 = ""
    st = st1 + st2 + st3 + st4

    # + ', pyiron_HDF_version=' + hdf_ver
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


def _add_vasp_contexts(method_dict):
    # TODO: check, expand on
    method_dict["@context"] = {}
    method_dict["@context"][
        "sample"
    ] = "http://purls.helmholtz-metadaten.de/cmso/AtomicScaleSample"
    method_dict["@context"]["path"] = "http://purls.helmholtz-metadaten.de/cmso/hasPath"
    method_dict["@context"][
        "dof"
    ] = "http://purls.helmholtz-metadaten.de/asmo/hasRelaxationDOF"
    method_dict["@context"][
        "inputs"
    ] = "http://purls.helmholtz-metadaten.de/asmo/hasInputParameter"
    method_dict["@context"]["label"] = "http://www.w3.org/2000/01/rdf-schema#label"
    method_dict["@context"]["unit"] = "http://purls.helmholtz-metadaten.de/asmo/hasUnit"
    method_dict["@context"][
        "value"
    ] = "http://purls.helmholtz-metadaten.de/asmo/hasValue"
    method_dict["@context"][
        "outputs"
    ] = "http://purls.helmholtz-metadaten.de/cmso/hasCalculatedProperty"
    method_dict["@context"][
        "workflow_manager"
    ] = "http://demo.fiz-karlsruhe.de/matwerk/E457491"
    # method_dict['@context']['software'] = '"@id": "https://www.vasp.at/","label": "VASP"'
    method_dict["@context"][
        "dft"
    ] = "http://purls.helmholtz-metadaten.de/asmo/DensityFunctionalTheory"
    method_dict["@context"][
        "xc_functional"
    ] = "https://w3id.org/mdo/calculation/hasXCFunctional"
    method_dict["@context"][
        "VASP"
    ] = "https://purls.helmholtz-metadaten.de/msekg/E425582"
    method_dict["@context"]["job_details"] = "http://id-from-pmdco-pending"
    method_dict["@context"][
        "periodic_boundary_condition"
    ] = "http://purls.helmholtz-metadaten.de/asmo/PeriodicBoundaryCondition"
    method_dict["@context"][
        "k_point_mesh"
    ] = "http://purls.helmholtz-metadaten.de/asmo/KPointMesh"
    method_dict["@context"][
        "k_point_generation"
    ] = "http://purls.helmholtz-metadaten.de/asmo/InputParameter"
    method_dict["@context"][
        "electronic_smearing"
    ] = "http://purls.helmholtz-metadaten.de/asmo/InputParameter"
    method_dict["@context"][
        "smearing_parameter_sigma"
    ] = "http://purls.helmholtz-metadaten.de/asmo/InputParameter"
    method_dict["@context"][
        "electronic_minimization_algorithm"
    ] = "http://purls.helmholtz-metadaten.de/asmo/InputParameter"
    method_dict["@context"][
        "ionic_minimization_algorithm"
    ] = "http://purls.helmholtz-metadaten.de/asmo/InputParameter"
    method_dict["@context"][
        "spin_polarization"
    ] = "http://purls.helmholtz-metadaten.de/asmo/InputParameter"
    method_dict["@context"][
        "input_magnetic_moments"
    ] = "http://purls.helmholtz-metadaten.de/asmo/InputParameter"
    method_dict["@context"][
        "electronic_energy_tolerance"
    ] = "http://purls.helmholtz-metadaten.de/asmo/InputParameter"
    method_dict["@context"][
        "ionic_energy_tolerance"
    ] = "http://purls.helmholtz-metadaten.de/asmo/InputParameter"
    method_dict["@context"][
        "force_tolerance"
    ] = "http://purls.helmholtz-metadaten.de/asmo/InputParameter"
    method_dict["@context"][
        "final_total_energy"
    ] = "http://purls.helmholtz-metadaten.de/asmo/TotalEnergy"
    method_dict["@context"][
        "final_potential_energy"
    ] = "http://purls.helmholtz-metadaten.de/asmo/PotentialEnergy"
    method_dict["@context"][
        "final_total_volume"
    ] = "http://purls.helmholtz-metadaten.de/asmo/Volume"
    method_dict["@context"][
        "final_maximum_force"
    ] = "http://purls.helmholtz-metadaten.de/asmo/Force"
    method_dict["@context"][
        "final_total_magnetic_moment"
    ] = "http://purls.helmholtz-metadaten.de/asmo/TotalMagneticMoment"
    method_dict["@context"][
        "number_ionic_steps"
    ] = "http://purls.helmholtz-metadaten.de/asmo/NumberOfIonicSteps"


def _identify_vasp_method(job, method_dict):
    # copy-pasted from pyiron-conceptual-dictionary
    indf = job.input.incar.to_dict()
    params = indf["data_dict"]["Parameter"]
    vals = indf["data_dict"]["Value"]
    mdict = {
        str(p).replace(" ", ""): str(v).replace(" ", "") for p, v in zip(params, vals)
    }

    method_dict["dft"] = {}
    method_dict["dft"]["inputs"] = []

    encut_dict = {}
    encut_dict["value"] = float(
        mdict.get("ENCUT", 0)
    )  # TODO: try to read OUTCAR if default (None now)
    encut_dict["label"] = "energy_cutoff"
    encut_dict["unit"] = "EV"
    method_dict["dft"]["inputs"].append(encut_dict)

    indf = job.input.to_dict()["kpoints/data_dict"]
    params = indf["Parameter"]
    vals = indf["Value"]

    kpoint_type = vals[2]
    kpoint_grid = vals[3]

    kpoint_dict = {}
    kpoint_dict["label"] = f"kpoint_{kpoint_type}"
    kpoint_dict["value"] = kpoint_grid
    method_dict["dft"]["inputs"].append(kpoint_dict)

    spin_dict = {}
    spin_dict["value"] = mdict.get("ISPIN") == "2"
    spin_dict["label"] = "spin_polarization"
    method_dict["dft"]["inputs"].append(spin_dict)

    if "MAGMOM" in mdict.keys():
        mg_val = mdict["MAGMOM"].split(" ")
        cleaned_mg_val = [x for x in mg_val if x]
        cleaned_mg_val = [float(i) for i in cleaned_mg_val]
        inp_magmom_dict = {}
        inp_magmom_dict["value"] = str(np.array(cleaned_mg_val, dtype="float32"))
        inp_magmom_dict["label"] = "input_magnetic_moments"
        inp_magmom_dict["unit"] = "BOHR-MAGNETON"
        method_dict["dft"]["inputs"].append(inp_magmom_dict)

    elec_min_algo_dict = {}
    elec_min_algo_dict["value"] = mdict["ALGO"]
    elec_min_algo_dict["label"] = "electronic_minimization_algorithm"
    method_dict["dft"]["inputs"].append(elec_min_algo_dict)

    ediff_dict = {}
    ediff_dict["value"] = float(mdict.get("EDIFF", 1e-4))
    ediff_dict["label"] = "electronic_energy_tolerance"
    ediff_dict["unit"] = "EV"
    method_dict["dft"]["inputs"].append(ediff_dict)

    dof = []
    if int(mdict.get("NSW", "0")):  # not a static calculation
        if mdict.get("ISIF", "0") in ["0", "1", "2", "3", "4", "8"]:
            dof.append(
                "http://purls.helmholtz-metadaten.de/asmo/AtomicPositionRelaxation"
            )
        if mdict.get("ISIF") in ["3", "4", "5", "6"]:
            dof.append("http://purls.helmholtz-metadaten.de/asmo/CellShapeRelaxation")
        if mdict.get("ISIF") in ["3", "6", "7", "8"]:
            dof.append("http://purls.helmholtz-metadaten.de/asmo/CellVolumeRelaxation")

        ion_min_algo_dict = {"label": "ionic_minimization_algorithm"}
        ibrion_algo_map = {"1": "rmm-diis", "2": "cg", "3": "damped_md"}
        ionic_min_algo_val = ibrion_algo_map.get(mdict.get("IBRION"))
        if ionic_min_algo_val:
            ion_min_algo_dict = {
                "label": "ionic_minimization_algorithm",
                "value": ionic_min_algo_val,
            }
            method_dict["dft"]["inputs"].append(ion_min_algo_dict)
        else:
            warnings.warn(
                f"Ionic minimization algorithm for IBRION='{mdict.get('IBRION')}' is not yet supported."
            )

        ediffg_val = float(mdict.get("EDIFFG", ediff_dict["value"] * 10))
        if ediffg_val < 0:
            ediffg_dict = {
                "value": abs(ediffg_val),
                "label": "force_tolerance",
                "unit": "EV-PER-ANGSTROM",
            }
        else:
            ediffg_dict = {
                "value": ediffg_val,
                "label": "ionic_energy_tolerance",
                "unit": "EV",
            }
        method_dict["dft"]["inputs"].append(ediffg_dict)
    method_dict["dof"] = dof

    smearing_dict = {"label": "electronic_smearing"}
    ismear_map = {
        "0": "Gaussian",
        "-1": "Fermi",
        "-4": "Tetrahedron",
        "-5": "Tetrahedron_Bloechl",
    }
    smearing_val = mdict.get("ISMEAR", "1")
    if int(smearing_val) > 0:
        smearing_dict["value"] = "Methfessel-Paxton"
        method_dict["dft"]["inputs"].append(smearing_dict)
    elif ismear_map.get(smearing_val):
        smearing_dict["value"] = ismear_map.get(smearing_val)
        method_dict["dft"]["inputs"].append(smearing_dict)
    else:
        warnings.warn(
            f"Smearing algorithm for ISMEAR='{mdict.get('ISMEAR')}' is not yet supported."
        )

    sigma_dict = {}
    if not mdict.get("SIGMA"):
        if smearing_val in ["0", "1"]:
            sigma_dict["value"] = 0.2
        else:
            sigma_dict["value"] = 0.0
    else:
        sigma_dict["value"] = float(mdict["SIGMA"])
    sigma_dict["label"] = "smearing_parameter_sigma"
    sigma_dict["unit"] = "EV"
    method_dict["dft"]["inputs"].append(sigma_dict)

    indf = job.input.to_dict()["potcar/data_dict"]
    xc = indf["Value"][0]
    method_dict["xc_functional"] = xc


def _extract_vasp_calculated_quantities(job, method_dict):
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

    indf = job.input.incar.to_dict()
    params = indf["data_dict"]["Parameter"]
    vals = indf["data_dict"]["Value"]
    mdict = {
        str(p).replace(" ", ""): str(v).replace(" ", "") for p, v in zip(params, vals)
    }

    outputs = []
    outputs.append(
        {
            "label": "final_total_energy",
            "value": np.round(job.output.energy_tot[-1], decimals=4),
            "unit": "EV",
        }
    )
    outputs.append(
        {
            "label": "final_potential_energy",
            "value": np.round(job.output.energy_pot[-1], decimals=4),
            "unit": "EV",
        }
    )
    outputs.append(
        {
            "label": "final_maximum_force",
            "value": np.round(job.output.force_max[-1], decimals=16),
            "unit": "EV-PER-ANGSTROM",
        }
    )
    outputs.append(
        {
            "label": "final_total_volume",
            "value": np.round(job.output.volume[-1], decimals=4),
            "unit": "ANGSTROM3",
        }
    )
    if job.input.incar["ISPIN"] == 2:
        outputs.append(
            {
                "label": "final_total_magnetic_moment",
                "value": np.round(
                    job["output"]["generic"]["dft"]["magnetization"][-1][-1], decimals=4
                ),
                "unit": "BOHR-MAGNETON",
            }
        )
    if "NSW" in mdict.keys() and mdict["NSW"] != "0":
        outputs.append(
            {
                "label": "number_ionic_steps",
                "value": len(job.output.steps) - 1,
            }
        )
        try:
            fpress = (
                1
                / 3
                * (
                    job.output.pressures[-1][0]
                    + job.output.pressures[-1][1]
                    + job.output.pressures[-1][2]
                )
            )
        except (
            IndexError
        ):  # sometimes the avg is already calculated an a single value stored
            fpress = job.output.pressures[-1]
        outputs.append(
            {
                "label": "final_pressure",
                "value": np.round(
                    fpress / units.GPa, decimals=4
                ),  # is in eV/A3, need GPa
                "unit": "GigaPA",
            }
        )
    method_dict["outputs"] = outputs


def export_env(path):
    """Exports to path+_environment.yml"""
    import os
    import platform

    if "Windows" in platform.system():
        os.system(
            'conda env export | findstr -v "^prefix: " > ' + path + "_environment.yml"
        )
    else:
        os.system(
            'conda env export | grep -v "^prefix: " > ' + path + "_environment.yml"
        )


def flatten_cdict(cdict):
    flat = {}
    for k, v in cdict.items():
        if k != "@context":
            if isinstance(v, dict):
                if "label" in v.keys():
                    flat[k] = v["label"]
                else:
                    flat = flat | flatten_cdict(v)
            elif k == "software":
                flat[k] = v[0]["label"]
            elif isinstance(v, list):
                for i in v:
                    if isinstance(i, dict):
                        try:
                            flat[i["label"]] = i["value"]
                        except (
                            KeyError
                        ):  # silently skips over terms that do not have label, value keys
                            pass
                    else:
                        flat[k] = v
            else:
                flat[k] = v
    return flat
