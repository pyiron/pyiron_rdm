def format_json_string(json_string):
    json_string = json_string.replace('\n', '<br>')
    result = []
    for index, char in enumerate(json_string):
        if char == " " and (index == 0 or json_string[index - 1] != ":"):
            result.append("&nbsp;&nbsp;")
        else:
            result.append(char)

    json_string = "".join(result)
    return json_string

def map_cdict_to_ob(o, cdict, concept_dict, add_missing=True):

    # concept_dict to be kicked out
    # cdict = flat concept_dict

    props = {}

    if add_missing: # TODO do this by default or leave as flag? i.e. are these properties *always* there?
        json_file = cdict['path'] + '_concept_dict.json'
        with open(json_file, 'r') as file:
            json_string = file.read()
        json_string = format_json_string(json_string)

        added = {
            'bam_username': o.get_session_info().userName, # TODO can we avoid o as input?
            'conceptual_dictionary': json_string,
            'description': '<p><span style="color:hsl(240,75%,60%);">' + \
                        '<strong>Scroll down below other properties to view conceptual dictionary with ontological ids of selected properties and values.</strong></span>' + \
                        '<br>The conceptual dictionary is in JSON-LD format. Learn more about it <a href="https://www.w3.org/ns/json-ld/">here</a></p>'
        }
        props |= added
    
    # TODO resolve whether we want to keep track of anything (from cdict) that didn't get used?
    
    # job, project, server
    if 'job_name' in cdict.keys():
        props |=  {'$name': cdict['job_name']}
    if 'workflow_manager' in cdict.keys():
        props |=  {'workflow_manager': cdict['workflow_manager']}
    if 'job_status' in cdict.keys():
        if cdict['job_status'] == 'finished': # TODO we also did True if status 'not_converged' or 'converged'??
            props |=  {'sim_job_finished': True}
        else:
            props |=  {'sim_job_finished': False}
    if 'job_starttime' in cdict.keys():
        props |=  {'start_date': cdict['job_starttime']}
    if 'job_starttime' in cdict.keys() and 'job_stoptime' in cdict.keys():
        from datetime import datetime
        import numpy as np
        delta = datetime.strptime(cdict['job_stoptime'], "%Y-%m-%d %H:%M:%S") - datetime.strptime(cdict['job_starttime'], "%Y-%m-%d %H:%M:%S")
        props |=  {'sim_walltime_in_hours': np.round(delta.total_seconds()/3600, 6)}
    if 'sim_coretime_hours' in cdict.keys():
        props |=  {'sim_coretime_in_hours': cdict['sim_coretime_hours']}
    if 'number_cores' in cdict.keys():
        props |=  {'ncores': cdict['number_cores']}

    # scientific
    if 'maximum iterations' in cdict.keys():
        props |=  {'max_iters': cdict['maximum iterations']}
    if 'ionic energy tolerance' in cdict.keys():
        props |=  {'atom_e_tol_ion_in_ev': cdict['ionic energy tolerance']}
    if 'force tolerance' in cdict.keys():
        props |=  {'atom_f_tol_in_ev_a': cdict['force tolerance']}
    if 'NumberIonicSteps' in cdict.keys():
        props |=  {'atom_ionic_steps': cdict['NumberIonicSteps']}
    if 'FinalMaximumForce' in cdict.keys():
        props |=  {'atom_force_max_in_ev_a': cdict['FinalMaximumForce']}
    if 'periodicity_in_x' in cdict.keys():
        props |=  {'periodic_boundary_x': cdict['periodicity_in_x']}
    if 'periodicity_in_y' in cdict.keys():
        props |=  {'periodic_boundary_y': cdict['periodicity_in_y']}
    if 'periodicity_in_z' in cdict.keys():
        props |=  {'periodic_boundary_z': cdict['periodicity_in_z']}
    if 'dof' in cdict.keys():
        props |=  {'atom_cell_vol_relax': True if 'http://purls.helmholtz-metadaten.de/asmo/CellVolumeRelaxation' in cdict['dof'] else False,
                   'atom_cell_shp_relax': True if 'CellShapeRelaxation' in cdict['dof'] else False,
                   'atom_pos_relax': True if 'http://purls.helmholtz-metadaten.de/asmo/AtomicPositionRelaxation' in cdict['dof'] else False
                   }
    if 'FinalTotalEnergy' in cdict.keys():
        props |=  {'atom_fin_tot_eng_in_ev': cdict['FinalTotalEnergy']}
    if 'FinalTotalVolume' in cdict.keys():
        props |=  {'atom_fin_vol_in_a3': cdict['FinalTotalVolume']}
    if 'FinalPotentialEnergy' in cdict.keys():
        props |=  {'atom_fin_pot_eng_in_ev': cdict['FinalPotentialEnergy']}
    if 'molecular_statics' in concept_dict.keys() and 'minimization_algorithm' in cdict.keys():
        description = f'{cdict["job_type"]} simulation using pyiron for energy minimization/structural optimization.' + props['description'] # TODO double check correctness
        if cdict['minimization_algorithm'] == 'fire':
            min_algo = 'MIN_ALGO_FIRE'
        elif cdict['minimization_algorithm'] == 'cg':
            min_algo = 'MIN_ALGO_CG'
        elif cdict['minimization_algorithm'] == 'hftn':
            min_algo = 'MIN_ALGO_HFTN'
        elif cdict['minimization_algorithm'] == 'lbfgs':
            min_algo = 'MIN_ALGO_LBFGS'
        elif cdict['minimization_algorithm'] == 'quickmin':
            min_algo = 'MIN_ALGO_QUICKMIN'
        elif cdict['minimization_algorithm'] == 'sd':
            min_algo = 'MIN_ALGO_STEEP_DESC'
        else:
            raise ValueError('Unknown minimization algorithm')
        props |=  {'atomistic_calc_type': 'atom_calc_struc_opt',
                   'atom_ionic_min_algo': min_algo, 
                   'description': description}
    
    return props

def dataset_job_h5(cdict):
    from datetime import datetime

    # TODO error handling if not job
    ds_props = {
        '$name': cdict['job_name'] + '.h5',
        'production_date': datetime.strptime(cdict['job_stoptime'], "%Y-%m-%d %H:%M:%S").date().strftime("%Y-%m-%d"),
        'file_format': 'HDF5',
        'reference': 'https://github.com/pyiron/pyiron_atomistics/blob/main/pyiron_atomistics/lammps/base.py',
    }
    ds_type = 'PYIRON_JOB'

    return ds_type, ds_props

def dataset_atom_struct_h5(cdict):
    # TODO error handling if not structure
    ds_props = {
        '$name': cdict['structure_name'] + '.h5',
        'multi_mat_scale': 'Electronic/Atomistic',
        'sw_compatibility': 'ASE',
        'file_format': 'HDF5',
    }
    ds_type = 'MAT_SIM_STRUCTURE'

    return ds_type, ds_props

def dataset_env_yml(cdict):
    # TODO error handling if not job
    ds_props = {'$name': cdict['job_name']+'_environment.yml', 
                'env_tool': 'conda'}
    ds_type = 'COMP_ENV'

    return ds_type, ds_props

def dataset_cdict_jsonld(cdict):
    if 'job_name' in cdict.keys(): # TODO could this just be 'name'?
        ds_props = {'$name': cdict['job_name'] + '_concept_dict.json'}
    elif 'structure_name' in cdict.keys():
        ds_props = {'$name': cdict['structure_name'] + '_concept_dict.json'}
    else:
        raise KeyError('Missing job_name or structure_name key missing in conceptual dictionary. Cannot upload.')
    ds_type = 'ATTACHMENT'

    return ds_type, ds_props

