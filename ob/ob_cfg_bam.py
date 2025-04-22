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

def map_cdict_to_ob(o, cdict, concept_dict):

    # concept_dict to be kicked out
    # cdict = flat concept_dict

    if 'structure_name' in cdict.keys():
        json_file = cdict['path'] + cdict['structure_name'] + '_concept_dict.json'
        props = {}
    else:
        json_file = cdict['path'] + '_concept_dict.json'
        props = {'bam_username': o.get_session_info().userName} # TODO can we avoid o as input?

    with open(json_file, 'r') as file:
        json_string = file.read()
    json_string = format_json_string(json_string)

    props |= {
        'conceptual_dictionary': json_string,
        'description': '<p><span style="color:hsl(240,75%,60%);">' + \
                    '<strong>Scroll down below other properties to view conceptual dictionary with ontological ids of selected properties and values.</strong></span>' + \
                    '<br>The conceptual dictionary is in JSON-LD format. Learn more about it <a href="https://www.w3.org/ns/json-ld/">here</a></p>'
            }
    description = props['description']
    
    # TODO resolve whether we want to keep track of anything (from cdict) that didn't get used?
    
    if 'workflow_manager' in cdict.keys():
        props['workflow_manager'] = cdict['workflow_manager']

    # structure
    if 'structure_name' in cdict.keys():
        props['$name'] = cdict['structure_name']
        props['description'] = 'Crystal structure generated using pyiron.' + props['description']
        map_struct_to_ob(props, cdict, concept_dict)
        return props

    # job, project, server
    elif 'job_name' in cdict.keys():
        props['$name'] = cdict['job_name']
        if 'job_status' in cdict.keys():
            if cdict['job_status'] == 'finished': # TODO we also did True if status 'not_converged' or 'converged'??
                props['sim_job_finished'] = True
            else:
                props['sim_job_finished'] = False
        if 'job_starttime' in cdict.keys():
            props['start_date'] = cdict['job_starttime']
        if 'job_starttime' in cdict.keys() and 'job_stoptime' in cdict.keys():
            from datetime import datetime
            import numpy as np
            delta = datetime.strptime(cdict['job_stoptime'], "%Y-%m-%d %H:%M:%S") - datetime.strptime(cdict['job_starttime'], "%Y-%m-%d %H:%M:%S")
            props['sim_walltime_in_hours'] = np.round(delta.total_seconds()/3600, 6)
        if 'sim_coretime_hours' in cdict.keys():
            props['sim_coretime_in_hours'] = cdict['sim_coretime_hours']
        if 'number_cores' in cdict.keys():
            props['ncores'] = cdict['number_cores']

        # scientific
        if 'maximum_iterations' in cdict.keys():
            props['max_iters'] = cdict['maximum_iterations']
        if 'ionic_energy_tolerance' in cdict.keys():
            props['atom_e_tol_ion_in_ev'] = cdict['ionic_energy_tolerance']
        if 'force_tolerance' in cdict.keys():
            props['atom_f_tol_in_ev_a'] = cdict['force_tolerance']
        if 'number_ionic_steps' in cdict.keys():
            props['atom_ionic_steps'] = cdict['number_ionic_steps']
        if 'final_maximum_force' in cdict.keys():
            props['atom_force_max_in_ev_a'] = cdict['final_maximum_force']
        if 'periodicity_in_x' in cdict.keys():
            props['periodic_boundary_x'] = cdict['periodicity_in_x']
        if 'periodicity_in_y' in cdict.keys():
            props['periodic_boundary_y'] = cdict['periodicity_in_y']
        if 'periodicity_in_z' in cdict.keys():
            props['periodic_boundary_z'] = cdict['periodicity_in_z']
        if 'dof' in cdict.keys():
            props |=  {'atom_cell_vol_relax': True if 'http://purls.helmholtz-metadaten.de/asmo/CellVolumeRelaxation' in cdict['dof'] else False,
                    'atom_cell_shp_relax': True if 'CellShapeRelaxation' in cdict['dof'] else False,
                    'atom_pos_relax': True if 'http://purls.helmholtz-metadaten.de/asmo/AtomicPositionRelaxation' in cdict['dof'] else False
                    }
        if 'final_total_energy' in cdict.keys():
            props['atom_fin_tot_eng_in_ev'] = cdict['final_total_energy']
        if 'final_total_volume' in cdict.keys():
            props['atom_fin_vol_in_a3'] = cdict['final_total_volume']
        if 'final_potential_energy' in cdict.keys():
            props['atom_fin_pot_eng_in_ev'] = cdict['final_potential_energy']
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
            if 'target_pressure' in cdict.keys() and cdict['target_pressure'] is not None:
                props['atom_targ_press_in_gpa'] = cdict['target_pressure']
        if 'molecular_dynamics' in concept_dict.keys():
            if 'http://purls.helmholtz-metadaten.de/asmo/MicrocanonicalEnsemble' in cdict['ensemble']:
                description = f'{cdict["job_type"]} simulation using pyiron for microcanonical ensemble.' + props['description']
                props['atom_md_ensemble'] = 'TD_ENSEMBLE_NVE'
            elif 'http://purls.helmholtz-metadaten.de/asmo/CanonicalEnsemble' in cdict['ensemble']:
                description = f'{cdict["job_type"]} simulation using pyiron for canonical ensemble.' + props['description']
                props['atom_md_ensemble'] = 'TD_ENSEMBLE_ATOM_ENS.NVT'
            elif 'http://purls.helmholtz-metadaten.de/asmo/IsothermalIsobaricEnsemble' in cdict['ensemble']:
                description = f'{cdict["job_type"]} simulation using pyiron for isothermal-isobaric ensemble.' + props['description'] # TODO double check correctness
                props['atom_md_ensemble'] = 'TD_ENSEMBLE_NPT'
            props |= {'atomistic_calc_type': 'Atom_calc_md',
                    'description': description}
            
            if 'timestep' in cdict.keys():
                props['atom_md_time_stp_in_ps'] = cdict['timestep']
            if 'simulation_time' in cdict.keys():
                props['atom_sim_time_ps_in_ps'] = cdict['simulation_time']
            if 'average_total_energy' in cdict.keys():
                props['atom_avg_tot_eng_in_ev'] = cdict['average_total_energy']
            if 'average_potential_energy' in cdict.keys():
                props['atom_avg_pot_eng_in_ev'] = cdict['average_potential_energy']
            if 'average_temperature' in cdict.keys():
                props['atom_md_avg_temp_in_k'] = cdict['average_temperature']
            if 'average_pressure' in cdict.keys():
                props['atom_avg_press_in_gpa'] = cdict['average_pressure']
            if 'average_total_volume' in cdict.keys():
                props['atom_avg_vol_in_a3'] = cdict['average_total_volume']
            if 'initial_temperature' in cdict.keys():
                props['atom_md_init_temp_in_k'] = cdict['initial_temperature']
            if 'target_temperature' in cdict.keys():
                props['atom_md_targ_temp_in_k'] = cdict['target_temperature']
            if 'initial_pressure' in cdict.keys():
                props['atom_md_init_press_in_gpa'] = cdict['initial_pressure']
            if 'target_pressure' in cdict.keys():
                props['atom_targ_press_in_gpa'] = cdict['target_pressure']

            
        if 'job_type' in cdict.keys() and 'Murn' in cdict['job_type']: # TODO general way to do this? Put together 
            props['description'] = 'Murnaghan job for structural optimization.' + props['description']
        if 'strain_axes' in cdict.keys():
            props['murn_strain_axes'] = cdict['strain_axes']
        if 'number_of_data_points' in cdict.keys():
            props['murn_n_data_points'] = cdict['number_of_data_points']
        if 'volume_range' in cdict.keys():
            props['murn_strainvol_range'] = cdict['volume_range']
        if 'equilibrium_bulk_modulus' in cdict.keys():
            props['atom_equil_k_mod_in_gpa'] = cdict['equilibrium_bulk_modulus']
        if 'equilibrium_total_energy' in cdict.keys():
            props['atom_equil_toteng_in_ev'] = cdict['equilibrium_total_energy']
        if 'equilibrium_volume' in cdict.keys():
            props['atom_equil_vol_in_a3'] = cdict['equilibrium_volume']
        if 'equation_of_state_fit' in cdict.keys():
            if cdict['equation_of_state_fit'] == 'http://purls.helmholtz-metadaten.de/asmo/BirchMurnaghan':
                props['murn_eqn_of_state'] = 'EOS_BIRCH_MURNAGHAN'
            elif cdict['equation_of_state_fit'] == 'http://purls.helmholtz-metadaten.de/asmo/Murnaghan':
                props['murn_eqn_of_state'] = 'EOS_MURNAGHAN'
            elif cdict['equation_of_state_fit'] == 'http://purls.helmholtz-metadaten.de/asmo/Vinet':
                props['murn_eqn_of_state'] = 'EOS_VINET'
            elif cdict['equation_of_state_fit'] == 'http://purls.helmholtz-metadaten.de/asmo/PolynomialFit':
                props['murn_eqn_of_state'] = 'EOS_POLYNOMIAL'
                props['murn_fit_eqn_order'] = cdict['fit_order']  # TODO test this
            else: 
                raise ValueError('Unknown equation of state')  # TODO is this necessary?
    
    else:
        print("Neither structure_name nor job_name in the object conceptual dictionary. \
              OpenBIS properties most likely incomplete.")

    return props

def map_struct_to_ob(props, cdict, concept_dict):
    if 'atoms' in concept_dict.keys():
        sorted_atoms = sorted(
            [atom for atom in concept_dict['atoms'] if atom['label'] != 'total_number_atoms'],
            key=lambda x: x['label']
        )
        species = {i['label']: i['value'] for i in sorted_atoms}
        props['chem_species_by_n_atoms'] = str(species)
    if 'total_number_atoms' in cdict.keys():
        props['n_atoms_total'] = cdict['total_number_atoms']
    if 'simulation_cell_lengths' in cdict.keys():
        props['sim_cell_lengths_in_a'] = cdict['simulation_cell_lengths']
    if 'simulation_cell_vectors' in cdict.keys():
        props['sim_cell_vectors'] = cdict['simulation_cell_vectors']
    if 'simulation_cell_angles' in cdict.keys():
        props['sim_cell_angles_in_deg'] = cdict['simulation_cell_angles']
    if 'simulation_cell_volume' in cdict.keys():
        props['sim_cell_volume_in_a3'] = cdict['simulation_cell_volume']  
    if 'crystal_orientation' in cdict.keys(): 
        props['crystal_orientation'] = cdict['crystal_orientation']
    if 'lattice_parameter_a' in cdict.keys(): 
        props['lattice_param_a_in_a'] = cdict['lattice_parameter_a']
    if 'lattice_parameter_b' in cdict.keys(): 
        props['lattice_param_b_in_a'] = cdict['lattice_parameter_b']
    if 'lattice_parameter_c' in cdict.keys(): 
        props['lattice_param_c_in_a'] = cdict['lattice_parameter_c']
    if 'lattice_parameter_c_over_a' in cdict.keys(): 
        props['lattice_c_over_a'] = cdict['lattice_parameter_c_over_a']
    if 'lattice_angle_alpha' in cdict.keys(): 
        props['lattice_angalpha_in_deg'] = cdict['lattice_angle_alpha']
    if 'lattice_angle_beta' in cdict.keys(): 
        props['lattice_angbeta_in_deg'] = cdict['lattice_angle_beta']
    if 'lattice_angle_gamma' in cdict.keys(): 
        props['lattice_anggamma_in_deg'] = cdict['lattice_angle_gamma']
    if 'lattice_volume' in cdict.keys(): 
        props['lattice_volume_in_a3'] = cdict['lattice_volume']
    if 'space_group' in cdict.keys():
        spg_map = get_space_group_mapping(cdict['space_group'])
        props['space_group'] = spg_map
    if 'bravais_lattice' in cdict.keys():
        bvl_map = get_bravais_lattice_mapping(cdict['bravais_lattice'])
        props['bravais_lattice'] = bvl_map

def dataset_job_h5(cdict):
    from datetime import datetime

    # TODO error handling if not job
    ds_props = {
        '$name': cdict['job_name'] + '.h5',
        'production_date': datetime.strptime(cdict['job_stoptime'], "%Y-%m-%d %H:%M:%S").date().strftime("%Y-%m-%d"),
        'file_format': 'HDF5'
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

def get_space_group_mapping(spg):
    if spg == 'Im-3m':
        return ('SPACE_GROUP.IM-3M').lower()
    elif spg == 'Fm-3m':
        return ('SPACE_GROUP.FM-3M').lower()
    elif spg == 'P6_3/mmc':
        return ('SPACE_GROUP.P63_MMC').lower()
    else:
        raise ValueError(f'Invalid Bravais lattice, maybe a formatting error?')

def get_bravais_lattice_mapping(bvl):
    if bvl == 'bcc':
        return ('BODY_CENTER_CUBIC').lower()
    elif bvl == 'fcc':
        return ('FACE_CENTER_CUBIC').lower()
    elif bvl == 'hcp':
        return ('HEX_CLOSE_PACK').lower()
    else:
        raise ValueError(f'Invalid Bravais lattice, maybe a formatting error?')