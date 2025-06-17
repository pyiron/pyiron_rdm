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

    # cdict = flat concept_dict

    if 'structure_name' in cdict.keys():
        json_file = cdict['path'] + cdict['structure_name'] + '_concept_dict.json'
        props = {}
    else:
        json_file = cdict['path'] + '_concept_dict.json'
        props = {'user_name': o.get_session_info().userName}

    with open(json_file, 'r') as file:
        json_string = file.read()
    json_string = format_json_string(json_string)

    props |= {
        'pyiron_conceptual_dictionary': json_string,
        'description_multiline': '<p><span style="color:hsl(240,75%,60%);">' + \
                    '<strong>Scroll down below other properties to view conceptual dictionary with ontological ids of selected properties and values.</strong></span>' + \
                    '<br>The conceptual dictionary is in JSON-LD format. Learn more about it <a href="https://www.w3.org/ns/json-ld/">here</a></p>'
            }
    description = props['description_multiline']
    
    if 'workflow_manager' in cdict.keys():
        props['workflow_manager'] = cdict['workflow_manager']

    # structure
    if 'structure_name' in cdict.keys():
        props['$name'] = cdict['structure_name']
        props['description_multiline'] = 'Crystal structure generated using pyiron.' + props['description_multiline']
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
                    'atom_cell_shp_relax': True if 'http://purls.helmholtz-metadaten.de/asmo/CellShapeRelaxation' in cdict['dof'] else False,
                    'atom_pos_relax': True if 'http://purls.helmholtz-metadaten.de/asmo/AtomicPositionRelaxation' in cdict['dof'] else False
                    }
        if 'final_total_energy' in cdict.keys():
            props['atom_fin_tot_eng_in_ev'] = cdict['final_total_energy']
        if 'final_total_volume' in cdict.keys():
            props['atom_fin_vol_in_a3'] = cdict['final_total_volume']
        if 'final_potential_energy' in cdict.keys():
            props['atom_fin_pot_eng_in_ev'] = cdict['final_potential_energy']
        if 'molecular_statics' in concept_dict.keys() and 'minimization_algorithm' in cdict.keys():
            description = f'{cdict["job_type"]} simulation using pyiron for energy minimization/structural optimization.' + props['description_multiline'] # TODO double check correctness
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
                    'description_multiline': description}
            if 'target_pressure' in cdict.keys() and cdict['target_pressure'] is not None:
                props['atom_targ_press_in_gpa'] = cdict['target_pressure']
        if 'molecular_dynamics' in concept_dict.keys():
            if 'http://purls.helmholtz-metadaten.de/asmo/MicrocanonicalEnsemble' in cdict['ensemble']:
                description = f'{cdict["job_type"]} simulation using pyiron for microcanonical ensemble.' + props['description_multiline']
                props['atom_md_ensemble'] = 'TD_ENSEMBLE_NVE'
            elif 'http://purls.helmholtz-metadaten.de/asmo/CanonicalEnsemble' in cdict['ensemble']:
                description = f'{cdict["job_type"]} simulation using pyiron for canonical ensemble.' + props['description_multiline']
                props['atom_md_ensemble'] = 'TD_ENSEMBLE_ATOM_ENS_NVT'
            elif 'http://purls.helmholtz-metadaten.de/asmo/IsothermalIsobaricEnsemble' in cdict['ensemble']:
                description = f'{cdict["job_type"]} simulation using pyiron for isothermal-isobaric ensemble.' + props['description_multiline'] # TODO double check correctness
                props['atom_md_ensemble'] = 'TD_ENSEMBLE_NPT'
            props |= {'atomistic_calc_type': 'Atom_calc_md',
                    'description_multiline': description}
            
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
            
        if 'job_type' in cdict.keys() and 'Murn' in cdict['job_type']: 
            props['description_multiline'] = 'Murnaghan job for structural optimization.' + props['description_multiline']
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
                props['murn_fit_eqn_order'] = cdict['fit_order'] 
            else: 
                import warnings
                warnings.warn('Unknown equation of state')

        if cdict.get('energy_cutoff'):
            props['atom_e_cutoff_in_ev'] = cdict['energy_cutoff']
        if 'xc_functional' in cdict.keys():
            if cdict['xc_functional'] == 'LDA':
                props['atom_xc_functional'] = 'XC_FUNC_LDA'
            elif cdict['xc_functional'] in ('PBE', 'GGA'):
                props['atom_xc_functional'] = 'XC_FUNC_PBE'
            else:
                import warnings
                warnings.warn(f"XC functional '{props['atom_xc_functional']}' is not yet mapped.")

        if 'spin_polarization' in cdict.keys():
            props['atom_spin_polarized'] = cdict['spin_polarization']
        if 'electronic_smearing' in cdict.keys():
            elsmear_map = {
                'Methfessel-Paxton': 'ELEC_SMEAR_MP',
                'Gaussian': 'ELEC_SMEAR_GAUSS',
                'Fermi': 'ELEC_SMEAR_FERMI',
                'Tetrahedron': 'ELEC_SMEAR_TET',
                'Tetrahedron_Bloechl': 'ELEC_SMEAR_TET_BL'
            }
            elsmear_val = elsmear_map.get(cdict.get('electronic_smearing'))
            if elsmear_val:
                props |= {'electronic_smearing': elsmear_val}
        if 'electronic_energy_tolerance' in cdict.keys():
            props['atom_el_e_tol_in_ev'] = cdict['electronic_energy_tolerance']
        if 'smearing_parameter_sigma' in cdict.keys():
            props['atom_sigma_in_ev'] = cdict['smearing_parameter_sigma']
        if 'final_pressure' in cdict.keys():
                props['atom_fin_press_in_gpa'] = cdict['final_pressure']
        if 'final_total_magnetic_moment' in cdict.keys():
                props['atom_fin_totmgmo_in_mub'] = str(cdict['final_total_magnetic_moment'])
        if 'dft' in concept_dict.keys():
            if len(cdict['dof']) == 0:
                pass
            else:
                description = f'{cdict["job_type"]} simulation using pyiron for energy minimization/structural optimization.' + props['description_multiline'] # TODO double check correctness
                min_algo = None
                if cdict['ionic_minimization_algorithm'] == 'rmm-diis':
                    min_algo = 'MIN_ALGO_RMM_DIIS'
                elif cdict['ionic_minimization_algorithm'] == 'cg':
                    min_algo = 'MIN_ALGO_CG'
                elif cdict['ionic_minimization_algorithm'] == 'damped_md':
                    min_algo = 'MIN_ALGO_DAMPED_MD'
                if min_algo:
                    props |=  {'atomistic_calc_type': 'atom_calc_struc_opt',
                        'atom_ionic_min_algo': min_algo, 
                        'description_multiline': description}
            if cdict.get('electronic_minimization_algorithm', '').lower() in ['normal', 'fast', 'veryfast']:
                elec_min_algo = 'MIN_ALGO_RMM_DIIS'
                props |=  {'atom_elec_min_algo': elec_min_algo}
            else:
                import warnings  
                warnings.warn(f"Electronic minimization algorithm for ALGO='{cdict.get('electronic_minimization_algorithm')}' not yet mapped.")
 
        if 'kpoint_Monkhorst_Pack' in cdict.keys():
            kpoint_algo = 'KPOINTS_MP'
            kpts_x = int(cdict['kpoint_Monkhorst_Pack'].split()[0])
            kpts_y = int(cdict['kpoint_Monkhorst_Pack'].split()[1])
            kpts_z = int(cdict['kpoint_Monkhorst_Pack'].split()[2])
            props |=  {'atom_kpoint_type': kpoint_algo,
            'atomistic_n_kpt_x': kpts_x,
            'atomistic_n_kpt_y': kpts_y,
            'atomistic_n_kpt_z': kpts_z}


    else:
        print("Neither structure_name nor job_name in the object conceptual dictionary. \
              OpenBIS properties most likely incomplete.")

    return props
    
def map_struct_to_ob(props, cdict, concept_dict):
    props['location'] = 'virtual'
    from datetime import date
    props['date'] = str(date.today())

    if 'atoms' in concept_dict.keys():
        import numpy as np
        sorted_atoms = sorted(
            [atom for atom in concept_dict['atoms'] if atom['label'] != 'total_number_atoms'],
            key=lambda x: x['label']
        )
        
        props['composition_desc'] = 'ATOMIC_FRACTION'
        i = 1
        for species in sorted_atoms: # TODO: include break for 9+ elements (?)
            prop_el = 'element_' + str(i)
            prop_el_pct = 'element_' + str(i) + '_at_percent'
            prop_el_num = 'element_' + str(i) + '_number'
            props[prop_el] = species['label']
            props[prop_el_pct] = np.round(species['value']*100/cdict['total_number_atoms'], 2)
            props[prop_el_num] = species['value']
            i += 1

    if 'simulation_cell_lengths' in cdict.keys():
        dim_list = [float(x) for x in cdict['simulation_cell_lengths'].strip('[]').split(',')]
        props['sample_dim'] = ' x '.join([f"{x} A" for x in dim_list])
        # props['unit_cell_lengths'] = cdict['simulation_cell_lengths'].strip('[]')

    if 'simulation_cell_vectors' in cdict.keys():
        props['sim_cell_vectors'] = cdict['simulation_cell_vectors']
    if 'space_group_number' in cdict.keys():
        props['space_group'] = cdict['space_group_number']

def dataset_job_h5(cdict, name_suffix='_0'):
    from datetime import datetime

    # TODO error handling if not job
    ds_props = {
        '$name': cdict['job_name'] + '.h5' + name_suffix
    }
    ds_type = 'PYIRON_HDF5'

    return ds_type, ds_props

def dataset_atom_struct_h5(cdict, name_suffix='_0'):
    # TODO error handling if not structure
    # !!! relies on concept_dict.py first adding chemical species info before anything else
    from datetime import date
    import numpy as np
    ds_props = {
        '$name': cdict['structure_name'] + '.h5' + name_suffix,
        'date': str(date.today()),
        'file_source': 'pyiron'
    }

    from itertools import takewhile
    sorted_atoms = dict(sorted(takewhile(lambda item: item[0] != 'total_number_atoms', cdict.items()),
                                key=lambda x: x[0]))
    i = 1
    for species, count in sorted_atoms.items(): # TODO: include break for 9+ elements (?)
        prop_el = 'element_' + str(i)
        prop_el_pct = 'element_' + str(i) + '_at_percent'
        prop_el_num = 'element_' + str(i) + '_number'
        ds_props[prop_el] = species
        ds_props[prop_el_pct] = np.round(count*100/cdict['total_number_atoms'], 2)
        ds_props[prop_el_num] = count
        i += 1
    ds_props['number_of_atoms'] = cdict['total_number_atoms']
    ds_props['number_of_species'] = len(sorted_atoms.keys())
    ds_props['list_of_species'] = ', '.join(sorted_atoms.keys())
    if 'simulation_cell_angles' in cdict.keys():
        angles = [float(x) for x in cdict['simulation_cell_angles'].strip('[]').split(',')]
        ds_props['angle_alpha'], ds_props['angle_beta'], ds_props['angle_gamma'] = angles
    if 'simulation_cell_lengths' in cdict.keys():
        lengths = [float(x) for x in cdict['simulation_cell_lengths'].strip('[]').split(',')]
        ds_props['box_length_a'], ds_props['box_length_b'], ds_props['box_length_c'] = lengths
    if 'space_group_number' in cdict.keys():
        ds_props['space_group'] = cdict['space_group_number']

    ds_type = 'CRYS-STRUCT_DATA'

    return ds_type, ds_props

def dataset_env_yml(cdict, name_suffix='_0'):
    # TODO error handling if not job
    ds_props = {'$name': cdict['job_name']+'_environment.yml' + name_suffix, 
                'env_tool': 'conda'}
    ds_type = 'COMP_ENV'

    return ds_type, ds_props

def dataset_cdict_jsonld(cdict, name_suffix='_0'):
    if 'job_name' in cdict.keys():
        ds_props = {'$name': cdict['job_name'] + '_concept_dict.json' + name_suffix}
    elif 'structure_name' in cdict.keys():
        ds_props = {'$name': cdict['structure_name'] + '_concept_dict.json' + name_suffix}
    else:
        raise KeyError('Missing job_name or structure_name key missing in conceptual dictionary. Cannot upload.')
    ds_type = 'PYIRON_CONCEPT_DICT_DATA'

    return ds_type, ds_props