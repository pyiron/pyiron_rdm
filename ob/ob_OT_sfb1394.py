# inventory parents _________________________________________

def material_par(props_dict):
    mat_dict_pct_str = species_by_num_to_pct(props_dict)
    # t = 'CRYSTALLINE_MATERIAL'
    t = 'MATERIAL'
    w = {'compo_atomic_percent': mat_dict_pct_str}
    a = ['compo_atomic_percent']
    return t, w, a

def intpot_par(cdict): 
    t = 'INTERATOMIC_POTENTIAL'
    w = {'external_identifier': cdict['potential']}
    a = ['external_identifier']
    return t, w, a

def workstation_par(cdict):
    t = 'COMPUTE_RESOURCE'
    c = '*' + cdict['host']
    return t, c

def pseudopot_par(cdict):     # TODO: complete
    return

def sw_par(cdict): 
    t = 'SOFTWARE'
    c = cdict['software']
    for ch in [' ', '-', '_', '.']:
        c = c.replace(ch, '')
    return t, c

# object types _______________________________________________

def struct_ot(): # TODO update material
    object_type = 'SAMPLE'
    datasets = ['structure_h5', 'cdict_json']
    parents = ['material']
    # parents = []
    return object_type, datasets, parents

def pyiron_job_ot():
    object_type = 'PYIRON_JOB_GENERIC'
    datasets = ['job_h5', 'env_yml', 'cdict_json']
    parents = ['workstation']
    return object_type, datasets, parents

def lammps_job_ot():
    object_type = 'PYIRON_JOB_LAMMPS'
    datasets = ['job_h5', 'env_yml', 'cdict_json']
    parents = ['interatomic_potential', 'software', 'workstation']
    return object_type, datasets, parents
    
def vasp_job_ot():
    object_type = 'PYIRON_JOB_VASP'
    datasets = ['job_h5', 'env_yml', 'cdict_json']
    parents = ['pseudopotential', 'software', 'workstation']
    return object_type, datasets, parents
    
def murn_job_ot():
    object_type = 'PYIRON_JOB_MURNAGHAN'
    datasets = ['job_h5', 'env_yml', 'cdict_json']
    parents = ['workstation']
    return object_type, datasets, parents

# calls ______________________________________________________
    
def get_ot_info(cdict):
    if 'structure_name' in cdict.keys():
        return struct_ot
    elif 'job_type' in cdict.keys():
        if 'lammps' in cdict['job_type'].lower():
            return lammps_job_ot
        elif 'vasp' in cdict['job_type'].lower():
            return vasp_job_ot
        elif 'murn' in cdict['job_type'].lower():
            return murn_job_ot
        else:
            return pyiron_job_ot
    else:
        raise ValueError('Neither structure_name nor job_type in conceptual dictionary. Cannot proceed.')
    
def get_inv_parent(parent_name, cdict, props_dict):
    t, w, a, c = '', {}, [], ''
    if parent_name == 'material':
        t, w, a = material_par(props_dict)
    elif parent_name == 'workstation':
        t, c = workstation_par(cdict)
    elif parent_name == 'software':
        t, c = sw_par(cdict)
    elif parent_name == 'interatomic_potential':
        t, w, a, = intpot_par(cdict)
    elif parent_name == 'pseudopotential':
        t, w, a = pseudopot_par(cdict)

    return t, w, a, c
    
# else ________________________________________________________

def species_by_num_to_pct(props, max_elements = 8):
    species_by_pct = {props[f'element_{i}']: props[f'element_{i}_at_percent'] for i in range(1, max_elements+1) 
          if f'element_{i}' in props and f'element_{i}_at_percent' in props}
    return str(species_by_pct)