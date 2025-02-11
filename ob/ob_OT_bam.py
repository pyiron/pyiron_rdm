# inventory parents _________________________________________

def material_par(props_dict):
    mat_dict_pct_str = species_by_num_to_pct(props_dict)
    t = 'MATERIAL_V1'
    w = {'chem_species_by_comp_in_pct': mat_dict_pct_str}
    a = ['chem_species_by_comp_in_pct']
    return t, w, a

def intpot_par(cdict):
    t = 'INTERATOMIC_POTENTIAL'
    w = {'$name': cdict['potential']}
    a = ['$name']
    return t, w, a

def pseudopot_par(cdict):     # TODO: complete
    return

def sw_par(cdict):
    t = 'SOFTWARE_CODE'
    c = cdict['software'].replace(' ', '_')
    return t, c

def localws_par(cdict):
    t = 'INSTRUMENT.LOCAL_WORKSTATION'
    c = cdict['host'].upper()  
    return t, c

def wfref_par(cdict):
    t = 'WORKFLOW_REFERENCE'
    if 'murn' in cdict['job_type'].lower():
        c = 'EQ_OF_STATE_FITTING'
    else:
        c = ''
    return t, c

# object types _______________________________________________

def crystal_struct_ot():
    object_type = 'MAT_SIM_STRUCTURE.CRYSTAL'
    datasets = ['structure_h5', 'cdict_json']
    parents = ['material']
    return object_type, datasets, parents

def pyiron_job_ot():
    object_type = 'PYIRON_JOB'
    datasets = ['job_h5', 'env_yml', 'cdict_json']
    parents = ['local_workstation']
    return object_type, datasets, parents

def lammps_job_ot():
    object_type = 'PYIRON_JOB.LAMMPS'
    datasets = ['job_h5', 'env_yml', 'cdict_json']
    parents = ['interatomic_potential', 'software', 'local_workstation']
    return object_type, datasets, parents
    
def vasp_job_ot():
    object_type = 'PYIRON_JOB.VASP'
    datasets = ['job_h5', 'env_yml', 'cdict_json']
    parents = ['pseudopotential', 'software', 'local_workstation']
    return object_type, datasets, parents
    
def murn_job_ot():
    object_type = 'PYIRON_JOB.MURNAGHAN'
    datasets = ['job_h5', 'env_yml', 'cdict_json']
    parents = ['wf_reference', 'local_workstation']
    return object_type, datasets, parents

# calls ______________________________________________________
    
def get_ot_info(cdict):
    if 'structure_name' in cdict.keys():
        return crystal_struct_ot
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
    elif parent_name == 'local_workstation':
        t, c = localws_par(cdict)
    elif parent_name == 'software':
        t, c = sw_par(cdict)
    elif parent_name == 'interatomic_potential':
        t, w, a, = intpot_par(cdict)
    elif parent_name == 'pseudopotential':
        t, w, a = pseudopot_par(cdict)
    elif parent_name == 'wf_reference':
        t, c, = wfref_par(cdict)

    return t, w, a, c
    
# else ________________________________________________________

def species_by_num_to_pct(props):
    import numpy as np
    species_by_num = eval(props['chem_species_by_n_atoms'])
    species_by_pct = {at: np.round(species_by_num[at]*100/props['n_atoms_total'], 2) for at in species_by_num.keys()}
    return str(species_by_pct)