# TODO better imports for multiple openbis instances

def openbis_login(url, username, s3_config_path=None):
    from getpass import getpass

    if s3_config_path:
        from ob.OpenbisAixTended import OpenbisWithS3
        o = OpenbisWithS3(url=url, verify_certificates=True, 
                          s3_config_path=s3_config_path, standardize_filenames=True)
    else:
        from pybis import Openbis
        o = Openbis(url)

    if not o.is_session_active():
        o.login(username, getpass('Enter openBIS password: '), save_token=True) # save the session token in ~/.pybis/example.com.token
    return o

# def format_json_string(json_string):
#     json_string = json_string.replace('\n', '<br>')
#     result = []
#     for index, char in enumerate(json_string):
#         if char == " " and (index == 0 or json_string[index - 1] != ":"):
#             result.append("&nbsp;&nbsp;")
#         else:
#             result.append(char)

#     json_string = "".join(result)
#     return json_string

def flatten_cdict(cdict):
        flat = {}
        for k, v in cdict.items():
            if k != '@context':
                if isinstance(v, dict):
                    if 'label' in v.keys():
                        flat[k] = v['label']
                    else:
                        flat = flat | flatten_cdict(v)
                elif k == 'software':
                    flat[k] = v[0]['label']
                elif isinstance(v, list):
                    for i in v: 
                        if isinstance(i, dict):
                            try:
                                flat[i['label']] = i['value']
                            except KeyError: # silently skips over terms that do not have label, value keys
                                pass
                        else:
                            flat[k] = v
                else:
                    flat[k] = v
        return flat

def openbis_upload(o, space, project, collection, concept_dict, ob_info):
    # TODO or allow to skip the last two a flag later?
    # TODO also type checking on values
    req_keys_ob_info = ['object_type', 'parents', 'datasets'] 
    if not all(k in ob_info.keys() for k in req_keys_ob_info):
        raise KeyError('The ob_info dictionary must contain all following keys: object_type, parents, datasets.')
    
    if 'S3' in str(type(o)):
        kind = 'LINK'
    else:
        kind = 'PHYSICAL'

    ob_coll = '/' + space + '/' + project + '/' + collection
    ob_project_obj = o.get_project('/' + space + '/' + project)

    ob_jobtype = ob_info['object_type']
    objects = ob_project_obj.get_objects(type = ob_jobtype)
    exists = False
    for object_ in objects:
        job_name = concept_dict['job_details'][0]['value'] # TODO do better !!!
        if object_.p.get('$name') == job_name:
            exists = True
            found_object = object_

    if exists:
        print("===================\n")
        print(f"Job already exists! Found job in: {found_object.identifier}\n")
        print("===================\n")
        print("Found job properties:\n")
        from IPython.display import display
        display(found_object.p)
        return found_object
    
    else:
        from ob.ob_cfg_bam import map_cdict_to_ob
        cdict = flatten_cdict(concept_dict)
        props_dict = map_cdict_to_ob(o, cdict, concept_dict)
        parents_to_link = ob_info['parents']
        ob_parents = []
        for k, v in parents_to_link.items(): # TODO refactor
            if 'MATERIAL' in k:
                mat_dict_pct_str = species_by_num_to_pct(props_dict)
                parent_material = o.get_objects(
                    type       = k,
                    where      = {'chem_species_by_comp_in_pct': mat_dict_pct_str},
                    attrs      = ['chem_species_by_comp_in_pct']
                )[0]       
                if parent_material:
                    ob_parents.append(parent_material)
                else:
                    print(f"No parents of the type {k} and chemical composition {mat_dict_pct_str} found, upload will not proceed.\
                        Please create them first and then try again.")
            elif 'WORKFLOW_REFERENCE' in k:
                parent = o.get_objects(
                    type       = k,
                    where      = {'$name': v}, # TODO map by code instead of name
                    attrs      = ['$name']
                )[0]         
                if parent:
                    ob_parents.append(parent)
                else:
                    print(f"No parents of the type {k} and name {v} found, upload will not proceed.\
                        Please create them first and then try again.")
            else:
                parent = o.get_objects(
                    type       = k,
                    where      = {'$name': cdict[v]}, # TODO map by code instead of name
                    attrs      = ['$name']
                )[0]         
                if parent:
                    ob_parents.append(parent)
                else:
                    print(f"No parents of the type {k} and name {cdict[v]} found, upload will not proceed.\
                        Please create them first and then try again.")
        if len(ob_parents) == len(parents_to_link.keys()): # Found all parents needed
            object_ = o.new_object(
                type       = ob_jobtype,
                space      = space,
                experiment = ob_coll,
                parents    = ob_parents,
                props      = props_dict
            )
            object_.save()

            ds_list = ob_info['datasets']
            for ds in ds_list:
                if ds == 'job_h5':
                    from ob.ob_cfg_bam import dataset_job_h5 as dataset_info
                    file_path = cdict['path'] + '.h5'
                elif ds == 'structure_h5':
                    from ob.ob_cfg_bam import dataset_atom_struct_h5 as dataset_info
                    file_path = cdict['path'] + cdict['structure_name'] + '.h5'
                elif ds == 'env_yml':
                    from ob.ob_cfg_bam import dataset_env_yml as dataset_info
                    file_path = cdict['path'] + '_environment.yml'
                elif ds == 'cdict_json':
                    from ob.ob_cfg_bam import dataset_cdict_jsonld as dataset_info
                    if 'structure_name' in cdict.keys():
                        file_path = cdict['path'] + cdict['structure_name'] + '_concept_dict.json'
                    else:
                        file_path = cdict['path'] + '_concept_dict.json'
                else:
                    raise ValueError(f'Dataset type {ds} not recognised. Supported datasets: job_h5, structure_h5, env_yml, cdict_json.')

                ds_type, ds_props = dataset_info(cdict)
                upload_dataset(o, object_, ob_coll, ds_type, ds_props, file_path, kind)

        # if show_object:
        #     from IPython.display import display
        #     display(object_.p)
        
        return object_
    
def upload_dataset(o, ob_object, collection, ds_type, ds_props, file_path, kind):
    # try:
    ds_hdf = o.new_dataset(
        type       = ds_type,
        collection = collection,
        object     = ob_object,
        files      = [file_path],
        kind       = kind,
        props      = ds_props
    )
    ds_hdf.save()
    # except ValueError: # TODO handling of other errors? or none
    #     print(f'Environment file not found in {file_path} and not uploaded.')
    # return ds_hdf
    
def species_by_num_to_pct(props):
    import numpy as np
    species_by_num = eval(props['chem_species_by_n_atoms'])
    species_by_pct = {at: np.round(species_by_num[at]*100/props['n_atoms_total'], 2) for at in species_by_num.keys()}
    return str(species_by_pct)




# flat_cdict:
#  {'temperature': None, 'pressure': None, 
#   'ensemble': None, 
#   'potential': '2022--Sun-Y--Fe--LAMMPS--ipr1',  
#   'AverageTotalEnergy': -8.0433, 
#   'AverageTotalVolume': 23.0275, 
#   'software': 'LAMMPS 2018.03.16', 'project_name': 'test_new', 
#   'job_type': 'Lammps', 'host': 'NB4267'}