# TODO better imports for multiple openbis instances

def openbis_login(url, username, s3_config_path=None, mapping_path=None, OT_path=None):
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
    
    o.mapping = mapping_path
    o.ot = OT_path

    return o

def openbis_upload(o, space, project, collection, concept_dict, parent_ids=None):
    # TODO or allow to skip the last two a flag later?
    # TODO also type checking on values
    import importlib
    module_spec = importlib.util.find_spec(o.mapping)
    if not module_spec:
        raise FileNotFoundError(f"There is no openBIS '{o.mapping}' mapping file in your pyiron resources. \
                                Please correct this before upload.")  # TODO correct the o.mapping here to just filename after pyiron_resources thing set up
    module_spec = importlib.util.find_spec(o.ot)
    if not module_spec:
        raise FileNotFoundError(f"There is no openBIS '{o.ot}' object type file in your pyiron resources. \
                                Please correct this before upload.")  # TODO correct the o.ot here to just filename after pyiron_resources thing set up

    # req_keys_ob_info = ['object_type', 'parents', 'datasets'] 
    # if not all(k in ob_info.keys() for k in req_keys_ob_info):
    #     raise KeyError('The ob_info dictionary must contain all following keys: object_type, parents, datasets.')
    
    if 'S3' in str(type(o)):
        kind = 'LINK'
    else:
        kind = 'PHYSICAL'

    ob_coll = '/' + space + '/' + project + '/' + collection
    ob_project_obj = o.get_project('/' + space + '/' + project)

    from ob.concept_dict import flatten_cdict
    cdict = flatten_cdict(concept_dict)
    ob_ot = importlib.import_module(o.ot).get_ot_info(cdict)

    object_type, ds_types, inv_parents = ob_ot()
    objects = ob_project_obj.get_objects(type = object_type)
    exists = False
    for object_ in objects:
        object_name = concept_dict['job_details'][0]['value'] # TODO do better !!!
        if object_.p.get('$name') == object_name:
            exists = True
            found_object = object_

    if exists:
        print("===================\n")
        print(f"Object already exists! Found object in: {found_object.identifier}\n")
        print("===================\n")
        print("Found job properties:\n")
        from IPython.display import display
        display(found_object.p)
        return found_object.identifier
    
    else:
        map_cdict_to_ob = importlib.import_module(o.mapping).map_cdict_to_ob
        get_inv_parent = importlib.import_module(o.ot).get_inv_parent
        props_dict = map_cdict_to_ob(o, cdict, concept_dict)
        ob_parents = []
        for inv_parent in inv_parents:
            t, w, a, c = get_inv_parent(inv_parent, cdict, props_dict)
            if c or w:
                parent = o.get_objects(
                    type = t,
                    code = c,
                    where = w,
                    attrs = a
                )[0]

                if parent:
                    ob_parents.append(parent)
                else:                              # TODO proper error or just a warning?
                    print(f"No objects of the type {t} and property {w} / code {c} found, \
                          upload will not proceed. Please create them first and then try again.")
                    return
            else:
                print(f"Not enough information to search for a parent object.\
                      Known information: type = {t}, code = {c}, attribute match: {w}")
                return

        if len(ob_parents) == len(inv_parents): # Found all parents needed
            object_ = o.new_object(
                type       = object_type,
                space      = space,
                experiment = ob_coll,
                parents    = ob_parents,
                props      = props_dict
            )
            object_.save()

            for ds in ds_types:
                if ds == 'job_h5':
                    # from ob.ob_cfg_bam import dataset_job_h5 as dataset_info
                    dataset_info = importlib.import_module(o.mapping).dataset_job_h5
                    file_path = cdict['path'] + '.h5'
                elif ds == 'structure_h5':
                    # from ob.ob_cfg_bam import dataset_atom_struct_h5 as dataset_info
                    dataset_info = importlib.import_module(o.mapping).dataset_atom_struct_h5
                    file_path = cdict['path'] + cdict['structure_name'] + '.h5'
                elif ds == 'env_yml':
                    # from ob.ob_cfg_bam import dataset_env_yml as dataset_info
                    dataset_info = importlib.import_module(o.mapping).dataset_env_yml
                    file_path = cdict['path'] + '_environment.yml'
                elif ds == 'cdict_json':
                    # from ob.ob_cfg_bam import dataset_cdict_jsonld as dataset_info
                    dataset_info = importlib.import_module(o.mapping).dataset_cdict_jsonld
                    if 'structure_name' in cdict.keys():
                        file_path = cdict['path'] + cdict['structure_name'] + '_concept_dict.json'
                    else:
                        file_path = cdict['path'] + '_concept_dict.json'
                else:
                    raise ValueError(f'Dataset type {ds} not recognised. Supported datasets: job_h5, structure_h5, env_yml, cdict_json.')

                ds_type, ds_props = dataset_info(cdict)
                upload_dataset(o, object_, ds_type, ds_props, file_path, kind)

            if parent_ids:
                link_parents(o, object_, parent_ids)

        # if show_object:
        #     from IPython.display import display
        #     display(object_.p)
        
        # return object_
            return object_.identifier # or object_.permId
    
def link_parents(o, ob_object, parent_ids):
    if type(ob_object) == str:
        ob_object = o.get_object(ob_object)
    if type(parent_ids) != list: # if string not list
        parent_ids = [parent_ids]
    for p in parent_ids:
        try:
            ob_object.add_parents(p)
        except ValueError:
            print(f'An object with the identifier {p} not found and hence not linked as parent.')
    ob_object.save()

def link_children(o, ob_object, children_ids):
    if type(ob_object) == str:
        ob_object = o.get_object(ob_object)
    if type(children_ids) != list: # if string not list
        children_ids = [children_ids]
    for ch in children_ids:
        try:
            ob_object.add_children(ch)
        except ValueError:
            print(f'An object with the identifier {ch} not found and hence not linked as child.')
    ob_object.save()

def upload_dataset(o, ob_object, ds_type, ds_props, file_path, kind):
    # try:
    ds_hdf = o.new_dataset(
        type       = ds_type,
        # collection = collection,
        object     = ob_object,
        files      = [file_path],
        kind       = kind,
        props      = ds_props
    )
    ds_hdf.save()
    # except ValueError: # TODO handling of other errors? or none
    #     print(f'Environment file not found in {file_path} and not uploaded.')
    # return ds_hdf