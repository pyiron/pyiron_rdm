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

    import importlib
    module_spec = importlib.util.find_spec(o.mapping)
    if not module_spec:
        raise FileNotFoundError(f"There is no openBIS '{o.mapping}' mapping file in your pyiron resources. "
                                "Please correct this before upload.")  # TODO correct the o.mapping here to just filename after pyiron_resources thing set up
    module_spec = importlib.util.find_spec(o.ot)
    if not module_spec:
        raise FileNotFoundError(f"There is no openBIS '{o.ot}' object type file in your pyiron resources. "
                                "Please correct this before upload.")  # TODO correct the o.ot here to just filename after pyiron_resources thing set up

    return o

def openbis_validate(o, space, project, collection, 
                     concept_dicts: dict|list):
    if isinstance(concept_dicts, dict):
        concept_dicts = [concept_dicts]
    all_issues = []
    outputs = []
    for concept_dict in concept_dicts:
        issues = validate_ob_destination(o, space, project, collection)
        from ob.concept_dict import flatten_cdict
        cdict = flatten_cdict(concept_dict)
        import importlib
        ob_ot = importlib.import_module(o.ot).get_ot_info(cdict)
        object_type, ds_types, inv_parents = ob_ot()
        map_cdict_to_ob = importlib.import_module(o.mapping).map_cdict_to_ob
        props_dict = map_cdict_to_ob(o, cdict, concept_dict)
        inv_issues, ob_parents = validate_inventory_parents(o, inv_parents, cdict, props_dict)
        issues += inv_issues
        object_name = concept_dict['job_details'][0]['value']

        if issues:
            issue_message = f'{object_name} object:\n' + '\n'.join(f'- {issue}' for issue in issues)
            all_issues.append(issue_message)
            outputs.append(())
        else:
            outputs.append((cdict, props_dict, object_type, ds_types, ob_parents, object_name))
        
    if all_issues:
        raise ValueError('The following issues were found and need to be fixed before upload:\n\n' +
                         '\n\n'.join(all_issues))
        # TODO perhaps more suitable exception? Custom ValidationError(Exception)?
    return outputs

def openbis_upload(o, space, project, collection, concept_dict: dict, parent_ids=None):
    '''Currently assumes a single concept_dict, not a list of them'''
    cdict, props_dict, object_type, ds_types, ob_parents, object_name = openbis_validate(o, space, project, collection, concept_dict)[0]
    object_id = openbis_upload_validated(o, space, project, collection, object_name, 
                                  object_type, ob_parents, props_dict, ds_types, cdict, parent_ids)
    return object_id

def openbis_upload_validated(o, space, project, collection, object_name, 
                             object_type, ob_parents, props_dict, ds_types, cdict,
                             parent_ids=None):
    # TODO or allow to skip the last two a flag later?
    # TODO also type checking on values

    if 'S3' in str(type(o)):
        kind = 'LINK'
    else:
        kind = 'PHYSICAL'

    ob_coll = o.get_collection(f'/{space}/{project}/{collection}')
    found_objects_ids = [obj.identifier for obj in ob_coll.get_objects() if obj.p['$name'] == object_name]

    if found_objects_ids:
        print("===================\n")
        print(f"Object with name {object_name} already exists! Found object(s) in: {found_objects_ids}\n")
        print("===================\n")

        return found_objects_ids   
            # TODO: should return a single id to match 'else' return format, however, multiple objects of the same name possible
    
    else:
        object_ = o.new_object(
            type       = object_type,
            space      = space,
            experiment = ob_coll,
            parents    = ob_parents,
            props      = props_dict
        )
        object_.save()

        import importlib
        for ds in ds_types:
            if ds == 'job_h5':
                dataset_info = importlib.import_module(o.mapping).dataset_job_h5
                file_path = cdict['path'] + '.h5'
            elif ds == 'structure_h5':
                dataset_info = importlib.import_module(o.mapping).dataset_atom_struct_h5
                file_path = cdict['path'] + cdict['structure_name'] + '.h5'
            elif ds == 'env_yml':
                dataset_info = importlib.import_module(o.mapping).dataset_env_yml
                file_path = cdict['path'] + '_environment.yml'
            elif ds == 'cdict_json':
                dataset_info = importlib.import_module(o.mapping).dataset_cdict_jsonld
                if 'structure_name' in cdict.keys():
                    file_path = cdict['path'] + cdict['structure_name'] + '_concept_dict.json'
                else:
                    file_path = cdict['path'] + '_concept_dict.json'
            else:
                raise ValueError(f'Dataset type {ds} not recognised. Supported datasets: job_h5, structure_h5, env_yml, cdict_json.')

            ds_type, ds_props = dataset_info(cdict)
            try:
                upload_dataset(o, object_, ds_type, ds_props, file_path, kind)
            except (ValueError, FileNotFoundError) as e: # pybis: ValueError; OpenbisAixTended - shutil: FileNotFoundError
                if ds == 'env_yml':
                    import warnings
                    warnings.warn('The environment file was not uploaded.')
                else:
                    raise e

        if parent_ids:
            link_parents(o, object_, parent_ids)

        return object_.identifier # or object_.permId
    
def link_parents(o, ob_object, parent_ids):
    if isinstance(ob_object, str):
        ob_object = o.get_object(ob_object)
    if not isinstance(parent_ids, list):      # if string not list
        parent_ids = [parent_ids]
    for p in parent_ids:
        try:
            ob_object.add_parents(p)
        except ValueError:
            print(f'An object with the identifier {p} not found and hence not linked as parent.')
    ob_object.save()

def link_children(o, ob_object, children_ids):
    if isinstance(ob_object, str):
        ob_object = o.get_object(ob_object)
    if not isinstance(children_ids, list):    # if string not list
        children_ids = [children_ids]
    for ch in children_ids:
        try:
            ob_object.add_children(ch)
        except ValueError:
            print(f'An object with the identifier {ch} not found and hence not linked as child.')
    ob_object.save()

def upload_dataset(o, ob_object, ds_type, ds_props, file_path, kind):
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

def validate_ob_destination(o, space, project, collection):
    try:
        o.get_space(space)
    except ValueError as e:
        return [str(e) + ' Project and collection will not be checked until space created/found.']
    try:
        o.get_project(f'/{space}/{project}')
    except ValueError as e:
        create_project = input(f"Project with the code {project} was not found in space {space}. Type 'yes' to create it.")
        if create_project.lower() in ['yes', 'y']:
            new_project = o.new_project(space=space, code=project)
            new_project.save()
        else:
            return([str(e) + ' and you chose not to create it. Collection will not be checked until project created/found.'])
    try:
        o.get_collection(f'/{space}/{project}/{collection}')
    except ValueError as e:
        create_coll = input(f"Collection with the code {collection} was not found in project {space}/{project}. Type 'yes' to create it.")
        if create_coll.lower() in ['yes', 'y']:
            new_coll = o.new_collection(project=f'/{space}/{project}', code=collection, type='DEFAULT_EXPERIMENT')
            new_coll.save()
        else:
            return [str(e) + ' and you chose not to create it.']
    return []

def validate_inventory_parents(o, inv_parents, cdict, props_dict):
    import importlib
    issues = []
    ob_parents = []
    get_inv_parent = importlib.import_module(o.ot).get_inv_parent
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
            else:
                issues.append(f'Parent object not found: No objects of the type {t} and property {w} / code "{c}" in inventory.')
        else:
            issues.append(f'Parent object not found: Not enough information to search. Known information: type = {t}, code = "{c}", attribute match: {w}')

    return issues, ob_parents