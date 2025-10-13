# TODO better imports for multiple openbis instances

import warnings


def openbis_login(
    url: str,
    mapping_path: str,
    OT_path: str,
    username: str | None = None,
    password: str | None = None,
    token: str | None = None,
    s3_config_path: str | None = None,
):
    if username is None and token is None:
        raise ValueError("Either username or token must be provided.")

    from getpass import getpass

    if s3_config_path:
        from pybis_aixtended.OpenbisAixTended import OpenbisWithS3

        o = OpenbisWithS3(
            url=url,
            verify_certificates=True,
            s3_config_path=s3_config_path,
            standardize_filenames=True,
        )
    else:
        from pybis import Openbis

        o = Openbis(url)

    if not o.is_session_active():
        if token is None:
            if password is None:
                password = getpass("Enter openBIS password: ")
            o.login(
                username, password, save_token=True
            )  # save the session token in ~/.pybis/example.com.token
        else:
            o.set_token(token)

    o.mapping = mapping_path
    o.ot = OT_path

    import importlib

    module_spec = importlib.util.find_spec(o.mapping)
    if not module_spec:
        raise FileNotFoundError(
            f"There is no openBIS '{o.mapping}' mapping file in your pyiron resources. "
            "Please correct this before upload."
        )  # TODO correct the o.mapping here to just filename after pyiron_resources thing set up
    module_spec = importlib.util.find_spec(o.ot)
    if not module_spec:
        raise FileNotFoundError(
            f"There is no openBIS '{o.ot}' object type file in your pyiron resources. "
            "Please correct this before upload."
        )  # TODO correct the o.ot here to just filename after pyiron_resources thing set up

    return o


def openbis_validate(
    o, concept_dicts: dict, options: dict, require_parents: bool = True
) -> list:
    outputs = {}
    for key, concept_dict in concept_dicts.items():
        from pyiron_rdm.concept_dict import flatten_cdict

        cdict = flatten_cdict(concept_dict)
        import importlib

        object_type, ds_types, inv_parents = importlib.import_module(o.ot).get_ot_info(
            cdict
        )
        map_cdict_to_ob = importlib.import_module(o.mapping).map_cdict_to_ob
        props_dict = map_cdict_to_ob(
            user_name=o.get_session_info().userName,
            cdict=cdict,
            concept_dict=concept_dict,
        )
        ob_parents = validate_inventory_parents(
            o, inv_parents, cdict, props_dict, options, require_parents=require_parents
        )
        object_name = concept_dict["job_details"][0]["value"]

        outputs[key] = {
            "cdict": cdict,
            "props_dict": props_dict,
            "object_type": object_type,
            "ds_types": ds_types,
            "ob_parents": ob_parents,
            "object_name": object_name,
        }

    return outputs


def openbis_upload_validated(
    o,
    space,
    project,
    collection,
    object_name,
    object_type,
    ob_parents,
    props_dict,
    ds_types,
    cdict,
    parent_ids=None,
):
    # TODO or allow to skip the last two a flag later?
    # TODO also type checking on values

    if "S3" in str(type(o)):
        kind = "LINK"
    else:
        kind = "PHYSICAL"

    ob_coll = o.get_collection(f"/{space}/{project}/{collection}")
    found_objects_ids = [
        obj.identifier for obj in ob_coll.get_objects() if obj.p["$name"] == object_name
    ]

    if found_objects_ids:
        print("===================\n")
        print(
            f"Object with name {object_name} already exists! Found object(s) in: {found_objects_ids}\n"
        )
        print("===================\n")

        return found_objects_ids
        # TODO: should return a single id to match 'else' return format, however, multiple objects of the same name possible

    object_ = o.new_object(
        type=object_type,
        space=space,
        experiment=ob_coll,
        parents=ob_parents,
        props=props_dict,
    )
    object_.save()

    from importlib import import_module

    for ds in ds_types:
        if ds == "job_h5":
            ds_type, ds_props = import_module(o.mapping).dataset_job_h5(cdict)
            file_path = cdict["path"] + ".h5"
        elif ds == "structure_h5":
            ds_type, ds_props = import_module(o.mapping).dataset_atom_struct_h5(cdict)
            file_path = cdict["path"] + cdict["structure_name"] + ".h5"
        elif ds == "env_yml":
            ds_type, ds_props = import_module(o.mapping).dataset_env_yml(cdict)
            file_path = cdict["path"] + "_environment.yml"
        elif ds == "cdict_json":
            ds_type, ds_props = import_module(o.mapping).dataset_cdict_jsonld(cdict)
            if "structure_name" in cdict.keys():
                file_path = (
                    cdict["path"] + cdict["structure_name"] + "_concept_dict.json"
                )
            else:
                file_path = cdict["path"] + "_concept_dict.json"
        else:
            raise ValueError(
                f"Dataset type {ds} not recognised. Supported datasets: job_h5,"
                " structure_h5, env_yml, cdict_json."
            )

        try:
            upload_dataset(
                o=o,
                object_=object_,
                ds_type=ds_type,
                ds_props=ds_props,
                file_path=file_path,
                kind=kind,
            )
        except (
            ValueError,
            FileNotFoundError,
        ) as e:  # pybis: ValueError; OpenbisAixTended - shutil: FileNotFoundError
            if ds == "env_yml":

                warnings.warn("The environment file was not uploaded.")
            else:
                raise e

    if parent_ids:
        link_parents(
            o=o,
            object_=object_,
            parent_ids=parent_ids,
        )

    return object_.identifier  # or object_.permId


def link_parents(o, ob_object, parent_ids):
    if isinstance(ob_object, str):
        ob_object = o.get_object(ob_object)
    if not isinstance(parent_ids, list):
        parent_ids = [parent_ids]
    for p in parent_ids:
        try:
            ob_object.add_parents(p)
        except ValueError:
            print(
                f"An object with the identifier {p} not found and hence not linked as parent."
            )
    ob_object.save()


def link_children(o, ob_object, children_ids):
    if isinstance(ob_object, str):
        ob_object = o.get_object(ob_object)
    if not isinstance(children_ids, list):
        children_ids = [children_ids]
    for ch in children_ids:
        try:
            ob_object.add_children(ch)
        except ValueError:
            print(
                f"An object with the identifier {ch} not found and hence not linked as child."
            )
    ob_object.save()


def upload_dataset(o, ob_object, ds_type, ds_props, file_path, kind):
    ds_hdf = o.new_dataset(
        type=ds_type,
        object=ob_object,
        files=[file_path],
        kind=kind,
        props=ds_props,
    )
    ds_hdf.save()


def validate_ob_destination(o, space: str, project: str, collection: str):
    try:
        o.get_space(space)
    except ValueError as e:
        return [f"{e}; available spaces: {[s.code for s in o.get_spaces()]}"]
    try:
        o.get_project(f"/{space}/{project}")
    except ValueError as e:
        create_project = input(
            f"Project with the code {project} was not found in space {space}."
            " Type 'yes' to create it."
        )
        if create_project.lower() in ["yes", "y"]:
            new_project = o.new_project(space=space, code=project)
            new_project.save()
        else:
            raise ValueError(
                f"{e}; available projects in space {space}:"
                f" {[p.code for p in o.get_projects(space=space)]}"
            )
    try:
        o.get_collection(f"/{space}/{project}/{collection}")
    except ValueError as e:
        create_coll = input(
            f"Collection with the code {collection} was not found in project"
            f" {space}/{project}. Type 'yes' to create it."
        )
        if create_coll.lower() in ["yes", "y"]:
            new_coll = o.new_collection(
                project=f"/{space}/{project}",
                code=collection,
                type="DEFAULT_EXPERIMENT",
            )
            new_coll.save()
        else:
            raise ValueError(
                f"{e}; available collections in project {space}/{project}:"
                f" {[c.code for c in o.get_collections(project=f'/{space}/{project}')]} "
            )


def validate_inventory_parents(
    o, inv_parents, cdict, props_dict, options, require_parents: bool = True
):
    import importlib

    ob_parents = []
    get_inv_parent = importlib.import_module(o.ot).get_inv_parent
    issues = []
    for inv_parent in inv_parents:
        ob_type, permids, where, attrs, code = get_inv_parent(
            inv_parent, cdict, props_dict, options
        )
        if permids is not None:
            permids = [o.get_object(p).permId for p in permids]
        if permids:  # multiple parents allowed when more permIds provided
            parents = o.get_objects(
                type=ob_type,
                permId=permids,
            )

            if parents:
                ob_parents += parents
            else:
                issues.append(
                    f"Parent object not found: No objects of the type {ob_type}"
                    f" and permId '{permids}' in inventory."
                )
        elif code or where:  # single parent allowed otherwise; taking the first
            parent = o.get_objects(
                type=ob_type, permId=permids, code=code, where=where, attrs=attrs
            )[0]
            if not parent:
                issues.append(
                    f"Parent object not found: No objects of the type {ob_type}"
                    f' and property {where} / code "{code}" in inventory.'
                )
            ob_parents.append(parent)
        elif ob_type == "PSEUDOPOTENTIAL":  # don't fail when pseudopotential missing
            warnings.warn(
                "Pseudopotential permId was not given. Please link job to"
                " appropriate pseudopotential on openBIS."
            )
        else:
            issues.append(
                "Parent object not found: Not enough information to search."
                f' Known information: type = {ob_type}, permId = "{permids}",'
                f' code = "{code}", attribute match: {where}'
            )
            ob_parents.append(None)
    if issues and require_parents:
        raise ValueError(" \n".join(issues))
    else:
        for issue in issues:
            warnings.warn(issue)
    return ob_parents
