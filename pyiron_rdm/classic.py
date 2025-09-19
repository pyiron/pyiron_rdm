def classic_structure(
    pr, structure, structure_name, options, is_init_struct: bool, init_structure=None
):
    # TODO rename is_init_struct and init_structure to better reflect that it needs to not be manipulated (e.g. repeated)
    structure_path = pr.path
    # structure_name = structure_name_prefix + '_input_structure'

    from pyiron_base.storage.hdfio import FileHDFio

    hdf = FileHDFio(structure_path + structure_name + ".h5")
    structure.to_hdf(hdf)

    from pyiron_rdm.concept_dict import (
        get_unit_cell_parameters,
        process_structure_crystal,
    )

    if is_init_struct:
        init_structure = structure
    if init_structure:
        try:
            struct_params = get_unit_cell_parameters(init_structure)
        except:
            struct_params = {}
    else:
        struct_params = {}
    struct_cdict = process_structure_crystal(
        pr, structure, structure_name, structure_path, struct_params, options
    )

    return struct_cdict


def classic_general_job(job, export_env_file: bool):
    from pyiron_rdm.concept_dict import process_general_job

    if export_env_file:
        from pyiron_rdm.concept_dict import export_env

        export_env(job.path)

    job_cdict = process_general_job(job)
    return job_cdict


def classic_lammps(lammps_job, export_env_file):
    """export_env_file: Bool"""
    from pyiron_rdm.concept_dict import process_lammps_job

    if export_env_file:
        from pyiron_rdm.concept_dict import export_env

        export_env(lammps_job.path)

    lammps_cdict = process_lammps_job(lammps_job)
    return lammps_cdict


def classic_vasp(vasp_job, export_env_file):
    """export_env_file: Bool"""
    from pyiron_rdm.concept_dict import process_vasp_job

    if export_env_file:
        from pyiron_rdm.concept_dict import export_env

        export_env(vasp_job.path)

    vasp_cdict = process_vasp_job(vasp_job)
    return vasp_cdict


def classic_murn(murn_job, export_env_file):
    # TODO this assumes only lammps child jobs - generalise!
    """export_env_file: Bool"""

    if export_env_file:
        import shutil

        from pyiron_rdm.concept_dict import export_env

        export_env(murn_job.path)

    from pyiron_rdm.concept_dict import (
        process_lammps_job,
        process_murnaghan_job,
        process_vasp_job,
    )

    child_jobs_cdict = []
    for job in murn_job.iter_jobs():
        if export_env_file:
            shutil.copy(
                murn_job.path + "_environment.yml", job.path + "_environment.yml"
            )
        if (
            "lammps" in job.to_dict()["TYPE"]
        ):  # if it's not possible to have multiple types, do this only once
            child_cdict = process_lammps_job(job)
        else:  # vasp
            child_cdict = process_vasp_job(job)
        child_jobs_cdict.append(child_cdict)

    job_cdict = process_murnaghan_job(murn_job)

    return job_cdict, child_jobs_cdict


def classic_murn_equil_structure(
    murn_job, options, is_init_struct: bool = True, init_structure=None
):

    if is_init_struct:
        init_structure = murn_job.structure
    equil_structure = murn_job.get_structure()
    structure_name = murn_job.name + "_equilibrium_structure"
    struct_cdict = classic_structure(
        murn_job.project,
        equil_structure,
        structure_name,
        options,
        is_init_struct,
        init_structure,
    )

    return struct_cdict


def get_datamodel(o):
    # TODO: adapt to also take url to use for openbis_login ?
    datamodels = {"bam": "bam", "imm.rwth": "sfb1394"}
    for key in datamodels:
        if key in o.hostname:
            return datamodels[key]
    raise KeyError(
        f"The {o.hostname} openBIS hostname is not paired with a supported data model yet ({datamodels.values()})."
    )


def validate_upload_options(o, options):
    import importlib

    options_cfg = importlib.import_module(o.ot)
    allowed_keys = options_cfg.allowed_keys
    allowed_defects = getattr(options_cfg, "allowed_defects", None)

    # currently not validating that materials/pseudopot is a string/list of strings
    invalid_keys = set(options) - allowed_keys
    if invalid_keys:
        raise KeyError(
            f'Unsupported key(s) in "options" dictionary: {sorted(invalid_keys)}. \
                       Allowed keys are: {sorted(allowed_keys)}'
        )
    defects = options.get("defects")
    if defects:
        if not isinstance(defects, list) or not all(
            isinstance(d, str) for d in defects
        ):
            raise TypeError('options["defects"] must be a list of strings.')
        if allowed_defects:
            invalid_defects = (
                set([d.lower().replace("_", " ") for d in defects]) - allowed_defects
            )
            if invalid_defects:
                raise ValueError(
                    f"Invalid defect(s) in \"options['defects']\": {sorted(invalid_defects)}. \
                    Allowed defects are: {sorted(allowed_defects)}"
                )
    materials = options.get("materials")
    if materials:
        if not isinstance(materials, list):
            options["materials"] = [options["materials"]]
    pseudopots = options.get("pseudopotentials")
    if pseudopots:
        if not isinstance(pseudopots, list):
            options["pseudopotentials"] = [options["pseudopotentials"]]
    return options


def openbis_login(url, username, instance="bam", s3_config_path=None):
    if instance != "bam" and instance != "sfb1394":
        raise ValueError(
            f"This script only supports upload to 'bam' and 'sfb1394' instances,\
                         {instance} not supported."
        )

    if instance == "bam":
        mapping_path = "pyiron_rdm.ob_cfg_bam"
        OT_path = "pyiron_rdm.ob_OT_bam"
        s3_config_path = None
    elif instance == "sfb1394":
        mapping_path = "pyiron_rdm.ob_cfg_sfb1394"
        OT_path = "pyiron_rdm.ob_OT_sfb1394"
        if not s3_config_path:
            s3_config_path = "test_sfb.cfg"

    from pyiron_rdm.ob_upload import openbis_login as ob_login

    o = ob_login(url, username, s3_config_path, mapping_path, OT_path)
    return o


def upload_classic_pyiron(
    job,
    o,
    space,
    project,
    collection=None,
    export_env_file=True,
    is_init_struct: bool = True,
    init_structure=None,
    options=None,
):
    # TODO should this return anything?

    # check options keys
    if options is not None:
        options = validate_upload_options(o, options)
    else:
        options = {}

    structure = job.structure
    if not structure:
        print(
            "This job does not contain a structure and will not be uploaded. \
                Please add structure before trying to upload."
        )
        return

    # Project env file - TODO what is this for??
    pr = job.project
    if export_env_file:
        from pyiron_rdm.concept_dict import export_env

        export_env(pr.path + pr.name)

    if not collection:
        collection = pr.name
    space = space.upper()
    project = project.upper()
    collection = collection.upper()

    # ------------------------------------VALIDATION----------------------------------------------
    cdicts_to_validate = []

    struct_dict = classic_structure(
        pr,
        structure,
        structure_name=job.name + "_structure",
        options=options,
        is_init_struct=is_init_struct,
        init_structure=init_structure,
    )
    cdicts_to_validate.append(struct_dict)

    job_type = job.to_dict()["TYPE"]
    proceed = True
    if "lammps" in job_type:
        job_cdict = classic_lammps(job, export_env_file=export_env_file)
        cdicts_to_validate.append(job_cdict)

    elif "vasp" in job_type:
        job_cdict = classic_vasp(job, export_env_file=export_env_file)
        cdicts_to_validate.append(job_cdict)

    elif "murn" in job_type:
        job_cdict, child_jobs_cdict = classic_murn(job, export_env_file=export_env_file)
        equil_struct_dict = classic_murn_equil_structure(
            job, options, is_init_struct, init_structure
        )
        cdicts_to_validate.append(job_cdict)
        cdicts_to_validate.append(equil_struct_dict)
        cdicts_to_validate += [child_cdict for child_cdict in child_jobs_cdict]

    else:
        print(f"The {job_type} job type is not implemented for OpenBIS upload yet.")
        proceed = input(
            "Type 'yes' to proceed with an upload to general pyiron job type."
        )
        if proceed.lower() == "yes" or proceed.lower() == "y":
            job_cdict = classic_general_job(job, export_env_file=export_env_file)
            cdicts_to_validate.append(job_cdict)
        else:
            print("Upload cancelled.")
            proceed = False

    datamodel = get_datamodel(o)
    upload_final_struct = datamodel == "sfb1394"

    if upload_final_struct and (not "murn" in job.to_dict()["TYPE"]):
        if is_init_struct:
            init_structure = structure
        final_structure = job.get_structure()
        final_struct_dict = classic_structure(
            pr,
            final_structure,
            structure_name=job.name + "_final_structure",
            options=options,
            is_init_struct=False,
            init_structure=init_structure,
        )
        cdicts_to_validate.append(final_struct_dict)

    from pyiron_rdm.ob_upload import openbis_validate

    validated_to_upload = openbis_validate(
        o, space, project, collection, cdicts_to_validate, options
    )
    if not proceed:
        raise ValueError("You asked to abort the upload.")
    # ---------------------------------------------------------------------------------------------

    # --------------------------------------UPLOAD-------------------------------------------------
    from pyiron_rdm.ob_upload import openbis_upload_validated

    # Structure
    cdict, props_dict, object_type, ds_types, ob_parents, object_name = (
        validated_to_upload[0]
    )
    ob_structure_id = openbis_upload_validated(
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
    )

    if datamodel == "sfb1394":
        job_parents = None  # job does not have init structure as parent
        str_parent = ob_structure_id  # equil structure has init as parent
    elif datamodel == "bam":
        job_parents = ob_structure_id  # job has init structure as parent
        str_parent = None  # equil structure does not have init as parent

    # (Main) job
    cdict, props_dict, object_type, ds_types, ob_parents, object_name = (
        validated_to_upload[1]
    )
    ob_job_id = openbis_upload_validated(
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
        parent_ids=job_parents,
    )

    if "murn" in job_type:
        ob_children_ids = []
        # Murn equilibrium structure
        cdict, props_dict, object_type, ds_types, ob_parents, object_name = (
            validated_to_upload[2]
        )
        ob_equil_struct_id = openbis_upload_validated(
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
            parent_ids=str_parent,
        )
        ob_children_ids.append(ob_equil_struct_id)
        # Murn children jobs
        for validated_child in validated_to_upload[3:]:
            cdict, props_dict, object_type, ds_types, ob_parents, object_name = (
                validated_child
            )
            ob_child_id = openbis_upload_validated(
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
                parent_ids=str_parent,
            )
            ob_children_ids.append(ob_child_id)
        from pyiron_rdm.ob_upload import link_children

        link_children(o, ob_job_id, ob_children_ids)

    # Final structure upload (already included as equilibrium for murn)
    elif upload_final_struct:
        cdict, props_dict, object_type, ds_types, ob_parents, object_name = (
            validated_to_upload[-1]
        )
        ob_final_structure_id = openbis_upload_validated(
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
            parent_ids=[ob_structure_id, ob_job_id],
        )


# ---------------------------------------------------------------------------------------------
