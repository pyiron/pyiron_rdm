def classic_structure(pr, structure, structure_name):
    structure_path = pr.name + '/'
    # structure_name = structure_name_prefix + '_input_structure'
    
    from pyiron_base.storage.hdfio import FileHDFio
    hdf = FileHDFio(structure_path + structure_name + '.h5')
    structure.to_hdf(hdf)

    # Cannot guarantee this will be before manipulations hence skipping get_unit_cell_parameters
    from ob.concept_dict import process_structure_crystal, get_unit_cell_parameters
    try:
        struct_params = get_unit_cell_parameters(structure)
    except:
        struct_params = {}
    struct_cdict = process_structure_crystal(pr, structure, structure_name, structure_path, struct_params)
    
    return struct_cdict

def classic_lammps(lammps_job):
    from ob.concept_dict import process_lammps_job, export_env

    export_env(lammps_job.path)

    lammps_cdict = process_lammps_job(lammps_job)
    return lammps_cdict

def classic_murn(murn_job):
    from ob.concept_dict import export_env
    export_env(murn_job.path)

    from ob.concept_dict import process_murnaghan_job, process_lammps_job
    child_jobs_cdict = []
    for jobs in murn_job.iter_jobs():
        # from pyiron_base.storage.hdfio import FileHDFio
        import platform
        if "Windows" in platform.system():
            import shutil
            shutil.copy(murn_job.path + '_environment.yml', jobs.path + '_environment.yml')
        else:
            import os
            os.system('cp ' + murn_job.path + '_environment.yml ' + jobs.path + '_environment.yml')
        child_cdict = process_lammps_job(jobs)
        child_jobs_cdict.append(child_cdict)

    job_cdict = process_murnaghan_job(murn_job)

    return job_cdict, child_jobs_cdict

def classic_murn_equil_structure(murn_job):

    from pyiron_base.storage.hdfio import FileHDFio
    from ob.concept_dict import process_structure_crystal, get_unit_cell_parameters

    equil_structure = murn_job.get_structure()
    structure_name = murn_job.name + '_equilibrium_structure'
    structure_path = murn_job.project.name + '/'
    hdf_equil = FileHDFio(structure_path + structure_name + '.h5')
    equil_structure.to_hdf(hdf_equil)

    # Must be conventional unit cell or primitive cell
    # Provide the structure before repetition or any other manipulations
    structure_parameters = get_unit_cell_parameters(equil_structure)
    
    struct_cdict = process_structure_crystal(murn_job.project, equil_structure, structure_name, structure_path, structure_parameters)
        
    return struct_cdict

def openbis_login(url, username, instance='bam'):
    # TODO this shouldn't be needed, default instance like a default queue?
    if instance != 'bam' and instance != 'sfb1394':    
        raise ValueError(f"This script only supports upload to 'bam' and 'sfb1394' instances,\
                         {instance} not supported.")
    
    if instance == 'bam':
        mapping_path = 'ob.ob_cfg_bam'
        OT_path = 'ob.ob_OT_bam'
        s3_config_path = None
    elif instance == 'sfb1394':
        mapping_path = 'ob.ob_cfg_sfb1394'
        OT_path = 'ob.ob_OT_sfb1394'
        s3_config_path = "test_sfb.cfg"

    from ob.ob_upload import openbis_login as ob_login

    o = ob_login(url, username, s3_config_path, mapping_path, OT_path)
    return o

def upload_classic_pyiron(job, o, space, project, collection=None):
    # TODO should this return anything?

    structure = job.structure
    if structure:

        # Project env file - TODO what is this for??
        pr = job.project
        from ob.concept_dict import export_env
        export_env(pr.path + pr.name)

        from ob.ob_upload import openbis_upload

        if not collection:
            collection = pr.name

        struct_dict = classic_structure(pr, structure, structure_name=job.name + '_structure')
        ob_structure_id = openbis_upload(o, space, project, collection, struct_dict)

        if 'lammps' in job.to_dict()['TYPE']:
            job_cdict = classic_lammps(job)
            ob_job_id = openbis_upload(o, space, project, collection, job_cdict, parent_ids=ob_structure_id)

        elif 'murn' in job.to_dict()['TYPE']: # TODO: Is it okay to upload the murn job last ?
            job_cdict, child_jobs_cdict = classic_murn(job)
            ob_children_ids = []
            for child_cdict in child_jobs_cdict:
                ob_child_id = openbis_upload(o, space, project, collection, child_cdict)
                ob_children_ids.append(ob_child_id)
            equil_struct_dict = classic_murn_equil_structure(job)
            ob_equil_struct_id = openbis_upload(o, space, project, collection, equil_struct_dict)
            ob_children_ids.append(ob_equil_struct_id)
            
            ob_murn_id = openbis_upload(o, space, project, collection, job_cdict, parent_ids=ob_structure_id)
            from ob.ob_upload import link_children
            link_children(o, ob_murn_id, ob_children_ids)

        else:
            job_type = job.to_dict()['TYPE']
            print(f'The {job_type} job type is not implemented for OpenBIS upload yet.')
        
    else:
        print('This job does not contain a structure and will not be uploaded. \
                Please add structure before trying to upload.')