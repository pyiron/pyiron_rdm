import os
import json
from pathlib import Path
import shutil

from ob.classic import openbis_login
from ob.ob_upload import openbis_validate
from ob.ob_upload import openbis_upload_validated

# Assumes just one file with the extension !
def rename_files_to_basename(path_in_cdict, dest_dir, extensions=('.h5', '.json', '.yml')):
    dest_dir = Path(dest_dir)
    base_name = Path(path_in_cdict).name

    rename_map = {
        '.json': f'{base_name}_concept_dict.json',
        '.yml': f'{base_name}_environment.yml'
    }

    for ext in extensions:
        matching_files = list(dest_dir.glob(f'*{ext}'))
        if matching_files:
            source = matching_files[0]
            target_name = rename_map.get(ext, f'{base_name}{ext}')
            target = dest_dir / target_name
            shutil.move(str(source), str(target))
    new_path_in_cdict = str(dest_dir / base_name)
    return new_path_in_cdict

def get_structures(o, job_to_copy):
    # Parent structure = structure to upload to BAM
    structure_child_id = [
    ch_id for ch_id in job_to_copy.children if 'SAMPLE' in ch_id
    ][0]
    structure_child = o.get_object(structure_child_id)
    structure_parent_id = [
        p_id for p_id in o.get_object(structure_child_id).parents if 'SAMPLE' in p_id
    ][0]
    structure_parent = o.get_object(structure_parent_id)
    return structure_child, structure_parent

def get_datasets(o, obj_to_copy, object_type='job', dest='temp'):
    '''
    object_type: 'job' or 'structure' -> different dataset types

    returns cdict_path = relative path to the concept_dict.json file
    '''
    cdict_path = ''
    
    allowed_object_types = ['job', 'structure']
    if object_type not in allowed_object_types:
        raise ValueError(f'Arg object_type must be one of {allowed_object_types}, got: {object_type}.')

    elif object_type == 'job':
        ds_types = ['PYIRON_HDF5', 'PYIRON_CONCEPT_DICT_DATA']
    elif object_type == 'structure':
        ds_types = ['CRYS-STRUCT_DATA', 'PYIRON_CONCEPT_DICT_DATA']

    for ds_type in ds_types:
        datasets = obj_to_copy.get_datasets(type=ds_type)
        for ds in datasets:
            ds = o.get_dataset(ds.permId)
            ds.download(
                destination=dest, 
                wait_until_finished=True,
                create_default_folders=False
            )
            if ds_type == 'PYIRON_CONCEPT_DICT_DATA':
                cdict_path = os.path.join(
                    dest,
                    ds.file_list[0].split("/")[-1]
                )

    return cdict_path

def get_pseudopotentials(o_bam, job_to_copy):
    pseudopot_ids_sfb = [p_id for p_id in job_to_copy.parents if 'POTENTIAL' in p_id]
    pseudopot_names = [p_id.split('/')[-1] for p_id in pseudopot_ids_sfb]
    pseudopots_bam = o_bam.get_objects(type='PSEUDOPOTENTIAL', code=pseudopot_names)
    pseudopot_ids_bam = [p.permId for p in pseudopots_bam]
    return pseudopot_ids_bam

def get_materials(o_sfb, o_bam, struct_to_copy):
    material_ids_sfb = [m_id for m_id in struct_to_copy.parents if 'MATERIAL' in m_id]
    materials_sfb = o_sfb.get_objects(material_ids_sfb)
    mats_sfb = [m.p.get('$name') for m in materials_sfb]
    material_ids_bam = [
        o_bam.get_objects(type='MATERIAL_V1', where={'$name': mat})[0].permId for mat in mats_sfb
    ]
    return material_ids_bam

def convert_and_upload(o, space, project, collection, cdict_path, 
                       dest='temp', parent_ids=None, options={}):
    with open(cdict_path, 'r') as f:
        loaded_cdict = json.load(f)

    # Returns a list of the following tuples for each cdict submitted
    # cdict, props_dict, object_type, ds_types, ob_parents, object_name
    validated_to_upload = openbis_validate(o, space, project, collection, loaded_cdict, options=options)
    cdict, props_dict, object_type, ds_types, ob_parents, object_name = validated_to_upload[0]
    
    # Rename datasets and change cdict path (locally) to match file location
    path_in_cdict = cdict["structure_name"] if 'structure_name' in cdict.keys() else cdict['path']
    cdict['path'] = rename_files_to_basename(path_in_cdict=path_in_cdict, dest_dir=dest)
    cdict['path'] = f'{dest}/' if 'structure_name' in cdict.keys() else cdict['path']
    object_id = openbis_upload_validated(
        o, space, project, collection, object_name, 
        object_type, ob_parents, props_dict, ds_types, 
        cdict, parent_ids = parent_ids
    )
    return object_id

# delete downloaded files
def cleanup(dest):
    dest_path = Path(dest)
    if dest_path.exists() and dest_path.is_dir():
        shutil.rmtree(dest_path)

def sfb_to_bam(
    username_bam, space_bam, project_bam, collection_bam,
    username_sfb, s3_config_path, job_ids=None, collection_sfb=None,
    include_structure=True, options={}
):
    """
    Copies job objects and its datasets from the SFB openBIS instance to the BAM instance.

    Parameters:
    -----------
    username_bam (str): The BAM instance username for authentication.
    space_bam (str): The BAM space to upload the job(s) to.
    project_bam (str): The BAM project to upload the job(s) to.
    collection_bam (str): The BAM collection to upload the job(s) to.

    username_sfb (str): The SFB instance username for authentication.
    s3_config_path (str): Path to the configuration file containing S3 credentials 
                            for accessing the SFB instance.
    
    job_ids (str or list, optional): The SFB job id(s) (string or list of strings) to copy 
                                      from the SFB instance. Default is None.
    
    collection_sfb (str, optional): The SFB collection (full string) whose all jobs will be copied. 
                                    If specified, it overwrites 'job_ids'. Default is None.
    
    include_structure (bool, optional): Whether to also copy the structure object and its datasets. 
                                        Defaults to True. Set to False for jobs like TableJob.
    
    options (dict, optional): Additional upload options, same function as ob.upload_classic_pyiron.
                                Specify 'materials' to avoid the automatic matching. 
                                Defaults to an empty dictionary.
    
    Notes:
    ------
    - If both 'job_ids' and 'collection_sfb' are specified, the 'collection_sfb' parameter takes 
    precedence, and all jobs in the specified collection will be copied.
    """

    if not (job_ids or collection_sfb):
        raise ValueError('One of the following arguments must be specified: job_ids, collection_sfb')

    # job_ids to be either a string or list of strings of job permIds or identifiers
    job_ids = [job_ids] if isinstance(job_ids, str) else job_ids

    print('SFB openBIS login:')
    o_sfb = openbis_login(url="https://openbis.imm.rwth-aachen.de/openbis/webapp/eln-lims/", 
                      username=username_sfb, instance="sfb1394", s3_config_path=s3_config_path)
    print('BAM openBIS login:')
    o_bam = openbis_login(url="https://test3.datastore.bam.de/", username=username_bam, instance="bam")

    if collection_sfb:
        job_ids = [job.permId for job in o_sfb.get_objects(type='PYIRON_JOB*', collection=collection_sfb)]

    dest = 'temp'
    job_dest = f'{dest}/job'
    child_struct_dest = f'{dest}/struct_final'
    par_struct_dest = f'{dest}/struct'

    for job_id in job_ids:
        job_to_copy = o_sfb.get_object(job_id)
        job_cdict_path = get_datasets(o_sfb, job_to_copy, object_type='job', dest=job_dest)

        material_ids = options.get('materials', [])
        material_ids = [material_ids] if isinstance(material_ids, str) else material_ids
        if include_structure:
            child_struct_to_copy, par_struct_to_copy = get_structures(o_sfb, job_to_copy)
            child_struct_cdict_path = get_datasets(o_sfb, child_struct_to_copy, object_type='structure', dest=child_struct_dest)
            par_struct_cdict_path = get_datasets(o_sfb, par_struct_to_copy, object_type='structure', dest=par_struct_dest)

            if not material_ids:
                material_ids = get_materials(o_sfb, o_bam, par_struct_to_copy)
                
        pseudopot_ids = options.get('pseudopotentials', [])
        pseudopot_ids = [pseudopot_ids] if isinstance(pseudopot_ids, str) else pseudopot_ids
        if not pseudopot_ids:
            pseudopot_ids = get_pseudopotentials(o_bam, job_to_copy)
        options = {
            'materials': material_ids,
            'pseudopotentials' : pseudopot_ids
        }

        par_struct_identifier = None
        if include_structure:
            par_struct_identifier = convert_and_upload(
                o_bam, space_bam, project_bam, collection_bam, par_struct_cdict_path, 
                dest = par_struct_dest, parent_ids = None, options = options
            )
        job_identifier = convert_and_upload(
            o_bam, space_bam, project_bam, collection_bam, job_cdict_path, 
            dest = job_dest, parent_ids = par_struct_identifier, options = options
        )
        if include_structure:
            child_struct_identifier = convert_and_upload(
                o_bam, space_bam, project_bam, collection_bam, child_struct_cdict_path, 
                dest = child_struct_dest, parent_ids = job_identifier, options = options
            )
        

        cleanup(dest)

    o_sfb.logout()
    o_bam.logout()

    print('Upload successful 🙂')