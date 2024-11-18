"""
Author: Khalil Rejiba, Ulrich Kerzel
Date: 2024-10-25
Description: Extension of pybis client to allow LinkedData stored in S3 
"""
from pybis import Openbis
from pybis.dataset import DataSet

import zlib
import xxhash

import boto3
import botocore

import logging
from datetime import datetime, timezone
from configparser import ConfigParser
from io import StringIO
import os
import shutil
import requests
import uuid


def get_bucket_from_client(s3_client: boto3.client):
    """
    Retrieves a bucket name from an S3 client.

    Args:
        s3_client (boto3.client): The S3 client to use for accessing buckets.

    Returns:
        str: The name of the first available bucket, or None if no buckets are found.

    Raises:
        RuntimeError: If no buckets are associated with the S3 client.
    """

    bucket = None

    try:
        response = s3_client.list_buckets()
        buckets = response["Buckets"]

        if len(buckets) == 0:
            # Coscine S3
            owner_display_name = response["Owner"]["DisplayName"]
            if owner_display_name.startswith("write_"):
                bucket = owner_display_name[6:]
        else:
            # AWS S3
            if len(buckets) > 1:
                logging.critical(
                    "Multiple buckets, choosing first bucket that contains objects.", e
                )
            for bucket_info in buckets:
                bucket_name = bucket_info["Name"]
                try:
                    s3_client.list_objects(Bucket=bucket_name)
                    bucket = bucket_name
                    break
                except botocore.exceptions.ClientError as e:
                    logging.error("Error when listing objects", e)

        if bucket is None:
            logging.error("No bucket associated with s3 client.")
            raise RuntimeError("No bucket associated with s3 client.")

    except botocore.exceptions.ClientError as e:
        logging.critical("User does not have access to buckets.")

    return bucket


def crc32(fileName: str):
    """since Python3 the zlib module returns unsigned integers (2.7: signed int)"""
    prev = 0
    for eachLine in open(fileName, "rb"):
        prev = zlib.crc32(eachLine, prev)
    # return as hex

    return "%x" % (prev & 0xFFFFFFFF)


def compute_xxhash64(path: str, block_size: int = 2 ** 32) -> str:
    """
    Computes the xxHash (64-bit) digest for a given file.

    Args:
        path (str): full path to the file to be hashed.
        block_size (int, optional): The block size for reading the file (default is 4096).

    Returns:
        str: The hexadecimal representation of the xxHash digest.

    Raises:
        IOError: If the file cannot be read or xxHash computation fails.
    """
    try:
        with open(path, "rb") as file:
            h = xxhash.xxh64()
            for chunk in iter(lambda: file.read(block_size), b""):
                h.update(chunk)
        return h.hexdigest()
    except IOError as e:
        logging.error("xxHash could not be computed: %s" % e)
        raise e


def get_dms_info(oBis: Openbis, filename: str, dms_code: str) -> tuple[str, dict]:
    """
    Generates the full file path on the external storage.

    The full path consists of the external storage, the relative path,
    and the name of the file.

    Args:
        oBis (Openbis): Openbis instance object, connecting to the openBIS server.
        filename (str): Filename to be stored.
        dms_code (str): Code of the openBIS external data storage service.

    Returns:
        tuple of str, dict: The full path name on the external storage.
    """
    try:
        dms = oBis.get_external_data_management_system(dms_code)
    except ValueError as e:
        logging.debug("Error: %s" % e)
        raise e
    logging.debug("Data will be stored in: %s" % dms.urlTemplate)

    dms_id = oBis.external_data_managment_system_to_dms_id(dms_code)
    logging.debug("openBIS DMS ID %s" % dms_id)

    file_path = dms.urlTemplate + "/" + filename
    return file_path, dms_id


def retrieve_data_store_url(dss_code: str, oBis: Openbis) -> str:
    """
    Retrieves the URL of the openBIS data storage server (DSS).

    Args:
        dss_code (str): Identifier of the openBIS DSS server to handle the request.
        oBis (openbis): Openbis instance object, connecting to the openBIS server.

    Returns:
        str: URL + json that is used by openBIS to process the request.
    """
    data_stores = oBis.get_datastores()
    data_store = data_stores[data_stores["code"] == dss_code]
    download_url = data_store["downloadUrl"][0]
    return "%s/datastore_server/rmi-data-store-server-v3.json" % download_url


def get_data_store_url(obis_link_metadata: dict, oBis: Openbis) -> str:
    """
    Extracts the URL and relevant JSON file of the openBIS data storage server (DSS).

    The output should resemble the following:
    https://openbis-t.imm.rwth-aachen.de:443/datastore_server/rmi-data-store-server-v3.json.

    Args:
        obis_link_metadata (dict): openBIS metadata and instructions for linked files.
        oBis (Openbis): Openbis instance object, connecting to the openBIS server.

    Returns:
        str: The server URL and JSON file combination.
    """

    server_url = retrieve_data_store_url(
        obis_link_metadata["dataStoreId"]["permId"], oBis=oBis
    )
    logging.debug("DSS URL: %s" % server_url)

    return server_url


def get_file_metadata(
    filename: str, dms_path: str, compute_crc32: bool = False
) -> list:
    """
    Generates metadata on the file to be linked itself: file_size, checksum.
    Note: Also needs to include the full qualified path and filename on the external datastore.

    Args:
        filename (str): Name of the file for which the metadata is required.
        dms_path (str): Path on the external storage where the file is to be stored.
        compute_crc32 (bool, optional): Whether or not a CRC32 checksum is to be computed (default is False).

    Raises:
        FileNotFoundError: If the file cannot be found.

    Returns:
        list: List of metadata for the file to be registered in openBIS.
    """

    # based on code from:
    # https://sissource.ethz.ch/sispub/openbis/-/blob/master/api-openbis-python3-pybis/src/python/pybis/data_set.py
    # param contents: A list of dicts...

    if not os.path.isfile(filename):
        logging.critical("File not found: %s" % filename)
        raise FileNotFoundError("File %s not found" % filename)
    file_size = os.path.getsize(filename)
    file_crc32 = 0
    # CRC32 can be slow for large files.

    if compute_crc32:
        file_crc32 = crc32(filename)
    try:
        file_xxhash = compute_xxhash64(filename)
    except IOError as e:
        logging.error("xxHash not computed: %s" % e)
        raise e
    logging.debug("file size: %s" % file_size)
    logging.debug("file_crc32: %s" % file_crc32)
    logging.debug("XXHash: %s" % file_xxhash)

    content = {
        "fileLength": file_size,
        "crc32": file_crc32,
        "crc32Checksum": file_crc32,
        "checksum": file_xxhash,
        "checksumType": "xxHash",
        "directory": False,
        "path": dms_path,
    }
    file_metadata = [content]

    logging.debug("File metadata: %s" % file_metadata)

    return file_metadata


def create_obis_link_metadata(
    oBis: Openbis,
    dms_path: str,
    data_set_type: str,
    dss_code: str,
    dms_id: dict,
    parent_ids: list,
    properties: dict,
    sample_name: str,
    experiment_name: str,
    data_set_code: str,
) -> dict:
    """
    Creates a metadata dictionary containing information on how openBIS should link the data,
    along with relevant metadata.

    Args:
        oBis (Openbis): Openbis instance object, connecting to the openBIS server.
        dms_path (str): Fully qualified path of the file on the external data store (including filename).
        data_set_type (str): openBIS data type for the linked file.
        dss_code (str): Identifier of the openBIS DSS server to handle the request.
        dms_id (dict): Details about the external storage (DMS).
        parent_ids (list): Parents of the file to be linked.
        properties (dict): Properties of the file in the ELN (e.g., name, metadata).
        sample_name (str): String identifier of the sample to link this file to (must specify either sample or experiment).
        experiment_name (str): String identifier of the experiment to link the file to (must specify either sample or experiment).
        data_set_code (str): Code of the dataset to attach the file to. If None, a new one will be created.

    Returns:
        dict: A metadata dictionary for the file to be registered in openBIS.
    """

    data_set_creation = {
        "linkedData": {
            "@type": "as.dto.dataset.create.LinkedDataCreation",
            "contentCopies": [
                {
                    "@type": "as.dto.dataset.create.ContentCopyCreation",
                    "path": dms_path,
                    # "gitCommitHash": None,
                    # "gitRepositoryId": None,
                    "externalDmsId": dms_id,
                }
            ],
        },
        "typeId": {
            "@type": "as.dto.entitytype.id.EntityTypePermId",
            "permId": data_set_type,
        },
        "dataStoreId": {
            "permId": dss_code,
            "@type": "as.dto.datastore.id.DataStorePermId",
        },
        "parentIds": parent_ids,
        "measured": False,
        "properties": properties,
        "@type": "as.dto.dataset.create.DataSetCreation",
    }

    if sample_name is not None:
        sample_id = oBis.sample_to_sample_id(sample_name)
        data_set_creation["sampleId"] = sample_id
    elif experiment_name is not None:
        experiment_id = oBis.experiment_to_experiment_id(experiment_name)
        data_set_creation["experimentId"] = experiment_id
    if data_set_code is not None:
        data_set_creation["code"] = data_set_code
        data_set_creation["autoGeneratedCode"] = False
    else:
        data_set_creation["autoGeneratedCode"] = True
    logging.debug("openBIS link metadata dict: %s" % data_set_creation)

    return data_set_creation


def create_obis_link_request(
    obis_link_metadata: dict, file_metadata: list, pat: str
) -> dict:
    """
    Combines openBIS metadata and file metadata into the service request
    that is then sent to the openBIS server.

    Args:
        obis_link_metadata (dict): openBIS metadata and instructions for linked files.
        file_metadata (list): Metadata about the file (CRC, path on DMS).
        pat (str): openBIS access token (session token or personal access token).

    Returns:
        dict: Request to send to the openBIS server.
    """

    data_set_creation = {
        "fileMetadata": file_metadata,
        "metadataCreation": obis_link_metadata,
        "@type": "dss.dto.dataset.create.FullDataSetCreation",
    }

    request = {
        "method": "createDataSets",
        "params": [pat, [data_set_creation]],
    }

    logging.debug("openBIS link request: %s" % request)

    return request


def register_file(
    oBis: Openbis,
    file_metadata: list,
    dms_path: str,
    dms_id: dict,
    dss_code: str,
    sample_name: str,
    experiment_name: str,
    properties: dict,
    data_set_code: str,
    data_set_type: str,
    parent_ids: list,
    token: str,
) -> str:
    """
    Registers the file in openBIS as a linked file on the external storage system.
    Assumes that the file is copied to this system manually or by other means:
    This function only registers the metadata in openBIS.

    Args:
        oBis (Openbis): Openbis instance object, connecting to the openBIS server.
        file_metadata (list): Metadata about the file (CRC, path on DMS).
        dms_path (str): Fully qualified path of the file on the external data store (including filename).
        dms_id (dict): Details about the external storage (DMS).
        dss_code (str): Identifier of the openBIS DSS server to handle the request.
        sample_name (str): String identifier of the sample to link this file to (must specify either sample or experiment).
        experiment_name (str): String identifier of the experiment to link the file to (must specify either sample or experiment).
        properties (dict): Properties of the file in the ELN (e.g., name, metadata).
        data_set_code (str): Code of the dataset to attach the file to. If None, a new one will be created.
        data_set_type (str): openBIS data type for the linked file.
        parent_ids (list): Parents of the file to be linked.
        token (str): openBIS access token (session token or personal access token).

    Raises:
        Exception: If the registration of the linked file in openBIS fails.

    Returns:
        str: openBIS permID of the newly linked file.
    """

    obis_link_metadata = create_obis_link_metadata(
        oBis=oBis,
        dms_path=dms_path,
        dms_id=dms_id,
        dss_code=dss_code,
        parent_ids=parent_ids,
        properties=properties,
        sample_name=sample_name,
        experiment_name=experiment_name,
        data_set_code=data_set_code,
        data_set_type=data_set_type,
    )
    logging.debug("openBIS link metadata: \n%s" % obis_link_metadata)

    obis_link_request = create_obis_link_request(
        obis_link_metadata=obis_link_metadata, file_metadata=file_metadata, pat=token
    )
    logging.debug("openBIS link request: \n%s" % obis_link_request)

    obis_dss_url = get_data_store_url(oBis=oBis, obis_link_metadata=obis_link_metadata)
    logging.debug("URL for the link request (DSS): %s" % obis_dss_url)

    # POST request to openBIS server

    try:
        response = oBis._post_request_full_url(obis_dss_url, obis_link_request)
        logging.debug("Server response to POST command: %s" % response)

        # openBIS permID of the newly linked dataset

        permid_linked_dataset = oBis.get_dataset(response[0]["permId"]).permId
        logging.debug("New openBIS permID %s" % permid_linked_dataset)
    except Exception as e:
        logging.error("Something went wrong in POST command: %s" % e)
        raise e
    return permid_linked_dataset


class ExtendedDataSet(
    DataSet,
    entity="dataSet",
    single_item_method_name="get_dataset",
):
    def __init__(
        self,
        openbis_client,
        type,
        data=None,
        files=None,
        zipfile=None,
        folder=None,
        kind=None,
        props=None,
        s3_client=None,
        s3_link_validity=604800,
        **kwargs,
    ):

        if kind == "PHYSICAL":
            if files is None and zipfile is None:
                raise ValueError("please provide at least one file")

            if files and zipfile:
                raise ValueError(
                    "please provide either a list of files or a single zipfile"
                )

            if zipfile:
                files = [zipfile]
                self.__dict__["isZipDirectoryUpload"] = True
            else:
                self.__dict__["isZipDirectoryUpload"] = False

            if files:
                if isinstance(files, str):
                    files = [files]

                for file in files:
                    if not os.path.exists(file):
                        raise ValueError(f"File {file} does not exist")

                self.__dict__["files"] = files

        # initialize the OpenBisObject
        super(ExtendedDataSet, self).__init__(
            openbis_client, type=type, data=data, props=props, **kwargs
        )

        self.__dict__["files_in_wsp"] = []

        # existing DataSet
        if data is not None:
            if data["physicalData"] is None:
                self.__dict__["shareId"] = None
                self.__dict__["location"] = None
            else:
                self.__dict__["shareId"] = data["physicalData"]["shareId"]
                self.__dict__["location"] = data["physicalData"]["location"]

        if kind is not None:
            kind = kind.upper()
            allowed_kinds = ["PHYSICAL", "CONTAINER", "LINK"]
            if kind not in allowed_kinds:
                raise ValueError(
                    f"only these values are allowed for kind: {allowed_kinds}"
                )
            self.a.__dict__["_kind"] = kind

        self.__dict__["folder"] = folder

        if getattr(self, "parents") is None:
            self.a.__dict__["_parents"] = []
        else:
            if not self.is_new:
                self.a.__dict__["_parents_orig"] = self.a.__dict__["_parents"]

        if getattr(self, "children") is None:
            self.a.__dict__["_children"] = []
        else:
            if not self.is_new:
                self.a.__dict__["_children_orig"] = self.a.__dict__["_children"]

        if getattr(self, "container") is None:
            self.a.__dict__["_container"] = []
        else:
            if not self.is_new:
                self.a.__dict__["_container_orig"] = self.a.__dict__["_container"]

        if getattr(self, "component") is None:
            self.a.__dict__["_component"] = []
        else:
            if not self.is_new:
                self.a.__dict__["_component_orig"] = self.a.__dict__["_component"]

        if getattr(self, "copy_path") is None:
            self.a.__dict__["copy_path"] = []

        filename = files[0]  # FIXME deal with multiple files

        time_stamp = datetime.now(timezone.utc)
        time_stamp = time_stamp.strftime("%Y-%m-%dT%H-%M-%S.%f")
        dataset_type = type.code
        username = openbis_client._get_username()
        prefix = time_stamp + "_" + dataset_type + "_" + username + "_"

        if openbis_client.standardize_filenames:
            destination = os.path.join(
                openbis_client.temp_dir, prefix + os.path.basename(filename)
            )
            shutil.copy(filename, destination)
            filename = destination
            self.a.__dict__["copy_path"].append(destination)

        if kind == "LINK":
            # upload_target=="s3" is a better choice if you have GitDataSet in addition to S3
            s3_client = openbis_client._s3_client

            # We need to uplaod the file to create the URL
            self._upload_file_to_s3(filename, s3_client)

            # Fill s3_download_link in self.props
            self._create_download_link(
                filename=filename, s3_client=s3_client, validity=s3_link_validity
            )

            props = self.props.all()
            dms_path, dms_id = get_dms_info(
                oBis=openbis_client, filename=filename, dms_code=openbis_client.dms_code
            )
            file_metadata = get_file_metadata(
                filename=filename,
                dms_path=dms_path,
                compute_crc32=False,
            )
            experiment_name = kwargs.get("experiment")
            sample_name = kwargs.get("sample")
            parent_ids = kwargs.get("parent_ids")
            dataset_code = None  # Generate a new one
            dss_code = openbis_client.get_datastores()["code"][0]

            self.a.__dict__["args_register_file"] = (
                openbis_client,
                file_metadata,
                dms_path,
                dms_id,
                dss_code,
                sample_name,
                experiment_name,
                props,
                dataset_code,
                dataset_type,
                parent_ids,
                openbis_client.token,
            )
            # Dataset is only registered in openBIS when ds.save() is called.

    def _upload_file_to_s3(self, filename, s3_client):
        """
        Uploads a file to an S3 bucket.

        This function uploads the specified file to an S3 bucket.

        Args:
            filename (str): The name of the file to upload.
            s3_client (boto3.client): The S3 client to use for uploading the file.

        """
        bucket = get_bucket_from_client(s3_client)

        try:
            s3_client.upload_file(
                Filename=filename,
                Bucket=bucket,
                Key=os.path.basename(filename),
            )
        except Exception as e:
            logging.error("Error when uploading to S3:", e)

    def _create_download_link(self, filename, s3_client, validity):
        """
        Creates a presigned URL for downloading a file from an S3 bucket.

        This function generates a presigned URL for the specified file in an S3 bucket,
        allowing the file to be downloaded within the specified validity period.
        The URL is property of the dataset that needs to be updated on a regular basis.

        Args:
            filename (str): The name of the file to generate the download link for.
            s3_client (boto3.client): The S3 client to use for generating the presigned URL.
            validity (int): The validity period of the presigned URL in seconds.
        """

        bucket = get_bucket_from_client(s3_client)

        try:
            s3_client.head_object(
                Bucket=bucket,
                Key=os.path.basename(filename),
            )
        except Exception as e:
            logging.error("Error when creating presigned URL:", e)

        url = s3_client.generate_presigned_url(
            "get_object",
            Params={"Bucket": bucket, "Key": os.path.basename(filename)},
            ExpiresIn=validity,
        )

        self.props["s3_download_link"] = url

        return url

    def __str__(self):
        if self.data is not None:
            return self.data["code"]
        else:
            return "Dataset permId will be created after calling save()"

    def save(self, permId=None):
        if self.kind == "PHYSICAL":
            super(ExtendedDataSet, self).save()
        elif self.kind == "LINK":
            permId = register_file(*self.a.__dict__["args_register_file"])
        if len(self.copy_path):
            if os.path.exists(self.copy_path[0]):
                os.unlink(self.copy_path[0])

    def download(
        self,
        files=None,
        destination=None,
        create_default_folders=True,
        wait_until_finished=True,
        workers=10,
        linked_dataset_fileservice_url=None,
        content_copy_index=0,
    ):
        # FIXME Not tested
        if self.kind == "PHYSICAL":
            super(ExtendedDataSet, self).download(
                files,
                destination,
                create_default_folders,
                wait_until_finished,
                workers,
                linked_dataset_fileservice_url,
                content_copy_index,
            )
        elif self.kind == "LINK":
            file_url = self.props["s3_download_link"]
            filename = file_url.split("/")[-1].split("?")[0]
            if destination is None:
                destination = "data/" + self.openbis.download_prefix

            if create_default_folders:
                filename_dest = os.path.join(
                    destination, self.permId, "original", filename
                )
            else:
                filename_dest = os.path.join(filename)
            try:
                with requests.get(file_url, stream=True) as response:
                    response.raise_for_status()
                    with open(filename_dest, "wb") as file:
                        for chunk in response.iter_content(chunk_size=8192):
                            if chunk:
                                file.write(chunk)
            except requests.exceptions.RequestException as e:
                logging.error(f"Error downloading {filename} from S3:", e)


class OpenbisWithS3(Openbis):
    """
    A class to extend Openbis functionality with S3 integration.

    This class reimplements some of the methods offered by pybis to enable users to interact with S3 storage,
    allowing for seamless integration of S3 with Openbis.
    Each S3 bucket is registered (by an instance admin) as an External Data Management System in openBIS.
    The files are not stored in the openBIS instance, therefore we use the LinkedData functionality of openBIS.
    """

    def __init__(
        self,
        url=None,
        verify_certificates=True,
        token=None,
        use_cache=True,
        allow_http_but_do_not_use_this_in_production_and_only_within_safe_networks=False,
        s3_config_path=None,
        standardize_filenames=True,
    ):
        super(OpenbisWithS3, self).__init__(
            url=url,
            verify_certificates=verify_certificates,
            token=token,
            use_cache=use_cache,
            allow_http_but_do_not_use_this_in_production_and_only_within_safe_networks=allow_http_but_do_not_use_this_in_production_and_only_within_safe_networks,
        )

        if s3_config_path is not None:
            self._configure_s3_client(s3_config_path)
        self.standardize_filenames = standardize_filenames  # FIXME only works for S3
        self.temp_dir = "pybisaixtented_tmp_" + str(uuid.uuid4())
        os.mkdir(self.temp_dir)

    def _configure_s3_client(self, config_path):
        """
        Configures the S3 client using the provided configuration file.

        This method initializes the S3 client used for file upload
        and fetches the bucket name and DMS code
        based on the configuration specified in the given file path.

        Args:
            config_path (str): The path to the S3 configuration file.

        Returns:
            None
        """
        s3_client, bucket_name, dms_code = self._get_s3_client(
            config_file=config_path,
            from_path=True,
        )
        self._s3_client = s3_client or None  # Upload
        self._bucket = bucket_name
        self.dms_code = dms_code

    def _get_s3_client(self, config_file=None, from_path=False):
        """Read in the config file and parse the configuration settings.

        Same settings can be provided in environment variables.

        Args:
            config (opened file): file-like object
        """
        (
            s3_region,
            s3_endpoint_url,
            s3_endpoint_port,
            s3_key,
            s3_secret,
            s3_bucket_name,
            obis_dms_code,
        ) = [None] * 7

        parser = ConfigParser(
            allow_no_value=True,
        )

        if config_file is not None:

            # Parse config file
            if from_path:
                if os.path.exists(config_file):
                    with open(config_file, "r") as fh:
                        data_str = fh.read()
                else:
                    logging.critical(f"Config file {config_file} was not found.\n")
            else:
                # Convert to string as the file is already opened
                stringio = StringIO(config_file.decode("utf-8"))
                data_str = stringio.read()
            parser.read_string(data_str)

            try:
                dms_codes = (
                    self.get_external_data_management_systems().df.code.to_list()
                )
            except:
                dms_codes = []

            try:

                # openBIS details
                obis_dms_code = parser.get(
                    "openbis", "dms_code", fallback=None
                ) or parser.get("openBIS", "dms_code", fallback=None)
                if obis_dms_code is None:
                    raise ValueError("dms_code missing in config file")
                if len(dms_codes) and obis_dms_code not in dms_codes:
                    raise ValueError(
                        "dms_code in config file is not recognized by openBIS"
                    )

                # S3 details
                s3_region = parser.get("s3", "s3_region", fallback="eu-central-1")
                s3_endpoint_url = parser.get("s3", "s3_endpoint_url", fallback=None)
                s3_endpoint_port = parser.get("s3", "s3_endpoint_port", fallback=None)
                s3_key = parser.get("s3", "s3_access_key")
                s3_secret = parser.get("s3", "s3_access_secret")
                s3_bucket_name = parser.get("s3", "s3_bucket")

            except Exception as e:
                logging.warning("Config file is NOT formatted correctly.\n")
                raise ValueError(f"Config file is NOT formatted correctly.\n>>> {e}")

        else:
            logging.warning(
                "Attempting to get settings from environment variables ...\n"
            )

            VARS = [
                "s3_region",
                "s3_endpoint_url",
                "s3_endpoint_port",
                "s3_access_key",
                "s3_access_secret",
                "s3_bucket",
                "dms_code",
            ]
            env_vars = []
            for var_name in VARS:
                try:
                    v = os.environ[var_name]
                except KeyError:
                    # logging.critical(f"Environment variable {var_name} was not found\n")
                    v = None
                env_vars.append(v)
            (
                s3_region,
                s3_endpoint_url,
                s3_endpoint_port,
                s3_key,
                s3_secret,
                s3_bucket_name,
                obis_dms_code,
            ) = env_vars

        if s3_endpoint_url is not None and s3_endpoint_port is not None:
            s3_url = s3_endpoint_url + ":" + str(s3_endpoint_port)
        else:
            s3_url = None

        s3_client = boto3.client(
            service_name="s3",
            endpoint_url=s3_url,
            region_name=s3_region,
            aws_access_key_id=s3_key,
            aws_secret_access_key=s3_secret,
            config=boto3.session.Config(
                signature_version="s3v4",
                connect_timeout=5,
                read_timeout=10,
            ),
        )
        try:
            # Test s3:GetObject
            s3_client.list_objects(Bucket=s3_bucket_name)

            # Create dummy file
            dummy_file = "dummy_file.txt"
            content = "OpenBISAixTended"
            n = 1024 ** 2 // len(content)
            with open(dummy_file, "w") as f:
                f.write(content * n)

            # Test s3:PutObject
            s3_client.upload_file(
                Filename=dummy_file, Bucket=s3_bucket_name, Key=dummy_file
            )

            # Test s3:DeleteObject
            s3_client.delete_object(Bucket=s3_bucket_name, Key=dummy_file)

            os.unlink(dummy_file)
            logging.info(f"S3 client configured successfully.\n")
            s3_upload_ok = True

        except botocore.exceptions.EndpointConnectionError as e:
            logging.error(f"Connection error: {e}")
            s3_upload_ok = False
        except botocore.exceptions.ReadTimeoutError as e:
            logging.error(f"Read timeout: {e}")
            s3_upload_ok = False
        except Exception as e:
            logging.error(f"Error when testing client: {e}")
            s3_upload_ok = False

        if not s3_upload_ok:
            raise RuntimeError(f"S3 client was NOT properly configured.")

        return s3_client, s3_bucket_name, obis_dms_code

    def __del__(self):
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)

    def new_dataset(
        self,
        type=None,
        kind="PHYSICAL",
        files=None,
        file=None,
        props=None,
        folder=None,
        **kwargs,
    ):
        try:
            dms_codes = self.get_external_data_management_systems().df.code.to_list()
        except:
            raise ValueError(
                "Please login using token or username and password before creating a dataset!"
            )

        if len(dms_codes) and kind == "LINK":
            if self.dms_code not in dms_codes:
                raise ValueError("dms_code in config file is not recognized by openBIS")

        if type is None:
            raise ValueError("Please provide a dataSet type")

        if file:
            files = [file]

        if isinstance(type, str):
            type_obj = self.get_dataset_type(type.upper())
        else:
            type_obj = type

        if "object" in kwargs:
            kwargs["sample"] = kwargs["object"]
            kwargs.pop("object", None)
        if "collection" in kwargs:
            kwargs["experiment"] = kwargs["collection"]
            kwargs.pop("collection", None)

        sample = kwargs.get("sample")
        experiment = kwargs.get("experiment")
        assert (sample is not None) or (
            experiment is not None
        ), "Either experiment or sample should be provided as an argument"

        if experiment is not None:
            try:
                kwargs["experiment"] = self.get_experiment(experiment).identifier
            except ValueError as e:
                raise ValueError(
                    f"Experiment / Collection {experiment} does NOT exist for user {self._get_username()}"
                )
        if sample is not None:
            try:
                kwargs["sample"] = self.get_sample(sample).identifier
                kwargs["experiment"] = self.get_sample(sample).experiment.identifier
            except ValueError as e:
                raise ValueError(
                    f"Sample / Object {sample} does NOT exist for user: {self._get_username()}"
                )

        return ExtendedDataSet(
            self,
            type=type_obj,
            kind=kind,
            files=files,
            folder=folder,
            props=props,
            s3_client=self._s3_client,
            **kwargs,
        )
