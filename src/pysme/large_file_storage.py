# -*- coding: utf-8 -*-
"""
System to store large data files on a server
Load them whem required by the user
Update the pointer file on github when new datafiles become available

Pro: Versioning is effectively done by Git
Con: Need to run server
"""

import gzip
import json
import logging
import os
import shutil
import tarfile
from os.path import basename
from pathlib import Path
from tempfile import NamedTemporaryFile

from astropy.utils.data import (
    clear_download_cache,
    download_file,
    import_file_to_cache,
)
from tqdm.auto import tqdm
from tqdm.utils import CallbackIOWrapper

from .config import Config
from .util import show_progress_bars

logger = logging.getLogger(__name__)

# We are lazy and want a simple check if a file is in the Path
Path.__contains__ = lambda self, key: (self / key).exists()


class LargeFileStorage:
    """
    Download large data files from data server when needed
    New versions of the datafiles are indicated in a 'pointer' file
    that includes the hash of the newest version of the files

    Raises
    ------
    FileNotFoundError
        If the datafiles can't be located anywhere
    """

    def __init__(self, server, pointers, storage):
        #:list[str]: ordered mirrors to try
        self.servers = self._normalize_servers(server)
        #:str: legacy single-server alias (first mirror)
        self.server = self.servers[0] if len(self.servers) > 0 else ""

        if isinstance(pointers, str):
            path = Path(__file__).parent / pointers
            pointers = LargeFileStorage.load_pointers_file(path)

        #:dict(fname:hash): points from a filename to the current newest object id, usually a hash
        self.pointers = pointers
        #:Directory: directory of the current data files
        cache_path = Path(storage).expanduser().resolve(strict=False)
        self.current = cache_path

        # set the folder to download the data file into
        # need to set environment variable because astropy will put things into home otherwise
        os.environ["XDG_CACHE_HOME"] = str(cache_path)
        # if someone is using astropy along with pysme, it might mess with their astropy file storage
        # not threadsafe, but multiprocessing safe, because threads shares environment variables
        self.PKGNAME = ""

        if not cache_path.exists():
            print("folder to store data file does not exist, creating")
        cache_path.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def _normalize_servers(server):
        if server is None:
            return []
        if isinstance(server, (list, tuple)):
            return [str(s).strip() for s in server if str(s).strip() != ""]
        value = str(server).strip()
        if value == "":
            return []
        return [value]

    @staticmethod
    def _is_uri(value):
        return value.startswith(("http://", "https://", "file://"))

    @staticmethod
    def _join_uri(base, path):
        return base.rstrip("/") + "/" + path.lstrip("/")

    @staticmethod
    def _unique_in_order(values):
        seen = set()
        unique = []
        for value in values:
            if value in seen:
                continue
            unique.append(value)
            seen.add(value)
        return unique

    @staticmethod
    def _get_nlte_element_from_key(key):
        key = str(key)
        prefix = "nlte_"
        suffix = "_pysme.grd"
        if not (key.startswith(prefix) and key.endswith(suffix)):
            return None
        return key[len(prefix) : -len(suffix)]

    def _store_processed_file(self, url, filename):
        import_file_to_cache(url, filename, pkgname=self.PKGNAME)
        return download_file(url, cache=True, pkgname=self.PKGNAME)

    def _detect_download_format(self, fname, url, compression):
        if compression is None:
            return "plain"
        if compression != "gzip":
            return compression

        url_lower = str(url).lower()
        if any(token in url_lower for token in (".tar.gz", ".tgz", ".tar")):
            if tarfile.is_tarfile(fname):
                return "tar.gz"

        if tarfile.is_tarfile(fname):
            return "tar.gz"
        return "gzip"

    @staticmethod
    def load_pointers_file(filename):
        try:
            with open(str(filename), "r") as f:
                pointers = json.load(f)
        except FileNotFoundError:
            logger.error("Could not find LargeFileStorage reference file %s", filename)
            pointers = {}
        return pointers

    def get(self, key):
        """
        Request a datafile from the LargeFileStorage
        Assures that tracked files are at the specified version
        And downloads data from the server if necessary

        Parameters
        ----------
        key : str
            Name of the requested datafile

        Raises
        ------
        FileNotFoundError
            If the requested datafile can not be found anywhere

        Returns
        -------
        fullpath : str
            Absolute path to the datafile
        """
        urls = self.get_urls(key)
        errors = []

        for i, url in enumerate(urls):
            try:
                # If its a direct file link, pass that directly to
                if url.startswith("file://"):
                    local_path = url[7:]
                    if os.path.exists(local_path):
                        return local_path
                    raise FileNotFoundError(
                        f"Local file URI does not exist: {local_path}"
                    )

                fname = download_file(url, cache=True, pkgname=self.PKGNAME)

                compression = self._test_compression(fname)
                file_format = self._detect_download_format(fname, url, compression)
                if file_format == "plain":
                    return fname
                if file_format == "gzip":
                    return self._unpack_gzip(fname, key, url)
                if file_format == "tar.gz":
                    return self._unpack_tar_gzip(fname, key, url)

                raise ValueError(
                    "The file is compressed using %s, which is not supported"
                    % compression
                )
            except Exception as exc:
                errors.append((url, exc))
                if i < len(urls) - 1:
                    logger.warning(
                        "Could not fetch %s from %s, trying fallback mirror",
                        key,
                        url,
                    )

        detail = "; ".join(f"{url}: {exc}" for url, exc in errors)
        raise FileNotFoundError(f"Could not fetch tracked file {key}. Attempts: {detail}")

    def get_urls(self, key):
        """
        Return ordered candidate URLs/URIs for a tracked key.

        For tracked files:
        - pointer value may be a string or a list of strings
        - each pointer string may be a full URI (http/https/file) or relative path
        - relative paths are combined with every configured mirror server in order
        """
        key = str(key)

        # Check if the file is tracked and/or exists in the storage directory
        if key not in self.pointers:
            if key not in self.current:
                if not os.path.exists(key):
                    raise FileNotFoundError(
                        f"File {key} does not exist and is not tracked by the Large File system"
                    )
                else:
                    return [Path(key).as_uri()]
            else:
                return [(self.current / key).as_uri()]

        newest = self.pointers[key]
        pointer_targets = newest if isinstance(newest, list) else [newest]

        urls = []
        for target in pointer_targets:
            target = str(target).strip()
            if target == "":
                continue
            if self._is_uri(target):
                urls.append(target)
            elif len(self.servers) > 0:
                for server in self.servers:
                    urls.append(self._join_uri(server, target))
            else:
                urls.append(target)

        return self._unique_in_order(urls)

    def _test_compression(self, fname):
        """Check filetype using the magic string"""
        with open(fname, "rb") as f:
            magic = f.read(6)
            if magic[:2] == b"\x1f\x8b":
                return "gzip"
            if magic[:6] == b"\x37\x7A\xBC\xAF\x27\x1C":
                return "7z"
            if magic[:5] == b"\x50\x4B\x03\x04":
                return "zip"
            return None

    def _unpack_gzip(self, fname, key, url):
        logger.debug("Unpacking data file %s", key)

        # We have to use a try except block, as this will crash with
        # permissions denied on windows, when trying to copy an open file
        # here the temporary file
        # Therefore we close the file, after copying and then delete it manually
        extracted_name = None
        try:
            with gzip.open(fname, "rb") as f_in:
                with NamedTemporaryFile("wb", delete=False) as f_out:
                    extracted_name = f_out.name
                    with tqdm(
                        # total=f_in.size,
                        desc="Unpack",
                        unit="B",
                        unit_scale=True,
                        unit_divisor=1024,
                        disable=~show_progress_bars,
                    ) as t:
                        fobj = CallbackIOWrapper(t.update, f_in, "read")
                        while True:
                            chunk = fobj.read(1024)
                            if not chunk:
                                break
                            f_out.write(chunk)
                        f_out.flush()
                        t.reset()
            return self._store_processed_file(url, extracted_name)
        finally:
            if extracted_name is not None:
                try:
                    os.remove(extracted_name)
                except OSError:
                    pass

    def _unpack_tar_gzip(self, fname, key, url):
        element = self._get_nlte_element_from_key(key)
        if element is None:
            raise ValueError(
                f"tar.gz archives are only supported for nlte_*_pysme.grd keys, got {key}"
            )

        target_name = f"nlte_{element.lower()}"
        extracted_name = None
        try:
            with tarfile.open(fname, "r:gz") as tf:
                members = []
                for member in tf.getmembers():
                    if not member.isfile():
                        continue
                    base = Path(member.name).name.lower()
                    if not base.endswith(".grd"):
                        continue
                    if "pysme" not in base:
                        continue
                    if not base.startswith(target_name):
                        continue
                    members.append(member)

                if len(members) != 1:
                    names = [m.name for m in members]
                    raise ValueError(
                        f"Expected exactly one pysme .grd for {key} in {url}, found {len(members)}: {names}"
                    )

                with tf.extractfile(members[0]) as f_in:
                    with NamedTemporaryFile("wb", suffix=".grd", delete=False) as f_out:
                        extracted_name = f_out.name
                        shutil.copyfileobj(f_in, f_out)

            return self._store_processed_file(url, extracted_name)
        finally:
            if extracted_name is not None:
                try:
                    os.remove(extracted_name)
                except OSError:
                    pass

    def get_url(self, key):
        urls = self.get_urls(key)
        if len(urls) == 0:
            raise FileNotFoundError(f"No download target configured for key {key}")
        return urls[0]

    def clean_cache(self):
        """Remove unused cache files (from old versions)"""
        clear_download_cache(pkgname=self.PKGNAME)

    def delete_file(self, fname):
        """Delete a file, including the cache file"""
        clear_download_cache(fname, pkgname=self.PKGNAME)

    def move_to_cache(self, fname, key=None):
        """Move currently used files into cache directory and use symlinks instead,
        just as if downloaded from a server"""
        if key is None:
            key = basename(fname)
        import_file_to_cache(key, fname, pkgname=self.PKGNAME)
        self.pointers[key] = key


def _get_file_servers(config):
    try:
        servers = config["data.file_servers"]
        if isinstance(servers, list) and len(servers) == 0:
            return config["data.file_server"]
        return servers
    except KeyError:
        return config["data.file_server"]


def setup_atmo(config=None):
    if config is None:
        config = Config()
    server = _get_file_servers(config)
    storage = config["data.atmospheres"]
    pointers = config["data.pointers.atmospheres"]
    lfs_atmo = LargeFileStorage(server, pointers, storage)
    return lfs_atmo


def setup_nlte(config=None):
    if config is None:
        config = Config()
    server = _get_file_servers(config)
    storage = config["data.nlte_grids"]
    pointers = config["data.pointers.nlte_grids"]
    lfs_nlte = LargeFileStorage(server, pointers, storage)
    return lfs_nlte


def setup_lfs(config=None, lfs_atmo=None, lfs_nlte=None):
    if config is None:
        config = Config()
    if lfs_atmo is None:
        lfs_atmo = setup_atmo(config)
    if lfs_nlte is None:
        lfs_nlte = setup_nlte(config)
    return config, lfs_atmo, lfs_nlte


def get_available_atmospheres(config=None):
    if config is None:
        config = Config()
    pointers = config["data.pointers.atmospheres"]
    storage = config["data.atmospheres"]
    data = get_available_files(pointers, storage)
    return data


def get_available_nlte_grids(config=None):
    if config is None:
        config = Config()
    pointers = config["data.pointers.nlte_grids"]
    storage = config["data.nlte_grids"]
    data = get_available_files(pointers, storage)
    return data


def get_available_files(pointers, storage):
    pointers = Path(__file__).parent / pointers
    storage = Path(storage).expanduser()
    data = LargeFileStorage.load_pointers_file(pointers)
    files = list(data.keys())
    files_non_lfs = [
        f
        for f in os.listdir(storage)
        if f not in data and not os.path.isdir(storage / f)
    ]
    files += files_non_lfs
    return files
