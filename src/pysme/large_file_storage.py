# -*- coding: utf-8 -*-
"""
System to store large data files on a server
Load them whem required by the user
Update the pointer file on github when new datafiles become available

Pro: Versioning is effectively done by Git
Con: Need to run server
"""

import gzip
import hashlib
import json
import logging
import os
import shutil
import tarfile
from dataclasses import dataclass
from html.parser import HTMLParser
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


class _HTMLDetectionParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.seen_html_tag = False

    def handle_starttag(self, tag, attrs):
        if tag.lower() == "html":
            self.seen_html_tag = True


@dataclass(frozen=True)
class _DownloadTarget:
    url: str
    md5: str | None = None
    sha256: str | None = None
    size: int | None = None


class LargeFileStorage:
    """
    Download large data files from data server when needed
    New versions of the datafiles are indicated in a pointer file.

    Pointer entries describe the downloadable object itself:
    - bare files are validated as bare files
    - `.gz` URLs are validated as the downloaded gzip payload
    - `.tar.gz` URLs are validated as the downloaded tarball

    Optional metadata in each pointer target can include `size`, `md5`,
    and/or `sha256`. Any validation fields that are present must match
    before the downloaded object is accepted.

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

        #:dict: maps tracked filenames to pointer targets for the downloadable object
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

    @staticmethod
    def _normalize_target_entry(target):
        if isinstance(target, str):
            return {"url": target}
        if isinstance(target, dict):
            if "url" not in target:
                raise ValueError(f"Pointer target dict must include 'url', got {target}")
            return dict(target)
        raise TypeError(f"Unsupported pointer target {target!r}")

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
    def _looks_like_html_error_page(fname):
        try:
            with open(fname, "rb") as f:
                sample = f.read(4096)
        except OSError:
            return False

        if not sample:
            return False

        stripped = sample.lstrip()
        if stripped[:1] != b"<":
            return False

        for encoding in ("utf-8", "latin-1"):
            try:
                text = stripped.decode(encoding, errors="ignore")
                break
            except Exception:
                text = None
        if not text:
            return False

        lowered = text.lower()
        if "<!doctype html" in lowered or "<html" in lowered:
            return True

        parser = _HTMLDetectionParser()
        try:
            parser.feed(text)
        except Exception:
            return False
        return parser.seen_html_tag

    @staticmethod
    def _compute_digest(fname, algorithm):
        digest = hashlib.new(algorithm)
        with open(fname, "rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()

    def _validate_download_payload(self, fname, target, file_format):
        if target.size is not None:
            actual_size = os.path.getsize(fname)
            if actual_size != target.size:
                raise ValueError(
                    f"Downloaded file size mismatch for {target.url}: expected {target.size}, got {actual_size}"
                )

        if target.md5 is not None:
            actual_md5 = self._compute_digest(fname, "md5")
            if actual_md5.lower() != str(target.md5).lower():
                raise ValueError(
                    f"Downloaded file MD5 mismatch for {target.url}: expected {target.md5}, got {actual_md5}"
                )

        if target.sha256 is not None:
            actual_sha256 = self._compute_digest(fname, "sha256")
            if actual_sha256.lower() != str(target.sha256).lower():
                raise ValueError(
                    f"Downloaded file checksum mismatch for {target.url}: expected {target.sha256}, got {actual_sha256}"
                )

        if file_format == "plain" and self._looks_like_html_error_page(fname):
            raise ValueError(
                f"Downloaded HTML page instead of data file from {target.url}"
            )

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
        targets = self.get_targets(key)
        errors = []

        for i, target in enumerate(targets):
            url = target.url
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
                self._validate_download_payload(fname, target, file_format)
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
                if i < len(targets) - 1:
                    logger.warning(
                        "Could not fetch %s from %s, trying fallback mirror",
                        key,
                        url,
                    )

        detail = "; ".join(f"{url}: {exc}" for url, exc in errors)
        raise FileNotFoundError(f"Could not fetch tracked file {key}. Attempts: {detail}")

    def get_targets(self, key):
        key = str(key)

        if key not in self.pointers:
            if key not in self.current:
                if not os.path.exists(key):
                    raise FileNotFoundError(
                        f"File {key} does not exist and is not tracked by the Large File system"
                    )
                return [_DownloadTarget(Path(key).as_uri())]
            return [_DownloadTarget((self.current / key).as_uri())]

        newest = self.pointers[key]
        raw_targets = newest if isinstance(newest, list) else [newest]

        targets = []
        seen = set()
        for raw_target in raw_targets:
            target = self._normalize_target_entry(raw_target)
            raw_url = str(target["url"]).strip()
            if raw_url == "":
                continue

            if self._is_uri(raw_url):
                candidate_urls = [raw_url]
            elif len(self.servers) > 0:
                candidate_urls = [self._join_uri(server, raw_url) for server in self.servers]
            else:
                candidate_urls = [raw_url]

            for url in candidate_urls:
                normalized = _DownloadTarget(
                    url=url,
                    md5=target.get("md5"),
                    sha256=target.get("sha256"),
                    size=target.get("size"),
                )
                if normalized in seen:
                    continue
                targets.append(normalized)
                seen.add(normalized)
        return targets

    def get_urls(self, key):
        """
        Return ordered candidate URLs/URIs for a tracked key.

        For tracked files:
        - pointer value may be a string or a list of strings
        - each pointer string may be a full URI (http/https/file) or relative path
        - relative paths are combined with every configured mirror server in order
        """
        return [target.url for target in self.get_targets(key)]

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
