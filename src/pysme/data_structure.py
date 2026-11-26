# -*- coding: utf-8 -*-
import logging
from copy import copy
from multiprocessing.sharedctypes import Value

import numpy as np
import pybtex
from flex.extensions.bindata import MultipleDataExtension
from flex.flex import FlexExtension

from . import __file_ending__, __version__, persistence
from .iliffe_vector import Iliffe_vector

logger = logging.getLogger(__name__)


def this(self, x):
    """This just returns the input"""
    return x


def notNone(func):
    def f(self, value):
        return func(self, value) if value is not None else None

    return f


def uppercase(func):
    def f(self, value):
        try:
            value = value.upper()
        except AttributeError:
            pass
        return func(self, value)

    return f


def lowercase(func):
    def f(self, value):
        try:
            value = value.casefold()
        except AttributeError:
            pass
        return func(self, value)

    return f


def oneof(*options, astype=None):
    def f(self, value):
        if value not in options:
            if astype is not None and astype(value) == value:
                return astype(value)
            raise ValueError(f"Received {value} but expected one of {options}")
        return value

    return f


def astype(func, allow_None=False):
    def f(self, value):
        if allow_None and value is None:
            return value
        if isinstance(value, func):
            return value
        return func(value)

    return f


asint = astype(int)
asfloat = astype(float)
asstr = astype(str)
asbool = astype(bool)


def absolute(self, value):
    try:
        return abs(value)
    except TypeError:
        logger.warning(
            "Tried setting number field with %s, trying to use the 0th element of it instead",
            type(value),
        )
        return abs(value[0])


def array(shape, dtype, allow_None=True):
    special_shape = False
    if isinstance(shape, str):
        shape = shape.split(",")
        for i, s in enumerate(shape):
            try:
                shape[i] = int(s)
            except ValueError:
                special_shape = True

    def f(self, value):
        if value is not None:
            value = np.asarray(value, dtype=dtype)
            value = np.atleast_1d(value)

            if shape is not None:
                if special_shape:
                    # Replace string placeholders with actual values
                    newshape = [
                        getattr(self, s) if isinstance(s, str) else s for s in shape
                    ]
                else:
                    newshape = shape
                value = value.reshape(newshape)
        elif not allow_None:
            raise ValueError(f"Received None but expected datatype {dtype}")
        return value

    return f


def vector(self, value):
    if value is None:
        return None
    elif np.isscalar(value):
        if self.wave is None:
            raise ValueError(
                "Can not set vector from single value without defining the wavelength grid"
            )
        wind = self.wave.sizes
        value = np.full(np.sum(wind), value)
        value = Iliffe_vector.from_indices(value, wind)
    elif isinstance(value, np.ndarray):
        value = np.require(value, requirements="W")
        if value.ndim == 1:
            wind = self.wave.sizes if self.wave is not None else None
            value = Iliffe_vector.from_indices(value, wind)
        else:
            value = Iliffe_vector(value)
    elif isinstance(value, list):
        value = Iliffe_vector(value)
    elif isinstance(value, Iliffe_vector):
        pass
    elif isinstance(value, FlexExtension):
        value = Iliffe_vector._load(value)
    elif isinstance(value, np.lib.npyio.NpzFile):
        value = Iliffe_vector._load_v1(value)
    else:
        raise TypeError("Input value is of the wrong type")

    return value


# Shorter versions, but less obvious (also "_name" is calculated each time, instead of once?)
# fget = lambda name, func: lambda self: func(getattr(self, f"_{name}"))
# fset = lambda name, func: lambda self, value: setattr(self, f"_{name}", func(self, value))


def fget(name, func):
    name = f"_{name}"

    def f(self):
        return func(self, getattr(self, name))

    return f


def fset(name, func):
    name = f"_{name}"

    def f(self, value):
        setattr(self, name, func(self, value))

    return f


def CollectionFactory(cls):
    """Decorator that turns Collection _fields into properties"""

    # Add properties to the class
    for name, _, setter, getter, doc in cls._fields:
        setattr(
            cls,
            name,
            property(fget(name, getter), fset(name, setter), None, doc),
        )
    cls._names = [f[0] for f in cls._fields]

    return cls


@CollectionFactory
class Collection(persistence.IPersist):
    _fields = [
        (
            "citation_info",
            "",
            asstr,
            this,
            "str: Bibtex representation of the citation",
        )
    ]  # [("name", "default", str, this, "doc")]

    def __init__(self, **kwargs):
        for name, default, *_ in self._fields:
            setattr(self, name, copy(default))

        for key, value in kwargs.items():
            if key in self._names:
                if isinstance(value, bytes):
                    value = value.decode()
                elif isinstance(value, np.ndarray) and value.dtype == np.dtype("O"):
                    value = value.astype(str)
                setattr(self, key, value)

    def __getitem__(self, key):
        key = key.casefold()
        return getattr(self, key)

    def __setitem__(self, key, value):
        key = key.casefold()
        setattr(self, key, value)

    def __contains__(self, key):
        return key in dir(self) and getattr(self, key) is not None

    def __str__(self):
        with np.printoptions(threshold=10):
            lines = []
            for name in self._names:
                lines += [f"{name}={self[name]!r}"]
            lines = ",\n".join(lines)
            output = f"{self.__class__.__name__}({lines})"
        return output

    def __repr__(self):
        return self.__str__()

    def create_citation(self, citation_info, output="string", sorting_style=None):
        if output == "bibtex":
            return citation_info
        elif output == "string":
            output = "plaintext"

        return pybtex.format_from_string(
            citation_info,
            style="plain",
            output_backend=output,
            sorting_style=sorting_style,
        )

    def citation(self, output="string"):
        return self.create_citation(self.citation_info, output=output)

    def to_dict(self):
        tree = {}
        for name in self._names:
            value = self[name]
            if hasattr(value, "to_dict"):
                value = value.to_dict()
            tree[name] = value
        return tree

    @classmethod
    def from_dict(cls, data):
        obj = cls(**data)
        return obj

    def _save(self):
        header = {}
        data = {}

        for name in self._names:
            value = self[name]
            if isinstance(value, np.ndarray) and value.size > 10:
                data[name] = value
            else:
                header[name] = value

        ext = MultipleDataExtension(header, data)
        return ext

    @classmethod
    def _load(cls, ext: MultipleDataExtension):
        header = ext.header
        if isinstance(ext, MultipleDataExtension):
            data = dict(ext.data)
            header.update(data)

        obj = cls(**header)
        return obj
