# -*- coding: utf-8 -*-
"""Base class for options classes
"""
from __future__ import annotations

import collections.abc
import re
import typing as _ty

_LINE_SEP = re.compile('\n {4,}')
# =============================================================================
# Fitter video options class
# =============================================================================


def _public(key: str) -> bool:
    """Is it a name of a public member?"""
    return not key.startswith('_')


def _fmt_sep(format_spec: str) -> _ty.Tuple[str, str, str]:
    """helper for Options.__format__: process `format_spec`."""
    if '#' not in format_spec:
        conv, next_spec = '', format_spec
    else:
        conv, next_spec = format_spec.split('#', maxsplit=1)
    sep = ',' + next_spec if next_spec else ', '
    conv = "!" + conv if conv else conv
    return sep, conv, next_spec


def _fmt_help(key: str, val: _ty.Any, conv: str, next_spec: str) -> str:
    """helper for Options.__format__: entry for one item"""
    if not conv == '!r' or _LINE_SEP.fullmatch(next_spec) is None:
        item = "{}={" + conv + "}"
        return item.format(key, val)
    val = repr(val).replace('\n', next_spec)
    return "{}={}".format(key, val)


# =============================================================================
# Base options class
# =============================================================================


# pylint: disable=too-many-ancestors
class Options(collections.abc.MutableMapping):
    """Base class for options classes

    The individual options can be accessed as object instance attributes
    (e.g. `obj.name`) or as dictionary items (e.g. `obj['name']`) for both
    getting and setting.

    If the name is not found, it will search the attributes whose names are
    listed in `mapping_attributes`. This does not apply to use as attributes
    (either `obj.name` or `getattr(obj, 'name')`). Iterating and unpacking do
    not recurse through these attributes, nor include properties and private
    attributes, unless their names are included in `property_attributes`.
    If an attribute's value is only set by a default value in a type hint, and
    not set in `__init__`, it will be omitted when iterating, unpacking or
    printing. If it is both a member of `self.__dict__` and listed in
    `prop_attributes`, it will appear twice.

    Setting attributes may be achieved via subscripting. If the name is not
    found, it will search the attributes listed in `mapping_attributes`. If
    a key is not found when recursing these attributes, a `KeyError` is raised.
    If an attribute's name is found in `mapping_attributes`, the attribute is
    updated when set rather than replaced like other attributes. These three
    statements do not apply to setting as an attribute. New keys may be added
    by setting as attributes, e.g. `obj.name=val` or `setattr(obj,'name',val)`.

    If a method `set_<name>(val)` exists, then `'<name>'` can be used as a key
    for setting but (unless it exists as a property or attribute) it cannot
    be used for getting. For such keys testing with `in` will return `False`,
    iteration will not include them and setting in a parent class will not
    propagate to this class.

    The attributes listed in `mapping_attributes` should be `MutableMapping`s.
    They, as well as the attributes in `property_attributes`, will raise a
    `TypeError` if you try to delete them.

    If the same item appears in more than one of the `mapping_attributes`, or
    in `self`, they can be partially synchronised by making it a property in
    the parent `Options` with a `set_<key>` method that updates the children.

    Parameters
    ----------
    All parameters are optional keywords. Any dictionary passed as positional
    parameters will be popped for the relevant items. Keyword parameters must
    be valid keys, otherwise a `KeyError` is raised.

    Raises
    ------
    KeyError
        If an invalid key is used when subscripting. This does not apply to use
        as attributes (either `obj.name` or `getattr(obj, 'name')`).
    """
    map_attributes: _ty.ClassVar[_ty.Tuple[str, ...]] = ()
    prop_attributes: _ty.ClassVar[_ty.Tuple[str, ...]] = ()

    def __init__(self, *args, **kwds) -> None:
        """The recommended approach to a subclass constructor is
        ```
        def __init__(self, *args, **kwds) -> None:
            self.my_attr = its_default
            self.other_attr = other_default
            ...
            self.last_attr = last_default
            order = ('my_attr', 'other_attr', ..., 'last_attr')
            args = sort_dicts(args, order, -1)
            kwds = sort_dict(kwds, order, -1)
            super().__init__(*args, **kwds)
        ```
        Mappings provided as positional arguments will be popped for the
        relevant items.
        """
        for mapping in args:
            self.pop_my_args(mapping)
        self.update(kwds)
        # for name in self.map_attributes:
        #     if isinstance(self[name], MasterOptions):
        #         raise TypeError(
        #             f"{name} is a {type(self[name]).__name__}, a subclass of"
        #             + "MasterOptions. It should not be used as a member of "
        #             + "another Options class, as it will block searches for "
        #             + "unknown keys when setting.")

    def __format__(self, format_spec: str) -> str:
        """formatted string representing object.

        Parameters
        ----------
        format_spec : str
            Formating choice. If it does not contain a `'#'` it is added to
            `","` as a separator and inserted before the first member.
            When it takes the form `'x#blah'`, any non `Options` members are
            converted as `"{}={!x}".format(key, val)`. `Options` members are
            converted as "{}={:x#blah   }".format(key, val)` if blah consists
            only of a newline followed by a minimum of four spaces, or
            "{}={!x:blah}".format(key, val)` otherwise.

        Returns
        -------
        str
            String reoesentation of object.
        """
        sep, conv, next_spec = _fmt_sep(format_spec)
        attrs = sep.join(_fmt_help(key, val, conv, next_spec)
                         for key, val in self.items())
        return type(self).__name__ + f"({sep[1:]}{attrs})"

    def __repr__(self) -> str:
        return self.__format__('r#\n    ')

    def __getitem__(self, key: str) -> _ty.Any:
        """Get an attribute"""
        try:
            return getattr(self, key)
        except AttributeError:
            for attr in self.map_attributes:
                try:
                    return getattr(self, attr)[key]
                except KeyError:
                    pass
        raise KeyError(f"Unknown key: {key}.")

    def __setitem__(self, key: str, value: _ty.Any) -> None:
        """Set an existing attribute"""
        if hasattr(self, 'set_' + key):
            getattr(self, 'set_' + key)(value)
        elif key in self.map_attributes:
            self[key].update(value)
        elif hasattr(self, key):
            setattr(self, key, value)
        else:
            for attr in self.map_attributes:
                if key in self[attr]:
                    getattr(self, attr).__setitem__(key, value)
                    return
                # try:
                #     getattr(self, attr).__setitem__(key, value)
                # except KeyError:
                #     pass
                # else:
                #     return
            raise KeyError(f"Unknown key: {key}.")

    def __delitem__(self, key: str) -> None:
        if key in self.map_attributes + self.prop_attributes:
            raise TypeError(f"`del {type(self).__name__}['{key}']` disallowed")
        try:
            delattr(self, key)
        except AttributeError:
            # raise KeyError(f"Unknown key: {key}.")
            for attr in self.map_attributes:
                try:
                    del getattr(self, attr)[key]
                except KeyError:
                    pass
                else:
                    return
            raise KeyError(f"Unknown key: {key}.")

    def __len__(self) -> int:
        # len(self.__dict__) + len(self.prop_attributes) includes privates.
        # tuple(self) appears to call len(self) -> infinite recursion.
        # return len(tuple(x for x in self))
        # barely any speed difference:
        count = 0
        for _ in self:
            count += 1
        return count

    def __iter__(self) -> _ty.Iterator[str]:
        yield from filter(_public, self.__dict__)
        yield from self.prop_attributes

    def pop_my_args(self, kwds: StrDict) -> None:
        """Pop any key from dict that can be set and use the value to set.
        """
        to_pop = []
        for key, val in kwds.items():
            try:
                self[key] = val
            except KeyError:
                pass
            else:
                to_pop.append(key)
        for key in to_pop:
            del kwds[key]


# =============================================================================
# Fallback options class
# =============================================================================


class AnyOptions(Options):
    """Same to `Options`, except it stores unknown keys as attributes.

    This can be used as a default place to store unknown items.
    """
    def __setitem__(self, key: str, val: _ty.Any) -> None:
        try:
            super().__setitem__(key, val)
        except KeyError:
            setattr(self, key, val)


# =============================================================================
# Master options class
# =============================================================================


class MasterOptions(Options):
    """Same to `Options`, except ot stores unknown keys in a mapping attribute.

    The name of the fallback mapping attribute is specified with the keyword
    `fallback` in the class definition.

    fallback : str|None
        Name of a mutable mapping attribute to store any completely unknown
        keys, i.e. if recursing through the mappings listed in `map_attributes`
        does not turn anything up. Any empty string indacates using `self` as
        the fallback storage, `None` indicates raising a `KeyError`.
    """
    def __init_subclass__(cls, fallback: _ty.Optional[str] = None) -> None:
        cls._fallback_mapping = fallback

    def __setitem__(self, key: str, val: _ty.Any) -> None:
        try:
            super().__setitem__(key, val)
        except KeyError:
            if self._fallback_mapping is None:
                raise
            if self._fallback_mapping:
                getattr(self, self._fallback_mapping).__setitem__(key, val)
            else:
                setattr(self, key, val)


# =============================================================================
# Helpers
# =============================================================================


def sort_dict(unordered: StrDict, order: _ty.Sequence[str],
              default: _ty.Optional[int] = None) -> StrDict:
    """Sort a dict by the order the keys appear in another list

    Parameters
    ----------
    unordered : Dict[str, Any]
        Dictionary whose entries we want to sort
    order : Sequence[str]
        Keys in order we want
    default : int|None, optional
        Sort keys for items that do not appear in `order`.
        By default `None` -> `len(order)`.

    Returns
    -------
    ordered : Dict[str, Any]
        Dictionary whose keys are in the same order as `order`
    """
    default = len(order) if default is None else default
    def key_fn(item: _ty.Tuple[str, _ty.Any]) -> int:
        """Key function for sorting"""
        return order.index(item[0]) if item[0] in order else default

    return dict(sorted(unordered.items(), key=key_fn))


def sort_dicts(unordered: _ty.Sequence[StrDict], order: _ty.Sequence[str],
               default: _ty.Optional[int] = None) -> _ty.List[StrDict]:
    """Sort a dict by the order the keys appear in another list

    Parameters
    ----------
    unordered : Sequence[Dict[str, Any]]
        Dictionary whose entries we want to sort
    order : Sequence[str]
        Keys in order we want
    default : int|None, optional
        Sort keys for items that do not appear in `order`.
        By default `None` -> `len(order)`.

    Returns
    -------
    ordered : Dict[str, Any]
        Dictionary whose keys are in the same order as `order`
    """
    default = len(order) if default is None else default
    return [sort_dict(arg, order, default) for arg in unordered]


# =============================================================================
# Hinting aliases
# =============================================================================
StrDict = _ty.Dict[str, _ty.Any]
Attrs = _ty.ClassVar[_ty.Tuple[str, ...]]
