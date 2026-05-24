#!/usr/bin/env python

import builtins
import importlib
import sys
from pathlib import Path

import numpy as np
import pytest
from anndata import AnnData

from bonesistools.sctools import _typing

_TYPING_PATH = Path(_typing.__file__)


def _load_typing_module(
    monkeypatch,
    module_name,
    *,
    force_no_mudata=False,
    force_legacy_importlib=False,
    force_typing_literal_fallback=False,
):
    spec = importlib.util.spec_from_file_location(module_name, _TYPING_PATH)
    module = importlib.util.module_from_spec(spec)

    if force_no_mudata and not force_legacy_importlib:
        find_spec = importlib.util.find_spec

        def find_spec_without_mudata(name, *args, **kwargs):
            if name == "mudata":
                return None
            return find_spec(name, *args, **kwargs)

        monkeypatch.setattr(importlib.util, "find_spec", find_spec_without_mudata)

    if force_legacy_importlib:
        monkeypatch.setattr(importlib, "find_loader", lambda name: None, raising=False)
        monkeypatch.delattr(importlib, "util", raising=False)

    if force_typing_literal_fallback:
        real_import = builtins.__import__
        raised = False

        def import_without_typing_literal(
            name,
            globals=None,
            locals=None,
            fromlist=(),
            level=0,
        ):
            nonlocal raised
            fromlist = fromlist or ()
            if name == "typing" and "Literal" in fromlist and not raised:
                raised = True
                raise ImportError("cannot import name 'Literal'")
            return real_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", import_without_typing_literal)

    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        sys.modules.pop(module_name, None)

    return module


def test_type_checker_accepts_direct_and_decorator_usage():
    @_typing.type_checker(value=_typing.UnionType(int, str))
    def stringify(value):
        return str(value)

    @_typing.type_checker(value=int)
    def double(value):
        return value * 2

    @_typing.type_checker
    def no_options(value):
        return value

    assert str(_typing.UnionType(int, str)) == "int or str"
    assert stringify(1) == "1"
    assert stringify(value="A") == "A"
    assert double(2) == 4

    with pytest.raises(TypeError, match="unsupported argument type for 'value'"):
        stringify(1.2)

    with pytest.raises(TypeError, match="unsupported argument type for 'value'"):
        double("2")

    with pytest.raises(Exception, match="Expected verification arguments"):
        no_options("A")


def test_anndata_checkers_validate_checked_arguments():
    adata = AnnData(np.ones((2, 2)))

    @_typing.anndata_checker(n=1)
    def n_obs(scdata):
        return scdata.n_obs

    assert n_obs(adata) == 2

    n_vars = _typing.anndata_or_mudata_checker(lambda scdata: scdata.n_vars)

    assert n_vars(adata) == 2

    @_typing.anndata_or_mudata_checker(n=1)
    def n_obs_from_factory(scdata):
        return scdata.n_obs

    assert n_obs_from_factory(adata) == 2

    with pytest.raises(TypeError, match="unsupported argument type"):
        n_obs(object())

    with pytest.raises(TypeError, match="unsupported argument type"):
        n_vars(object())


def test_mudata_checker_matches_optional_dependency_state():
    @_typing.mudata_checker(n=1)
    def identity(scdata):
        return scdata

    if not _typing._mudata_is_available:
        with pytest.raises(ModuleNotFoundError, match="mudata"):
            identity(object())
        return

    from mudata import MuData

    with pytest.warns(FutureWarning):
        mudata = MuData({"rna": AnnData(np.ones((2, 2)))})

    assert identity(mudata) is mudata

    with pytest.raises(TypeError, match="unsupported argument type"):
        identity(object())


def test_anndata_or_mudata_checker_falls_back_to_anndata_when_mudata_is_missing(
    monkeypatch,
):
    typing_module = _load_typing_module(
        monkeypatch,
        "bonesistools_sctools_typing_without_mudata",
        force_no_mudata=True,
    )
    adata = AnnData(np.ones((2, 2)))

    assert typing_module._mudata_is_available is False
    assert typing_module.ScData is typing_module.AnnData

    @_typing.anndata_checker(n=1)
    def expected_behavior(scdata):
        return scdata.n_obs

    @typing_module.anndata_or_mudata_checker(n=1)
    def n_obs(scdata):
        return scdata.n_obs

    n_vars = typing_module.anndata_or_mudata_checker(lambda scdata: scdata.n_vars)

    assert n_obs(adata) == expected_behavior(adata)
    assert n_vars(adata) == 2

    with pytest.raises(TypeError, match="unsupported argument type"):
        n_obs(object())

    with pytest.raises(TypeError, match="unsupported argument type"):
        n_vars(object())

    @typing_module.mudata_checker(n=1)
    def unavailable_mudata(scdata):
        return scdata

    with pytest.raises(ModuleNotFoundError, match="mudata"):
        unavailable_mudata(adata)

    with pytest.raises(ModuleNotFoundError, match="mudata"):
        typing_module.mudata_checker(lambda scdata: scdata)(adata)


def test_typing_module_supports_import_compatibility_fallbacks(monkeypatch):
    literal_module = _load_typing_module(
        monkeypatch,
        "bonesistools_sctools_typing_literal_fallback",
        force_typing_literal_fallback=True,
    )
    legacy_module = _load_typing_module(
        monkeypatch,
        "bonesistools_sctools_typing_legacy_importlib",
        force_legacy_importlib=True,
    )

    assert legacy_module._mudata_is_available is False
    assert literal_module.Literal is not None
