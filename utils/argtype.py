#!/usr/bin/env python

from typing import Union

import argparse

class Range(argparse.Action):
    
    def __init__(
        self,
        min=None,
        max=None,
        *args,
        **kwargs
    ):
        self.min = min
        self.max = max
        kwargs.update({
            "nargs":None,
            "metavar":f"[{self.min}-{self.max}]"
        })
        super(Range, self).__init__(*args, **kwargs)

    def __call__(
        self,
        parser,
        namespace,
        value,
        option_string=None
    ):
        if not (self.min <= value <= self.max):
            raise argparse.ArgumentError(self, f"value {value} not in range [{self.min}-{self.max}].")
        setattr(namespace, self.dest, value)

class Store_boolean(argparse.Action):

    def __init__(
        self,
        *args,
        **kwargs
    ):
        kwargs.update({
            "type":str,
            "metavar":"BOOL"
        })
        super(Store_boolean, self).__init__(*args, **kwargs)

    def __call__(
        self,
        parser,
        namespace,
        value,
        option_string=None
    ):
        if value.lower() in ("0", "n", "no", "false"):
            value = False
        elif value.lower() in ("1", "y", "yes", "true"):
            value = True
        else:
            raise argparse.ArgumentTypeError("Boolean value expected.")
        setattr(namespace, self.dest, value)

class Store_prefix(argparse.Action):

    def __init__(
        self,
        *args,
        **kwargs
    ):
        kwargs.update({
            "type": str,
            "metavar": "LITERAL"
        })
        super(Store_prefix, self).__init__(*args, **kwargs)

    def __call__(
        self,
        parser,
        namespace,
        value,
        option_string=None
    ):
        if value:
            value = value if value[-1] in ["-","_"] else value + "_"
        setattr(namespace, self.dest, value)

class Store_dict(argparse.Action):

    def __init__(
        self,
        type_key: type=str,
        type_value: type=str,
        *args,
        **kwargs
    ):

        if isinstance(type_key, type):
            metavar_key = "LITERAL" if type_key == str else type_key.__name__.upper()
            self.type_key = type_key
        else:
            raise TypeError(f"`type_key` is of type {type(type_key)} instead of {type}") 

        if isinstance(type_value, type):
            metavar_value = "LITERAL" if type_value == str else type_value.__name__.upper()
            self.type_value = type_value
        else:
            raise TypeError(f"`type_key` is of type {type(type_value)} instead of {type}")
        
        if "nargs" not in kwargs:
            kwargs["nargs"] = "+"
        
        kwargs["metavar"] = f"{metavar_key}={metavar_value}"
        super(Store_dict, self).__init__(*args, **kwargs)

    def __call__(
        self,
        parser,
        namespace,
        values,
        option_string=None
    ):

        setattr(namespace, self.dest, dict())
        for element in values:
            key, value = element.split("=")
            key = self.type_key(key)
            value = self.type_value(value)
            getattr(namespace, self.dest)[key] = value

class Store_organism(argparse.Action):

    def __init__(
        self,
        *args,
        **kwargs
    ):
        default = kwargs["default"] if "default" in kwargs else "human"
        kwargs.update({
            "type": str,
            "metavar": "ORGANISM",
            "default": default,
            "help": kwargs["help"] if "help" in kwargs else f"common name or identifier of the organism of interest (default: {default})"
        })
        super(Store_organism, self).__init__(*args, **kwargs)

    def __call__(
        self,
        parser,
        namespace,
        value,
        option_string=None
    ):
        if value.is_integer():
            value = int(value)
        else:
            pass
        setattr(namespace, self.dest, value)