#!/usr/bin/env python

import os, io
import contextlib

import datetime

@contextlib.contextmanager
def disable_print(disable: bool=True):
    if disable is True:
        with open(os.devnull, "w") as f, contextlib.redirect_stdout(f):
            yield
    else:
        with io.StringIO() as f:
            yield

class Section(object):

    def __init__(
        self,
        init: int = 1,
        verbose: bool = True
    ):
        self.init = init
        self._i = init
        self._verbose = verbose
    
    def __call__(
        self,
        v: str,
        reset: bool = False
    ):
        self._i = self.init if reset else self._i
        if self._verbose is True:
            print(f"{self._i}) {v}")
        self._i+=1
        return None
    
    def reset(self):
        self._i = self.init
        return None
    
    def quiet(self):
        self._verbose = False
    
    def verbose(self):
        self._verbose = True

def print_task(message:str=None) -> None:
    print(f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]} - TASK - {message}")
    return None

def print_info(message:str=None) -> None:
    print(f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]} - INFO - {message}")
    return None

def print_warning(message:str=None) -> None:
    print(f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]} - WARNING - {message}")
    return None
