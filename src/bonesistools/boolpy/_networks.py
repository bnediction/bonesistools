#!/usr/bin/env python

import importlib

from typing import List, Dict
from ._typing import MPBooleanNetwork
    
try:
    _mpbn_is_available = importlib.util.find_spec("mpbn") is not None
except:
    _mpbn_is_available = importlib.find_loader("mpbn") is not None

if _mpbn_is_available:
    from mpbn.minibn import struct_of_dnf
else:
    def struct_of_dnf(ba, f, container=frozenset, sort=False):
        import boolean.boolean as bpy
        def make_lit(l):
            if isinstance(l, ba.NOT):
                return (l.args[0].obj, False)
            else:
                return (l.obj, True)
        def make_clause(c):
            if isinstance(c, ba.AND):
                lits = c.args
            else:
                lits = [c]
            lits = map(make_lit, lits)
            return container(sorted(lits) if sort else lits)
        if isinstance(f, bpy._TRUE):
            return True
        elif isinstance(f, bpy._FALSE):
            return False
        if not isinstance(f, ba.OR):
            clauses = [f]
        else:
            clauses = f.args
        clauses = map(make_clause, clauses)
        return container(sorted(clauses) if sort else clauses)

class BooleanNetworkEnsemble(list):
    """
    Class storing ensemble of Boolean networks.

    Parameters
    ----------
    bns
        List of Boolean networks.
    """

    def __init__(
        self,
        components: List[str] = None,
        bns: List[MPBooleanNetwork] = None
    ) -> None:

        if components is None and bns is None:
            raise TypeError(f"unsupported type for instancing BooleanNetworkEnsemble: require either 'components' or 'bns'")
        elif components is not None and bns is None:
            super().__init__()
            self.__components = set(components)
        elif components is None and bns is not None:
            if not isinstance(bns, list):
                raise TypeError(f"unsupported type for instancing BooleanNetworkEnsemble: '{type(bns)}', not {list}")
            elif not all(isinstance(bn, MPBooleanNetwork) for bn in bns):
                raise TypeError(f"unsupported argument type: expected all elements being {MPBooleanNetwork}")
            else:
                components = set(bns[0])
                if not all(components == set(bn) for bn in bns[1:]):
                    raise ValueError("invalid value: different components between elements")
            super().__init__(bns)
            self.__components = components
        else:
            if not all(set(components) == set(bn) for bn in bns[1:]):
                raise ValueError("invalid value: different components between elements")
            super().__init__(bns)
            self.__components = set(components)

    def __setitem__(self, key, value) -> None:
        
        if isinstance(value, MPBooleanNetwork):
            super().__setitem__(key, value)
        else:
            raise TypeError(f"unsupported argument type: expected {MPBooleanNetwork}, but received {type(value)}")

    def append(self, other) -> None:
        
        if isinstance(other, MPBooleanNetwork):
            if set(other) == self.__components:
                super().append(other)
            else:
                raise ValueError("invalid value: missing or overflow components")
        else:
            raise TypeError(f"unsupported argument type: expected {MPBooleanNetwork}, but received {type(other)}")

    def insert(self, index, other) -> None:
        
        if isinstance(other, MPBooleanNetwork):
            if set(other) == self.__components:
                super().insert(index, other)
            else:
                raise ValueError("invalid value: missing or overflow components")
        else:
            raise TypeError(f"unsupported argument type: expected {MPBooleanNetwork}, but received {type(other)}")

    def get_components(self) -> set:
        
        return self.__components.copy()
    
    def get_clauses(self) -> Dict[str, set]:
        
        clauses = {node: [] for node in self[0]}
        for bn in self:
            for node in bn.keys():
                clauses[node].append(bn[node] if (bn[node] is True or bn[node] is False) else struct_of_dnf(bn.ba, bn[node]))

        return clauses
    
    def get_transcription_factors(self):
        # target => source

        def get_transcription_factors_from_clause(clause: frozenset) -> dict:
            transcriptions = {}
            for conjunction in clause:
                for factor, sign in conjunction:
                    if factor in transcriptions:
                        pass
                    else:
                        transcriptions[factor] = sign
            return transcriptions

        clauses = self.get_clauses()
        transcription_factors = {components: {} for components in self.__components}
        for component, clause_set in clauses.items():
            for clause in clause_set:
                if clause is True or clause is False:
                    pass
                else:
                    transcription_factors_for_one_node = get_transcription_factors_from_clause(clause)
                    for factor, sign in transcription_factors_for_one_node.items():
                        if factor in transcription_factors[component]:
                                if sign in transcription_factors[component][factor]:
                                    transcription_factors[component][factor][sign] += 1
                                else:
                                    transcription_factors[component][factor] = {sign: 1}
                        else:
                            transcription_factors[component][factor] = {sign: 1}

        return transcription_factors
    
    def get_influences(self):

        # source => target
        influences = {components: {} for components in self.__components}
        transcription_factors = self.get_transcription_factors()        
        for target, sources in transcription_factors.items():
            for source, influence in sources.items():
                influences[source][target] = influence
        return influences
