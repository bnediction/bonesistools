from typing import (
    Union,
    List
)

try:
    from typing import Self
except:
    from typing_extensions import Self # type: ignore

class Hypercube(dict):

    def __init__(self, mapping: dict) -> Self:

        for v in mapping.values():
            if v not in [0,1,"*"]:
                raise ValueError("unsupported initialization: value not equal to 0, 1 or '*'")
        super().__init__(mapping)

    @property
    def update(self, other: Self) -> None:

        if not isinstance(other, dict):
            raise TypeError(f"unsupported method type for 'update': '{type(self)}' and '{type(other)}'")
        elif not isinstance(other, Hypercube):
            return super().update(Hypercube(other))
        else:
            return super().update(other)
    
    @property
    def drop(
        self,
        keys: List[str],
        inplace: bool=False
    ) -> Union[Self, None]:

        if inplace is False:
            self = self.copy()

        for k in keys:
            del self[k]
        
        return self if inplace is False else None
    
    def is_fixed_point(self) -> bool:

        return False if "*" in self.values() else True

    def __eq__(self, other: Self) -> bool:
        
        if not isinstance(other, dict):
            raise TypeError(f"'==' not supported between instances of '{type(self)}' and '{type(other)}'")
        elif self.keys() != other.keys():
            raise ValueError("'==' not supported for different components")
        elif not isinstance(other, Hypercube):
            other = Hypercube(other)
        
        return super().__eq__(other)

    def __ne__(self, other: Self) -> bool:
        
        if not isinstance(other, dict):
            raise TypeError(f"'!=' not supported between instances of '{type(self)}' and '{type(other)}'")
        elif self.keys() != other.keys():
            raise ValueError("'==' not supported for different components")
        elif not isinstance(other, Hypercube):
            other = Hypercube(other)
        
        return super().__ne__(other)

    def __le__(self, other: Self) -> bool:

        if not isinstance(other, dict):
            raise TypeError(f"'<=' not supported between instances of '{type(self)}' and '{type(other)}'")
        elif self.keys() != other.keys():
            raise ValueError("'<=' not supported for different components")
        elif not isinstance(other, Hypercube):
            other = Hypercube(other)

        for c, v1 in self.items():
            v2 = other[c]
            if v2 != "*" and v1 != v2:
                return False
        return True
    
    def __lt__(self, other: Self) -> bool:

        if not isinstance(other, dict):
            raise TypeError(f"'<' not supported between instances of '{type(self)}' and '{type(other)}'")
        elif self.keys() != other.keys():
            raise ValueError("'<' not supported for different components")
        elif not isinstance(other, Hypercube):
            other = Hypercube(other)

        smaller = False
        for c, v1 in self.items():
            v2 = other[c]
            if v2 == "*" and v1 in [0,1]:
                smaller = True
            if v2 != "*" and v1 != v2:
                return False
        return smaller
    
    def __ge__(self, other: Self) -> bool:

        if not isinstance(other, dict):
            raise TypeError(f"'>=' not supported between instances of '{type(self)}' and '{type(other)}'")
        elif self.keys() != other.keys():
            raise ValueError("'>=' not supported for different components")
        elif not isinstance(other, Hypercube):
            other = Hypercube(other)

        for c, v1 in self.items():
            v2 = other[c]
            if v1 != "*" and v1 != v2:
                return False
        return True

    def __gt__(self, other: Self) -> bool:

        if not isinstance(other, dict):
            raise TypeError(f"'>' not supported between instances of '{type(self)}' and '{type(other)}'")
        elif self.keys() != other.keys():
            raise ValueError("'>' not supported for different components")
        elif not isinstance(other, Hypercube):
            other = Hypercube(other)

        larger = False
        for c, v1 in self.items():
            v2 = other[c]
            if v1 == "*" and v2 in [0,1]:
                larger = True
            if v1 != "*" and v1 != v2:
                return False
        return larger

    is_smaller_than = __le__
    is_larger_than = __ge__

    def identical(self, other: Self) -> set:
        
        if not isinstance(other, dict):
            raise TypeError(f"unsupported method type for 'intersection': '{type(self)}' and '{type(other)}'")
        elif self.keys() != other.keys():
            raise ValueError("invalid method value: different components")
        elif not isinstance(other, Hypercube):
            other = Hypercube(other)
        
        subset = set()
        for c, v1 in self.items():
            v2 = other[c]
            if v1 == v2:
                subset.add(c)
        return subset

    def different(self, other: Self) -> set:
        
        if not isinstance(other, dict):
            raise TypeError(f"unsupported method type for 'different': '{type(self)}' and '{type(other)}'")
        elif self.keys() != other.keys():
            raise ValueError("invalid method value: different components")
        elif not isinstance(other, Hypercube):
            other = Hypercube(other)
        
        subset = set()
        for c, v1 in self.items():
            v2 = other[c]
            if v1 != v2:
                subset.add(c)
        return subset

class HypercubeCollection(list):

    def __init__(
        self,
        hypercubes: List[Hypercube] = None
    ) -> None:

        if hypercubes is None:
            super().__init__()
        elif not isinstance(hypercubes, list):
            raise TypeError(f"unsupported type for instancing HypercubeCollection: '{type(hypercubes)}', not {list}")
        else:
            super().__init__([hypercube if isinstance(hypercube, Hypercube) else Hypercube(hypercube) for hypercube in hypercubes])
#        for hypercube in hypercubes:
#            if isinstance(hypercube, Hypercube):
#                super().append(hypercube)
#            else:
#                super().append(Hypercube(hypercube))

    def append(self, other):

        if isinstance(other, dict):
            super().append(other if isinstance(other, Hypercube) else Hypercube(other))
        elif isinstance(other, list):
            for hypercube in other:
                super().append(hypercube if isinstance(other, Hypercube) else Hypercube(other))
            super().append([hypercube if isinstance(hypercube, Hypercube) else Hypercube(hypercube) for hypercube in self])

    def are_fixed_points(self):

        fixed_points = HypercubeCollection()
        for hypercube in self:
            if hypercube.is_fixed_point():
                fixed_points.append(hypercube.copy())
        return fixed_points

    def are_smaller_than(self, other):

        if not isinstance(other, dict):
            raise TypeError(f"unsupported method types for are_smaller_than: '{type(self)}' and '{type(other)}'")
        elif not isinstance(other, Hypercube):
            other = Hypercube(other)

        smaller_hypercubes = []
        for hypercube in self:
            if hypercube.is_smaller_than(other):
                smaller_hypercubes.append(hypercube.copy())
        return smaller_hypercubes

    def are_larger_than(self, other):

        if not isinstance(other, dict):
            raise TypeError(f"unsupported method types for are_larger_than: '{type(self)}' and '{type(other)}'")
        elif not isinstance(other, Hypercube):
            other = Hypercube(other)

        larger_hypercubes = []
        for hypercube in self:
            if hypercube.is_larger_than(other):
                larger_hypercubes.append(hypercube.copy())
        return larger_hypercubes