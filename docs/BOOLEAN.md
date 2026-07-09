# Boolean Modelling

Boolean modelling utilities are exposed through `bt.bpy`.

```python
import bonesistools as bt
```

## Boolean Algebra

Boolean-algebra utilities are exposed through `bt.bpy.ba`.

### Partial Boolean Values

`PartialBoolean` represents a set-theoretic partial Boolean value:

```python
bt.bpy.ba.PartialBoolean(0)    # admissible values: {0}
bt.bpy.ba.PartialBoolean(1)    # admissible values: {1}
bt.bpy.ba.PartialBoolean("*")  # admissible values: {0, 1}
```

`"*"` means an unspecified or free Boolean value. It is not a Kleene truth
value.

### Hypercubes

`Hypercube` represents a partial Boolean configuration.

```python
cube = bt.bpy.ba.Hypercube({"A": 0, "B": "*"})

cube.contains({"A": 0, "B": 1})
# True
```

Missing components are interpreted as free:

```python
bt.bpy.ba.Hypercube({"A": 0}) == bt.bpy.ba.Hypercube({"A": 0, "B": "*"})
# True
```

### Configuration Sets

`ConfigurationSet` represents an exact finite set of complete Boolean
configurations over a fixed list of components.

The implementation is compact and hidden. Users interact with configurations,
not with the internal representation.

```python
states = bt.bpy.ba.ConfigurationSet(["A", "B"])
states.add({"A": 0, "B": 0})
states.add({"A": "*", "B": 0})
```

`len(...)` returns the number of represented complete configurations:

```python
len(states)
# 2
```

Membership accepts complete or partial configurations:

```python
{"A": 0, "B": 0} in states
# True

{"A": 1, "B": 0} in states
# True

{"B": 0} in states
# True

{"A": 1, "B": 1} in states
# False
```

Iteration yields complete configurations:

```python
list(states)
# [
#     {"A": 0, "B": 0},
#     {"A": 1, "B": 0},
# ]
```

Use `enumerate()` to materialize all configurations:

```python
states.enumerate()
# (
#     {"A": 0, "B": 0},
#     {"A": 1, "B": 0},
# )
```

Use `sample()` when the represented set is too large to enumerate:

```python
states.sample(seed=0)
# {"A": 0, "B": 0}

states.sample(3, seed=0)
# (
#     {"A": 0, "B": 0},
#     {"A": 1, "B": 0},
#     {"A": 0, "B": 0},
# )
```

### Redundancy Policy

`ConfigurationSet` keeps an exact, non-redundant internal representation, but
does not guarantee a globally minimal representation.

Equality is semantic: `a == b` compares the complete configurations
represented by each object, not the internal hypercubes used to store them.

For example:

```python
states = bt.bpy.ba.ConfigurationSet(["A", "B"])
states.add({"A": 0, "B": 0})
states.add({"A": "*", "B": 0})
```

The first configuration is covered by the second addition, so it does not stay
stored as a separate redundant subspace. The public result is:

```python
states.enumerate()
# (
#     {"A": 0, "B": 0},
#     {"A": 1, "B": 0},
# )
```

`compress()` can be called to opportunistically reduce the internal
representation without changing the represented configurations:

```python
states = bt.bpy.ba.ConfigurationSet(["A", "B"])
states.add({"A": 0, "B": 0})
states.add({"A": 1, "B": 0})
states.compress()

states.enumerate()
# (
#     {"A": 0, "B": 0},
#     {"A": 1, "B": 0},
# )
```

### Boolean Implicants

`dnf_implicants(...)` extracts implicants from an expression already written in
DNF-like form. `prime_implicants(...)` computes prime implicants.

```python
from boolean import BooleanAlgebra

ba = BooleanAlgebra()
rule = ba.parse("A & B | A & C")

bt.bpy.ba.prime_implicants(rule)
```

## Reachable Attractors

`BooleanNetwork.reachable_attractors(...)` returns exact reachable attractors as
`ConfigurationSet` objects.

```python
bn = bt.bpy.bn.BooleanNetwork({"A": "B", "B": "A"})

attractors = bn.reachable_attractors(
    {"A": 0, "B": 1},
    update="asynchronous",
)
```

Each attractor behaves as an exact set of complete configurations:

```python
for attractor in attractors:
    print(len(attractor))
    print(attractor.enumerate())
```

The initial state may be partial. Missing or free components define the set of
all compatible initial configurations, and the explicit backend explores the
union of reachable states only once:

```python
attractors = bn.reachable_attractors(
    {"A": 0},
    update="asynchronous",
)
```

Supported update semantics are:

- `"synchronous"`: update every unstable component at once;
- `"asynchronous"`: update one unstable component at a time;
- `"general"`: update any non-empty subset of unstable components;
- `"most-permissive"`: return reachable minimal trap spaces under
  most-permissive reachability.

For `"synchronous"`, `"asynchronous"` and `"general"` dynamics, the
`"explicit"` backend explores the reachable state space explicitly. The
optional `"bdd"` backend computes the reachable state set symbolically using
binary decision diagrams for these finite-state update semantics. Install it
with:

```bash
pip install "bonesistools[bdd]"
```

For `"most-permissive"` dynamics, the `"asp"` backend computes reachable
minimal trap spaces. The `"bdd"` backend is not used for most-permissive
reachability.

## Boolean Model I/O

Boolean-model readers are exposed through `bt.bpy.io`.

```python
model = bt.bpy.io.read_zginml("model.zginml")
graph = model.get("influence_graph")
bn = model.get("boolean_network")
```

### Readers

- `read_bnet(path)` reads a `.bnet` Boolean network.
- `read_bnet_directory(path)` reads a directory of `.bnet` files into a
  Boolean-network ensemble.
- `read_ginml(path)` reads a GINML logical model.
- `read_zginml(path)` reads a ZGINML archive and preserves companion files.

### Imported Logical Models

GINML and ZGINML files may contain more than a Boolean network. They can also
store an influence graph, initial states, perturbations, simulation settings,
layout data, annotations and other GINsim-specific sections.

`read_ginml(...)` and `read_zginml(...)` therefore return a model container
with:

- `boolean_network`: a `BooleanNetwork` when rules can be safely converted;
- `influence_graph`: an `InfluenceGraph` when signed interactions are present;
- `initial_states`: named `Hypercube` objects when states are Boolean;
- `perturbations`: parsed perturbation records when recognized;
- `metadata`: supported and unsupported model metadata.

Use:

```python
model.get("boolean_network")
model.get("influence_graph")
```

to request a required object explicitly. `get(...)` raises `ValueError` with a
clear reason when the requested object is unavailable.

### GINML/ZGINML Policy

GINML/ZGINML import is conservative. BoNesisTools does not try to fully
reimplement GINsim.

Signed interactions are converted to an `InfluenceGraph` when possible.
Logical rules are converted to a `BooleanNetwork` only when they are directly
representable or when a multi-valued component can be safely encoded with
monotone Boolean threshold variables:

```text
RA with maxvalue=2
    -> RA_b1
    -> RA_b2
```

For example, `RA:2` is converted to `RA_b1 & RA_b2`.

Unsupported logical parameters or companion files are not silently discarded:
they are reported in `metadata`. In that case, `boolean_network` may be `None`
while `influence_graph`, initial states and metadata remain available.
