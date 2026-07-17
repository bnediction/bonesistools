# Logic Implementation

This document describes the internal representations and algorithms used by
`bonesistools.logic`. It is intended for developers changing the logic engine.
Public call signatures and usage examples belong in docstrings, not here.

When an implementation strategy changes, update this document together with
the relevant regression and golden tests.

## Boolean Rules

`BooleanNetwork` stores rules as `boolean.py` expressions associated with one
shared `BooleanAlgebra`. Construction validates closure but deliberately keeps
the supplied expression structure. Simplification is explicit so importing or
copying a model does not silently rewrite its rules.

### Equivalence

Logical equivalence uses progressively more expensive exact checks:

1. simplify both expressions and compare their structures;
2. evaluate all assignments simultaneously with integer bitsets when the union
   of symbols contains at most 15 variables;
3. ask Clingo for a counterexample for larger expressions.

The ASP path encodes expression trees directly. It avoids constructing an
exponential truth table while preserving exact semantics.

### Prime Implicants

Prime-implicant computation also selects its strategy from the rule size:

- constant rules are handled directly;
- rules with at most four symbols use truth-table minterms followed by exact
  cube merging;
- larger rules are compiled to an ROBDD, then Clingo enumerates minimal partial
  assignments that force the selected value.

`dnf_implicants` has a different contract: it extracts clauses from a DNF-like
expression and does not complete them with consensus implicants. Dynamics code
must request prime implicants only when that minimality and completeness are
required.

## Boolean Subspaces

### Hypercubes

A `Hypercube` is the public mapping representation of a partial Boolean
configuration. Missing components and `PartialBoolean("*")` denote free
dimensions. It is convenient at API boundaries but expensive for large sets of
configurations.

### Configuration Sets

`ConfigurationSet` stores each internal hypercube as two Python integers:

```text
(fixed_mask, value_mask)
```

For component `i`, bit `i` in `fixed_mask` indicates whether the component is
fixed. The corresponding bit in `value_mask` stores its value. Bits outside
`fixed_mask` are ignored.

The representation maintains these invariants:

- it represents an exact finite set of complete configurations;
- stored hypercubes do not overlap;
- no stored hypercube is covered by another stored hypercube;
- component order is fixed for the lifetime of the object;
- the number of stored hypercubes is not guaranteed to be globally minimal.

Insertion subtracts the new hypercube from overlapping stored hypercubes, then
adds the new region. This makes counting, sampling and semantic equality exact
without materializing every configuration. `compress()` opportunistically
merges compatible adjacent regions but does not solve a global minimization
problem.

Dynamics code uses complete configurations as integer bitsets. Conversions to
mapping objects should remain at public boundaries or final result decoding.

### Symbolic Configuration Sets

`SymbolicConfigurationSet` stores one private BDD node and a strong reference
to its `SymbolicTransitionSystem`. It is an immutable working value for
composable symbolic operations, not a replacement for `ConfigurationSet`.

All symbolic sets produced by one transition system share the same manager,
state-variable order and transition partitions. Sets from different systems
must never be combined or converted implicitly, including systems compiled
from the same network. This identity constraint prevents cross-manager BDD
operations and gives the manager a clear lifetime.

Set algebra, inclusion, counting, containment and dynamic navigation operate
directly on BDD nodes. Iteration is lazy but can still enumerate exponentially
many configurations. Materialization is therefore explicit through
`configurations()`.

The materialization path calls `pick_iter` without care variables to obtain an
exact disjoint cube cover, translates those cubes directly to
`(fixed_mask, value_mask)` pairs, then calls
`ConfigurationSet._from_encoded_hypercubes(...)`. It must not use repeated
`add()`, decode intermediate `Hypercube` objects or call `compress()`.

### ROBDD

The internal ROBDD implementation represents one scalar Boolean function with
a fixed variable order, terminal nodes `0` and `1`, a unique-node table and
memoized Boolean operations. Reduction merges identical decision nodes and
eliminates nodes whose low and high branches are equal.

It is used where an exact compact representation of one rule is needed,
notably for large prime-implicant computations and non-unate most-permissive
rules. Finite-state symbolic dynamics instead use the required `dd` package,
which provides efficient relational operations over state sets.

## Finite-State Dynamics

Let `x` be a complete configuration, `f(x)` its synchronous image and
`U(x) = {i | f_i(x) != x_i}` its unstable components.

- synchronous semantics has the single successor `f(x)`;
- asynchronous semantics updates exactly one component in `U(x)`;
- general semantics updates any non-empty subset of `U(x)`.

Identity transitions are omitted when they do not change reachability or
terminal strongly connected components.

### Explicit Representation

Rules are compiled once into callables operating directly on integer bitsets.
For a configuration `x`, `x XOR f(x)` is the unstable-component mask.

- asynchronous successors are generated by repeatedly extracting the least
  significant set bit;
- general successors are generated by enumerating non-empty submasks;
- synchronous trajectories repeatedly apply the compiled image function until
  a configuration repeats.

Explicit asynchronous attractors use an iterative, on-the-fly Tarjan search.
The algorithm discovers terminal SCCs without first materializing the complete
reachable graph. General dynamics currently materialize the reachable
successor mapping and extract terminal SCCs with iterative Kosaraju. The
synchronous implementation memoizes the terminal cycle reached by each visited
configuration so trajectories with a shared suffix are traversed once.

### Symbolic BDD Representation

`SymbolicTransitionSystem` is the public name of the BDD transition engine. It
uses `dd.cudd` when available and falls back to `dd.autoref`. Current-state and
next-state variables are interleaved:

```text
x0, y0, x1, y1, ...
```

The transition relation is partitioned to avoid constructing one large BDD:

- synchronous and general semantics use conjunctive per-component clusters;
- asynchronous semantics use one disjunctive partition per updated component;
- current variables are existentially quantified after the last partition that
  depends on them;
- `and_exists` is used when the selected `dd` backend provides it.

General clusters encode either an unchanged component or its updated rule
value. The direct general transition relation additionally requires at least
one changed component. Asynchronous partitions use precomputed unchanged
prefixes and suffixes around the updated component.

Public symbolic `post()` and `pre()` preserve this non-reflexive general
relation without constructing a monolithic BDD. A lazily built changed-state
partition is composed with the existing per-component clusters;
`and_exists` and early quantification remain active. Forward and backward
closures may use the reflexive cluster relation because identity edges do not
change a fixed point.

Synchronous and general forward reachability, together with their backward
reachability, use frontier fixed points over BDD state sets. Asynchronous
forward and backward reachability use partition chaining with per-transition
worksets as described below.
The implementation never materializes an STG unless it intentionally switches
to the explicit path for a small symbolic result.

#### Asynchronous Reachability Chaining

The asynchronous relation is a disjunction of transition partitions
`T_0, ..., T_n`, one for each updated component. Applying every partition to a
global breadth-first frontier introduces a synchronization barrier: states
discovered by one partition cannot be consumed by another partition until the
next breadth-first round. Each round also constructs the union of every
partition's successors before removing configurations that were already
reached. On large symbolic regions, both the extra rounds and this transient
union can create expensive intermediate BDDs.

Forward closure therefore uses the transition-partition chaining described in
Section 2.5 of [4], refined with one workset per transition partition. Let `R`
contain every configuration reached so far and let `P_i` contain the
configurations already propagated through `T_i`:

```text
R = initial
P_i = false for every transition partition i

repeat
    for each transition partition T_i
        fresh = R - P_i
        P_i = P_i union fresh
        R = R union (post_i(fresh) intersect within)
until R no longer changes
```

The `within` intersection is omitted when no region is supplied. The
`P_i` sets ensure that each pair consisting of a reached configuration and a
transition partition is propagated at most once. States discovered by a later
partition remain fresh for earlier partitions and are handled on the next
pass. At termination, every partition maps `R` into `R`, so the result is the
same least forward-reachability closure as the global-frontier algorithm.

The performance gain comes from the order of symbolic operations, not from
discarding transitions or configurations. A partition updates `R` immediately,
so every later partition in the same pass can consume the newly reached states.
This immediate reuse is the defining property of chaining [4]. The `P_i`
worksets additionally avoid propagating a configuration through the same
partition twice. Subtracting `R` and merging new states after each partition
also avoids a monolithic union of all partition successors. The per-partition
`fresh` sets often retain more BDD sharing than separate graph-distance
frontiers, which reduces both relational-product work and peak live nodes.

This citation applies specifically to the chaining strategy in Section 2.5 of
[4]. It does not claim that this implementation reproduces the
decision-diagram saturation algorithm from Section 3.5. That algorithm
recursively saturates local MDD nodes according to event locality, whereas this
implementation remains a top-level worklist iteration over BDD transition
partitions.

The implementation is `_bdd_asynchronous_forward_chaining()` in
`_asynchronous.py`. Backward closure uses the symmetric recurrence with
`pre_i(fresh)` through `_bdd_asynchronous_backward_chaining()`. Its sweep
direction alternates after every productive pass, while worksets remain tied
to their transition labels. This reduces dependence on component order without
processing the same configuration-transition pair twice.

`_forward_reachable_states()` and `_backward_reachable_states()` dispatch to
these helpers solely for asynchronous semantics; synchronous and general
semantics retain their conjunctive-partition frontier algorithm. There is
deliberately no size threshold or public backend option.

This choice was validated on the 36 local ZGINML models. The benchmark covered
the all-zero configuration, all-one and alternating configurations, a
region-restricted closure, and up to two named partial initial conditions per
model:

| Workload | Paired cases | Frontier | Chaining | Speedup | Extra completions |
| --- | ---: | ---: | ---: | ---: | ---: |
| all-zero | 33 | 5.90 s | 0.99 s | 5.96-fold | 2 |
| varied and region-restricted | 82 | 17.34 s | 1.84 s | 9.42-fold | 12 |
| named partial conditions | 27 | 15.46 s | 0.83 s | 18.73-fold | 2 |
| **total** | **142** | **38.70 s** | **3.66 s** | **10.58-fold** | **16** |

An extra completion is a case where chaining terminated but the frontier
implementation exceeded its isolated timeout. The 35 measured regressions
were confined to already trivial closures and the largest absolute difference
was 0.139 milliseconds. For all 27 cases taking at least 10 milliseconds with
the frontier algorithm, chaining reduced peak live BDD nodes; the aggregate
reduction was 4.48-fold.

Result equality was checked directly at the BDD-node level on 34 ZGINML models
and on 700 randomly generated networks with two to eight components, including
partial initial configurations and 349 restricted regions. The two remaining
ZGINML equality checks exceeded the isolated verification timeout; every case
where both implementations completed had identical state counts.

Backward chaining was benchmarked separately on 181 target-region pairs from
the same ZGINML corpus. Among the 121 pairs completed by both implementations,
frontier iteration took 15.17 seconds and alternating chaining took 0.86
seconds, a 17.64-fold aggregate speedup. Every one of the 28 frontier cases
taking at least 10 milliseconds became faster, with a 13.41-fold median
speedup. Chaining completed 32 additional cases within the isolated timeout,
introduced no new timeout and reduced aggregate peak live BDD nodes 4.20-fold.
Exact BDD equality was additionally checked on 700 randomized networks with
partial target sets and regions.

`BooleanNetwork.symbolic()` captures a network copy and compiles one reusable
transition system. Subsequent changes to the source network do not alter the
compiled dynamics. The public system exposes no manager, BDD node, variable
name, renaming map or transition partition. Most-permissive dynamics remain in
their ASP/hypercube implementation and are deliberately excluded from this
finite-state symbolic API.

For asynchronous attractor detection, interleaved transition-guided reduction
(ITGR) treats each component update as a transition label and removes states
that cannot belong to a terminal SCC [1]. The reduced set is only guaranteed to
contain every attractor; it need not be forward-closed. Terminal SCC extraction
therefore uses the Xie-Beerel procedure against the original transition
relation rather than treating the reduced set as an induced closed graph [2].

### Reachability Quantifiers

For partial source and target hypercubes, the three exact contracts are:

- `exists`: at least one source configuration reaches at least one target
  configuration;
- `robust`: every source configuration reaches at least one target
  configuration;
- `universal`: every source configuration reaches every target configuration.

The universal implementation keeps source-target pairs symbolic by introducing
a third BDD variable family. It must not be replaced by independent set-level
reachability tests, which would lose the pairing information.

Concrete synchronous sources bypass BDD construction and follow their unique
trajectory directly. Partial synchronous sources and asynchronous/general
subspace queries use the symbolic relation.

### Attractors

Finite-state attractors are terminal SCCs of the reachable transition system.
All BDD terminal-SCC requests enter one strategy selector in
`SymbolicTransitionSystem`, whether they originate from
`BooleanNetwork.attractors()` or the composable symbolic API. If the analyzed
region contains at most 128 configurations, it is decoded and its induced SCCs
are computed explicitly. Larger regions remain symbolic:

- synchronous terminal components use either recurrent-state cycle extraction
  or the general symbolic SCC search;
- asynchronous candidates are reduced with ITGR and then enumerated with
  Xie-Beerel [1, 2];
- general terminal SCCs are isolated with alternating symbolic forward and
  backward closures;
- the backward basin of each terminal SCC is removed before searching for the
  next one.

The complete state space and every forward-reachability closure are transition
closed. This permits the asynchronous ITGR and synchronous cycle
specializations. For an arbitrary symbolic region that is not transition
closed, the selector retains the general induced-subgraph SCC algorithm; an
edge leaving the supplied region must not affect terminality inside that
region.

The synchronous choice is deliberately conservative. Networks with at least
90 components use the general symbolic SCC search. For smaller networks, the
same search is selected only when the conjunctive transition partitions have
an aggregate DAG size above 256 and BDD cofactor tests show that more than one
third of regulatory influences are non-unate. Other cases use deterministic
cycle extraction. These thresholds were selected from isolated benchmarks on
the 36 local ZGINML models and the golden models; changing them requires the
same time, peak-memory and result-equivalence checks.

Final state sets are decoded directly into exact disjoint bitset hypercubes and
used to construct `ConfigurationSet` objects in one operation. Avoid rebuilding
large results through repeated public `add()` and `compress()` calls.

The historical dynamics API and the composable symbolic API converge on the
same engine but not on each other's public wrappers:

```text
BooleanNetwork.transition/reachability/... ----┐
                                               ├─ SymbolicTransitionSystem
BooleanNetwork.symbolic() + symbolic sets -----┘
```

`BooleanNetwork.attractors()` selects only between the explicit and BDD engines
for finite-state semantics; most-permissive dynamics remain a separate
dispatch. Each engine owns its internal optimizations. Once the BDD path is
selected, every terminal-SCC choice belongs to `SymbolicTransitionSystem` and
is therefore shared with the symbolic API. Do not route historical calls
through `SymbolicConfigurationSet`; wrapper allocation would weaken their
specialized performance contracts.

### Automatic Strategy Selection

The public API does not expose finite-state backends. Selection is based on a
conservative upper bound: the number of free dimensions in the smallest closed
hypercube containing the initial subspace.

Current attractor selection is:

| Semantics | Explicit path | BDD path |
| --- | --- | --- |
| synchronous | at most 1,024 initial configurations | larger initial sets |
| asynchronous | at least 2,048 initial configurations always select BDD; otherwise at most 13 bounded free dimensions select explicit | remaining cases |
| general | at most 10 bounded free dimensions | remaining cases |

For reachable-configuration enumeration, asynchronous dynamics use the
explicit traversal up to 12 bounded free dimensions and general dynamics up to
8.

These values are conservative engineering heuristics, not semantic constants.
Any change requires benchmarks on small, medium and pathological models, plus
result-equivalence checks between explicit and symbolic implementations.

## Most-Permissive Dynamics

Most-permissive transition and reachability use the same relation. They are not
implemented by constructing a finite STG. The semantics and its reachability
characterization follow Pauleve et al. [3].

### Direct Reachability

One reachability query is encoded as an ASP decision problem. The network
program is grounded once and source/target values are supplied as assumptions,
which allows repeated robust or universal checks to reuse the compiled solver.

Unate rules are encoded through DNF implicants. Non-unate rules are encoded as
exact ROBDD decision facts so mixed positive and negative dependencies are not
approximated by an influence-graph monotonicity assumption.

### Reachable Configuration Enumeration

Enumeration uses a compact decomposition into transition regions. Each region
records:

- the components allowed to relax;
- the smallest hypercube closed by their Boolean functions;
- the free components in that hypercube;
- components whose change is irreversible within the region.

Regions are refined by removing irreversible components from subsequent
closure sets. ROBDD rule queries determine whether a component function can
take a value inside a partial hypercube. This representation avoids storing all
reachable configurations before iteration.

The smallest closed hypercube is also used as an inexpensive upper bound for
finite-state backend selection. It is a closure bound, not by itself the exact
set of most-permissive reachable configurations: irreversible-component
conditions are required for exact reachability.

### Trap Spaces and Attractors

Minimal trap spaces are enumerated with Clingo. Reachable most-permissive
attractors combine the trap-space and reachability constraints in one ASP
program and enumerate inclusion-minimal solutions with domain recursion
(`domRec`) and domain heuristics.

The principal trap space containing one configuration is computed directly by
iteratively relaxing values that are not preserved in the current hypercube;
it does not require global trap-space enumeration.

## Logical-Model Import

GINML and SBML Level 3 Qual parsers share `_LogicalModel`, an intermediate
representation containing maximum levels, inputs, desired threshold rules and
normalized component names. Format-specific parsing should end at this
boundary; threshold Booleanization belongs in `_booleanization.py`.

### Threshold Booleanization

A component `C` with maximum level `m > 1` is represented by:

```text
C_b1, C_b2, ..., C_bm
```

`C_bk = 1` means `C >= k`. Consequently, admissible encodings satisfy:

```text
C_b(k+1) -> C_bk
```

Desired threshold rules are regularized so lower thresholds are established
before a higher threshold can activate and existing higher activity remains
consistent while decreasing. Exact levels and level intervals are translated
to conjunctions/disjunctions over these monotone threshold variables.

Initial states use the same encoding: a source value `v` sets `C_bk` to
`int(v >= k)`. Component names are normalized for Boolean expressions;
normalization collisions are rejected rather than silently merged.

### Influence-Graph Consistency

When Booleanization succeeds, the imported influence graph has exactly the
Boolean-network components. Its edges are extracted syntactically without
resimplifying potentially large rules, so it can conservatively retain signed
influences that disappear from `BooleanNetwork.to_influence_graph()` after
simplification. Source layout, style and annotation attributes are attached to
the corresponding Boolean nodes and edges. Attributes of a multi-valued
component are copied to every threshold node, with origin and threshold
metadata retained.

If Booleanization cannot be justified, the Boolean network remains unavailable
and the parser falls back to signed source interactions when possible. It must
not invent logical rules from edge signs alone.

### Conservative Parsing

Unsupported GINsim companion sections, logical parameters and SBML Qual/MathML
constructs must not be interpreted heuristically. Recognized information is
stored in the structured fields of `ExecutableModel`; unsupported information
and failure reasons are retained in metadata whenever practical.

`ExecutableModel` stores protected copies of its network, graph, named initial
conditions, parameters, perturbations and metadata. Accessors return fresh
copies, while dynamic analyses delegate to the protected Boolean network. This
prevents external mutations from making the stored representations diverge.

`ExecutableModel.save()` writes a ZGINML archive. For an imported Boolean
GINML model, the writer reconstructs the retained logical parameters, graph
attributes, annotations, visual settings and recognized companion files. A
multi-valued source is instead exported as the Boolean threshold network held
by the executable model; initial conditions and fixed or interval
perturbations are translated to the same threshold encoding. Unknown companion
files cannot be reproduced because their byte content is intentionally not
retained, so their omission is reported explicitly.

Threshold nodes retain the source component's visual attributes but are placed
horizontally from the original coordinates using the node width plus a fixed
gap. This prevents all thresholds from overlapping in GINsim. Legacy `point`
visuals are converted to GINsim's default ellipse only along this Booleanization
path; visual XML reconstructed from an already Boolean source is left intact.

Only SBML Level 3 documents using the Qual package enter the SBML logical-model
pipeline. Other SBML documents fail before Booleanization.

## Validation

Algorithm changes should be validated at three levels:

1. regression tests for local invariants and edge cases;
2. cross-strategy checks that explicit, BDD and ASP formulations agree where
   their semantics overlap;
3. golden tests on fixed real models for end-to-end stability.

Performance changes should report both elapsed time and peak memory. A faster
implementation is accepted only after exact result comparison; output ordering
is not a correctness requirement unless a public contract explicitly says so.

## References

[1] Benes, N., Brim, L., Pastva, S., and Safranek, D. (2021). Computing bottom
SCCs symbolically using transition guided reduction. International Conference
on Computer Aided Verification, 505-528.

[2] Xie, A., and Beerel, P. A. (2000). Implicit enumeration of strongly
connected components and an application to formal verification. IEEE
Transactions on Computer-Aided Design of Integrated Circuits and Systems,
19(10), 1225-1230.

[3] Pauleve, L., Kolcak, J., Chatain, T., and Haar, S. (2020). Reconciling
qualitative, abstract, and scalable modeling of biological networks. Nature
Communications, 11(1), 4256.

[4] Ciardo, G., Marmorstein, R., & Siminiceanu, R. (2006). The saturation
algorithm for symbolic state-space exploration. International Journal on
Software Tools for Technology Transfer, 8(1), 4-25.
doi:10.1007/s10009-005-0188-7.
