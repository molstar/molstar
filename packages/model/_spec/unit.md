# Units: static props, dynamic data, the parent

`Unit` (`src/model/structure/structure/unit.ts`): a sorted set of atoms (or coarse elements) of one model, an operator
that places them, and data derived from them. Kinds: `Unit.Atomic`, `Unit.Spheres`, `Unit.Gaussians`.

## 1. Static and dynamic

**Static props** (`unit.props`) are derived from the unit's elements and model alone. They use model (invariant)
coordinates, so they do not depend on the unit's id or operator; the one exception is index-pair bonds listed per
operator (§4). One object is shared by the copies of one unit: symmetry copies (`applyOperator`), `getCopy`, the mates
of a symmetry group in `Structure.remapModel`; units with equal elements made otherwise (`Unit.create`, two equal
selections) have their own. It is filled lazily by the getters of the units that share it; everyone else reads it
(`Readonly`). Why: a 60-mer computes its lookups, bonds and boundaries once, not 60 times.

| Static prop                                                                               | Depends on                   |
| ----------------------------------------------------------------------------------------- | ---------------------------- |
| `polymerElements`, `gapElements`, `nucleotideElements`, `proteinElements`, `residueCount` | elements, hierarchy          |
| `boundary`, `lookup3d`, `principalAxes`                                                   | elements, coordinates        |
| `bonds`                                                                                   | elements, model, coordinates |
| `rings`, `resonance`                                                                      | `bonds`                      |

**The bond provider** (`unit.bondProvider`, atomic units; `unit/bonds/provider.ts`) says where the unit's intra-unit
bonds come from. It is a constructor argument, not a static prop: fixed when the unit is made, passed on by `getChild`,
`getCopy`, `applyOperator` and `remapModel`, and given by `Unit.create` (the built-in `BondProvider.DefaultInstance`
when not). Why: two reads of `unit.bonds` must return the same object, so the provider cannot be looked up when bonds
are read; and it is model-independent, so a remapped frame keeps it.

**Dynamic data** belongs to one unit instance: `conformation` (the operator mapping) and `transientCache` (valence
model, interaction features; keyed by the caller). Never inherited by children; `getCopy({ propagateTransientCache })`
is the one way to share it.

## 2. The parent

`unit.parent`: the unit `getChild` was called on, flattened to that unit's own parent, as `Structure.parent` is
flattened to the root. A child holds at most one other unit; no chains of ancestors build up.

- Same kind and model as the child (by construction: `getChild` keeps the model, `remapModel` drops a parent of another
  model). A superset of its elements is expected, not checked: `getChild` takes the caller's word, and a child whose
  elements are not all in the parent's computes its bonds.
- The operator may differ: the subset builders give a symmetry mate `pivot.applyOperator(...)`, which keeps the pivot's
  parent. So a child uses only the parent's static props, and not its bonds when those depend on the operator (§4).
- The root is the unit's origin, not always a unit of the structure's root: a symmetry builder that selects part of a
  unit (an assembly's chains out of merged waters) applies operators to a child, which keeps the asymmetric unit's unit
  as its parent.
- A child peeks at `parent.props` and never calls the parent's getters: a one-residue child must not compute the bonds
  of a 186k-atom unit. The one exception is a provider that asks for it (`childBonds: 'root'`, §4).
- Flattened, not immediate: a grandchild cannot use bonds only the intermediate unit computed, but holding the immediate
  parent would keep every ancestor's props alive. Bonds are usually computed on the root (full-unit components share the
  root's `Unit` objects).

## 3. Operations

| Operation                     | Static props                   | `parent`                                 | `transientCache`        |
| ----------------------------- | ------------------------------ | ---------------------------------------- | ----------------------- |
| `getChild(elements)`          | new, empty                     | `this.parent ?? this`                    | new                     |
| `getChild`, as many elements  | returns `this`                 | –                                        | –                       |
| `getCopy(id, …)`              | shared                         | same                                     | new, or shared if asked |
| `applyOperator(id, op)`       | shared                         | same                                     | new                     |
| `remapModel(m, dynamicBonds)` | new, filled by the rules below | kept if `m === this.model`, else dropped | new                     |
| `Structure.remapModel` mates  | the first mate's new object    | as `remapModel`                          | new                     |

`remapModel` keeps:

- topology props always (frames share the hierarchy);
- `boundary`, `lookup3d`, `principalAxes` if `Unit.isSameConformation` (through `Model.CoordinatesHistory` if the root
  model has one);
- `bonds` by `tryRemapBonds`: index-pair bonds while the index pairs are the same object, `canRemap` bonds, any bonds
  without `dynamicBonds`, else only on the same coordinates;
- `rings`, `resonance` exactly when it keeps `bonds`, as they are computed from them.

A remapped unit loses its parent because the parent still belongs to the old model. `Structure.remapModel` can receive
the remapped root (structure.md §3); passing the root's unit with the same id as the unit's parent is still open.

## 4. What a child derives from its parent

- **Bonds**: the induced subgraph of `parent.props.bonds` when they exist and the parent has the same bond provider
  (`subsetIntraUnitBonds`, `unit/bonds/subset.ts`); else `ElementSetIntraBondCache` under the provider's `cacheKey`,
  else the provider's `intra`, or `computeIntraUnitBonds` without one. Exact for the built-in bonds: every decision of
  the computation is per atom pair (struct_conn, the component dictionary, the hydrogen fallback, alt-locs, thresholds,
  index pairs). The exception is index-pair bonds listed per operator (`molstar_bond_site` with operator ids:
  `findIndexPairBonds`, `getIndexPairBonds`), which depend on the unit's operator: with index pairs, a child restricts
  only a parent of the same operator key. Edge order is the parent's; `canRemap` and `cacheable` are the parent's; bonds
  the parent carried over from another frame are carried over to the child. Restricted bonds are not put in
  `ElementSetIntraBondCache`: cheap to redo, and that cache keeps only the 64 sets used last.
- **Bonds of a provider** (bond-providers plan §2.3): a provider whose decisions are per atom pair keeps the rule above
  (`childBonds: 'self'`, the default). One whose decisions span rings or conjugated systems sets `childBonds: 'root'`: a
  child without restricted bonds computes the root unit's (`parent.bonds`) and restricts them, so a residue never
  differs from its chain, at the cost of the root unit once. `intra(unit)` depends on the model, the elements and its
  props only, which is what makes a `cacheKey` sound. The cache is one per model and key, the key under the provider's
  name (`intraCacheKey`), so two providers never share bonds; the built-in bonds keep the model's original cache for its
  life, the other keys only the last `MaxProviderKeys` (an LRU: a key per param value, a slider, would keep every
  value's bonds). No key caches nothing. A provider reads its unit's links to other units from `ctx.modelStructure()`
  (the model as a structure with the built-in bonds, once per model), never from a structure of its own (a unit is
  shared by structures), and finds rings over the graph it is building with `UnitRings.create(unit, bonds)`, as
  `unit.bonds` is not set while it runs.
- **Cost O(child)**: the parent index of each child element by binary search in the parent's elements, each neighbour's
  child index by binary search among those; a child with `n·log2(p) ≥ p` uses a merged walk and a parent-sized map
  instead, O(p). A parent-sized pass per small child is the hazard: 2000 one-residue children of a 186k-atom unit took
  1601 ms that way, 77 ms recomputed, 4 ms restricted.
- **Not derived**:
  - rings: found per residue by a depth-limited BFS over bonds. The parent's rings whose atoms are all in the child are
    likely exact for children of whole residues, not shown for partial ones;
  - polymer, gap, nucleotide and protein elements and the residue count are O(child) from the hierarchy already. A
    child's `polymerElements` holds trace atoms it may lack, so it is not `parent ∩ child`;
  - boundary, lookup3d and principal axes have no exact derivation cheaper than O(child).

## 5. Memory

One pointer per unit. A child keeps its parent unit alive, static props and transient cache included. Child structures
already keep their root structure (`Structure.parent`), whose units are these parents, so usually nothing more is
retained. Units whose root is not their structure's root (§2: symmetry built from part of a unit) keep one unit of the
discarded asymmetric unit alive each: at most one per unit, never a chain. `Loci.toStructure` of a root's loci makes a
structure without `Structure.parent` whose units still keep the root's units (and their grids and bonds) alive; keeping
such structures (e.g. per trajectory frame) keeps those. Coarse units keep a parent they never use.
