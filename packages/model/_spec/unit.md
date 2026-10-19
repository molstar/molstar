# Units: static props, dynamic data, the parent

`Unit` (`src/model/structure/structure/unit.ts`): a sorted set of atoms (or coarse elements) of one model, an operator
that places them, and data derived from them. Kinds: `Unit.Atomic`, `Unit.Spheres`, `Unit.Gaussians`.

## 1. Static and dynamic

**Static props** (`unit.props`) are derived from the unit's elements and model alone. They use model (invariant)
coordinates, so they do not depend on the unit's id or operator. One object is shared by every unit with the same
elements and model: symmetry copies (`applyOperator`), `getCopy`, the mates of a symmetry group in
`Structure.remapModel`. It is filled lazily by the getters of the units that share it; everyone else reads it
(`Readonly`). Why: a 60-mer computes its lookups, bonds and boundaries once, not 60 times.

| Static prop                                                                               | Depends on                   |
| ----------------------------------------------------------------------------------------- | ---------------------------- |
| `polymerElements`, `gapElements`, `nucleotideElements`, `proteinElements`, `residueCount` | elements, hierarchy          |
| `boundary`, `lookup3d`, `principalAxes`                                                   | elements, coordinates        |
| `bonds`                                                                                   | elements, model, coordinates |
| `rings`, `resonance`                                                                      | `bonds`                      |

**Dynamic data** belongs to one unit instance: `conformation` (the operator mapping) and `transientCache` (valence
model, interaction features; keyed by the caller). Never inherited by children; `getCopy({ propagateTransientCache })`
is the one way to share it.

## 2. The parent

`unit.parent`: the unit `getChild` was called on, flattened to that unit's own parent, as `Structure.parent` is
flattened to the root. A child holds at most one other unit; no chains of ancestors build up.

- Same kind and model as the child, and a superset of its elements (set by construction: `getChild` keeps the model,
  `remapModel` drops a parent of another model).
- The operator may differ: the subset builders give a symmetry mate `pivot.applyOperator(...)`, which keeps the pivot's
  parent. So a child uses only the parent's static props.
- A child peeks at `parent.props` and never calls the parent's getters: a one-residue child must not compute the bonds
  of a 186k-atom unit.
- Flattened, not immediate: a grandchild cannot use bonds only the intermediate unit computed, but holding the immediate
  parent would keep every ancestor's props alive. Bonds are usually computed on the root (full-unit components share the
  root's `Unit` objects).

## 3. Operations

| Operation                     | Static props                   | `parent`                                 | `transientCache`        |
| ----------------------------- | ------------------------------ | ---------------------------------------- | ----------------------- |
| `getChild(elements)`          | new, empty                     | `this.parent ?? this`                    | new                     |
| `getChild` with every element | returns `this`                 | –                                        | –                       |
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

A remapped unit loses its parent because the parent still belongs to the old model. Once `Structure.remapModel` receives
the remapped root (05-model-structure B1), it can pass the root's unit with the same id instead.

## 4. What a child derives from its parent

- **Bonds**: the induced subgraph of `parent.props.bonds` when they exist (`subsetIntraUnitBonds`,
  `unit/bonds/subset.ts`); else `ElementSetIntraBondCache`, else `computeIntraUnitBonds`. Exact: every decision of the
  computation is per atom pair (struct_conn, the component dictionary, the hydrogen fallback, alt-locs, thresholds,
  index pairs). Edge order is the parent's; `canRemap` and `cacheable` are the parent's; bonds the parent carried over
  from another frame are carried over to the child. Restricted bonds are not put in `ElementSetIntraBondCache`: cheap to
  redo, and that cache never evicts.
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
retained. Units of structures made without `parent` (`plugin/ui/src/sequence.tsx`) keep one root unit each: at most one
per child, never a chain.
