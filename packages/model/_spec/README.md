# Mol\* model: the specification

`@molstar/model` (`packages/model/`): molecular data without graphics. It reads files into models (atoms, coarse
elements, coordinates), builds structures from them (units, operators, assemblies, bonds), selects their parts (loci,
queries, MolScript), and holds volumes, particle lists, sequences and the custom-property framework. It depends on
`@molstar/core`, `@molstar/io` and `@molstar/query-language`, never on graphics or the plugin. This spec is for working
on the code without reading all of it: the model, the invariants the code relies on, and the decisions (the rule, then
why). Code comments cite sections (`unit.md §2`). The code is the truth for how. Open items: [../TODO.md](../TODO.md).

| File                           | What it covers                                                                          |
| ------------------------------ | --------------------------------------------------------------------------------------- |
| [model.md](model.md)           | `Model`: hierarchy, conformation, frames and what they share, property bags, trajectory |
| [formats.md](formats.md)       | files → `BasicData` → models; format properties; volumes, particles; writers, export    |
| [structure.md](structure.md)   | `Structure`: building, symmetry, children, lazy data, elements and loci, bonds          |
| [unit.md](unit.md)             | `Unit`: static props, dynamic data, the parent, what children derive                    |
| [properties.md](properties.md) | custom properties (model, structure, volume) and the computed ones                      |
| [query.md](query.md)           | the query engine, selections, MolScript compilation, `Script`                           |

## The model

| Piece              | Where                                        | What it is                                                                                    |
| ------------------ | -------------------------------------------- | --------------------------------------------------------------------------------------------- |
| `Model`            | `src/model/structure/model/model.ts`         | one model of a file: hierarchy (who), conformation (where), entities, sequence, property bags |
| `Trajectory`       | `src/model/structure/trajectory.ts`          | frames of models; `ArrayTrajectory`                                                           |
| `Unit`             | `src/model/structure/structure/unit.ts`      | a sorted set of a model's elements with an operator and derived static props                  |
| `Structure`        | `src/model/structure/structure/structure.ts` | units sorted by id, its root (`parent`), lazy derived data (bonds, boundary, symmetry groups) |
| `StructureElement` | `src/model/structure/structure/element/`     | `Location`, `Loci` (UnitIndex sets per unit), `Bundle`, `Stats`, `Schema`                     |
| query              | `src/model/structure/query/`, `src/script/`  | `StructureQuery` over a `QueryContext`, `StructureSelection`, MolScript → `QueryFn`           |
| properties         | `src/props/`                                 | `CustomModelProperty`, `CustomStructureProperty`, `CustomVolumeProperty`; computed providers  |
| formats            | `src/formats/`                               | readers to `Trajectory`, `Coordinates`, `Topology`, `Volume`, `ParticleList`; writers         |
| `Volume`           | `src/model/volume/`                          | a `Grid` (tensor, transform, stats) and instances                                             |
| `ParticleList`     | `src/model/particles/`                       | particles placing target structures, shapes or volumes                                        |

Index spaces: `ElementIndex` (an atom or coarse element of a model), `ResidueIndex`, `ChainIndex`, `EntityIndex`
(`indexing.ts`), and `StructureElement.UnitIndex` (a position in `unit.elements`). Loci, bonds and Locations of a unit
use `UnitIndex`; hierarchies and conformations use `ElementIndex`.

## Invariants

1. A model's atoms are sorted so that chains and residues are contiguous segments (`residueAtomSegments`,
   `chainAtomSegments`); `atomSourceIndex` maps back to file rows (model.md §2).
2. `unit.elements` is a sorted, duplicate-free `SortedArray<ElementIndex>`; `structure.units` is sorted by `unit.id`.
3. A unit keeps its `id` from the root structure to every child; loci, bundles, `remapModel` and the bond subsets key on
   it. A child's unit is a subset of the root's unit with the same id. `Structure.instances` breaks this (TODO).
4. Static unit props are in invariant (model) coordinates and shared by the copies of one unit (unit.md §1).
5. `Structure` and `Model` are immutable: lazy data is computed once and never invalidated; a change is a new object
   (`remapModel`, a new frame).
6. `structure.parent` is always the root (flattened at creation); `unit.parent` is the unit's root, flattened likewise
   (unit.md §2).

## Decisions

- **Hierarchy and conformation apart.** Frames of a trajectory and same-topology models share the hierarchy, static
  property data and `customProperties`; each has its own conformation and dynamic property data (model.md §3). A
  trajectory costs coordinates per frame, not a model.
- **Units share their atoms and their derived data.** An assembly or crystal is the asymmetric unit's units under more
  operators; symmetry copies share static props (unit.md §1). A 60-mer costs one set of lookups and bonds.
- **Children restrict, not recompute.** Child structures restrict the root's inter-unit bonds and child units the parent
  unit's intra-unit bonds when those exist (structure.md §5, unit.md §4). Selections, components and focus are children
  of a root that has usually computed them.
- **Loci in unit indices.** A selection is per-unit sets of `UnitIndex`, cheap to intersect and to map onto the same
  units of another structure; `Bundle` is the structure-independent form, `Schema` the file-identifier form.
- **Readers converge on mmCIF.** Every structure format builds `BasicData` (a subset of mmCIF tables) and one parser
  makes models from it (formats.md §1): one place for hierarchy, entities and sorting.
- **No graphics.** `Shape` here is structural; rendering is `@molstar/graphics`.
