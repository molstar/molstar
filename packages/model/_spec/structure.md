# Structures

`Structure` (`src/model/structure/structure/structure.ts`): units sorted by id, a root, and lazy derived data. Units:
[unit.md](unit.md).

## 1. Building

- `Structure.create(units, props)`: sorts `units` by id (in place), builds `unitMap` and `unitIndexMap`. Props:
  `parent`, `interUnitBonds`, `dynamicBonds`, `bondProvider`, `interBondsValidUnit(Pair)`, `coordinateSystem`, `label`,
  `masterModel`, `representativeModel`; what is not given is inherited from `parent`.
- `bondProvider` (`unit/bonds/provider.ts`; when neither given nor inherited, the first atomic unit's, so a root made of
  another's units (`Loci.toStructure` of a root's loci) keeps their bonds' source; else `BondProvider.DefaultInstance`):
  the structure's inter-unit bonds come from it, and `StructureBuilder.addUnit` gives it to the units it makes (unit.md
  §1). Paths that make a new root without a parent pass it on: `remapModel`, `ofTrajectory` (its props),
  `mergeUnitsByOperator` and `mergeUnitsWithSameOperator` (the merged unit takes its first unit's), and the
  `StructureSymmetry` builders. A structure over units of several providers is allowed: each unit's decides its
  intra-unit bonds.
- `ofModel(model)`: one atomic unit per chain (operator from `atomicChainOperatorMappinng`); consecutive water chains
  with one operator merge (`Trait.Water`), as do consecutive single-atom chains of one entity and auth asym id
  (`Trait.MultiChain`); one coarse unit per coarse chain.
- `ofTrajectory`: every frame's units, `invariantId` and `chainGroupId` offset per frame, so frames never share symmetry
  groups.
- `StructureBuilder`: `addUnit`, `copyUnit`, `addWithOperator`, chain groups; single-element units at the same position
  with one `invariantId` are dropped (degenerate spacegroups).
- `transform(s, m)` (rigid only; same unit ids, `parent: s`), `instances(s, ms)`, `mergeUnitsByOperator(s)` (one unit
  per model, kind and operator; no parent).

## 2. Symmetry

`StructureSymmetry` (`symmetry.ts`), tasks over a structure of one model: `buildAssembly` (the assembly's operator
groups select by query), `buildSymmetryAssembly`, `builderSymmetryMates` (radius), `buildSymmetryRange` (cells),
`buildNcs`. Results carry the label, `dynamicBonds` and `bondProvider` and have no parent: an assembly is a new root
(its units may still have a unit parent, unit.md §2).

Units that differ only by operator form a `Unit.SymmetryGroup` (`unitSymmetryGroups`: equal `invariantId`, model and
elements). Builders keep a group's units sharing static props: `StructureSymmetry` adds each source unit under each
operator (`addWithOperator`, composed), the subset builders apply the operator of a group's other units to its first
child (`applyOperator(id, op, true)`). Graphics draws a group as one geometry with instances.

## 3. Children and the root

- `state.parent` is flattened to the root (`props.parent.parent ?? props.parent`); `root` is `parent ?? this`.
  Root-scoped properties live on the root (properties.md §2).
- Children come from `subsetBuilder`/`StructureSubsetBuilder` (reuses the parent unit when a unit is whole), selections,
  `Loci.toStructure`, `Bundle.toStructure`, `Schema.toStructure`, `transform`. Their units are `getChild` of the root's
  (unit.md §2).
- `asParent()`: a proxy of the parent in the child's coordinate system, for representations that need the context of the
  whole.
- `remapModel(m, parent?)`: remaps each symmetry group's first unit and passes its props to the others; keeps
  `interUnitBonds` unless `dynamicBonds`. `parent` is the root already remapped to `m`: the components of a trajectory
  frame share it, and with it the root-scoped properties and `areRootsEquivalent`; a unit with the elements of the
  parent's unit of its id and its operator is that unit (as `getChild` gives it), so a whole-unit component shares the
  root's units and their boundaries and lookups; a unit moved from the root's (a component of `Structure.transform`'s
  result) is remapped on its own. Without it the parent is remapped too, a root of the component's own.

## 4. Lazy data

Memoized in `state` (sentinel -1 or `undefined`), never invalidated: counts (`elementCount` eager; residues, polymer
units, gaps, unique elements), hashes (`hashCode`: unit ids and element sets; `transformHash`: unit ids), `boundary`
(from the unit spheres, `util/boundary.ts`), `lookup3d` (a grid of unit spheres; queries go to each unit's grid in its
frame), `unitSymmetryGroups`, `interUnitBonds`, `carbohydrates`, `coordination`, `models`, serial and bond mappings,
property bags (`currentPropertyData`; `inheritedPropertyData` is the root's).

`model` is the representative, the master or the single model, and throws for several models without either.

## 5. Bonds

- Intra-unit: unit.md §4 and `unit/bonds/intra-compute.ts`: index pairs if the model has them; else struct_conn, the
  component dictionary per residue (not at microheterogeneity; pairs within `maxRadius`, from the residue's atoms or,
  above 256 atoms, the grid), then distance thresholds (`common.ts`) over the unit's grid, queried at 3.1 Å, above every
  threshold (coarse-grained models: `maxRadius`). Edge order follows the grid. Edge props `order`, `flags`, `key` (the
  struct_conn row, index-pair key or dictionary key; -1 computed). `OrderUnknown`: a covalent struct_conn bond whose
  order neither the file nor the inter-residue table gives (order 1), for providers that perceive orders.
- Both come from the bond provider when it has `intra` (a unit's) or `inter` (the structure's), which may call the
  built-in computation through its context (`ctx.builtIn(props)`); else from the built-in computation below. The generic
  providers are in `unit/bonds/providers.ts`: `'none'` (the built-in computation with `noDistance`) and `'per-unit'`
  (other providers by rules on a unit's first element: entity type, chain, entity; bond-providers plan §5).
- `noDistance` (`BondComputationProps`): only the bonds the format states (index pairs, struct_conn, the component
  dictionary); the distance pass over the grid is skipped, intra and inter, as is the distance fallback for hydrogens
  the dictionary lacks. The rest of the computation is the same pass, so the result is the default bonds without the
  `Computed` ones.
- Inter-unit (`inter-compute.ts`): unit pairs from `eachUnitPair` (structure grid, `validUnit`, `validUnitPair`: same
  model or the master model; waters and ions skipped unless `dynamicBonds` or a common struct_conn record). Per pair,
  A's atoms near B come from A's grid (P4). Index pairs, else struct_conn, else distance.
- A child's inter-unit bonds: empty if the root's are empty and `dynamicBonds` and `bondProvider` are the root's; else
  restricted from the root's (`subsetInterUnitBonds`, `unit/bonds/subset.ts`) when the bond props and the provider are
  the root's and the units' operators are the root's after one common rigid motion (the identity included; same operator
  keys); else computed. A child moved that way (`Structure.transform`, its components) computes the root's first if the
  root has none: such children are made again at every move, the root stays. The restriction keeps what the root shows
  where a recompute could differ (struct_conn partners outside the child, which unit of a pair is iterated, water and
  ion pairing by a unit's first element) and for bonds the root carried over from another frame: a child shows a subset
  of its root. A provider's `inter` must allow both: its bonds of a subset are among those of the whole, and a common
  rigid motion of all units keeps them (no absolute coordinates: a box, a slab, `structure.boundary`). The shipped
  providers do.
- `IntraUnitBonds` store both directions (`edgeCount` undirected); `InterUnitBonds` list a pair under both units
  (`edgeCount` directed; `Structure.bondCount` halves it).

## 6. Elements and selections

- `StructureElement.Location`: `{ structure, unit, element: ElementIndex }`, mutable and reused by iterators; clone to
  keep one. `StructureProperties` (`properties.ts`) are Location → value accessors.
- `StructureElement.Loci`: `{ structure, elements: { unit, indices: OrderedSet<UnitIndex> }[] }`. Set operations and
  `remap` key on unit ids. `extendTo*` grow to residues, chains, entities, models, operators, a radius. `toExpression`
  writes MolScript by source index and operator name.
- `Bundle`: a loci without a structure (unit ids grouped by `invariantId`, ranges and a set), checked against
  `structure.hashCode`; `Schema`: rows or columns of file identifiers, to an expression or a substructure.
- `Bond.Location`, `Bond.Loci`: bonds by unit and `UnitIndex` on both sides.

## 7. Other derived data

- Carbohydrates (`carbohydrates/`): saccharide elements from rings, links from bonds; SNFG shapes and colors
  (`SaccharideCompIdMap` is global: `setSaccharideCompIdMapType`).
- Coordination (`coordination/`): sites with four or more covalent or metallic bonds.
- `util/`: polymer element sets, principal axes, subset builders, unit merging, superposition (RMSD, SIFTS, ligand MCCS,
  TM-align).
