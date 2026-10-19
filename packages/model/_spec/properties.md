# Custom properties

`src/props/`: data attached to models, structures and volumes on demand, and the computed properties.

## 1. The framework

- `CustomPropertyDescriptor` (`model/custom-property.ts`): a `name`, optional `cifExport` and MolScript `symbols`.
  `CustomProperties` (on models, structures, volumes): the attached descriptors, reference counts, and assets per
  descriptor and owner (disposed when replaced).
- `CustomProperty.Provider` (`props/common/custom-property.ts`): `descriptor`, `getParams`, `defaultParams`,
  `isApplicable`, `attach(ctx, data, props?, addRef?)`, `ref`, `get(data)` (a `ValueBox`), `set`, `props`. `Registry`
  lists providers per data kind; the plugin owns three (models, structures, volumes) and builds the auto-attach
  parameters from the applicable ones.
- `attach` computes only when there is no value or the merged props differ from the stored ones; it never notices a
  change of the data: a new `Model` or `Structure` is new data. `set` keeps the box and bumps its version, which
  graphics reads as a change hash.

| Factory                   | `type`    | Stored in                                                         |
| ------------------------- | --------- | ----------------------------------------------------------------- |
| `CustomModelProperty`     | `static`  | `model._staticPropertyData`: shared by frames (model.md §3)       |
|                           | `dynamic` | `model._dynamicPropertyData`: per frame; assets owned by model id |
| `CustomStructureProperty` | `root`    | the root's `currentPropertyData`; computed once for all children  |
|                           | `local`   | the structure's own `currentPropertyData`                         |
| `CustomVolumeProperty`    | –         | `volume._propertyData`                                            |

- Asset owners: frames share `customProperties` but not dynamic data, so dynamic assets are owned per model id and one
  frame's attach does not dispose another's.
- Root-scoped values index through the root: readers map a child's location to the root's unit (ASA by
  `root.unitIndexMap`; secondary structure by `invariantId`).
- `createSimple(name, type)`: a hidden value holder (`Model.Index`, `Structure.Index`, ...).
- Format properties (`FormatPropertyProvider`, formats.md §2) are a separate, read-only mechanism over the same bags.

## 2. Computed properties

| Provider                        | Scope            | What, cached where                                                                                                 |
| ------------------------------- | ---------------- | ------------------------------------------------------------------------------------------------------------------ |
| `AccessibleSurfaceAreaProvider` | structure, root  | Shrake-Rupley per atom, summed per residue; neighbours prefiltered with `structure.lookup3d`                       |
| `SecondaryStructureProvider`    | structure, root  | per symmetry group (`invariantId`): the model's (archive or present), else DSSP, else Zhang-Skolnick               |
| `ValenceModelProvider`          | structure, local | charges, implicit H, geometry per unit; computed on the root's units for the whole topology; `unit.transientCache` |
| `InteractionsProvider`          | structure, local | features per symmetry group (`unit.transientCache`), intra- and inter-unit contacts (§2.1), refinement             |
| `HelixOrientationProvider`      | model, dynamic   | helix centres per residue (GROMACS helixorient); used by cartoons                                                  |
| `CrossLinkRestraintProvider`    | structure, local | pairs from `ihm_cross_link_restraint`                                                                              |
| `SIFTSMapping.Provider`         | model, static    | UniProt mapping per residue from `atom_site.pdbx_sifts_xref_db_*`                                                  |
| `StreamlinesProvider`           | volume           | streamlines through a vector field                                                                                 |

`interactions/params.ts` has no computation, so state and UI code can describe the parameters without loading the
providers.

### 2.1 Interactions and symmetry

- Features and intra-unit contacts are computed once per symmetry group, from its first unit, and shared.
- Inter-unit contacts: `symmetry: 'all'` (the default) searches every close unit pair. `'per-pair-class'` searches once
  per class (group of A, group of B, `opB⁻¹·opA` within 1e-7) and replays the `(i, j, type)` candidates on the others,
  checking line of sight per pair (it sees third units). Edges keep the order of 'all'.
- A pair is replayed only if both units have a rigid operator (testers measure world distances and angles), no covalent
  inter-unit bond (`validPair` and the testers' bond walks read them), no operator-scoped `IndexPairBonds`, and the same
  `idealGeometry` as their group's first such unit (the valence model is per unit, from the root's topology). Then every
  input of the search is a function of the class, and the result equals 'all' except a contact on a distance or angle
  threshold within the operators' precision.
- The default stays 'all' until the maintainers decide; `computeInteractions(…, { stats })` counts pairs, classes and
  replays.
