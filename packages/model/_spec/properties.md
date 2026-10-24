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
  change of the data: a new `Model` or `Structure` is new data. A structure property's attach with the props of one in
  flight awaits it (callers outside the state queue overlap with builds) and computes itself if that one fails; while it
  waits, its own context's abort and progress do not apply. `set` keeps the box and bumps its version, which graphics
  reads as a change hash.

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
| `InteractionsProvider`          | structure, local | features per symmetry group (`unit.transientCache`), intra- and inter-unit contacts, refinement                    |
| `HelixOrientationProvider`      | model, dynamic   | helix centres per residue (GROMACS helixorient); used by cartoons                                                  |
| `CrossLinkRestraintProvider`    | structure, local | pairs from `ihm_cross_link_restraint`                                                                              |
| `SIFTSMapping.Provider`         | model, static    | UniProt mapping per residue from `atom_site.pdbx_sifts_xref_db_*`                                                  |
| `StreamlinesProvider`           | volume           | streamlines through a vector field                                                                                 |

`interactions/params.ts` has no computation, so state and UI code can describe the parameters without loading the
providers.
