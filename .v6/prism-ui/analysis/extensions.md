# Extensions and apps: UI inventory and prism-ui ports

Every UI piece contributed outside `@molstar/plugin-ui`, how it plugs in today, and its proposed prism-ui form. The
prism-ui contribution kinds (command, panel section, panel, dialog, viewport overlay) are defined in
[../design/shell.md](../design/shell.md) and [../design/commands.md](../design/commands.md).

## 1. Findings

- Only 8 extensions ship React (12 `.tsx` files). All panels extend `CollapsableControls` and build on
  `ParameterControls`, `ApplyActionControl`, and `UpdateTransformControl`.
- Their parameters are already `PD.Params` (the backend-generated animation params, `GeometryParams`, model-export
  `Params`, `AssemblySymmetryProvider.getParams()`, `ZenodoImportParams`), so a generic PD form plus command arguments
  covers most of them.
- Most extensions register only data providers (presets, themes, selection queries, loci labels, actions, formats).
  These need no ported UI; they surface through generic UI and, in prism-ui, through commands generated from the
  registry (§3 of [../design/commands.md](../design/commands.md)).
- Several extensions hold their logic in UI-free controllers already (`GeometryControls`, `VideoControls`, volseg entry
  data, `MAPairwiseScorePlot`'s interactivity subject). Ports wrap these in view models.
- The PAE plot (`MAPairwiseScorePlot`) is used standalone by `examples/alphafolddb-pae`, so it must stay embeddable.
- `kinemage` probes a non-existent `ctx.customControls`; `zenodo`'s behavior display name says "Export".

## 2. Extensions with UI

| Extension                 | Today (LOC)                                                                                          | Integration                                                 | prism-ui form                                                                                                                                                                                                                                     |
| ------------------------- | ---------------------------------------------------------------------------------------------------- | ----------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| geo-export                | `GeometryExporterUI` (116): format select, Save, View in AR                                          | `customStructureControls`                                   | A section of the Export panel (`export` slot), plus commands (opening the panel at the section, or exporting directly): `geo-export.export` (format arg: glb/stl/obj/usdz), `geo-export.view-ar` (when supported). Blocked when nothing is drawn. |
| model-export              | `ModelExportUI` (75): format select, Export (ZIP for several)                                        | `customStructureControls`                                   | Commands: `model-export.export` (format arg: cif/bcif).                                                                                                                                                                                           |
| video-export              | `VideoEncoderUI` (138): animation, its params, quantization, Render, Save                            | `customStructureControls`; embedded by mesoscale-explorer   | A section of the Export panel (animation params, render progress, save), plus command `video-export.export` opening the panel at it; backed by `VideoExportModel` over `VideoControls`.                                                           |
| zenodo                    | `ZenodoImportUI` (345): record ID → file list → typed import                                         | `customImportControls` (its only user)                      | Command `zenodo.import` with a custom UI (record ID → file list → typed import params), inline in the command menu or as a dialog.                                                                                                                |
| assembly-symmetry         | `AssemblySymmetryControls` (196): Enable, symmetry params, C1 notice                                 | `customStructureControls` + `genericRepresentationControls` | Section in the Structure tab for structure nodes (shown when applicable) + commands: enable, choose symmetry, color by cluster.                                                                                                                   |
| kinemage                  | `KinemageControls` (453): views, animate, group/subgroup/master toggles                              | `customStructureControls`                                   | Section in the Structure tab (toggle tree, shown when a kinemage exists) + commands: apply view N, next animate frame.                                                                                                                            |
| volumes-and-segmentations | `VolsegUI` + sub-controls (380): entry, global options, fitted models, segments, annotations, volume | `customStructureControls`; `LoadVolseg` action              | Tab "Segmentation" (contributed tab) + commands: load entry, show/hide all segments, set opacity, select segment.                                                                                                                                 |
| model-archive (QA)        | `MAPairwiseScorePlotPanel` (719): PAE heatmap, hover, box-select                                     | `customStructureControls` (config-gated)                    | Tab "PAE plot" (right, movable to the bottom dock) + commands: toggle PAE plot, color by pLDDT/QMEAN. Plot core stays embeddable (custom element).                                                                                                |

## 3. Extensions without UI

anvil, dnatco, rcsb, pdbe, wwpdb, sb-ncbr, g3d, backgrounds, debug-helpers, interactions, json-cif, meshes, volume-tools
(controllers only), volseg-api, alpha-orbitals. They contribute through `ctx.register(...)` and behavior params.
prism-ui gains commands for free from the registry: representation and hierarchy presets ("Apply preset…"), color themes
("Color by…"), selection queries ("Select…"), state actions ("Download tunnels…", "Load G3D…"), background styles ("Set
background…"). Behavior params (debug-helpers) stay editable through the Settings view's behavior tree.

`genericRepresentationControls` (anvil "Membrane Orientation", assembly-symmetry "Symmetries", kinemage) becomes a typed
"component entries" contribution on the UI spec (see [../design/shell.md](../design/shell.md) §5).

`leftPanelTabName.next('data')` (sb-ncbr, volseg, g3d) becomes `shell.reveal('state-tree')` or, better, nothing: the
command menu and toasts make the jump unnecessary.

## 4. Apps

| App                | Today                                                                                                     | prism-ui                                                                                                                                                                                 |
| ------------------ | --------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| viewer             | `DefaultPluginUISpec()` + layout flags, region state, `controls.*: 'none'`, ~20 `PluginConfig.Viewport.*` | Stays on the classic UI, unchanged. Prism is a new app (`apps/prism`), not a mode of the viewer; options it shares map to the shell's initial layout and viewport-toolbar contributions. |
| docking-viewer     | Replaces `viewport.view` to add 4 preset buttons                                                          | A viewport-overlay contribution (button group) + 4 commands; no viewport replacement.                                                                                                    |
| mesoscale-explorer | Own left/right panels (~3k LOC): entity tree, selection info, quick styles, animation, snapshots          | Its own views contributed to the shell; loaders, sessions, styles and animation presets become commands. The best test of shell modularity; port last.                                   |
| mvs-stories        | `<mvs-stories-viewer>` + `<mvs-stories-snapshot-markdown>` custom elements                                | Embeddable `msp-viewport` + `msp-snapshot-stepper` elements over shared models; commands: snapshot next/previous, load story.                                                            |

## 5. Examples

`volume-tools` (bodies and mask panels, viewport drawing canvas; ~1.7k), `ligand-editor` (edit panel; ~500),
`interactions`, `alphafolddb-pae`, `alpha-orbitals`, `react`, `basic-wrapper`, `proteopedia-wrapper` build their own
controls around the plugin. Ports are optional and serve as API checks: `react` (React adapter over view models),
`alpha-orbitals` (PD form element outside the shell), `ligand-editor` (commands + a contributed panel).
