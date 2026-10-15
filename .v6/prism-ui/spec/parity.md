# Parity with the v5 Viewer

Every v5 Viewer feature and where it lives in the Prism viewer: its place in the UI, its API area (README §3.3; each
area's operations have commands, tools, or both), and its state today. The inventory is
[../analysis/current-ui.md](../analysis/current-ui.md) §5, [../analysis/extensions.md](../analysis/extensions.md) and
the Viewer's API and URL parameters; the state of the plugin behind each row is [state.md](state.md).

Status: **accepted, first pass** (2026-10-15). State: **mocked** (a view over a mock view model in Storybook), **built** (over a real
plugin), **missing** (not mocked), **dropped** (decided not to carry over), **later**.

## 1. Loading and data

| v5                                                                    | Prism                                                                 | API area                     | State  |
| --------------------------------------------------------------------- | --------------------------------------------------------------------- | ---------------------------- | ------ |
| Home: Download Structure, Open Files, Download Density, Load Trajectory | Home tab, action blocks (`HomeAction`, `triggerLabel`)                | `load.*`                     | mocked |
| custom import controls (Zenodo)                                       | Home's `home.import` slot (Zenodo's extension)                        | `zenodo.*`                   | missing |
| drag and drop (files, sessions)                                       | the shell's drop target → open files                                  | `load.files`                 | mocked |
| `Viewer.load*` (structure, PDB, PDB-IHM, EMDB, AFDB, ModelArchive, volume, trajectory, MVS, files, URL, snapshot) | `@molstar/plugin-extension/loaders`, imported as needed; `viewer.run` for a history step; the loaders feature for the API | `load.*` (loaders feature) | missing |
| URL parameters (`pdb`, `afdb`, `emdb`, `url`, `mvs-url`, `snapshot-url`, …) | run at start by the loaders feature                             | `load.*`, options            | missing |
| remote states (list, upload, fetch, remove)                           | Home's list and the Remote states card, behind a setting; the listing and removal now live in the classic UI's view (`state/snapshots.tsx`, `RemoteStateSnapshots`) and move into the feature's API | `session.remote.*`           | mocked |
| sessions `.molx` / `.molj` (open, save)                               | Snapshots tab; drop; Export › Session                                 | `session.*`                  | open and save built (`export-session`: `session.save`; the Snapshots tab's Save and Open) |

## 2. The scene

| v5                                                                    | Prism                                                                 | API area                     | State  |
| --------------------------------------------------------------------- | --------------------------------------------------------------------- | ---------------------------- | ------ |
| Structure source (hierarchy picker, model index, structure type, conformation) | the Structure tab's scene tree and the structure inspector   | `structure.*`                | mocked |
| (new) the structure's contents                                        | "Structure" group, first under each structure in the scene tree: its entities, then their chains, built on first expand (mesoscale-safe); each row Select, Focus, Color…, Create component, and visibility | `structure.contents.*` (targets by entity and chain) | mocked |
| Components (preset, add, options, groups, visibility, delete, per-representation editors) | scene tree rows, node menus, inspectors          | `component.*`, `representation.*` | mocked |
| Quick styles                                                          | Structure tab's styles; command menu Styles group                     | `styles.*`, `preset.*`       | mocked |
| generic representation controls (membrane orientation, symmetry, kinemage) | structure rows and inspector sections (extensions)              | the extension's              | mocked (rows) |
| Volume (isovalue, representation), volume streaming                   | volume nodes and inspectors; the streaming structure row              | `volume.*`                   | built (`get`, `isovalue`, `representation`, `visibility`, `remove`, `streaming`) |
| Particles                                                             | particle nodes and inspectors                                         | `particles.*`                | built |
| Measurements (distance, angle, dihedral, orientation, plane; options) | Selection › Measure; measurement nodes                                | `measure.*`                  | mocked |
| Superposition (chains, TM-align, ligands, atoms, SIFTS)               | Selection › Superpose                                                 | `superpose.*`                | mocked |
| State tree (nodes, current, delete, visibility, decorators, applicable actions) | State tree tab and its inspector                            | `state.*`                    | mocked |
| Undo (components, state tree)                                         | the viewer's history (one per step and per agent turn)                | `history.*`                  | missing |

## 3. Selection and focus

| v5                                                                    | Prism                                                                 | API area                     | State  |
| --------------------------------------------------------------------- | --------------------------------------------------------------------- | ---------------------------- | ------ |
| selection mode, granularity                                           | the canvas mode switch                                                | `selection.mode`             | mocked |
| selection toolbar (union, subtract, intersect, set; by query or residue list; theme layers: color, transparency, emissive, material, clipping, resets; component and Modify by Selection; help) | the selection bar and Selection tab; `select` commands with targets ([targets.md](targets.md)) | `select.*`, `selection.*` | mocked |
| focus (click, history; the focus picker's per-ligand and per-residue entries) | focus from the selection bar, the sequence, the scene tree; the picker's entries as `focus.set` targets in the command menu (`focus` group, ligands and residues as completions); the history as recent targets | `focus.*` | mocked (partly) |
| procedural animation                                                  | Animation card › Procedural                                           | `animation.*`                | mocked |

## 4. The canvas

| v5                                                                    | Prism                                                                 | API area                     | State  |
| --------------------------------------------------------------------- | --------------------------------------------------------------------- | ---------------------------- | ------ |
| camera reset, zoom, orient and reset axes                             | canvas toolbar Reset menu                                             | `camera.*`                   | mocked |
| screenshot (preview, crop, copy, download; state import/export; camera info) | Export › Image                                                 | `export.image`               | built (`export-image`; WebGL paths checked in the browser) |
| canvas settings (background, lighting, camera, clipping)              | Canvas settings card                                                  | `canvas.*`                   | built (`canvas-settings`: `canvas.settings`, `.lighting`, `.camera`, `.clipping`) |
| illumination, XR                                                      | toolbar Rendering menu                                                | `canvas.illumination`, `canvas.xr` | illumination built (`canvas-settings`); XR mocked |
| expand, fullscreen                                                    | toolbar Expand / Fullscreen                                           | `view.expand`, `view.fullscreen` | mocked |
| toggle controls                                                       | the rail (compact rail)                                               | `navigation.*`               | mocked |
| animation (snapshots, trajectory, spin, rock), snapshot stepper, description | Animation card; playback card top left; description markdown   | `animation.*`, `snapshot.*`  | mocked |
| trajectory controls                                                   | trajectory bar                                                        | `trajectory.*`               | mocked |
| hover labels                                                          | hover label at the cursor                                             | none (state only)            | built  |
| toasts, task progress (abortable), log                                | status column, busy overlay, History tab (journal.md §3)              | none (state only)            | mocked |
| drag-and-drop overlay                                                 | the shell's drop target                                               | `load.files`                 | mocked |

## 5. Panels and settings

| v5                                                                    | Prism                                                                 | API area                     | State  |
| --------------------------------------------------------------------- | --------------------------------------------------------------------- | ---------------------------- | ------ |
| sequence viewer (structure, mode, entity, chain selectors)           | the sequence band, Text and Tracks ([sequence.md](sequence.md)); no pickers | `sequence.*`           | built (first slice) |
| Snapshots (add, edit, move, replace, remove, options; save/open)      | Snapshots tab                                                         | `snapshot.*`                 | built (`snapshots`) |
| Settings (Canvas3D params, behavior tree)                             | Settings tab (groups; the behaviors and the plugin's settings)        | `settings.*`                 | mocked |
| Help (prose, bindings)                                                | Help tab (bindings); the command menu's `?`                           | none                         | mocked |
| layout (regions, controls display, expanded)                          | the shell (docks, compact, pinned)                                    | `navigation.*`, `viewState`  | mocked |
| keyboard: canvas and behavior keys (trackball W/A/S/D/R/F/Q/E, focus and camera keys, snapshot Ctrl+arrows), Enter/Escape in forms | the same canvas keys while the pointer is over the canvas; command keybindings for the rest, scoped to the viewer (README decision 29, design/commands.md §6) | — | mocked (command shortcuts) |

## 6. Extensions

| Extension                                                             | Prism                                                                 | State  |
| --------------------------------------------------------------------- | --------------------------------------------------------------------- | ------ |
| model-export, geo-export                                              | Export › Models, Geometry (sections); commands and tools               | mocked (sections) |
| video-export                                                          | Export › Animation; lazy                                              | mocked (section) |
| zenodo                                                                | Home's import slot; lazy                                              | missing |
| assembly-symmetry                                                     | structure row, inspector section                                      | mocked (row) |
| kinemage                                                              | inspector section                                                     | missing |
| model-archive quality assessment (pLDDT, QMEAN), PAE plot             | designed anew as an extension with commands first (README decision 27): color by pLDDT or QMEAN, show PAE; sequence behaviors for the pLDDT track and the PAE dot plot; a new PAE view later. The viewer core no longer depends on it (`@molstar/prism-viewer`'s dependency and the first slice's pLDDT coloring move into the extension) | built (pLDDT track, to move) |
| volumes-and-segmentations                                             | dropped for now: to be redesigned                                     | dropped |
| anvil (membrane orientation)                                          | an extension with commands (show, hide and fit the membrane) and a structure row for its options (v5's generic representation control) | missing |
| dnatco, rcsb, pdbe, wwpdb, sb-ncbr, g3d, backgrounds, debug-helpers  | data-only (`pluginExtension`): presets, themes, formats through the generic commands and tools; their annotations through one `annotations.apply` command (below) | missing |
| MolViewSpec                                                           | the viewer's base (load, snapshots)                                   | built (first slice) |

## 6a. Not copied verbatim

- **"Automatic (w/ Annotation)"** (v5's default preset: pLDDT or QMEAN coloring for predicted models, SB-NCBR charges)
  is not ported as a preset (README decision 27): loading a predicted model applies the default preset, and an
  `annotations.apply` command (and tool) applies what the loaded structure carries (quality scores, validation,
  charges), offered in Home's result line after a load ("This model has pLDDT scores: color by them?"). The
  model-archive and sb-ncbr extensions register what they can apply.

## 7. New in Prism (no v5 counterpart)

The "Structure" contents group in the scene tree (§2), the command menu, the agent and its Inspector, the history, the sequence's Text view, the canvas ink, compact mode on
phones, runtime themes (light, dark, blue, auto; colorblind-safe and high-contrast to add).

## 8. Hiding part of a structure

Two ways, both API operations over targets (README §2.18): at the entity and chain level (a contents row's eye: what
is hidden is left out of the structure's representations) and as transparency overlays (a layer over any target, from
opaque to hidden, keeping the geometry). Select, focus, color and create component map onto existing operations through
targets ([targets.md](targets.md)).
