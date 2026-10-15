# Parity with the v5 Viewer

Every v5 Viewer feature and where it lives in the Prism viewer: its place in the UI, its API area (README §3.3; each
area's operations have commands, tools, or both), and its state today. The inventory is
[../analysis/current-ui.md](../analysis/current-ui.md) §5, [../analysis/extensions.md](../analysis/extensions.md) and
the Viewer's API and URL parameters; the state of the plugin behind each row is [state.md](state.md).

Status: **draft** (2026-10-15). State: **mocked** (a view over a mock view model in Storybook), **built** (over a real
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
| sessions `.molx` / `.molj` (open, save)                               | Snapshots tab; drop; Export › Session                                 | `session.*`                  | mocked |

## 2. The scene

| v5                                                                    | Prism                                                                 | API area                     | State  |
| --------------------------------------------------------------------- | --------------------------------------------------------------------- | ---------------------------- | ------ |
| Structure source (hierarchy picker, model index, structure type, conformation) | the Structure tab's scene tree and the structure inspector   | `structure.*`                | mocked |
| (new) the structure's contents                                        | "Structure" group, first under each structure in the scene tree: its entities, then their chains, built on first expand (mesoscale-safe); each row Select, Focus, Color…, Create component, and visibility | `structure.contents.*` (targets by entity and chain) | mocked |
| Components (preset, add, options, groups, visibility, delete, per-representation editors) | scene tree rows, node menus, inspectors          | `component.*`, `representation.*` | mocked |
| Quick styles                                                          | Structure tab's styles; command menu Styles group                     | `styles.*`, `preset.*`       | mocked |
| generic representation controls (membrane orientation, symmetry, kinemage) | structure rows and inspector sections (extensions)              | the extension's              | mocked (rows) |
| Volume (isovalue, representation), volume streaming                   | volume nodes and inspectors; the streaming structure row              | `volume.*`                   | mocked |
| Particles                                                             | particle nodes (to add) and inspectors                                | `particles.*`                | missing |
| Measurements (distance, angle, dihedral, orientation, plane; options) | Selection › Measure; measurement nodes                                | `measure.*`                  | mocked |
| Superposition (chains, TM-align, ligands, atoms, SIFTS)               | Selection › Superpose                                                 | `superpose.*`                | mocked |
| State tree (nodes, current, delete, visibility, decorators, applicable actions) | State tree tab and its inspector                            | `state.*`                    | mocked |
| Undo (components, state tree)                                         | the viewer's history (one per step and per agent turn)                | `history.*`                  | missing |

## 3. Selection and focus

| v5                                                                    | Prism                                                                 | API area                     | State  |
| --------------------------------------------------------------------- | --------------------------------------------------------------------- | ---------------------------- | ------ |
| selection mode, granularity                                           | the canvas mode switch                                                | `selection.mode`             | mocked |
| selection toolbar (union, subtract, intersect, set; by query or residue list; theme; component; undo; help) | the selection bar and Selection tab; `select` commands with targets ([targets.md](targets.md)) | `select.*`, `selection.*` | mocked |
| focus (click, history)                                                | focus from the selection bar, the sequence, the scene tree            | `focus.*`                    | mocked |
| procedural animation                                                  | Animation card › Procedural                                           | `animation.*`                | mocked |

## 4. The canvas

| v5                                                                    | Prism                                                                 | API area                     | State  |
| --------------------------------------------------------------------- | --------------------------------------------------------------------- | ---------------------------- | ------ |
| camera reset, zoom, orient and reset axes                             | canvas toolbar Reset menu                                             | `camera.*`                   | mocked |
| screenshot (preview, crop, copy, download; state import/export; camera info) | Export › Image                                                 | `export.image`               | mocked |
| canvas settings (background, lighting, camera, clipping)              | Canvas settings card                                                  | `canvas.*`                   | mocked (background built) |
| illumination, XR                                                      | toolbar Rendering menu                                                | `canvas.illumination`, `canvas.xr` | mocked |
| expand, fullscreen                                                    | toolbar Expand / Fullscreen                                           | `view.expand`, `view.fullscreen` | mocked |
| toggle controls                                                       | the rail (compact rail)                                               | `navigation.*`               | mocked |
| animation (snapshots, trajectory, spin, rock), snapshot stepper, description | Animation card; playback card top left; description markdown   | `animation.*`, `snapshot.*`  | mocked |
| trajectory controls                                                   | trajectory bar                                                        | `trajectory.*`               | mocked |
| hover labels                                                          | hover label at the cursor                                             | none (state only)            | built  |
| toasts, task progress (abortable), log                                | status column, busy overlay, Log tab                                  | none (state only)            | mocked |
| drag-and-drop overlay                                                 | the shell's drop target                                               | `load.files`                 | mocked |

## 5. Panels and settings

| v5                                                                    | Prism                                                                 | API area                     | State  |
| --------------------------------------------------------------------- | --------------------------------------------------------------------- | ---------------------------- | ------ |
| sequence viewer (structure, mode, entity, chain selectors)           | the sequence band, Text and Tracks ([sequence.md](sequence.md)); no pickers | `sequence.*`           | built (first slice) |
| Snapshots (add, edit, move, replace, remove, options; save/open)      | Snapshots tab                                                         | `snapshot.*`                 | mocked |
| Settings (Canvas3D params, behavior tree)                             | Settings tab (groups; the behaviors and the plugin's settings)        | `settings.*`                 | mocked |
| Help (prose, bindings)                                                | Help tab (bindings); the command menu's `?`                           | none                         | mocked |
| layout (regions, controls display, expanded)                          | the shell (docks, compact, pinned)                                    | `navigation.*`, `viewState`  | mocked |
| keyboard (none beyond Enter/Escape)                                   | keybindings over commands                                             | —                            | mocked (command shortcuts) |

## 6. Extensions

| Extension                                                             | Prism                                                                 | State  |
| --------------------------------------------------------------------- | --------------------------------------------------------------------- | ------ |
| model-export, geo-export                                              | Export › Models, Geometry (sections); commands and tools               | mocked (sections) |
| video-export                                                          | Export › Animation; lazy                                              | mocked (section) |
| zenodo                                                                | Home's import slot; lazy                                              | missing |
| assembly-symmetry                                                     | structure row, inspector section                                      | mocked (row) |
| kinemage                                                              | inspector section                                                     | missing |
| model-archive quality assessment (pLDDT, QMEAN), PAE plot             | presets and themes; pLDDT and PAE in the sequence (Seq\*'s dot plot); a refreshed PAE view later | built (pLDDT track) |
| volumes-and-segmentations                                             | dropped for now: to be redesigned                                     | dropped |
| anvil, dnatco, rcsb, pdbe, wwpdb, sb-ncbr, g3d, backgrounds, debug-helpers | data-only (`pluginExtension`): presets, themes, formats through the generic commands and tools | missing |
| MolViewSpec                                                           | the viewer's base (load, snapshots)                                   | built (first slice) |

## 7. New in Prism (no v5 counterpart)

The "Structure" contents group in the scene tree (§2), the command menu, the agent and its Inspector, the history, the sequence's Text view, the canvas ink, compact mode on
phones, runtime themes (light, dark, blue, auto; colorblind-safe and high-contrast to add).

## 8. Hiding part of a structure

Two ways, both API operations over targets (README §2.18): at the entity and chain level (a contents row's eye: what
is hidden is left out of the structure's representations) and as transparency overlays (a layer over any target, from
opaque to hidden, keeping the geometry). Select, focus, color and create component map onto existing operations through
targets ([targets.md](targets.md)).
