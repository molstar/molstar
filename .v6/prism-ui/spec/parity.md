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
| Home: Download Structure, Open Files, Download Density, Load Trajectory | Home tab, action blocks (`HomeAction`, `triggerLabel`): the registry's actions (Mol*'s Download Structure among them), Open files | `load.action`, `load.entry`, `load.files` | built (`load.action` runs a registry action; `load.entry` for commands and the agent; `load.files` is Open files; the Download Density and Load Trajectory entries are registry actions) |
| custom import controls (Zenodo) | Home's `home.import` slot (Zenodo's extension) | `zenodo.*` | later (phase 5) |
| drag and drop (files, sessions) | the shell's drop target → open files; the drop overlay shows the formats hint | `load.files` | built (`load.files`; the overlay's hint lists the formats the registry reads) |
| `Viewer.load*` (structure, PDB, PDB-IHM, EMDB, AFDB, ModelArchive, volume, trajectory, MVS, files, URL, snapshot) | `@molstar/plugin-extension/loaders`, imported as needed; `viewer.run` for a history step; the loaders feature for the API | `load.*` (loaders feature) | built (`loaders`, ported over `@molstar/plugin`: `plugin-extension` pulls in React; `load.svs` later (with SeqViewSpec documents); `load.data`, raw data, later (phase 4, for the agent)) |
| URL parameters (`pdb`, `afdb`, `emdb`, `url`, `mvs-url`, `snapshot-url`, …) | run at start by the loaders feature                             | `load.*`, options            | built (`loadsFromUrl`, `runUrlLoads`; the app runs them with origin `url`) |
| remote states (list, upload, fetch, remove) | Home's list and the Remote states card, behind a setting; the listing and removal now live in the classic UI's view (`state/snapshots.tsx`, `RemoteStateSnapshots`) and move into the feature's API | `session.remote.*` | built (`remote-states`, lazy, behind `ui.showRemoteStates`; the card from the toolbar); the list on Home later (phase 6) |
| sessions `.molx` / `.molj` (open, save)                               | Snapshots tab; drop; Export › Session                                 | `session.*`                  | open and save built (`export-session`: `session.save`; the Snapshots tab's Save and Open) |
| loading raw data (a string or bytes with a format); the all-models and assembly URL option | the loaders feature (`load.data`); URL parameters | `load.data`, options | later (phase 4, for the agent) |

## 2. The scene

| v5                                                                    | Prism                                                                 | API area                     | State  |
| --------------------------------------------------------------------- | --------------------------------------------------------------------- | ---------------------------- | ------ |
| Structure source (hierarchy picker, model index, structure type, conformation) | the Structure tab's scene tree and the structure inspector | `structure.*` | built (`structure.model`, `structure.assembly`, `structure.preset`, `structure.symmetry`: mates, cells, off; the symmetry-assembly type through `state.update` only); the Conformation Transform (move and rotate by hand) is missing: later (phase 6) |
| crystal symmetry and the unit cell (Symmetry Mates, Symmetry indices, unit cell representation) | the structure inspector (symmetry); the unit cell as model-level representation rows | `structure.symmetry`, `structure.*` | symmetry built (`structure.symmetry`: mates, cells, off); the unit cell's model-level representation rows missing: later (phase 6) |
| Conformation Transform (move, rotate a structure by hand) | the structure inspector | `structure.*` | missing: later (phase 6) |
| (new) the structure's contents | "Structure" group, first under each structure in the scene tree: its entities, then their chains, built on first expand (mesoscale-safe); each row Select, Focus, Color…, Create component, and visibility | `contents.get`, `contents.visibility` (targets by entity and chain) | built (`contents.get`, `contents.visibility`) |
| Components (preset, add, options, groups, visibility, delete, per-representation editors) | scene tree rows, node menus, inspectors | `component.*`, `representation.*` | built (`component.add`, `.list`, `.modify`, `.options`, `.remove`, `.rename`, `.visibility`; `representation.*`) |
| Quick styles | Structure tab's styles; command menu Styles group; presets are `structure.preset` | `styles.apply`, `structure.preset` | built (`styles.apply`, `structure.preset`) |
| generic representation controls (membrane orientation, symmetry, kinemage) | structure rows and inspector sections (extensions) | the extension's | later (phase 5) |
| Volume (isovalue, representation), volume streaming | volume nodes and inspectors; the streaming structure row | `volume.*` | built (`get`, `isovalue`, `representation`, `visibility`, `remove`, `streaming`: on and off only; the view type, detail and channels go through `state.update`, later (phase 6)) |
| Particles                                                             | particle nodes and inspectors                                         | `particles.*`                | built |
| Measurements (distance, angle, dihedral, orientation, plane; options) | Selection › Measure; measurement nodes                                | `measure.*`                  | built (`measurements`; rows in the scene tree) |
| Superposition (chains, TM-align, ligands, atoms, SIFTS)               | Selection › Superpose                                                 | `superpose.*`                | built (`superposition`, lazy) |
| State tree (nodes, current, delete, visibility, decorators, applicable actions) | State tree tab and its inspector                            | `state.*`                    | built (`state-tree`, lazy) |
| Undo (components, state tree) | the viewer's history (one per step and per agent turn) | `history.*` | built (`history.undo`, `history.undo-turn`, `history.list`; the History tab) |

## 3. Selection and focus

| v5                                                                    | Prism                                                                 | API area                     | State  |
| --------------------------------------------------------------------- | --------------------------------------------------------------------- | ---------------------------- | ------ |
| selection mode, granularity | the canvas mode switch | `select.mode`, `select.granularity` | built (`select.mode`, `select.granularity`) |
| selection toolbar (union, subtract, intersect, set; by query or residue list; theme layers: color, transparency, emissive, material, clipping, resets; component and Modify by Selection; help) | the selection bar and Selection tab; `select` commands with targets ([targets.md](targets.md)) | `select.*`, `representation.*` | built (API: `select.*`, `representation.color`, `.transparency`, `.emissive`, `.material`, `.clipping`, `.reset`); the Selection tab has color, theme, transparency, glow and Make a component; material, clipping, layer reset and Modify by selection are commands only |
| focus (click, history; the focus picker's per-ligand and per-residue entries) | focus from the selection bar, the sequence, the scene tree; the picker's entries as `focus.set` targets in the command menu (`focus` group, ligands and residues as completions); the history as recent targets | `focus.*` | built (`focus.set`, `.get`, `.clear`, with target completions); the focus history is missing: later (phase 4) |
| procedural animation                                                  | Animation card › Procedural                                           | `animation.*`                | built (`animation.wiggle`) |

## 4. The canvas

| v5                                                                    | Prism                                                                 | API area                     | State  |
| --------------------------------------------------------------------- | --------------------------------------------------------------------- | ---------------------------- | ------ |
| camera reset, zoom, orient and reset axes | canvas toolbar Reset menu (Reset, Orient, Reset axes) | `camera.*` | built (`camera.reset`, `.zoom`, `.orient`, `.frame`, `.reset-axes`; the toolbar's Reset with Orient and Reset axes) |
| screenshot (preview, crop, copy, download; state import/export; camera info) | Export › Image | `export.image` | built (`export-image`; WebGL paths checked in the browser); state import/export and camera info: see below |
| viewport state import/export, camera info (v5's screenshot panel) | Export › Image | `export.*`, `camera.get` | later (phase 6) |
| canvas settings (background, lighting, camera, clipping)              | Canvas settings card                                                  | `canvas.*`                   | built (`canvas-settings`: `canvas.settings`, `.lighting`, `.camera`, `.clipping`) |
| illumination, XR | Canvas settings card and a command; XR: none yet | `canvas.illumination`, `canvas.xr` | illumination built (`canvas-settings`: `canvas.illumination`, in the Canvas card and as a command; a toolbar toggle later (phase 6)); `canvas.xr`: later (phase 6) |
| illumination toggle on the toolbar | toolbar Rendering menu | `canvas.illumination` | later (phase 6) (the Canvas card and a command have it now) |
| expand, fullscreen | toolbar item `expand` (fullscreen) | none (the toolbar item) | built (fullscreen only: the toolbar item `expand`); no expanded-viewport mode and no `view.expand` or `view.fullscreen`: later (phase 6) |
| toggle controls | the rail (always there; compact rail); `ui.collapsed` starts it collapsed | `ui.*` options | built (no `navigation.*` operations: the rail is always there; `ui.collapsed` starts collapsed); a toggle later (phase 6) |
| animation (snapshots, trajectory, spin, rock), snapshot stepper, description | Animation card; playback card top left; description markdown   | `animation.*`, `snapshot.*`  | built (`animation`; the snapshot card is an overlay for now) |
| trajectory controls                                                   | trajectory bar                                                        | `trajectory.*`               | built (the trajectory bar over `trajectory.*`) |
| hover labels                                                          | hover label at the cursor                                             | none (state only)            | built  |
| toasts, task progress (abortable), log | status column, busy overlay, History tab (journal.md §3) | none (state only) | built (status column, busy overlay, History tab); the plugin's toasts (`plugin.managers.toast`) are not wired: later (phase 5) |
| drag-and-drop overlay | the shell's drop target | `load.files` | built (`load.files`) |

## 5. Panels and settings

| v5                                                                    | Prism                                                                 | API area                     | State  |
| --------------------------------------------------------------------- | --------------------------------------------------------------------- | ---------------------------- | ------ |
| sequence viewer (structure, mode, entity, chain selectors)           | the sequence band, Text and Tracks ([sequence.md](sequence.md)); no pickers | `sequence.*`           | built (first slice) |
| Snapshots (add, edit, move, replace, remove, options; save/open)      | Snapshots tab                                                         | `snapshot.*`                 | built (`snapshots`) |
| Settings (Canvas3D params, behavior tree)                             | Settings tab (groups; the behaviors and the plugin's settings)        | `settings.*`                 | built (Performance and the picking level are not undone) |
| behaviors tree (plugin behaviors and their params) | folded into the Settings tab's Behaviors group | `settings.*` | built (`settings.get`, `settings.set`) |
| Help (prose, bindings) | Help tab (bindings); the command menu's `?` | none | built (the bindings and the `?` command); the help prose later (phase 6) |
| layout (regions, controls display, expanded) | the shell (docks, compact, pinned) | `options.ui.*` (`buttons`, `showSequence`, `showHistory`, `collapsed`), `view.*` | built (docks, compact rail, `options.ui.buttons`, `showSequence`, `showHistory`, `collapsed`) |
| keyboard: canvas and behavior keys (trackball W/A/S/D/R/F/Q/E, focus and camera keys, snapshot Ctrl+arrows), Enter/Escape in forms | the same canvas keys while the pointer is over the canvas; command keybindings for the rest, scoped to the viewer (README decision 29, design/commands.md §6) | — | built (canvas keys through the plugin's behaviors, Ctrl+arrows for snapshots; command shortcuts ⌘K, ⌘Z, `?`) |

## 6. Extensions

| Extension                                                             | Prism                                                                 | State  |
| --------------------------------------------------------------------- | --------------------------------------------------------------------- | ------ |
| model-export, geo-export                                              | Export › Models, Geometry (sections); commands and tools               | mocked (sections) |
| video-export                                                          | Export › Animation; lazy                                              | mocked (section) |
| zenodo | Home's import slot; lazy | later (phase 5) |
| assembly-symmetry                                                     | structure row, inspector section                                      | mocked (row) |
| kinemage                                                              | inspector section                                                     | missing |
| model-archive quality assessment (pLDDT, QMEAN), PAE plot             | designed anew as an extension with commands first (README decision 27): color by pLDDT or QMEAN, show PAE; sequence behaviors for the pLDDT track and the PAE dot plot; a new PAE view later. The viewer core no longer depends on it: `@molstar/model-archive-extension/prism-viewer` brings the themes, `annotations.apply`'s `plddt` and `qmean`, the commands and the pLDDT track | built (pLDDT, QMEAN; PAE later) |
| volumes-and-segmentations                                             | dropped for now: to be redesigned                                     | dropped |
| anvil (membrane orientation) | an extension with commands (show, hide and fit the membrane) and a structure row for its options (v5's generic representation control) | later (phase 5) |
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
