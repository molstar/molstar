# Branch review, 2026-10-16

A review of `v6-prism-impl-proto` against `master` (base `f67e419d9`), after the smaller core, and the fixes that
followed. For David: **what was fixed** is below for reference; **what needs your call** is the
[decisions](#decisions-for-david) section, each with a recommendation; nothing in it is implemented.

## How

- **Scope**: the branch's diff, 76k lines of source over 609 files (tests and generated files apart), mostly
  `packages/prism/viewer` (41k), the extensions (5.8k), prism-ui (4.7k) and the assistant (1.7k).
- **Review** (one workflow, 36 agents): 15 slices of 3–4k changed lines, one reviewer each (Opus) on four lenses at
  once: bugs, performance (pointer-rate and per-frame work, allocations, leaks), API surface and simplicity of use,
  package size; three cross-cutting reviewers (bundle size with esbuild metafiles, Sonnet; the public API; the hot paths).
  Each reviewer's findings were checked against the code by one verifier (Opus), who set the severity, the fix's scope
  and an owner.
- **Result**: 191 findings, 188 confirmed (9 high, 70 medium, 109 low): 90 bugs, 51 API, 30 performance, 17 size.
- **Fixes** (one workflow, 16 agents, three lanes in parallel, an Opus pass then a Sonnet pass per area, each fix with a
  test, committed per area): 116 fixed, 9 skipped (5 already fixed by another lane, 4 deliberate). The shared core's
  findings (the runner, `viewer.ts`, the spec, the shell) were fixed by the lead: 23. Every check passes: lint, format,
  styles, typecheck, the workspace checks and 9,239 unit tests (128 more than before).

## The smaller core (before the review)

The core is now what everything stands on; the rest are built-in features in `DefaultPrismViewerSpec` (viewer-spec.md
§3.1): the scene's API (`structure`, `volume`, `trajectory`, `annotations`), the panels (`home`, `structure-panel`,
`selection-panel`, `help`), and the query languages from the spec (§2.7). The minimal example
(`examples/minimal-prism-viewer`, `structure` only) went from 3,253.8 kB to 2,892.9 kB as one file (899.9 kB to
803.9 kB gzipped); the Prism layer above the engine from about 0.29 MB to about 0.20 MB gzipped.

## Fixed (highlights)

The full list is the branch's commits after `5dc085f7e`; the most consequential:

- **`api.group` is all or nothing for every call in it**: a call refused before it ran (invalid values, blocked,
  unavailable), an unknown operation, and a read or view call that failed now fail the group and undo it (they never
  reached the history before, so the group committed).
- **The scene's tree**: a representation row's Remove removes that one (it removed every one of its type, with the type
  it had when the row was made); inspections follow what they show (facts, contents, a volume's levels and transfer
  functions went stale).
- **Coarse models**: targets, hover and selection work on integrative (PDB-IHM) models.
- **The sequence band** draws a loaded SeqViewSpec, follows a replaced document, and keeps its band when the followed
  structure changes; marks and hover of another structure no longer land on its residues.
- **Marks** survive a removed representation, and repaint only what changed.
- **Loading**: a file that does not parse fails (it reported success); a session with one unreachable part opens the
  rest; a duplicate load in flight says it is already loaded.
- **Pointer rate**: the 3D hover, quick-settings sliders, the transfer function, the snapshot list, the Image preview,
  the animation source, the settings and selection sources, the journal (no `describe` per run, `record: false` for a
  scrubber's steps and link hovers).
- **Typing**: the eager built-in features' areas and their operations in core areas (`load.entry`, `canvas.lighting`,
  `session.save`, `measure`, `snapshot`, `export`…) are typed as present like the scene's areas, and a missing one names
  its feature; `api.area(id)` is typed; the spec's options are partial inside `ui.buttons`.
- **The spec check** runs once every command is in (overriding a page-only command failed creation), says an unknown
  keybinding, and `viewer.events.error` reports an extension that did not load.
- **History** shows a redone step as done again; the shell keeps a theme only once the person picks one, retries an
  assistant that failed to start, and cleans up when it fails partway; a plugin task no longer wipes an extension's.
- **Size**: the assistant loads a model provider's SDK on first use (~350 kB not downloaded), the app builds carry
  zod's English locale only (~260 kB), extensions import `SceneRow` without the other panels.
- **PDB export** keeps modified residues in their chain; the PAE row no longer reads the whole matrix to show a size.

## Decisions for David

Each is confirmed, but a public API shape, a recorded design decision or a package boundary: yours to judge.

### A. The API, simply used

1. **Absent vs given.** `checkValues` fills every default, so operations use sentinels for "not given" (`-1`, `0`, an
   identity matrix, `''`): `structure.unitCell({ color: 'red' })` also makes the cell visible, `structure.transform({
   axis })` resets a moved structure. **Recommend**: optional params left `undefined`, no sentinels (one sweep).
2. **On/off.** `on: 'toggle'|'on'|'off'` defaulting to toggle (`view.expand`, `view.fullscreen`, `canvas.xr`), `on:
   boolean`, `'keep'|'on'|'off'` (`canvas.lighting`), `visible`, `enabled`, `play`; `canvas.clipping` overwrites both
   fields and `canvas.get` cannot read them. **Recommend**: desired state everywhere (`on: boolean`, `visible` for
   visibility), never a toggle by default; with (1), an absent field keeps its value; `canvas.get` reports everything.
3. **A representation's address.** By ref (`visibility`, `update`), component+type (`set`, `remove`), component
   (`material`), owner+ref-or-type (`volume.*`, `particles.*`). **Recommend**: one `representation` field everywhere,
   a ref or `{ component, type }` (the Structure panel's Remove now passes a ref, as an addition).
4. **An empty `structure`.** Three rules (`oneStructure`, the target resolver, `EveryStructureParam`); the param says
   "the one in view", which no code does. **Recommend**: one documented rule; destructive operations need a structure
   when several are loaded.
5. **Words.** `get`/`list`/`layers`; `add/remove/intersect` vs `union/subtract/intersect`; `load.entry.representation`
   is a preset; `structure.list` a state while `component.list` is an operation; a structure's label in some reads,
   `{ ref, label }` in others; `remove()` with nothing means "all", "the only one" or fails. **Recommend**: `list` for
   collections, `get` for one; one mode vocabulary; `preset`; `{ ref, label }`; empty `remove` = the only one, `all:
   true` to clear; old names as aliases until release.
6. **Loading.** `load.url` defaults to mmCIF though `'auto'` works; `load.volume` has no auto and defaults to binary
   CCP4; `load.mvs` ignores `.mvsx`; content is `data`, `bytes`, `json` or `file` by operation; a session opens three
   ways in two areas. **Recommend**: `format: 'auto'` by default with binary inferred; `{ url } | { data } | { file }`
   everywhere; one `session.open`.
7. **One state, two operations.** `camera.set` (no step) and `canvas.camera` (a step) set the camera mode; `structure.model`
   and `trajectory.frame` the model index (the first is now view state too). **Recommend**: one operation per state.
8. **Smaller ones.** `animation.play` needs full ids and plays the first one without a name; video's size vs width and
   height; extension areas not named for their extension (`anvil` → `membrane`, `video-export` → `video`); a sequence
   rail change bypasses the API; `measure.*` points cannot each name a structure; `export.camera`, `measure.model`,
   `state.tree` are view plumbing on the public API; the animation presets' commands are made once (an animation
   registered later gets none).

### B. What is public

9. **The package root and `./*`**: the root re-exports modes-era view models and every panel, not `createPrismViewer`
   or the spec; `./*` makes every file public; `@molstar/prism-viewer/api` exports the runner; `viewer.sources` has
   writers (`marks.set`); `targets/index` exports unused ported code. **Recommend**: the root as the embedder's entry, a
   curated export list before release, `sources` read-only, the unused targets code removed.
10. **The agent's spec and API**: extending skills imports the toolkit (zod) into the page; `ViewerAgent.resume` has its
    arguments in the reverse order of the session it exposes; configuring the assistant needs a hand-written loader.
    **Recommend**: default skill lists in a dependency-free module, additive spec fields, options objects,
    `assistant(options)`.
11. **Duplicates**: `StatusModel` = `NotificationsManager`; `ShellModel.layout` beside `viewState.save()`; `persist.key`
    namespaces nothing (each store repeats it); `StructureTab` has three ways to run Preset and Styles.
12. **The query languages' registry is page-wide**: two viewers on a page read and offer each other's languages.
    **Recommend**: per viewer if pages with several viewers matter; page-wide is simpler and costs nothing today.

### C. Size

13. **Identifier minification is off** in the app, distribution and slim builds (`minifyIdentifiers: false`, from
    #1542): about a third more bytes, ~280 kB gzipped on the full app. **Recommend**: on, unless something reads
    function or class names. **Resolved**: the override was a no-op (esbuild ignores it beside `minify: true`), so
    identifiers were minified all along and the estimate above was wrong; the override is removed (the output is
    byte-identical), and nothing reads names at runtime but diagnostics (architecture.md §8.1).
14. **MolViewSpec in the startup bundle** (~197 kB with io-ts and fflate; viewer-spec §2.1 lists it among the default
    behaviors): a lazy MVS feature with stub format providers.
15. **The Seq\* chunk** (1.23 MB, 446 kB gzipped, loaded with the first structure) is built from Seq\*'s default spec
    and registry, with dot plots, graphs and pairs the Prism documents do not use: a Prism Seq\* composition.
16. **Eager extensions in the app** (~190 kB: kinemage, anvil, assembly symmetry, model archive, the data-only ones;
    plan.md records which are lazy): lazy stubs owning their formats.
17. **Panels that show on demand are static imports**: the command menu, History, Settings (with the params form, 50
    kB) and the Structure panel's volume inspectors are bundled when their tabs are overridden away: features or `import()`.
18. **`sideEffects`** is declared by almost no package (~48 kB off the minimal build), but every module that registers
    on import (prism-ui's styles, Mol\*'s registries) must be listed.

### D. Choices made during the night (keep or change)

- Built-in areas, operations and view models are **typed as present and throw a clear error** when the spec left the
  feature out (simplicity for scripts), rather than optional types.
- The panels **keep their tab ids** (`home`, `structure`, `selection`, `help`; `core.tab.*` parts) though features bring
  them, so `layout.start` and `view.open` keep their simple values.
- **History and Settings stay in the core**; the volume inspectors stay with the Structure panel.
- **The runner re-checks `exists` only for loads** (a duplicate `load.*` in the history's queue says it is loaded); a
  generic re-check in the runner would cover every operation with `exists`.
- **The kept-records index grows** (journal.md §1 defers its compaction): compact it or bound it.

## Not done

- A fast path in the runner for pointer-rate view operations (measure first: the hover is per element, not per frame).
- `sources.formats` is read when asked, but the drop overlay and Home read it once when the shell mounts.
- The PAE row's size before the matrix is read is the model's polymer residues: right for a full ModelCIF PAE, not for a
  partial metric.
