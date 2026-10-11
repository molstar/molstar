# ui-next plan

Phases are buildable and reviewable on their own. Storybook comes first so the UI is mocked before it is wired.

Where it stands: phases 0, 1 and 3 are done in Storybook (phase 3 ran ahead of phase 2: every mock is built on the `dom`
runtime). Next, in order: the designs describe what the mocks settled, the renderer check (phase 2), then the review
gate. Phases 4–8 have not started, apart from the core `triggerLabel` (phase 4) and the shell's own models the mocks run
on (phase 6).

## Phase 0: package and Storybook skeleton

**Done.**

- `packages/plugin/ui-next` (`@molstar/plugin-ui-next`, private, `6.0.0-dev.0`): `package.json` with the explicit
  exports map (`molstar-src` to `src/`, `_stories` and `_test` not exported), `tsconfig.json` (a reference in the root
  `tsconfig.json`), workspace registration, Biome, and jest through the root config (`moduleNameMapper` to `src/`; tests
  in `_test` folders).
- Storybook as a **devDependency of ui-next only** (`.storybook/`):
  - Framework `@storybook/web-components-vite` (ui-next renders real DOM nodes; stories return them).
  - Vite is used only by Storybook. Its config resolves workspace packages to source (`molstar-src` before the default
    conditions), dedupes `rxjs`, allows the monorepo root, and runs the CSS-module generator as a plugin.
  - Addons: docs (the MDX introduction: where to look, the principles), a11y. A toolbar global for the theme (light,
    dark, blue, auto) sets `data-msp-theme` on the document root.
  - `storybook` / `build-storybook` scripts; not part of CI builds. Also `styles` / `styles:check` / `styles:watch` (CSS
    modules) and `icons` / `icons:check` (the icon set).
- Not done: React stories for `ui-next/react` (a second config or composition); the adapter is covered by jest tests.

## Phase 1: primitives

**Done.**

- any-ui's state, lib, and command code as new Mol\* code with the Mol\* header
  ([design/architecture.md](design/architecture.md) §4): `reactive/` (`observableState`, `ReadonlyState`,
  `ReactiveModel`, transactions, `combineBehaviors`, `DisposableGroup`, shallow equality) and `commands/`
  (`CommandManager`, `CommandInvocation`, fuzzy matching, the argument tail reader, key chords, the command types with
  `CommandAssistant`); jest tests.
- `dom/` runtime: `jsx`, reactive props/children, scopes, `Show`/`For`/`Switch`/`Match`/`Dynamic`, `MspElement`,
  `defineElements`, deferred disconnect, and the view contract (`domView`, `elementView`).
- `react/` adapter: `useObservable`, `useObservableProp`, `useObservableSelection`, `useReactiveModel`, `reactView`,
  `MspElementHost`.
- Tokens, themes, reset; the CSS Modules build step (lightningcss into `_generated/*.css.ts`,
  [design/styling.md](design/styling.md) §3); the icon set (used Lucide icons extracted by `scripts/icons.mjs`, plus
  hand-written molecular glyphs).
- `forms/`: `ParamsFormModel` and the flat-field helpers that the command menu's inline sentence uses.

Storybook: Overview › Foundations (color tokens, surfaces, themes, type and space), Components › Icon, Internals ›
Runtime (bindings, keyed lists, conditionals).

## Phase 2: renderer check (decision gate)

**Pending; next after the designs are synced.** The mocks already build the PD form (`components/params`) and the
command menu (`views/command-menu`) on the `dom` runtime, so two of the spike's three pieces exist on one side; the
Preact side was not built. What the mocks do not show is cost at scale: the State tree tab's mock holds one structure
and one map. A story for the `dom` side exists (`_stories/performance.stories.tsx`, Internals › Performance: the State
tree tab at 5,025 cells under a burst of relabels, select, collapse and expand, next to plain rows that follow their
record by id; wall time per case, layout forced, no mutation counts); no result is recorded yet.

The check: the state tree with 5,000 cells and bursts of updates on the `dom` runtime, measured in DOM mutations and
time per burst against the budget under Cross-cutting (a single-node change in under one frame), with the same tree in
Preact for comparison ([design/rendering.md](design/rendering.md) §3). Record LOC, update cost, and friction in
`.v6/ui-next/spike.md`. Decide; if Preact wins, `dom/` shrinks to the element base and the view contract stays.

## Phase 3: mocked UI in Storybook

**Done** (the mocks the review will look at). Dumb components (`components/`) and views over mocked view models,
iterated on until the layout and the command UX settled ([design/shell.md](design/shell.md),
[design/commands.md](design/commands.md), [design/agent.md](design/agent.md)):

- **Shell › Shell**: the whole app over a placeholder canvas: Empty (the wordmark), Loaded, a session with playback,
  many structures, an EMDB map, a dark canvas background, the command menu open, Export at a section, the pinned panel,
  remote states (few and many), compact mode, Expand or Fullscreen by context (real frames), the assistant (with a
  model, and before setup).
- **Views**: Home tab, Structure tab (representation, structure, streaming density, many structures, EMDB map and model,
  the tree alone). The other tabs and cards (Selection with measurements and superposition, Snapshots, State tree,
  Settings, Log and History, Assistant, Animation, Canvas settings, Export, Remote states) are seen in the Shell
  stories.
- **Canvas**: Playback (snapshots and trajectory, many snapshots, trajectory alone).
- **Components**: action footer, button, choices, drag to reorder, icon, markdown, params form (every PD type in its
  states), scroll area, select, text input, tooltip, wordmark.
- Interaction tests with `play` in the stories, and jest tests next to the components and views.

Not built in the mocks: dialogs (open as form, Zenodo import), a Help view (help is the command menu's `?` list: the
prefixes and every command with a shortcut), content in the left and bottom docks (`ShellView` renders them; no story
fills them).

Story levels:

1. **Components**: props only, inline data. Built.
2. **Views over fake data**: built over **mocked view models**: `model/*.ts` holds the view-model interfaces and
   `_stories/mocks/*` implements them with in-memory data and the shell's own models (`ShellModel`, `StatusModel`,
   `ParamsFormModel` are real). The design's form, real view models over in-memory fake sources
   (`createFakeSources({ structures: […] })`), comes with phase 4, when the view models exist.
3. **Integration**: a real `PluginContext` (no canvas, or with canvas in one story), local fixtures from `data/`, real
   sources; used for the scene tree, inspector, commands. Not built; comes with phases 4–6.

## Review gate: spec review before the full build

**Next after the renderer check.** The designs (`PluginUIContext`, `PluginUISpec`, contributions, command and view
contracts) are updated to describe what the mocks settled, each with its "State (mocks)" and "Not yet" notes; then a
deep multi-agent review (ultracode) of the spec and mocks runs before phases 4–8 build it out. Findings are folded into
the designs; the gate passes when the contracts are frozen for the first implementation.

Carry into the proper implementation (from the mock review):

- **Every PD control in Storybook.** Done (`Components/Params form`): every type alone and in its states, on a panel and
  in a well, including multi-select, vec3, mat4, interval, line graph, color list, object list, url, converted,
  conditioned and script (More types; plain values have no row), and help with legends (Help and legends).
- **The view-model interfaces are the contracts the build implements.** `model/*.ts` holds them (scene, selection,
  measurements, superposition, snapshots, playback, home, settings, state tree, image export, agent); the mocks in
  `_stories/mocks/*` implement them over in-memory data and stay as the stories' data until real view models over fake
  sources replace them. What the review freezes there is what phase 4 builds against.
- **The agent stays opt-in.** `views/agent/_test/boundary.test.ts` keeps the rest of the package (stories aside) from
  importing `views/agent` and `model/agent`; the wiring (phase 5) goes through the same ports (the command menu's
  `assistant`, `ShellView`'s `agent` slot in the status column, contributed tabs and commands), not new imports from the
  core.
- **Theming and override docs after the review** (data-msp-part names, component tokens, a `ui.theme` spec option, a
  Storybook "Theming & overrides" page): once the contracts are frozen.
- **A UX pass on every panel, not a port.** The mocks prove the layout; the content is still a translation of the old UI
  in places. Export's format sections (rows of GLB / STL / OBJ / USDZ / mmCIF buttons) are the clearest case: design
  each panel for its task (choose what, then one action) rather than listing every variant as a button.

## Phase 4: sources, view models, actions

**Done: an action's trigger label.** `display` of state actions and transformers (`packages/core/src/state/action.ts`,
`transformer.ts`) has `triggerLabel?: string`: what running it does, for its button. Set on the root actions:
`DownloadStructure`, `DownloadDensity`, `DownloadFile` → Download; `OpenFiles` → Open; `LoadTrajectory` → Load;
`AddTrajectory` → Add. The old UI ignores it; ui-next's action blocks are to take it as `HomeAction.applyLabel` (the
mocks set the same labels by hand; nothing reads `triggerLabel` until the view models exist) and fall back to "Apply"
(transformer updates stay "Update"). Icons stay a UI concern: ui-next maps action ids to icons in its contributions,
core does not learn about them.

**Transform controls store their params correctly.** Forms for a cell's transform (updates, the state tree, inspectors)
emit the complete new params, each of its PD type (Vec3, Mat4 as arrays of their size; colors as numbers; assets kept as
`Asset` objects; object lists as new lists), never a partial object or a shared reference, and skip the update when
nothing changed (`PD.areEqual`), so `state.build().to(ref).update(params)` records exactly the new version.

**Fewer re-renders: combine, then debounce.** Where a view depends on several streams (a row's value, its visibility,
its help; a tree row's label, description, visibility), combine them (`combineLatest`) and add a short debounce (a
frame, `debounceTime(0)` or `auditTime(16)` with a scheduler) where a burst of updates is expected (state tree rebuilds,
selection changes), so one change renders once.

`model/sources/*` with tests against a real plugin and fixtures; view models in the order Scene → ParamsForm/Inspector →
Selection/Focus → Measurements → Snapshots/Animation/Trajectory → Viewport/Screenshot → Sequence → Tasks/
Notifications/Log → StateTree → Help. Actions extracted from the old UI where logic lives in views.

## Phase 5: commands

Registry wiring, `CommandContext`, built-in commands (map in [design/commands.md](design/commands.md) §8), generated
commands from registries, contextual commands, formats, keybinding service, command menu. Command tests run headless.
The four ways a command runs ([design/commands.md](design/commands.md) §2.1): `opens` as a declared field (commands that
open a panel or section; the mocks call the shell from `run`), color as an inline argument (the parser in
`components/color`), and "open as form" with the generated form in a dialog for params that are not flat (today they run
at their defaults).

Agents ([design/agent.md](design/agent.md)): the PD → JSON Schema generator, the state and focus accounts, the UI client
API over the command registry, checkpoints per turn (one at most: the next action releases it), the any-harness wiring
and the MCP server, and the Assistant tab's real key store and model lists (the mocks show them: keys kept until reload,
for the session, or in this browser, kept ones stored scrambled; models from the providers' `/models`). Then, as a
separate slice after the UI is done, Stories untold's MolQL selection skills.

## Phase 6: shell and integration

`ShellModel`, contributions, `PluginUISpec` (ui-next), `createPluginUI`, `<msp-plugin>`, persistence, Layout bridge.
Default contributions reach functional parity with [analysis/current-ui.md](analysis/current-ui.md) §5. The mocks
already run on the real `ShellModel` (docks, tabs, overlays and their sections, toolbar items, pinned, compact),
`StatusModel`, `ThemeModel` and `FullscreenModel` (`shell/`); this phase connects them to a plugin and a spec.

**`apps/viewer-next`: a temporary copy of the viewer on the new UI**, so it can be built, deployed and shown to people
(and used day to day) before anything moves over. Same extensions, presets and URL options as `apps/viewer`, the plugin
created with ui-next's `createPluginUI`; built and served next to the viewer (its own entry in the app build). The
classic viewer stays untouched until viewer-next reaches parity; then the viewer switches over and viewer-next goes
away. Started as soon as the shell and the first tabs run on a real plugin, so feedback comes early.

## Phase 7: extension ports

In order of effort: model-export, geo-export (Export sections + commands) → assembly-symmetry, kinemage (inspector
sections + commands) → mp4-export (an Export section, Animation, whose render holds the plugin as a busy-overlay task) →
zenodo (dialog) → model-archive PAE (view + embeddable element) → volumes-and-segmentations (view). The Export sections
are mocked (Shell › Export at section: Image, Session & state, Animation, Geometry, Models). Generated commands cover
the data-only extensions. Each port is a `ui-next` module in the extension package, next to its existing `ui.tsx`.

## Phase 8: apps

docking-viewer (viewport items), mvs-stories (embeddable elements), mesoscale-explorer (contributed views; the
modularity test), examples as API checks (`react`, `alpha-orbitals`, `ligand-editor`).

## Cross-cutting

- API ledger and migration-map entries for new public APIs and for deprecations (core `layout`, `leftPanelTabName`, UI
  maps on `PluginContext`).
- Accessibility pass per view (keyboard-only use of every view; command menu as combobox).
- Performance budget: scene tree with 1,000 components and state tree with 5,000 cells update in under one frame for a
  single-node change.
