# ui-next plan

Phases are buildable and reviewable on their own. Storybook comes first so the UI is mocked before it is wired.

## Phase 0: package and Storybook skeleton

- Create `packages/plugin/ui-next` (`@molstar/plugin-ui-next`, private): `package.json`, `tsconfig.json` (references
  core, plugin, model, graphics), exports map, Biome and jest coverage, workspace registration.
- Storybook as a **devDependency of ui-next only**:
  - Framework `@storybook/web-components-vite` (ui-next renders real DOM nodes; stories return them). React stories for
    `ui-next/react` use a second, small config or the React renderer in the same Storybook via composition.
  - Vite is used only by Storybook. Its config sets `resolve.conditions: ['molstar-src', …]` so workspace packages
    resolve to source, dedupes `rxjs`, and allows the monorepo root for fixtures.
  - Addons: docs (MDX for the design principles and tokens), a11y. Toolbar global for theme (`data-theme`).
  - `pnpm --filter @molstar/plugin-ui-next storybook` / `build-storybook`; not part of CI builds at first.
- Delegation: this phase is mechanical; a smaller model (Sonnet or Haiku) can do it from this section and the existing
  package templates.

Acceptance: an empty story renders; `pnpm typecheck` and `pnpm lint` pass with the new package.

## Phase 1: primitives

- Copy any-ui's state, lib, and command code into `reactive/` and `commands/` as new Mol\* code with the Mol\* header,
  then adapt ([design/architecture.md](design/architecture.md) §4); jest tests.
- `dom/` runtime: `jsx`, reactive props/children, scopes, `Show`/`For`/`Switch`/`Dynamic`, `MspElement`,
  `defineElements`, deferred disconnect.
- `react/` adapter.
- Tokens, themes, reset; the CSS Modules build step.

Acceptance: unit tests for runtime and primitives; a Storybook "Foundations" page (tokens, themes, icons).

## Phase 2: renderer spike (decision gate)

Build the PD form, the state tree (5,000 cells, update bursts), and the command menu in the `dom` runtime and in Preact,
in Storybook ([design/rendering.md](design/rendering.md) §3). Record LOC, update cost, and friction in
`.v6/ui-next/spike.md`. Decide; if Preact wins, `dom/` shrinks to the element base and the view contract stays.

## Phase 3: mocked UI in Storybook

Dumb components (`components/`) and views over **fake sources**: the shell with docks and activity bars, scene tree,
inspector, command menu (with inline arguments), dialogs, notifications, status bar, viewport overlays with a
placeholder canvas, compact mode. Interaction tests with `play`. This is where the layout and command UX are iterated on
before any wiring.

Story levels:

1. **Components**: props only, inline data.
2. **Views over fake sources**: view models are real, sources are in-memory fakes
   (`createFakeSources({ structures: […] })`), so stories exercise real model logic without a plugin.
3. **Integration**: a real `PluginContext` (no canvas, or with canvas in one story), local fixtures from `data/`, real
   sources; used for the scene tree, inspector, commands.

## Review gate: spec review before the full build

Once the mocks settle the UX, update the designs (`PluginUIContext`, `PluginUISpec`, contributions, command and view
contracts) to match what the mocks taught, then run a deep multi-agent review (ultracode) of the spec and mocks before
phases 4–8 build it out. Findings are folded into the designs; the gate passes when the contracts are frozen for the
first implementation.

Carry into the proper implementation (from the mock review):

- **Every PD control in Storybook.** One story per param type (text, numeric/slider, select, boolean, color, mapped,
  flat and nested groups, on/off mapped, lists, …) in each state, so the renderer's spacing and alignment are checked in
  isolation, not found in whichever panel happens to use a control.
- **A UX pass on every panel, not a port.** The mocks prove the layout; the content is still a translation of the old UI
  in places. Export's format sections (rows of GLB / STL / OBJ / USDZ / mmCIF buttons) are the clearest case: design
  each panel for its task (choose what, then one action) rather than listing every variant as a button.

## Phase 4: sources, view models, actions

**Core tweak first: an action's apply label.** `StateAction.build`'s `display` (`packages/core/src/state/action.ts`,
`{ name, description? }`) gains `applyLabel?: string`, so an action names what applying does ("Download", "Open",
"Load") instead of a generic "Apply": `DownloadStructure` and `DownloadDensity` → Download, `OpenFiles` and the session
action → Open, `LoadTrajectory` → Load. The old UI ignores it; ui-next's action blocks use it and fall back to "Apply"
(transformer updates stay "Update"). Icons stay a UI concern: ui-next maps action ids to icons in its contributions,
core does not learn about them.

`model/sources/*` with tests against a real plugin and fixtures; view models in the order Scene → ParamsForm/Inspector →
Selection/Focus → Measurements → Snapshots/Animation/Trajectory → Viewport/Screenshot → Sequence → Tasks/
Notifications/Log → StateTree → Help. Actions extracted from the old UI where logic lives in views.

## Phase 5: commands

Registry wiring, `CommandContext`, built-in commands (map in [design/commands.md](design/commands.md) §8), generated
commands from registries, contextual commands, formats, keybinding service, command menu. Command tests run headless.

## Phase 6: shell and integration

`ShellModel`, contributions, `PluginUISpec` (ui-next), `createPluginUI`, `<msp-plugin>`, persistence, Layout bridge,
viewer `?ui=next`. Default contributions reach functional parity with [analysis/current-ui.md](analysis/current-ui.md)
§5.

## Phase 7: extension ports

In order of effort: model-export, geo-export (commands only) → assembly-symmetry, kinemage (inspector sections +
commands) → zenodo, mp4-export (dialogs) → model-archive PAE (view + embeddable element) → volumes-and-segmentations
(view). Generated commands cover the data-only extensions. Each port is a `ui-next` module in the extension package,
next to its existing `ui.tsx`.

## Phase 8: apps

docking-viewer (viewport items), mvs-stories (embeddable elements), mesoscale-explorer (contributed views; the
modularity test), examples as API checks (`react`, `alpha-orbitals`, `ligand-editor`).

## Cross-cutting

- API ledger and migration-map entries for new public APIs and for deprecations (core `layout`, `leftPanelTabName`, UI
  maps on `PluginContext`).
- Accessibility pass per view (keyboard-only use of every view; command menu as combobox).
- Performance budget: scene tree with 1,000 components and state tree with 5,000 cells update in under one frame for a
  single-node change.
