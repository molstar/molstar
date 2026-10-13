# Prism UI architecture

## 1. Rewrite or refactor?

**Rewrite, alongside the old UI.** Reasons:

- The current UI is 17.9k LOC TS/TSX and 3.7k LOC SCSS, and almost no logic lives only there: the plugin, its managers,
  builders, and state tree do the work. The UI's job (projecting state, invoking operations) is exactly what the
  view-model layer replaces, and its structure (class components reading managers ad hoc, cards with embedded
  algorithms, fixed Sass layout) does not survive the change incrementally.
- The pieces that do hold logic move into actions or view models during the rewrite (superposition modes, lazy volume
  loading, quick styles, focus targets, snapshot key handling).
- Commands can be added to the old UI only by building the same registry and view models anyway; doing that first and
  then re-rendering with a new shell is the same work without the old React constraints.
- Rewriting next to the old UI keeps risk low: `@molstar/plugin-ui` stays the default until prism-ui reaches parity;
  apps opt in per page.

What is reused rather than rewritten: PD semantics (the params form in `components/params` renders `PD.Params` as they
are, one control per PD type; `value-ref` and `data-ref` not yet, they show their type name), and, when the views are
wired to a plugin, sequence wrappers (`sequence/*.ts`, UI-free), screenshot helpers, markdown extension handling,
controllers in extensions (`GeometryControls`, `VideoControls`, volseg entry data). Icons are not reused: the icon set
is Lucide (extracted by `scripts/icons.mjs`) plus a few molecular glyphs drawn in the same style
(`components/icon/molecular`).

## 2. Packages

Prism is a family of packages under `packages/prism/`, layered as stories-untold is (an application layer without UI,
the screens on it, a headless engine beside them; [../analysis/stories-untold.md](../analysis/stories-untold.md)):

| Path                       | Package                                       | Holds                                                                                                                                                                                                                                          | State                                                      |
| -------------------------- | --------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------- |
| `packages/prism/core`      | `@molstar/prism`                              | `PrismContext` (what every mode shares: the plugin behind sources and actions, the command registry, status, history and runs, settings, the module registry; [contexts.md](contexts.md)); `reactive/`, `forms/`, `commands/`, a core icon key | Planned (Phase 4): carved out of prism-ui with the context |
| `packages/prism/ui`        | `@molstar/prism-ui`                           | the common base every mode builds on: the DOM runtime, components, the shell frame and its views, the command menu, Settings, `PrismUIContext`, the React/Preact adapter, `testing/` (the mock frame)                                          | Built (mocks)                                              |
| `packages/prism/viewer`    | `@molstar/prism-viewer`                       | the viewer mode (Mol\*'s own UI): Home, the Structure tab, the state tree, snapshots, export, the sequence strip, the canvas's mode switch and selection bar, their `model/` contracts and mocks; opt-in                                       | Built (mocks)                                              |
| `packages/prism/stories`   | `@molstar/prism-stories`                      | story mode: the story model, the Story tab, its mocks; opt-in ([modes.md](modes.md))                                                                                                                                                           | Built (sketch)                                             |
| `packages/prism/assistant` | `@molstar/prism-assistant`                    | the assistant: its model contracts, the turn card, History, the Assistant tab, its mocks; opt-in (§6)                                                                                                                                          | Built (mocks)                                              |
| later                      | `@molstar/prism-headless`, `@molstar/prism-…` | the headless engine (Node: tests, CLI, evals, the agent's tools without a page; on `@molstar/prism` as `@molstar/plugin-headless` is on the plugin), the new sequence viewer, the narrative features                                           | Planned                                                    |

`apps/prism` is the app built on them (§5), composed of the modes it includes: the viewer, story mode, or both
([modes.md](modes.md)). The viewer keeps its files' paths (`views/structure/…`, `model/scene`), so the designs' paths
hold relative to whichever package has them. All are private until stable; public names are decided at release. Same
conventions as the rest of v6: ESM, `molstar-src` source condition, `tsc -b` library build.

- **Per-file subpath exports.** `exports` maps `./*` to `src/*.ts` (`molstar-src`) and `lib/*.js`, so a module is
  imported by its path: `@molstar/prism-ui/reactive/observable-state`; the JSX import source is `@molstar/prism-ui/dom`
  (`dom/jsx-runtime`). The `molstar-src` target names `.ts`, so `.tsx` modules (components, views) do not resolve by
  path from source; code inside the package imports relatively, and packages built on prism-ui import `.` instead.
- **One entry, `.` (`src/index.ts`)**, for packages built on prism-ui: the runtime, reactive state, the shell's models,
  commands, the components, and the shell views a contribution sits in (`ShellView`, `StatusPill`, `LogView`,
  `SectionList`, `CommandMenu`). `@molstar/prism-ui/testing/index` is the same for stories and tests. The top-level
  `./_stories/*` and `./_test/*` are `null` (nested `_test/` paths are not); `tsc -b` leaves tests and stories out of
  `lib/`, and `files` leaves them and `*.stories.*` out of the published package. Every Prism package keeps its stories
  in `_stories/` folders next to the code they show.
- **Dependencies.** `@molstar/core`, `rxjs`, `tslib`; `react >= 18` is an optional peer used only by `react/`. The
  `tsconfig.json` references core only: the mocks need nothing from `@molstar/plugin`, `model` or `graphics` yet, and
  the references grow when sources are wired to a plugin. Storybook's types, lit, lightningcss and lucide-static are
  devDependencies (Storybook itself and Vite are the top-level `storybook/` package's, §8); React, React DOM and Preact
  are devDependencies too (the `react/` tests and stories).

```
packages/prism/ui/
  src/
    reactive/        any-ui state: observable-state, readonly-state, reactive-model, combine-behaviors, transaction,
                     disposable-group, shallow-equal, same-shallow                (no Mol*, no DOM, rxjs only)
    forms/           PD as forms: flat-fields (params read as inline fields), params-form (ParamsFormModel, the PD
                     form state of inspectors, settings, command forms)           (PD from @molstar/core, no DOM)
    commands/        any-ui command core: manager, invocation, match, tail, keys, types (no DOM)
    dom/             JSX runtime (jsx-runtime, render, scope, control-flow), MspElement + defineElements (element),
                     the view contract ViewRenderer/domView/elementView (view)   (DOM, no Mol*)
    react/           hooks, reactView, MspElementHost                             (optional React/Preact peer)
    model/           view-model contracts (SceneModel, HomeModel, SnapshotsModel, …): interfaces the views read,
                     implemented by mocks for now; canvas-theme (canvas color → chrome theme)
    shell/           ShellModel (docks, tabs, overlays, toolbar, section slots), StatusModel, ThemeModel,
                     FullscreenModel
    components/      dumb components, one folder each: button, icon, choices, select, text-input, color, params
                     (PD controls), scroll-area, tooltip, markdown, reorder, action-footer
    views/           built-in views: shell (shell view, docks, canvas chrome, log), home, structure, state,
                     state-tree, settings, export, sequence, canvas, command-menu
    styles/          registry (CSS registered by generated modules, adopted per document), base, reset.css,
                     theme.css (tokens, themes)
    testing/         in-memory mocks of the model/ contracts, the mock shell (ShellStory, taking a module's
                     contributions as an add-on), story() (a story's JSX in an element)
    _stories/        introduction.mdx, foundations, runtime and performance stories
    _test/           package-level tests (the agent boundary)
    index.ts         the public entry (above); createPrismUI later
  scripts/           css.mjs (CSS → generated modules, for prism-ui and packages built on it), icons.mjs (Lucide)
  package.json, tsconfig.json, tsconfig.storybook.json

packages/prism/assistant/
  src/
    model.ts         AgentModel, AISettingsModel and the turn, run and history types
    views/           agent-turn (turn card, run rows, Undo turn), history-tab, assistant-tab, ai-setup
    testing/         MockAgentModel, MockAISettings, mockAssistant, mockAgentAddOn (the agent in the mock shell)
    _stories/        assistant.stories.tsx: Shell › Assistant (set up, not set up)
    index.ts         the public entry
  package.json       depends on @molstar/prism-ui; `styles` runs prism-ui's css.mjs with `--prefix agent`
```

Planned, not built (Phases 4–6 in [../plan.md](../plan.md)): `@molstar/prism` with `PrismContext`, `model/sources/`
(plugin → `ReadonlyState` adapters) and `model/actions/` (operations shared by view models and commands), view-model
implementations over a real plugin, the built-in, generated and contextual commands; in prism-ui the keybinding service,
the built-in custom elements (`<prism-viewer>`, embeddable views), `context.ts` (`PrismUIContext`), `spec.ts`
(`PrismUISpec`), `default-spec.ts` (`DefaultUIContributions`, `DefaultPrismUISpec`) and `createPrismUI`.

Each `.css` file compiles into a committed module next to it (`x.css` → `_generated/x.css.ts`) that exports the class
names of a CSS module and registers its CSS with `styles/registry`; `MspElement` adopts the registered styles into its
document on connect. So tsc, Jest, bundlers and browsers need no CSS support, and a page gets the styles of exactly the
components it imports ([styling.md](styling.md)). `_generated/` folders are excluded from Biome.

## 3. Dependency rules

As the code holds them now (rxjs is allowed everywhere, and so is `import type` of `IconName` from `components/icon`:
icon names are data in models, commands and tabs). With the core carve-out (Phase 4) that type would make core import
prism-ui, so core gets its own icon key (a string-literal type) that prism-ui maps to glyphs (review M4):

| Folder                       | May import                                                                                         |
| ---------------------------- | -------------------------------------------------------------------------------------------------- |
| `reactive/`                  | rxjs                                                                                               |
| `forms/`                     | `reactive/`, `@molstar/core` (`util/param-definition`, `util/string`)                              |
| `commands/`                  | `reactive/`, `forms/`, `@molstar/core/util/param-definition`                                       |
| `styles/`                    | nothing (generated modules import the registry)                                                    |
| `dom/`                       | `reactive/`, `styles/` (an element adopts the registered styles)                                   |
| `react/`                     | `reactive/`, `dom/` (view contract types), `react` (peer, optional)                                |
| `model/`                     | rxjs, `reactive/`, `forms/` (types); no DOM, no components, no plugin yet                          |
| `shell/`                     | `reactive/`, `dom/` (view contract types)                                                          |
| `components/`                | `reactive/`, `dom/`, `forms/`, `styles/`, `@molstar/core` utils; **not** `model/`, `shell/`, views |
| `views/`                     | everything above                                                                                   |
| `testing/`                   | everything above; imported by stories and tests only (here and in packages built on prism-ui)      |
| `_stories/`, `*.stories.tsx` | anything; imported by nothing                                                                      |

Components are dumb: props and states in, callbacks out; views bind them to a model. Views read models through the
`model/` contracts and the shell through `ShellModel`, never the plugin directly (PD from `@molstar/core` is fine), so
the same views run over mocks now and over plugin-backed models later. Once sources land, `model/sources/` and
`model/actions/` add `@molstar/plugin`, `@molstar/model` and `@molstar/graphics`; views and components still do not
import them.

Nothing in prism-ui imports `@molstar/plugin-ui`, and nothing outside `react/` imports React. `reactive/`, `forms/`,
`commands/` and `model/` move to `@molstar/prism` with `PrismContext` (Phase 4); keeping them free of DOM and plugin
imports now makes that a file move. `@molstar/prism-assistant` imports prism-ui (its entry, and `.ts` modules by path);
prism-ui never depends on it or names it.

**The plugin stays behind one layer** ([contexts.md](contexts.md) §7), by folder role in any Prism package: only
`sources/` and `actions/` folders import `@molstar/plugin` (and `@molstar/model`/`@molstar/graphics` for the plugin's
data), in core, in a mode, or in an extension's prism entry. Views and view models reach the plugin through sources,
commands through actions.

Enforcement: between packages, the workspace check (`scripts/workspace/check.mjs`: declared dependencies, no
cross-package relative imports, no cycles) and `src/_test/opt-in-boundary.test.ts` (prism-ui's manifest and files never
name the agent or a mode). Planned with `@molstar/prism` (Phase 4): one reusable plugin-boundary check run on every
package with Prism modules, failing an import of `@molstar/plugin` outside `sources/` and `actions/`. Inside prism-ui
the folder rules hold by review.

## 4. Bringing in any-ui

The any-ui code (`src/state`, `src/lib`, `src/commands`, and their tests, as found in stories-untold's
`vendor/anystar/any-ui`) is copied into `src/reactive/`, `src/forms/` and `src/commands/` (hooks into `src/react/`) as
**new Mol\* code**: each file gets the standard Mol\* header, and from then on it is simply part of the new UI. There is
no provenance file, no change ledger, and no `STORIES-UNTOLD` markers; stories-untold is an inspiration, not an upstream
to track. Code adapted from a third party under a compatible license keeps the Mol\* header and says where it came from
in it (`components/scroll-area`, from Base UI's ScrollArea, MIT).

```ts
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */
```

Adaptations made while copying or soon after:

| Area              | Change                                                                                                                                                                                                                                                                                                                                                                                                                            |
| ----------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Style             | Done: Biome format (2 spaces, repo quote style), `.js` relative imports instead of `.ts`, Mol\* headers, upstream issue markers removed                                                                                                                                                                                                                                                                                           |
| Compiler          | Done: the repo targets ES2022 with `useDefineForClassFields: false` (the Jest transform sets the flag too); `#private` fields and field initialisers are fine; `ReactiveModel` keeps its `Object.defineProperties` lifecycle guard, and its lifecycle tests run under that setting                                                                                                                                                |
| Tests             | Done: vitest tests ported to Jest under `_test/**/*.test.ts(x)`; the root `testRegex` takes `.tsx`, so hook and DOM tests stay JSX (§9)                                                                                                                                                                                                                                                                                           |
| `ReactiveModel`   | Not yet: a mount-scope `signal: AbortSignal`                                                                                                                                                                                                                                                                                                                                                                                      |
| New helpers       | Not as `reactive/` helpers: keyed lists are the DOM runtime's keyed `For` (`dom/control-flow`); `fromEvent` not added                                                                                                                                                                                                                                                                                                             |
| Hooks             | Done: in `react/` (`useObservable`, `useObservableProp`, `useObservableSelection`, `useReactiveModel`); React is an optional peer of the package                                                                                                                                                                                                                                                                                  |
| Commands: forms   | Done: `params(ctx, values) → PD.Params` replaces the declarative definition; `forms/flat-fields.ts` reads text, string selects, booleans and numbers as inline fields (`fields` hints add before-words and formats); other params are form-only and get their PD defaults; the invocation keeps its own value store; `forms/params-form.ts` (`ParamsFormModel`) holds the full PD form for inspectors, settings and command forms |
| Commands: values  | Done differently: `CommandValues` is `FieldValues` (flat `string \| number \| boolean` by param name), not `PD.Values`; `fits()` and `readValue()` check a value against its field read from PD (select options, numeric bounds) instead of value contracts; `argsFor(params, values)` builds the run's arguments                                                                                                                 |
| Commands: formats | Done: the `CommandFormat` port (`CommandFormats` in the manager's options), attached per field by name via `fields: { x: { format } }`                                                                                                                                                                                                                                                                                            |
| Commands: API     | Done: owner-scoped registration (each `register` call is one owner and returns its disposer), command groups (`registerGroups`, a browsable tree), ranking by `category`; prefix filters are the command menu's (`CommandMenuPrefix`, provided by the host); `assist` became a real port (`CommandAssistant.propose` in the manager's options, `assist` on a command; [agent.md](agent.md)). Not yet: `execute(id, values)`       |
| Dropped           | `src/declarative/**`, `@anystar/ui-declarative` (zod, yaml), generated schema                                                                                                                                                                                                                                                                                                                                                     |

## 5. Coexistence and migration

- prism-ui and `@molstar/plugin-ui` can run in the same repo and the same page (different class names, no shared
  globals). They do not share a `PluginContext` at the same time.
- **`PrismContext` and `PrismUIContext`** ([contexts.md](contexts.md); not built yet). `PrismContext` embeds the plugin
  (`{ plugin: PluginContext }`; it does not extend it) and holds what every mode shares, with no UI, so it runs
  headless; a mode's, an extension's and the agent's models are modules; `PrismUIContext` extends it with the UI. An app
  runs one mode, chosen before it exists; `<prism-viewer>` is its element, with a shadow root. Neither is
  `@molstar/plugin-ui`'s `PluginUIContext`, so code that uses both has no clash. The UI-only maps on core
  `PluginContext` (`customStructureControls`, …) are not read by Prism; contributions replace them. Removing them from
  core waits for the classic UI's retirement (recorded in
  [../../plans/breaking-v6-changes.md](../../plans/breaking-v6-changes.md)).
- **Commands.** `PluginCommands` stay in the plugin core, the classic UI's bus. Prism's command registry is part of
  `PrismContext`, so commands run headless and are what the agent calls.
- **Layout.** `PrismUIContext`'s shell is Prism's layout. Prism neither reads nor writes core `PluginLayout` or
  `behaviors.layout.leftPanelTabName`; they stay as they are for the classic UI and are not deprecated in v6.
- **`apps/prism` is a new app**, not a copy of the viewer: the base Mol\* viewer plus Stories Untold's features for
  building narratives, with the new sequence viewer as a contribution. It need not keep `apps/viewer`'s URL options or
  globals; the classic viewer stays as it is ([../plan.md](../plan.md) Phase 6; not started).
- Old React building blocks used by extensions (`CollapsableControls`, `ParameterControls`, `ApplyActionControl`,
  `UpdateTransformControl`) keep working in the old UI; ports provide prism-ui contributions next to them
  ([../analysis/extensions.md](../analysis/extensions.md)).
- When Prism covers what the classic UI does, the classic UI is retired (then the deferred core removals land) or kept
  for a major version. Decided at that point.

## 6. How an app includes features

Everything beyond the frame is registered, and every registration returns its disposer, so a feature is added and
removed as a unit:

- `ShellModel`: `registerTabs` (a tab lands in its dock by `order`, a restored layout keeps its place),
  `registerOverlays`, `registerToolbar` (rail and canvas-stack buttons, by `group` and `order`), `registerSections`
  (sections in a slot: the Export panel's `export` slot, `inspector:<kind>` slots in the Structure tab). Extensions add
  sections to existing panels rather than new buttons or tabs, so the chrome does not grow with every extension.
- `CommandManager`: `register` (commands, one owner per call) and `registerGroups`.
- What is shown is chosen by the spec, not by unregistering: `ShellOptions.buttons` (`layout.buttons`, `{ xr: false }`)
  hides toolbar buttons by id; optional models (`HomeModel.remoteStates`) are absent unless the spec turns them on.
- Views a contribution renders follow the view contract `ViewRenderer` (`dom/view`): `domView` for the built-in JSX
  runtime, `elementView` for a custom element, `reactView` (`react/react-view`) for React or Preact
  ([rendering.md](rendering.md)).

Opt-in modules are packages the UI does not depend on. The agent (`@molstar/prism-assistant`) is the first: an app gets
it only by depending on it, prism-ui never imports it (§3), and the UI reaches it only through ports: the command menu's
`assistant` prop (`CommandMenuAssistant`), the shell view's `agent` slot after the status line's tasks, and tabs and
commands registered like any other. An app that does not import it ships none of it ([agent.md](agent.md) §6).

State (mocks): the registrations above are real `ShellModel`/`CommandManager` APIs; the mock frame (`testing/shell.tsx`,
`ShellStory`, with the modes it is given: `mockViewerMode` from the viewer, `mockStoryMode` from story mode,
`mockBlankMode` for the frame alone) assembles the UI from them by hand and takes a module's contributions as an add-on
(`MockShellAddOn`: tabs, one replacing a tab of the same id; the command menu's assistant; the agent slot; finished
commands). The agent's `mockAgentAddOn` (`@molstar/prism-assistant`, `testing/`) is the first. Not yet: the declarative
`UIContribution` object listed in the spec and activated after `plugin.init()` ([shell.md](shell.md) §5),
`DefaultUIContributions`, `agentContribution`, extension packages' `…/prism-ui` entries.

## 7. Changes to core

prism-ui asks core for data a UI needs, never for UI concepts (icons, layout stay in prism-ui):

- **`display.triggerLabel`** on state actions and transformers (`packages/core/src/state/action.ts`, `transformer.ts`):
  what running it does, for the button that runs it. Set on the root actions (`DownloadStructure`, `DownloadDensity`,
  `DownloadFile` → Download; `OpenFiles` → Open; `LoadTrajectory` → Load; `AddTrajectory` → Add). The old UI ignores it;
  prism-ui's action blocks take it as `HomeAction.applyLabel` and fall back to "Apply" (transform updates stay
  "Update"). The mocks set `applyLabel` directly; reading it from `display` comes with the sources.
- Planned, with the old UI's retirement: removing the UI-only maps from core `PluginContext` (§5).

## 8. Storybook: the mock harness

The UI is built as Storybook mocks first: views are real, the models they read are fakes, so layout and interaction
settle before anything talks to a plugin.

- Storybook 10 with `@storybook/web-components-vite` (stories return real DOM nodes), addons docs and a11y, in the
  top-level `storybook/` package (`@molstar/storybook`, private). `storybook/packages.ts` registers the packages it
  shows: every `*.stories.ts(x)` and `*.mdx` under a registered package's `src/` is collected (by convention in
  `_stories/` folders), and an entry names its sidebar group and the Vite plugins its stories need. Vite runs only for
  Storybook: `molstar-src` condition (workspace packages from source), `rxjs` deduped, the monorepo root allowed; for
  each Prism package `cssModulesPlugin()` from prism-ui's `scripts/css.mjs` keeps its generated CSS modules current, and
  for prism-ui the Preact plugin serves Internals › React and Preact ([rendering.md](rendering.md)).
  `pnpm -C storybook storybook` / `build-storybook` (or `mise run storybook`); not in CI.
- A toolbar global sets the theme (`data-msp-theme` on the document root: light, dark, blue, auto).
- Story order, all under one Prism group (the registry's `group`, Storybook's `titlePrefix`): Overview (Introduction,
  Foundations), Shell, Views, Canvas, Components, Internals (the DOM runtime, and Performance: the renderer check).
  Every component but the color picker has a story; `Components/Params form` shows every PD control in its states.
- `testing/story-root` (`story()`) hosts a story's JSX in an `MspElement`, so its bindings are disposed when Storybook
  removes it.
- `testing/` implements the `model/` contracts with in-memory data (scene, home, snapshots, selection, …), and
  `testing/shell` is the whole mock shell; stories, a few view tests (command menu, Structure tab) and the agent's
  stories import them. The agent's mocks (an agent with scripted turns) are in its package.
- The Prism packages built on prism-ui are registered too (viewer: Shell › Viewer mode, Views; stories: Shell › Story
  mode; assistant: Shell › Assistant). Stories over a real `PluginContext` (plan's level 3) are not built yet.

## 9. Tests and checks

- Jest from the repo root (`pnpm jest`), tests in `_test/` folders next to the code as `*.test.ts` or `*.test.tsx`. The
  root `moduleNameMapper` maps `@molstar/prism-ui/*` and `@molstar/prism-assistant/*` to their `src/`. DOM tests opt
  into jsdom per file (`/** @jest-environment jsdom */`) and name their JSX runtime with the pragma
  ([rendering.md](rendering.md) §4.1). Interaction is tested in Jest (keyboard, focus, events on real nodes), not with
  Storybook `play`.
- `tsc -b` builds `lib/` without tests and stories; `typecheck:storybook` (`tsconfig.storybook.json`) checks everything,
  stories and tests included.
- `styles:check` fails when a generated CSS module is missing, stale or orphaned, `icons:check` when the Lucide module
  is missing or stale; generated modules are committed. They are package scripts, not yet part of the root `test`.
