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
controllers in extensions (`GeometryControls`, `Mp4Controls`, volseg entry data). Icons are not reused: the icon set is
Lucide (extracted by `scripts/icons.mjs`) plus a few molecular glyphs drawn in the same style
(`components/icon/molecular`).

## 2. Package

`packages/plugin/prism-ui` → `@molstar/plugin-prism-ui` (private until stable; public name decided at release). Same
package conventions as the rest of v6: ESM, `molstar-src` source condition, no barrels, `tsc -b` library build.

- **Per-file subpath exports.** `exports` maps `./*` to `src/*.ts` (`molstar-src`) and `lib/*.js`, so a module is
  imported by its path: `@molstar/plugin-prism-ui/reactive/observable-state`; the JSX import source is
  `@molstar/plugin-prism-ui/dom` (`dom/jsx-runtime`). The `molstar-src` target names `.ts`, so `.tsx` modules
  (components, views) resolve only from `lib/` for now; code inside the package imports relatively. `.` is
  `src/index.ts`, a placeholder until `createPluginUI` lands. The top-level `./_stories/*` and `./_test/*` are `null`
  (nested `_test/` paths are not); `tsc -b` leaves tests and stories out of `lib/`, and `files` leaves them and
  `*.stories.*` out of the published package.
- **Dependencies.** `@molstar/core`, `rxjs`, `tslib`; `react >= 18` is an optional peer used only by `react/`. The
  `tsconfig.json` references core only: the mocks need nothing from `@molstar/plugin`, `model` or `graphics` yet, and
  the references grow when sources are wired to a plugin. Storybook, Vite, lit, lightningcss and lucide-static are
  devDependencies pinned in the workspace catalog as prism-ui only; React and React DOM are devDependencies too (the
  `react/` tests), from the shared catalog.

```
packages/plugin/prism-ui/
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
                     state-tree, settings, export, sequence, canvas, command-menu; agent (opt-in, §6)
    styles/          registry (CSS registered by generated modules, adopted per document), base, reset.css,
                     theme.css (tokens, themes)
    _stories/        introduction.mdx, foundations, runtime and performance stories, story-root, mocks/
                     (view-model fakes)
    index.ts         placeholder; createPluginUI later
  scripts/           css.mjs (CSS → generated modules), icons.mjs (Lucide extraction)
  .storybook/        (dev only)
  package.json, tsconfig.json, tsconfig.storybook.json
```

Planned, not built (Phases 4–6 in [../plan.md](../plan.md)): `model/sources/` (plugin → `ReadonlyState` adapters) and
`model/actions/` (operations shared by view models and commands), view-model implementations over a real plugin,
`ui-commands/` (built-in, generated and contextual commands, formats, keybinding service), the built-in custom elements
(`<msp-plugin>`, embeddable views), `context.ts` (`PluginUIContext`), `spec.ts` (`PluginUISpec`), `default-spec.ts`
(`DefaultUIContributions`, `DefaultPluginUISpec`) and `createPluginUI` in `index.ts`.

Each `.css` file compiles into a committed module next to it (`x.css` → `_generated/x.css.ts`) that exports the class
names of a CSS module and registers its CSS with `styles/registry`; `MspElement` adopts the registered styles into its
document on connect. So tsc, Jest, bundlers and browsers need no CSS support, and a page gets the styles of exactly the
components it imports ([styling.md](styling.md)). `_generated/` folders are excluded from Biome.

## 3. Dependency rules

As the code holds them now (rxjs is allowed everywhere, and so is `import type` of `IconName` from `components/icon`:
icon names are data in models, commands and tabs):

| Folder                        | May import                                                                                         |
| ----------------------------- | -------------------------------------------------------------------------------------------------- |
| `reactive/`                   | rxjs                                                                                               |
| `forms/`                      | `reactive/`, `@molstar/core` (`util/param-definition`, `util/string`)                              |
| `commands/`                   | `reactive/`, `forms/`, `@molstar/core/util/param-definition`                                       |
| `styles/`                     | nothing (generated modules import the registry)                                                    |
| `dom/`                        | `reactive/`, `styles/` (an element adopts the registered styles)                                   |
| `react/`                      | `reactive/`, `dom/` (view contract types), `react` (peer, optional)                                |
| `model/`                      | rxjs, `reactive/`, `forms/` (types); no DOM, no components, no plugin yet                          |
| `shell/`                      | `reactive/`, `dom/` (view contract types)                                                          |
| `components/`                 | `reactive/`, `dom/`, `forms/`, `styles/`, `@molstar/core` utils; **not** `model/`, `shell/`, views |
| `views/`                      | everything above                                                                                   |
| `views/agent/`, `model/agent` | everything above; imported by nothing outside the agent module but stories (§6)                    |
| `_stories/`, `*.stories.tsx`  | anything; imported only by stories and tests (view tests use the mocks)                            |

Components are dumb: props and states in, callbacks out; views bind them to a model. Views read models through the
`model/` contracts and the shell through `ShellModel`, never the plugin directly (PD from `@molstar/core` is fine), so
the same views run over mocks now and over plugin-backed models later. Once sources land, `model/sources/` and
`model/actions/` add `@molstar/plugin`, `@molstar/model` and `@molstar/graphics`; views and components still do not
import them.

Nothing in prism-ui imports `@molstar/plugin-ui`, and nothing outside `react/` imports React. `reactive/`, `forms/` and
`commands/` are candidates to move to `@molstar/core` or `@molstar/plugin` later (headless commands); keeping them free
of DOM and plugin imports now makes that a file move.

Enforcement: only the agent rule is checked today (`views/agent/_test/boundary.test.ts`, a Jest test that scans the
sources for relative imports of the agent module). The other rules hold by review; import checks like the other v6
boundaries are not yet in CI.

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
- prism-ui has its **own `PluginUIContext`** (`@molstar/plugin-prism-ui/context`, [shell.md](shell.md) §6; not built
  yet) with the new UI fields (shell, commands, contributions, keybindings, sources, view models). It is a different
  class from `@molstar/plugin-ui`'s `PluginUIContext`; code that needs both imports them under different names. The
  UI-only maps on core `PluginContext` (`customStructureControls`, …) are not read by prism-ui; their typed replacements
  live on the new context, and removing them from core is the plugin-composition follow-up and an API-ledger entry once
  the old UI is retired.
- The viewer is not switched by a flag: `apps/prism`, a temporary copy of the viewer on prism-ui (same extensions,
  presets and URL options, its own entry in the app build), is built and served next to `apps/viewer` so the new UI can
  be shown and used early; the classic viewer stays untouched until `apps/prism` reaches parity, then switches over and
  `apps/prism` goes away ([../plan.md](../plan.md) Phase 6; not started).
- Old React building blocks used by extensions (`CollapsableControls`, `ParameterControls`, `ApplyActionControl`,
  `UpdateTransformControl`) keep working in the old UI; ports provide prism-ui contributions next to them
  ([../analysis/extensions.md](../analysis/extensions.md)).
- When parity is reached: prism-ui becomes `@molstar/plugin-ui`, the old package becomes `@molstar/plugin-ui-legacy` for
  one major version, or is removed. Decided at that point.

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

Opt-in modules are not imported by the core at all. The agent (`views/agent/`, `model/agent`) is the first: nothing else
in the package imports it (a test enforces it, §3), and the core reaches it only through ports: the command menu's
`assistant` prop (`CommandMenuAssistant`), the shell view's `agent` slot after the status line's tasks, and tabs and
commands registered like any other. An app that does not import it ships none of it ([agent.md](agent.md) §6).

State (mocks): the registrations above are real `ShellModel`/`CommandManager` APIs; the shell story
(`views/shell/shell.stories.tsx`) assembles the default UI and the agent (with `ai` on) from them by hand. Not yet: the
declarative `UIContribution` object listed in the spec and activated after `plugin.init()` ([shell.md](shell.md) §5),
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

- Storybook 10 with `@storybook/web-components-vite` (stories return real DOM nodes), addons docs and a11y; a
  devDependency of prism-ui only. Vite runs only for Storybook: `molstar-src` condition (workspace packages from
  source), `rxjs` deduped, the monorepo root allowed, and `cssModulesPlugin()` from `scripts/css.mjs` keeps generated
  CSS modules current while it runs. `pnpm --filter @molstar/plugin-prism-ui storybook` / `build-storybook`; not in CI.
- A toolbar global sets the theme (`data-msp-theme` on the document root: light, dark, blue, auto).
- Story order: Overview (Introduction, Foundations), Shell, Views, Canvas, Components, Internals (the DOM runtime, and
  Performance: the renderer check). Every component but the color picker has a story; `Components/Params form` shows
  every PD control in its states.
- `_stories/story-root` hosts a story's JSX in an `MspElement`, so its bindings are disposed when Storybook removes it.
- `_stories/mocks/` implements the `model/` contracts with in-memory data (scene, home, snapshots, selection, agent with
  scripted turns, …); only stories and a few view tests (command menu, Structure tab) import them. Stories over a real
  `PluginContext` (plan's level 3) are not built yet.

## 9. Tests and checks

- Jest from the repo root (`pnpm jest`), tests in `_test/` folders next to the code as `*.test.ts` or `*.test.tsx`. The
  root `moduleNameMapper` maps `@molstar/plugin-prism-ui/*` to `src/`. DOM tests opt into jsdom per file
  (`/** @jest-environment jsdom */`) and name their JSX runtime with the pragma ([rendering.md](rendering.md) §4.1).
  Interaction is tested in Jest (keyboard, focus, events on real nodes), not with Storybook `play`.
- `tsc -b` builds `lib/` without tests and stories; `typecheck:storybook` (`tsconfig.storybook.json`) checks everything,
  stories and tests included.
- `styles:check` fails when a generated CSS module is missing, stale or orphaned, `icons:check` when the Lucide module
  is missing or stale; generated modules are committed. They are package scripts, not yet part of the root `test`.
