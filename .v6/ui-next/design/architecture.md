# ui-next architecture

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
- Rewriting next to the old UI keeps risk low: `@molstar/plugin-ui` stays the default until ui-next reaches parity; apps
  opt in per page.

What is reused rather than rewritten: PD semantics (the param-editor contract `controlFor(type)`), icon SVG data,
sequence wrappers (`sequence/*.ts`, UI-free), screenshot helpers, markdown extension handling, controllers in extensions
(`GeometryControls`, `Mp4Controls`, volseg entry data).

## 2. Package

`packages/plugin/ui-next` → `@molstar/plugin-ui-next` (private until stable; public name decided at release). Same
package conventions as the rest of v6: ESM, `molstar-src` source condition, explicit subpath exports, no barrels,
`tsc -b` library build.

```
packages/plugin/ui-next/
  src/
    reactive/        any-ui state + lib: observable-state, reactive-model, combine, transaction, disposable-group,
                     keyed-list, from-event                                     (no Mol*, no DOM, rxjs only)
    commands/        any-ui command core: registry, invocation, match, tail, keys + PD form adapter
                     (no Mol*, no DOM)
    dom/             JSX runtime, MspElement, control flow, behaviours (focus, popover, roving)   (DOM, no Mol*)
    react/           hooks + reactView + React 18 element wrappers                (optional React/Preact peer)
    model/
      sources/       plugin → ReadonlyState adapters                              (Mol*, no DOM)
      actions/       operations shared by view models and commands               (Mol*, no DOM)
      *.ts           view models (PluginUIModel, SceneModel, ParamsFormModel, …)  (Mol*, no DOM)
    ui-commands/     built-in, generated, contextual commands; formats; keybinding service
    shell/           ShellModel, contributions, spec types, layout persistence
    components/      dumb components: button, icon, menu, tree, tabs, slider, color, PD field controls
    views/           built-in views and viewport items (scene, inspector, state tree, snapshots, sequence, …)
    elements/        custom element classes and `defineElements`
    styles/          tokens.css, themes, reset
    context.ts       PluginUIContext (ui-next): extends PluginContext with the UI fields
    spec.ts          PluginUISpec (ui-next) type
    default-spec.ts  DefaultUIContributions, DefaultPluginUISpec
    index.ts         createPluginUI
  .storybook/        (dev only)
  package.json, tsconfig.json
```

## 3. Dependency rules

Enforced like the other v6 boundaries (import checks in CI):

| Folder                   | May import                                                                                          |
| ------------------------ | --------------------------------------------------------------------------------------------------- |
| `reactive/`              | rxjs                                                                                                |
| `commands/`              | `reactive/`, `@molstar/core/util/param-definition`                                                  |
| `dom/`                   | `reactive/`                                                                                         |
| `react/`                 | `reactive/`, `dom/` types, `react` (peer, optional)                                                 |
| `model/`                 | `reactive/`, `commands/`, `@molstar/plugin`, `@molstar/core`, `@molstar/model`, `@molstar/graphics` |
| `ui-commands/`, `shell/` | the above                                                                                           |
| `components/`            | `reactive/`, `dom/`, styles; **not** `model/` (dumb, props only)                                    |
| `views/`, `elements/`    | everything above                                                                                    |

Nothing in ui-next imports `@molstar/plugin-ui` or React outside `react/`. `reactive/` and `commands/` are candidates to
move to `@molstar/core` or `@molstar/plugin` later (headless commands); keeping them dependency-free now makes that a
file move.

## 4. Bringing in any-ui

The any-ui code (`src/state`, `src/lib`, `src/commands`, and their tests, as found in stories-untold's
`vendor/anystar/any-ui`) is copied into `src/reactive/` and `src/commands/` as **new Mol\* code**: each file gets the
standard Mol\* header, and from then on it is simply part of the new UI. There is no provenance file, no change ledger,
and no `STORIES-UNTOLD` markers; stories-untold is an inspiration, not an upstream to track.

```ts
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */
```

Adaptations made while copying or soon after:

| Area              | Change                                                                                                                                                                                                                                                                                                              |
| ----------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Style             | Biome format (2 spaces, repo quote style), `.js` relative imports instead of `.ts`, Mol\* headers, upstream issue markers removed                                                                                                                                                                                   |
| Compiler          | Repo sets `useDefineForClassFields: false` and ES2022; `#private` fields and field initialisers are fine, verify `ReactiveModel`'s `Object.defineProperties` lifecycle guard under it                                                                                                                               |
| Tests             | Port vitest tests to jest under `_test/**/*.test.ts`; hook tests (`.tsx`) need the jest regex widened or moving to `.ts` with `createElement`                                                                                                                                                                       |
| `ReactiveModel`   | Add mount-scope `signal: AbortSignal`                                                                                                                                                                                                                                                                               |
| New helpers       | `keyedList`, `fromEvent`                                                                                                                                                                                                                                                                                            |
| Hooks             | Move to `react/`; React becomes an optional peer of the package                                                                                                                                                                                                                                                     |
| Commands: forms   | Done: `params(ctx, values) → PD.Params` replaces the declarative definition; `forms/flat-fields.ts` reads text, string selects, booleans and numbers as inline fields (`fields` hints add before-words and formats); other params are form-only and get their PD defaults; the invocation keeps its own value store |
| Commands: values  | `CommandValues` → `PD.Values`; `fits()` checks against PD (select options, numeric bounds) instead of value contracts                                                                                                                                                                                               |
| Commands: formats | Keep the `CommandFormat` port; attach per field via `formats`                                                                                                                                                                                                                                                       |
| Commands: API     | Add `execute(id, values)`, owner-scoped registration, category/prefix filters; keep `assist` as a type-only port                                                                                                                                                                                                    |
| Dropped           | `src/declarative/**`, `@anystar/ui-declarative` (zod, yaml), generated schema                                                                                                                                                                                                                                       |

## 5. Coexistence and migration

- ui-next and `@molstar/plugin-ui` can run in the same repo and the same page (different class names, no shared
  globals). They do not share a `PluginContext` at the same time.
- ui-next has its **own `PluginUIContext`** (`@molstar/plugin-ui-next/context`, [shell.md](shell.md) §6) with the new UI
  fields (shell, commands, contributions, keybindings, sources, view models). It is a different class from
  `@molstar/plugin-ui`'s `PluginUIContext`; code that needs both imports them under different names. The UI-only maps on
  core `PluginContext` (`customStructureControls`, …) are not read by ui-next; their typed replacements live on the new
  context, and removing them from core is the plugin-composition follow-up and an API-ledger entry once the old UI is
  retired.
- The viewer gets `?ui=next` (and a `Viewer.create(…, { ui: 'next' })` option) during development.
- Old React building blocks used by extensions (`CollapsableControls`, `ParameterControls`, `ApplyActionControl`,
  `UpdateTransformControl`) keep working in the old UI; ports provide ui-next contributions next to them
  ([../analysis/extensions.md](../analysis/extensions.md)).
- When parity is reached: ui-next becomes `@molstar/plugin-ui`, the old package becomes `@molstar/plugin-ui-legacy` for
  one major version, or is removed. Decided at that point.
