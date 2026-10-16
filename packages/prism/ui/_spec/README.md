# Prism UI: the specification

`@molstar/prism-ui`: the common base of Prism, Mol\*'s new UI. Its own DOM runtime (JSX without React), reactive state,
the shell and its managers, the command menu, the components, the styles. It knows nothing of Mol\*'s plugin: the viewer
(`@molstar/prism-viewer`), the stories and the assistant build on it, and it imports none of them
(`src/_test/opt-in-boundary.test.ts`). It is a library, not a framework: no modules, contributions registry or modes;
each app composes it in its own code. The code is the truth for how.

| File                                 | What it covers                                                                                    |
| ------------------------------------ | ------------------------------------------------------------------------------------------------- |
| [library.md](library.md)             | what the package gives an app and leaves to it; its entries, shadow roots and keys, package rules |
| [rendering.md](rendering.md)         | the DOM runtime: JSX, scopes (active and paused), control flow, elements, views, React/Preact     |
| [view-models.md](view-models.md)     | reactive primitives, model conventions, sources and records, lifetimes, results and state         |
| [shell.md](shell.md)                 | the shell: rail, panel, cards, pages, docks, canvas chrome, status column, compact, the model     |
| [commands.md](commands.md)           | commands, the manager, the command menu, keybindings                                              |
| [styling.md](styling.md)             | CSS modules, tokens, themes, spacing, the build, the look's rules, performance (§5)               |
| [accessibility.md](accessibility.md) | keyboard, focus, ARIA, the live region, touch                                                     |

## Folders

| Folder        | Holds                                                                                                                                                                                                   |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `dom/`        | the runtime: `render`, JSX (`jsxImportSource: @molstar/prism-ui/dom`), `Show`, `For`, scopes (`bind`, `onCleanup`, `onActive`), events, motion, views (`domView`), the trees' keyboard (`tree-keys.ts`) |
| `reactive/`   | `observableState`, `ReadonlyState`, `mapState`, keyed selectors, transactions: state as rxjs observables                                                                                                |
| `shell/`      | `ShellModel` (tabs, docks, cards, toolbar, compact mode), the status column, the theme, fullscreen                                                                                                      |
| `managers/`   | the history (checkpoints, steps, groups, undo, redo, undo turn), navigation, view state, notifications                                                                                                  |
| `commands/`   | the command manager (registration, matching, invocation, shortcuts), the command types                                                                                                                  |
| `forms/`      | `ParamsForm` over Mol\*'s PD params (no plugin: PD is `@molstar/core`'s)                                                                                                                                |
| `components/` | buttons, choices, legends, menus, icons (a generated Lucide subset plus molecular icons), markdown, …                                                                                                   |
| `views/`      | the shell's views (`ShellView`, the canvas chrome, the status line), the command menu, Settings, Help, …                                                                                                |
| `model/`      | DOM-free models the views share: the canvas ink, playback, settings                                                                                                                                     |
| `styles/`     | the tokens (`theme.css`: light, dark, blue, high contrast; `prefers-contrast`), the reset, the style registry that adopts sheets into shadow roots                                                      |
| `react/`      | an adapter for React hosts (`reactView`, hooks over observables)                                                                                                                                        |
| `testing/`    | stories' fixtures                                                                                                                                                                                       |

## Contracts

- **Views** are functions over models: they bind to observables, hold no state of their own beyond the view, and release
  what they subscribe (`bind`, or a `const` unsubscribed in `onCleanup`: `src/_test/subscriptions.test.ts`). Hidden UI
  is inert and paused (`onActive`, `Scope.setActive`; [styling.md](styling.md) §5).
- **Styles** are CSS modules compiled to `_generated/` (`pnpm -C packages/prism/ui styles`; `styles:check` verifies the
  files and the spacing scale `--msp-space-1..6`); every color and size from `--msp-*` tokens, which cross a shadow
  root, so a host restyles through them. In compact mode, fixed-size targets grow to `--msp-touch-target` (44px).
- **The history** (`managers/history.ts`) takes a checkpoint through a `HistoryPort` (the viewer's: the plugin's state
  without view state) before each step; three slow checkpoints in a row turn checkpoints off (undo off) until the scene
  is cleared.
- **Commands** are `{ id, title, category, icon?, alias?, shortcut?, params?, run, … }`; shortcuts go to the viewer that
  has focus (`ownsKey`), never the window ([commands.md](commands.md)).
- **DOM-free parts** (`reactive/`, `commands/`, `forms/`' model, `managers/`, `model/`, `shell/`) run in Node:
  `src/_test/dom-free.test.ts` ([library.md](library.md) §1).
