# prism-ui: the library

What `@molstar/prism-ui` gives an app built on it (the Prism viewer now, Prism stories later), and what it leaves to the
app (README §2.1). The designs behind it are [../design/rendering.md](../design/rendering.md),
[../design/styling.md](../design/styling.md), [../design/commands.md](../design/commands.md) and
[../design/shell.md](../design/shell.md) §1–4. Status: **accepted** (2026-10-15).

## 1. What it holds

| Folder        | What                                                                                                              | DOM |
| ------------- | ----------------------------------------------------------------------------------------------------------------- | --- |
| `reactive/`   | `observableState`, `ReadonlyState`, `batch`, `ReactiveModel`, `combineBehaviors`, `keyedSelector`, transactions (the optimistic coordinator), `DisposableGroup` | no  |
| `commands/`   | `CommandManager`, `CommandInvocation`, matching, tails, key chords, command types, the `CommandAssistant` port   | no  |
| `forms/`      | `ParamsFormModel` (PD as form state), flat fields (PD as inline sentences)                                        | no  |
| `managers/`   | the generic managers an app composes (README §3.3): `transactions`, `history` (steps, checkpoints through a port, runs), `navigation`, `viewState`, `notifications` (today's `StatusModel`) | no  |
| `model/`      | the canvas ink (`CanvasInk`, `canvasInk`, `fixedInk`), playback contracts (`SnapshotStepperModel`, `TrajectoryControlsModel`, `HoverLabelModel`), `SettingsModel` | no  |
| `dom/`        | the JSX runtime (`Show`, `For`, `Switch`, `Dynamic`, scopes, `activity`), `MspElement`, `defineElements`, the view contract (`domView`, `elementView`) | yes |
| `react/`      | hooks and `reactView` for React and Preact hosts (optional peer)                                                  | yes |
| `styles/`     | the registry of compiled CSS modules, tokens, themes (light, dark, blue, auto; colorblind-safe and high-contrast to add) | yes |
| `components/` | dumb components: button, choices, color, context menu, legend, icon and wordmark, markdown, params form, scroll area, select, text input, tooltip, reorder, virtual list, action footer | yes |
| `shell/`      | `ShellModel` (docks, tabs, cards, sections, toolbar, pages), `ThemeModel`, `FullscreenModel`                      | no  |
| `views/`      | `ShellView`, the canvas chrome (status line, hover label, playback card, trajectory bar), the command menu, Log, Help, Settings, section lists | yes |
| `testing/`    | mocks and the story harness, for stories and tests only                                                           | yes |

- **The DOM-free folders** (`reactive/`, `commands/`, `forms/`, `managers/`, `model/`, `shell/`) are what a headless
  viewer uses; a check keeps them free of DOM imports, so no separate core package is needed (README §3.3).
- **`managers/` is new**: the generic halves of today's `ShellModel` (`navigation`: what is open; `viewState`: docks and
  sizes, serializable) and `StatusModel` move there; `history` takes its checkpoints through a port the app gives
  (the viewer's: plugin snapshot plus the sequence's state, sequence.md §7).

## 2. What it leaves to the app

- **Composition**: no modules, contributions, modes or spec; the app creates the managers and views it needs and
  registers into the shell (`registerTabs`, `registerSections`, `registerToolbar`, `registerOverlays`), each returning
  its undo.
- **The plugin**: prism-ui never imports `@molstar/plugin`, nor a mode, the assistant or the AI packages
  (`ui/src/_test/opt-in-boundary.test.ts`).
- **What a feature means**: no scene, selection or snapshots model; those contracts move to the viewer's `model/`
  (several are in `ui/src/model` only because the mocks started there).

## 3. Rules it keeps

- Snappy and quiet (packages/prism/AGENTS.md): no layout reads after writes in a frame, pointer-rate paths that update
  only what changed, `transform` and `opacity` only for animation; hidden UI inert and paused.
- Styling only through `--msp-*` tokens and the spacing scale (`styles:check`); hosts theme through the same tokens
  (README §2.20).
- Every view checked light and dark, wide and narrow, in Storybook.
