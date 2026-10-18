# prism-ui as a library

What `@molstar/prism-ui` gives an app built on it (the Prism viewer, the assistant, later Prism stories), what it leaves
to the app, and its package rules. The code is the truth for how.

## 1. What it holds

- **Two entries.** `.` (`src/index.ts`): everything, the DOM runtime, components and views included.
  `@molstar/prism-ui/core` (`src/core.ts`): the DOM-free part a headless app imports: reactive state, commands, the
  managers, the shell's model, the theme model, the canvas ink, plus types (`IconName`, `ViewRenderer`). Modules are
  also importable by path (`./*` → `src/*.ts` under `molstar-src`, `lib/*.js` built); `.tsx` modules only through `.`.
- **The DOM-free folders** (`reactive/`, `commands/`, `forms/`, `managers/`, `model/`, `shell/`) take no value import
  from the DOM folders (types are fine): `src/_test/dom-free.test.ts`. The viewer checks that its headless module graph
  reaches none (`viewer/src/_test/headless-graph.test.ts`), so no separate core package is needed.
- **`managers/`**: the generic managers an app composes: `history` (`HistoryManager`: steps and groups, checkpoints
  through the app's `HistoryPort`, undo, redo, undo turn), `navigation` (what is open), `viewState` (the layout,
  serializable), `notifications` (`NotificationsManager`: messages, tasks, chips; `announce` for screen readers only).
  The app also composes a `CommandManager` and the transaction coordinator ([view-models.md](view-models.md) §5).
- **`model/`**: DOM-free models views share: the canvas ink (`canvasInk`, `fixedInk`), playback contracts
  (`SnapshotStepperModel`, `TrajectoryControlsModel`, `HoverLabelModel`), `SettingsModel`.

## 2. What it leaves to the app

- **Composition**: no modules, contributions registry, modes or app spec. The app makes the managers and views it needs
  and registers into the shell (`registerTabs`, `registerSections`, `registerToolbar`, `registerOverlays`,
  `registerPages`), each returning its undo; what shows is the app's options (`ShellOptions.buttons`), not
  unregistering.
- **The plugin**: prism-ui never imports `@molstar/plugin`, the viewer, the stories, the assistant or the AI packages,
  nor names them, stories and tests included (`src/_test/opt-in-boundary.test.ts`). It reaches the opt-in parts only
  through ports: the command menu's `assistant`, `ShellView`'s `agent` slot, tabs and commands registered like any
  other.
- **What a feature means**: no scene, selection or snapshots model; those are the app's.

## 3. Shadow roots and keys

- **Styles adopt into shadow roots**: `adoptStyles(document | shadowRoot)` and `releaseStyles` (`styles/registry.ts`);
  every adopted root gets sheets registered later (a lazy extension's); tokens are declared on `:host` beside `:root`.
  Several viewers on a page each adopt into their own root.
- **Events across a shadow root**: `eventTarget`, `isInside`, `activeElement` (`dom/events.ts`) read the composed path,
  since a document listener sees an event from inside retargeted to the host.
- **Keys go to one viewer**: listeners check `ownsKey(node, event)`: the event is the viewer's when focus is inside it,
  or nothing on the page has focus and the pointer is over it (`viewerHost`: the shadow root's host, else the
  `.msp-root`). Outside a viewer (a story, a test) every key is.

## 4. Package rules

- **Dependencies**: `@molstar/core`, `rxjs`, `tslib`; `react >= 18` an optional peer used only by `react/`. Nothing
  imports `@molstar/plugin-ui`.
- **Folder rules**: `reactive/` imports only rxjs; `forms/` and `commands/` add `@molstar/core`'s PD; `components/` are
  dumb (props and states in, callbacks out) and never import `model/`, `shell/` or `views/`; views bind components to
  models. `IconName` may be imported as a type anywhere (icon names are data).
- **Tests and stories stay out of the package** (`_test/`, `_stories/`, `*.stories.*` excluded from `lib/` and `files`);
  `testing/` is for stories and tests only.
- **Code from elsewhere** (any-ui's reactive state, commands and hooks; Base UI's scroll area) is Mol\* code with the
  Mol\* header, which names a third-party source where there is one; no upstream to track.
- **Checked**: styles and icons (`styles:check`, `icons:check`), `typecheck:storybook` for stories and tests, every view
  in Storybook light and dark, wide and narrow.
