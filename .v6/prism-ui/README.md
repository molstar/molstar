# Prism UI

A new Mol* plugin UI built on a view-model layer (rxjs, from any-ui), a VS Code–style command system, a dynamic shell,
and web components authored in TSX. It lives in `packages/plugin/prism-ui` next to the current `@molstar/plugin-ui`
until it reaches parity. Status: **mocked in Storybook**; the designs describe what the mocks settled. Next: a renderer
bench, then the spec review ([plan.md](plan.md)). Nothing is wired to a real plugin yet.

## Summary

- **Rewrite, not refactor**, side by side with the old UI. The current UI (17.9k LOC TS/TSX + 3.7k SCSS) holds little
  logic of its own; its structure does not survive moving to view models and commands.
  ([design/architecture.md](design/architecture.md) §1)
- **any-ui is copied in as new Mol\* code** (Mol\* header, no provenance or change ledger; stories-untold is only an
  inspiration): state primitives, `ReactiveModel`, transactions, and the command manager. Its declarative form stack
  (zod/yaml) is dropped; command arguments are `ParamDefinition` params (for now).
  ([design/architecture.md](design/architecture.md) §4)
- **Two-layer state**: _sources_ turn plugin managers/events into clean `ReadonlyState`s and hide the plugin's
  observability gaps; _view models_ are UI-shaped keyed trees with granular fields and action methods. Storybook and
  tests swap sources, not the plugin. In the mocks, the view-model interfaces (`model/*.ts`) are implemented over
  in-memory data (`_stories/mocks/*`); sources come with the build. ([design/view-models.md](design/view-models.md))
- **Rendering**: a small fine-grained JSX runtime (`dom/`; no VDOM, states bind straight to DOM) with light-DOM custom
  elements at island granularity; a framework-neutral view contract so contributions can use React or Preact;
  `prism-ui/react` hooks are built for both (`useSyncExternalStore`; tested with React, Preact not yet). Every mock is
  built on the runtime; the choice is still to be confirmed by a check at scale (the state tree with 5,000 cells and
  update bursts) before the review. Preact inside the elements is the fallback.
  ([design/rendering.md](design/rendering.md))
- **Commands everywhere**: built-in, contributed, generated from registries (presets, themes, queries, state actions),
  and contextual (applicable actions for the current cell); a command menu with groups, aliases, prefixes and inline
  arguments; four ways a command runs (right away, an inline sentence, a generated form in a dialog, or custom UI: the
  panel or section it opens, `opens`); scoped keybindings; menus and buttons reference commands. Mocked: the command
  menu, inline sentences, command shortcuts (`CommandManager.handleKey`); not yet the dialog form, `opens`, keybinding
  scopes, or menus and buttons that reference command ids. ([design/commands.md](design/commands.md))
- **Agents through the commands**, opt-in: commands are the tools an agent (in-app or over MCP via any-harness) calls;
  each prompt carries a state account and the person's focus; agent turns are attributed and undo in one step (one
  checkpoint, held only until the next action); the command menu's Ask bar hands typed words to the assistant. The agent
  is a module an app includes (`views/agent`, `model/agent`); the core knows it only through ports, and a test keeps the
  core from importing it. Mocked: the Ask bar, the turn card, History, the Assistant tab; not yet the client API or the
  harness. ([design/agent.md](design/agent.md))
- **New `PluginUIContext`** (prism-ui's own, extending `PluginContext`) carries the UI fields: `shell`, `commands`,
  `keybindings`, `contributions`, `sources`, `model`. Not built yet: the mocks run on the `shell/` models alone.
  ([design/shell.md](design/shell.md) §6)
- **Shell** (`ctx.shell`, `ShellModel`) replaces `PluginContext.layout`, canvas first: a right panel floating over the
  canvas (pinnable beside it) with a rail at the right edge: Home (the Mol\* mark), the tabs (Structure, Selection,
  Snapshots, State tree, extension tabs), the scene's and the app's overlay cards (Animation, Canvas settings, Export
  with contributed sections, Remote states when on), and at its foot Log (History, and the Assistant, with the agent),
  Settings, the theme and the command menu; left/top/bottom docks optional (the sequence floats on top); a canvas stack
  (mode, camera, illumination & XR, expand or fullscreen, sequence); playback top left (snapshots card, trajectory bar);
  hover labels at the cursor; a status column of fading messages and tasks bottom left, selection and focus chips bottom
  right; a busy overlay for tasks that hold the plugin; compact mode as a bottom sheet.
  ([design/shell.md](design/shell.md))
- **PluginUISpec refresh** (a starting point, iterated as views are built): `ui.contributions` (explicit list, like
  `registry`), `ui.layout`, `ui.theme`, `ui.keybindings`. Extensions export UI contributions instead of registering
  React components from behaviors. ([design/shell.md](design/shell.md) §5–6)
- **Styling**: CSS Modules compiled by lightningcss into generated TypeScript modules (`_generated/*.css.ts`), `--msp-*`
  tokens, runtime themes (light, dark, blue, auto), `@layer msp`; flat surfaces, few lines, glass over the canvas that
  follows the canvas background; every scrolling area a `ScrollArea`. ([design/styling.md](design/styling.md))
- **Storybook first** (devDependency of prism-ui): Overview, Shell, Views, Canvas, Components, Internals, with `play`
  interaction tests. Built so far: components and views over mocked view models; integration with a real plugin comes
  with the build. Then the renderer check and a deep multi-agent (ultracode) review of the spec before the full build.
  ([plan.md](plan.md))

## Documents

| Document                                                 | Content                                                                 |
| -------------------------------------------------------- | ----------------------------------------------------------------------- |
| [analysis/current-ui.md](analysis/current-ui.md)         | Functional inventory of the current UI and its extension points         |
| [analysis/plugin-state.md](analysis/plugin-state.md)     | Observable plugin state, gaps, command/input surface                    |
| [analysis/stories-untold.md](analysis/stories-untold.md) | How stories-untold uses any-ui; take/adapt/leave                        |
| [analysis/extensions.md](analysis/extensions.md)         | Extension and app UI inventory with prism-ui port forms                 |
| [design/architecture.md](design/architecture.md)         | Rewrite decision, package layout, dependency rules, any-ui import       |
| [design/view-models.md](design/view-models.md)           | Sources, view models, optimistic edits, actions                         |
| [design/rendering.md](design/rendering.md)               | JSX runtime, custom elements, React/Preact                              |
| [design/commands.md](design/commands.md)                 | Command definitions, registry, command menu, keybindings                |
| [design/agent.md](design/agent.md)                       | Agents: commands as tools, state and focus, checkpoints, the client API |
| [design/shell.md](design/shell.md)                       | Shell, default views, contributions, UI spec, embeddables               |
| [design/styling.md](design/styling.md)                   | CSS Modules, tokens, themes, build                                      |
| [plan.md](plan.md)                                       | Phases and acceptance                                                   |

## Open questions

1. **Renderer**: decided, the own `dom` runtime ([spike.md](spike.md)): time slicing, per-node records and virtualized
   long lists carry the state tree at 5,000 cells.
2. **Where the command registry lives long-term**: prism-ui only, or `@molstar/plugin` for headless/scripting/agent use.
3. **Core follow-ups** ([analysis/plugin-state.md](analysis/plugin-state.md) §5): registry-changed event, richer
   `cell.stateUpdated`, task list, stable `PluginCommands` ids. Do them in v6 alongside prism-ui, or work around them?
   (`display.triggerLabel` on state actions and transformers is done.)
4. **Deprecating `PluginLayout` and `leftPanelTabName`** in v6, or only when the old UI is retired.
5. **Package name and release**: ship `@molstar/plugin-prism-ui` as `dev`/experimental in 6.0, or keep it private until
   it replaces `@molstar/plugin-ui`.

Settled by the mocks:

- **Command menu UX**: the modal command menu only (⌘K / ⇧⌘P, the rail's button; Preset and Styles open it in their
  groups), with the assistant as a bar under its list; no inline command line in the status area.
- **`assist` hook**: kept. `CommandAssistant.propose` (`commands/types.ts`), given to `CommandManager` as `assistant`,
  reads words into an open command's values (`CommandInvocation.assist`, on Tab; jest-tested, no story wires it). The
  command menu's Ask bar has its own port, `CommandMenuAssistant.propose` (`views/command-menu`), which picks one
  command and fills its sentence from the words. Either way it proposes; nothing runs until Enter.
- **CSS tool**: lightningcss (`scripts/css.mjs`), which names module classes deterministically
  (`msp-<component>-<class>`).
