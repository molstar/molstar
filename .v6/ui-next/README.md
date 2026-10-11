# ui-next

A new Mol* plugin UI built on a view-model layer (rxjs, from any-ui), a VS Code–style command system, a dynamic shell,
and web components authored in TSX. It lives in `packages/plugin/ui-next` next to the current `@molstar/plugin-ui` until
it reaches parity. Status: **analysis / proposal**.

## Summary

- **Rewrite, not refactor**, side by side with the old UI. The current UI (17.9k LOC TSX + 3.7k SCSS) holds little logic
  of its own; its structure does not survive moving to view models and commands.
  ([design/architecture.md](design/architecture.md) §1)
- **any-ui is copied in as new Mol\* code** (Mol\* header, no provenance or change ledger; stories-untold is only an
  inspiration): state primitives, `ReactiveModel`, transactions, and the command manager. Its declarative form stack
  (zod/yaml) is dropped; command arguments are `ParamDefinition` params (for now).
  ([design/architecture.md](design/architecture.md) §4)
- **Two-layer state**: _sources_ turn plugin managers/events into clean `ReadonlyState`s and hide the plugin's
  observability gaps; _view models_ are UI-shaped keyed trees with granular fields and action methods. Storybook and
  tests swap sources, not the plugin. ([design/view-models.md](design/view-models.md))
- **Rendering**: a small fine-grained JSX runtime (no VDOM; states bind straight to DOM) with light-DOM custom elements
  at island granularity; a framework-neutral view contract so contributions can use React or Preact; `ui-next/react`
  hooks work with both. Gated by a spike against Preact. ([design/rendering.md](design/rendering.md))
- **Commands everywhere**: built-in, contributed, generated from registries (presets, themes, queries, state actions),
  and contextual (applicable actions for the current cell); command menu with prefixes and inline arguments; four ways a
  command runs (right away, an inline sentence, a generated form in a dialog, or custom UI that is the panel or section
  it opens); scoped keybindings; menus and buttons reference commands. ([design/commands.md](design/commands.md))
- **Agents through the commands**: commands are the tools an agent (in-app or over MCP via any-harness) calls; each
  prompt carries a state account and the person's focus; agent turns are attributed, checkpointed, and undo in one step;
  the command menu can ask the assistant. ([design/agent.md](design/agent.md))
- **New `PluginUIContext`** (ui-next's own, extending `PluginContext`) carries the UI fields: `shell`, `commands`,
  `keybindings`, `contributions`, `sources`, `model`. ([design/shell.md](design/shell.md) §6)
- **Shell** (`ctx.shell`) replaces `PluginContext.layout`, canvas first: a tabbed right panel (Home when empty,
  Structure, State tree, Settings, extension tabs); left/top/bottom docks optional and empty by default (sequence on
  top); short tasks as canvas overlays (save state, screenshot, log, help); a playback cluster (snapshot stepper and
  markdown description, trajectory, animation); hover labels at the cursor; a fading status line instead of the log
  panel. ([design/shell.md](design/shell.md))
- **PluginUISpec refresh** (a starting point, iterated as views are built): `ui.contributions` (explicit list, like
  `registry`), `ui.layout`, `ui.theme`, `ui.keybindings`. Extensions export UI contributions instead of registering
  React components from behaviors. ([design/shell.md](design/shell.md) §5–6)
- **Styling**: CSS Modules, `--msp-*` tokens, runtime themes, `@layer msp`. ([design/styling.md](design/styling.md))
- **Storybook first** (devDependency of ui-next), with three story levels: components, views over fake sources,
  integration with a real plugin; then a deep multi-agent (ultracode) review of the spec before the full build.
  ([plan.md](plan.md))

## Documents

| Document                                                 | Content                                                                 |
| -------------------------------------------------------- | ----------------------------------------------------------------------- |
| [analysis/current-ui.md](analysis/current-ui.md)         | Functional inventory of the current UI and its extension points         |
| [analysis/plugin-state.md](analysis/plugin-state.md)     | Observable plugin state, gaps, command/input surface                    |
| [analysis/stories-untold.md](analysis/stories-untold.md) | How stories-untold uses any-ui; take/adapt/leave                        |
| [analysis/extensions.md](analysis/extensions.md)         | Extension and app UI inventory with ui-next port forms                  |
| [design/architecture.md](design/architecture.md)         | Rewrite decision, package layout, dependency rules, any-ui import       |
| [design/view-models.md](design/view-models.md)           | Sources, view models, optimistic edits, actions                         |
| [design/rendering.md](design/rendering.md)               | JSX runtime, custom elements, React/Preact                              |
| [design/commands.md](design/commands.md)                 | Command definitions, registry, command menu, keybindings                |
| [design/agent.md](design/agent.md)                       | Agents: commands as tools, state and focus, checkpoints, the client API |
| [design/shell.md](design/shell.md)                       | Shell, default views, contributions, UI spec, embeddables               |
| [design/styling.md](design/styling.md)                   | CSS Modules, tokens, themes, build                                      |
| [plan.md](plan.md)                                       | Phases and acceptance                                                   |

## Open questions

1. **Renderer**: own `dom` runtime vs Preact inside elements. Proposal: own runtime, decided by the phase-2 spike.
2. **Where the command registry lives long-term**: ui-next only, or `@molstar/plugin` for headless/scripting/agent use.
3. **Core follow-ups** ([analysis/plugin-state.md](analysis/plugin-state.md) §5): registry-changed event, richer
   `cell.stateUpdated`, task list, stable `PluginCommands` ids. Do them in v6 alongside ui-next, or work around them?
4. **Deprecating `PluginLayout` and `leftPanelTabName`** in v6, or only when the old UI is retired.
5. **Package name and release**: ship `@molstar/plugin-ui-next` as `dev`/experimental in 6.0, or keep it private until
   it replaces `@molstar/plugin-ui`.
6. **Command menu UX**: modal command menu only, or also an inline command line (stories-untold style) in the status
   bar.
7. **`assist` hook** (LLM-filled command arguments, from any-ui): keep it as a type-only hook, or drop it.
8. **CSS tool**: lightningcss for deterministic CSS-module names vs esbuild's `local-css` (verify naming control).
