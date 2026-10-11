# Commands, command menu, and keybindings

A VS Code–style command layer: every user-facing operation is a named, searchable, context-aware command; the command
menu, keybindings, menus, toolbar buttons, and scripts all go through it. This is what lets the new UI surface far less
nested UI: what used to be a card with a form becomes a command with arguments.

## 1. Relationship to `PluginCommands`

`PluginCommands` stay as they are: an internal, typed message bus with UUID ids that behaviors subscribe to, and the
classic UI's. Prism's commands are a separate, user-facing layer: string ids, titles, categories, arguments, enablement,
keybindings. A command's `run` calls actions ([view-models.md](view-models.md) §6), never a `PluginCommand` (the plugin
boundary, [contexts.md](contexts.md) §7).

The command core (`commands/`: `CommandManager<Ctx>`, `CommandInvocation<Ctx>`, the types, tail reading, matching, key
chords) is framework-free and generic over its context; it moves to `@molstar/prism` with `PrismContext`, which owns the
one registry (contexts.md §4), so commands run headless (tests, the agent's tools, MCP). Its icons become a core icon
key that prism-ui maps to glyphs (review M4). The command menu (`views/command-menu`) is the only view over it so far.

The skeleton's command contract is in [contexts.md](contexts.md) §4–5: the composed command context, `execute` as the
one run path, results with their effect, the metadata that makes a command a tool, origins and history. This document
covers definitions, the menu, keybindings and the built-in map.

## 2. Definition

Adapted from any-ui's `CommandDefinition<Ctx>`; the main change is that arguments are `ParamDefinition` params instead
of declarative definitions. As built (`commands/types.ts`):

```ts
interface CommandDefinition<Ctx> {
  id: string;                         // 'data.download-structure', 'select.pymol'
  title: string;                      // 'Download structure…'
  category: string;                   // 'Data', 'Select', 'Selection', 'View', 'Export', …
  icon?: IconName;                    // before the title in the command menu and on buttons that run it
  group?: string;                     // the group it sits in, by id ('select.expression'); top level without one
  alias?: string | readonly string[]; // other names it is found by, and words that start it with arguments: 'color'
  shortcut?: KeyChord;                // { key, mod?, alt?, shift? }; display derived per platform
  when?(ctx: Ctx): true | string;     // a reason string: listed, greyed, explained
  params?(ctx: Ctx, values: CommandValues): PD.Params;  // the arguments; without them, picking runs the command
  fields?: Record<string, FieldHint>; // how params read inline, by name: { before?: 'from', format?: 'color' }
  action?: string;                    // the run button's label: 'Download' (default 'Run')
  reshapes?: readonly string[];       // params whose value changes the params (built again, keeping what fits)
  presets?(ctx: Ctx): readonly CommandPreset[];   // search-only ways in: 'from AlphaFold DB' with its values
  defaults?(ctx: Ctx): CommandValues; // initial values (e.g. from the selection)
  preview?(args: CommandArgs, ctx: Ctx): unknown; // handed to the PreviewSink while the sentence is complete
  check?(args: CommandArgs, ctx: Ctx): string | undefined;  // a last check across fields; a message blocks the run
  run(args: CommandArgs, ctx: Ctx, signal: AbortSignal): CommandResult | Promise<CommandResult>;
  assist?: { instruction: string; answers?: boolean };  // an assistant may fill it in from words (§5)
}

/** A node of the command tree for browsing (Select › Expression); groups nest by `parent`. */
interface CommandGroup {
  id: string;                         // 'select', 'select.query', 'select.expression'
  title: string;                      // 'Select', 'By query', 'Expression'
  icon?: IconName;
  parent?: string;                    // top level without one
  order?: number;                     // among its siblings, lower first
}
```

`CommandValues` are the inline values, primitive as inline fields hold them (`string | number | boolean`); `CommandArgs`
are what `run`, `check` and `preview` get: the inline values plus the params' PD defaults for everything else (`argsFor`
in `forms/flat-fields`), so params that are not inline run at their defaults. A command without `params` runs with its
`defaults` as they are.

Planned, not built: `opens?: CommandTarget` for custom UI (§2.1; the mock's panel commands call the shell in `run`);
`order` on commands; several chords per command; and on groups, `when`, `shared` arguments (one value for every command
in the group, shown once in its header: §10's selection mode) and `children` computed from context or registries
(entities, chains, ligands of the loaded structures).

```ts
/** Where a command with custom UI takes the user: a tab or an overlay panel, optionally at one of its sections. */
type CommandTarget =
  | { tab: string; section?: string }        // { tab: 'selection', section: 'measure' }
  | { overlay: string; section?: string };   // { overlay: 'export', section: 'animation' }
```

- **Arguments are PD params** (for now; the argument schema can be revisited). `formFields` (`forms/flat-fields`) reads
  the params that are inline: text, selects with string options, booleans (on/off) and numbers (bounded ones as
  `slider`); each with the quiet word before it and its format from `fields`. Everything else (colors, groups, mapped,
  object lists…) is `formOnlyParams`: edited in the full form only, and run at its default from the inline sentence.
  Hidden params never show. Existing transformer and action params, extension params (`Mp4AnimationParams`,
  `GeometryParams`, …) become command arguments without rewriting.
- **Formats** are named argument types (`CommandFormat<Ctx>`, looked up through the `CommandFormats` port by name):
  `hint(ctx)` says in a few words what can be typed, `check(text, ctx)` gives `{ ok, value }` or a message (with a hint,
  and `unreadable` when the text is words rather than a broken value), `complete(text, ctx)` gives completions; a
  completion with `fill` (an example expression) fills the text for editing instead of choosing it. The mocks have
  `color` (names or hex), `residues` ("A 10-20, B 5"), `ligand` (the ligands of 1TQN) and `expression.<language>` (an
  example list and a light syntax check standing in for the transpilers). Planned built-ins: `structure` (pick from the
  scene), `component`, `cell`, `selection` (residue list or MolScript, compiled with the existing
  `compileIdListSelection`/script support), `color` (names, hex, palette), `query` (registered selection queries),
  `snapshot`, `file`, `url`, `pdb-id`.
- **Custom UI is a panel or a section, not something inside the command menu** (§2.1). A command whose arguments are not
  a form (`opens`) takes the user to the panel or section that already holds that UI, focused on it (the other sections
  closed), so the command still has one id, appears in the command menu, takes keybindings and can sit in a menu; the
  panel's own actions report outcomes like any command. The command menu itself only ever shows an inline sentence or
  the generated form in a dialog.
- **Results** (`CommandResult`: `message`, one line for the person; `detail`; `ops`, the equivalent operation in mono;
  `undoable`; `data` for the host) reach the host through the manager's `outcomes`. In the mocks the message goes to the
  status line and, with the agent module, the run to History ([agent.md](agent.md) §5); a failure or a refused
  keybinding is a warning there. Not yet: `undoable` offering Undo through `state.data.undo()`.
- **Groups** nest commands where there are too many to list flat (presets, quick styles, selection queries, expression
  languages). Nesting is for browsing only: a search from the top level still matches every command at any depth and
  shows its path ("Select › Expression"), so nothing is harder to find. In the command menu, Enter or → opens a group,
  Backspace or ← on an empty input goes up, and the head shows the path as a breadcrumb (each part clickable). The
  command menu can open inside a group (its `group` prop), which is how the Structure tab's Preset and Styles buttons
  open theirs, teaching the command menu as they are used.

### 2.1 How a command runs

Every command runs one of four ways, decided by what it needs from the user:

| How it runs                    | When                                                                                                    | Commands                                                                                                                                                                                                                                                                           |
| ------------------------------ | ------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Right away**                 | Nothing to ask (no `params`)                                                                            | camera reset / orient, toggles (sequence, illumination, expand), panels, presets, quick styles, next / previous snapshot, the fixed selection queries (ligands, water, helices…), clear selection, wiggle, undo                                                                    |
| **Inline sentence**            | A few flat arguments: text, select, number, boolean, and formatted text                                 | download structure / volume, color selection, transparency, focus on, select residues or by an expression, set isovalue, model index, rename, add representation, color by, copy / download screenshot                                                                             |
| **Generated form in a dialog** | Rich arguments that are still a plain form (group, mapped, object list…): "open as form", or by default | component options, measurement options, add component, a hierarchy preset with its params, extension actions (Load G3D, tunnels, assembly symmetry, volume streaming)                                                                                                              |
| **Custom UI** (`opens`)        | Several steps, picking in the 3D view, or a preview                                                     | open files (a list, a format per file, as the original), Zenodo import (record, then a file list), image export (preview, size), animation / MP4 export (render, then save, with progress), measure (picking atoms), superpose (structures and chains), PAE plot, keybindings help |

The custom UIs are the panels and sections that hold them: `Measure…` opens the selection tab at its `measure` section,
`Superpose…` at `superpose`, `Export animation…` the export overlay at `animation`, `Open files…` the home tab's open
section. A command that needs the user to act in the 3D view (pick atoms) is always custom UI, never a form.

**State (mocks):** right away and the inline sentence are built (`CommandManager.open` runs a command without `params`
at once, with its `defaults`, and opens a `CommandInvocation` for one with them). A color is inline as a text field with
the `color` format, since `PD.Color` is form-only. Not yet: the generated form in a dialog (the inline sentence notes "+
N more options at their defaults" and runs with them), and `opens`: the mock's custom-UI commands (`measure.open`,
`superpose.open`, `export.image`, `export.animation`, `export.geometry`, `canvas.background`) call the shell in `run`
(`shell.show(tab, section)`, `shell.openOverlay(overlay, …, section)`) and report "Showing …" or "Opened …".

`CommandContext` (planned; the mocks run with an empty context, and `when` reads closures such as the structure count)
is a `ReadonlyState` composed of a frame base and slices that contributions add (contexts.md §4):

```ts
interface CommandContext {
  busy: boolean;
  mode: string;
  activeView?: string;                     // shell focus, for scoped keybindings (UI only)
  slice<S>(key: CommandContextKey<S>): S | undefined;
}

// The viewer's slice (in @molstar/prism-viewer): what its commands' `when` and `params` read.
interface SceneCommandContext {
  structures: readonly SceneNode[];        // the structures in the hierarchy selection
  target?: SceneNode;                      // the node acted on (context menu, inspector)
  selection: { elementCount: number; structureCount: number };
  focus?: { label: string };
  selectionMode: boolean;
}
```

## 3. Where commands come from

| Source                    | Examples                                                                                                                                                                                                                              | Lifetime                                                                |
| ------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------- |
| Built-in                  | Load file/URL/PDB/AlphaFold, undo, reset camera, screenshot, toggle selection mode, clear selection, add measurement, superpose, quick styles, snapshot add/next/previous/play, toggle view, change theme                             | Registered by `PrismContext` (the frame's and the mode's contributions) |
| Contributed               | Extension contributions (geo-export, model-export, mp4 dialog, zenodo, assembly symmetry, …)                                                                                                                                          | While the contribution is active                                        |
| Generated from registries | "Apply preset: …" (hierarchy and representation presets), "Color by: …" (themes), "Select: …" (queries), "Add representation: …", state actions ("Download structure", "Load G3D", "Download tunnels"), animations, background styles | Recomputed on `registry.version`                                        |
| Contextual                | `actions.fromCell(currentCell)` as "Apply to <cell>: …"; node actions in scene/tree context menus                                                                                                                                     | While the cell/node is current                                          |
| View-scoped               | Tree navigation, command menu navigation, dialog confirm                                                                                                                                                                              | While the view is mounted                                               |

Generated commands are the largest win: today these choices are buried in nested menus inside cards. The command menu
lists them flat, fuzzy-searchable, with `when` explaining why one is unavailable ("Select a structure first").

**State (mocks, `testing/commands.ts`):** `mockCommands` registers a hand-written built-in set against the mock shell
(§8 lists its ids); representation presets, quick styles, selection queries and expression languages are listed by hand
where they would be generated. Not yet: generation from registries, contextual and view-scoped commands, and
contributions registering their own.

## 4. Registry

`CommandManager<Ctx>` (`commands/manager.ts`) is any-ui's `CommandManager`, adapted. Its options are the `context` (a
`ReadonlyState<Ctx>`), the `formats` port, and optionally a `PreviewSink`, a `CommandAssistant`, a `RecentsStore`
(`maxRecents`, default 5) and `previewDeps` (what in the context an open preview reads, so it is redrawn only when that
changes).

- `register(defs)` returns an unregister function; contributions and views register by mount scope. Unregistering
  cancels a run of theirs and closes their open form, and an outcome that settles after that reports to no one. A
  duplicate id throws unless the newcomer declares `replaces` (a shadow stack: unregistering restores the previous one;
  contexts.md §6). Today the mock replaces silently.
- `registerGroups(groups)` likewise; `groups`, `groupPath(id)`, `children(id)` (a group's subgroups by `order` and its
  commands) and `inGroup(command, id)` serve browsing.
- `query(text)` gives `recent` (on an empty text), `matches`, and `direct`: the text starts with a command's alias and
  goes on, so that command opens with the rest as its tail (§5). Matching is a fuzzy subsequence score (`fuzzyScore` in
  `commands/match.ts`) over the title and the aliases (a little lower); equal matches keep their category together, in
  the order categories were registered, then go by title. A preset is listed only when the text says more about it than
  the command alone (its title, its keywords, or both titles). Prefixes are the command menu's, not the manager's (§5).
- `open(id, tail?, preset?)`: a command with `params` opens a `CommandInvocation` (the inline sentence; `active` holds
  it); one without runs at once. Refused with a reason (`{ kind: 'blocked' }`, never a throw) while a run is in flight
  and when `when` says no. `close()`, `cancel()` (an explicit act: the run's signal aborts), `running`.
- `handleKey(event)` runs a command's chord (§6).
- `outcomes`: a non-replaying event stream (`done`, `failed`, `cancelled`, and `blocked` for a refused keybinding, which
  has no caller to tell).
- Recents are remembered on `done` and saved through the injected store (localStorage in apps, none in embeds); the
  mocks keep them in memory.
- Not yet: `execute(id, args, { source, signal })` as the one run path (contexts.md §4: writes queued, reads alongside,
  arguments validated, a `RunHandle`; today one run at a time); contextual commands ranked first.

`CommandInvocation` (`commands/invocation.ts`) is one open command: its state is `editing`, `resolving` (the assistant
is reading words), `running`, `cancelling`, `done` or `closed`; a failed run goes back to `editing` with the error, and
the form stays open (`InvocationStatus` also has `failed`, never set). Its value store is the flat inline values
(`keptValues`, `readValue` in `forms/flat-fields`), not a `ParamsFormModel`; `args()` adds the PD defaults. It starts
from `defaults`, then a preset's values, then the typed tail (typed words win over defaults); a field in `reshapes`
builds the params again when set or cleared, keeping the values that still fit. The preview is drawn once every required
field has a value and `when` allows the run, with the highlighted completion in place of the active field's value, so it
shows what Enter would run; it stays while the run runs and clears when it is done or the form closes. Before a run,
`when`, each value's format and the command's `check` are asked again.

## 5. Command menu

`CommandMenu` (`views/command-menu`) over a `CommandManager`: props `manager`, `onClose`, and optionally `initial` text,
`group` (open inside it), `prefixes` and `assistant`.

- Opens with ⌘K / Ctrl+K and ⇧⌘P / Ctrl+Shift+P (both; the first is common on the web, the second familiar from VS
  Code), from the rail's command button (its tooltip shows ⌘K; for touch), from Home, and inside a group from the
  Structure tab's Preset and Styles buttons. In the mocks the host (the shell story) listens for the two chords; they
  are not commands.
- Modal over the viewport (a scrim; a click outside closes it); the canvas stays live (previews draw in 3D).
- **Browsing** (empty input): at the top level, Recent (the recently used commands), then Browse (the top-level groups),
  then the ungrouped commands under their category captions; inside a group, its subgroups, then its commands, under
  category captions when there is more than one. **Searching**: every command at any depth (inside a group, only its
  subtree), with a direct alias match first, then the groups whose title matches, then the commands, best first.
- Rows: icon, title with the matched letters marked, then at the right the reason it cannot run, or its group path below
  the group being browsed, or (at the top level) its category; its shortcut (`formatChord`). A group row ends in a
  chevron. A preset is its own row ("Download structure from AlphaFold DB"). Picking a blocked row shows the reason in a
  line above the list.
- **Prefixes** turn the menu into another list. The host supplies them as `CommandMenuPrefix` (`char`, `title`,
  `placeholder`, `items(text)` returning `CommandMenuItem`s whose `run` may return true to keep the menu open); `?` is
  built in: the prefixes (picking one types it), then every command with a shortcut, with its chord. The mocks
  (`mockPrefixes`): `@` goes to a scene node (structure, component, volume…: shows the Structure tab and picks it
  there); `#` goes to a snapshot by title; `:` selects by residues (`:A 10-20`) or by an expression after a language's
  name (`:pymol resn HEM`), opening that command's sentence with the rest read into it.
- Commands carry the **other names people would type** as aliases, searched like their titles (a little lower):
  Superpose: align, superimpose, RMSD, TM-align; Export image: screenshot, PNG; Clear selection: deselect; Download
  structure: fetch, load. Aliases of a command without arguments are only searched.
- **Aliases with arguments.** Typing an alias of a command with `params` followed by more text opens that command with
  the rest read into its arguments ("fetch 1tqn from alphafold", "color red on cartoon", "pymol resn HEM"); the longest
  alias that fits wins, and among commands sharing it, one that can run. Text that is still the start of a title or an
  alias ("focus on" of "Focus on…") stays a search. The row shows the tail after the title.
- **Reading the tail** (`readTail`, `commands/tail.ts`), in three passes: the words after a field's before-word, up to
  the next one, are its value when they read as an option or in its format, or as they are for a plain text field
  ("onto", "by", "on"); then words (up to four at a time) that name an option, or a value its format lists ("sticks",
  "blue", "by element"), fill the first open field that takes them; then the words left, when they are one run, may read
  as a still-open field's value in its format. What remains is typed into the first empty field, and taken when that
  field's format reads it.
- **The inline sentence**: the command's title, then each inline field as a quiet before-word and a value token (its
  label until set; a required one still empty is marked), then the run button with the command's `action` label and ↵.
  The active field is typed into the input (its label as placeholder); under it, completions from the field's options or
  format, an exact name first, and text that already reads as a value first "as typed". ↑↓ move through them, Tab takes
  one (or moves to the next field, Shift+Tab back), Enter takes what is typed and runs when that leaves no field empty
  (else moves on to the next empty one), and with nothing typed runs once no required field is missing (else goes
  there), Backspace on an empty input clears the field's value, then goes back a field, then to the list; Esc goes back
  to the list, and Esc there closes. A number takes "50%" as 0.5. A field's error shows under the sentence with its
  hint. Params that are not inline show as "+ N more options at their defaults".
- A command that ran closes the command menu (its message goes to the status line); one that fails keeps its sentence
  open with the error. Commands' own shortcuts run them while the menu is closed.
- **The assistant** (`assistant?: CommandMenuAssistant`, only when the app includes the agent module;
  [agent.md](agent.md) §5). The port: `ready()`, `propose(text)` giving `{ command, tail? }` when the words fit one
  command, `ask(text)` handing them to the agent as a turn of its own, and `setup()`. Under the list, always there, an
  _Ask the assistant_ bar shows what is typed and takes it with ⌘↵ (Ctrl+↵), or with ↵ when no command matches (it is
  then the default, and the empty list says "No command matches: ask the assistant."). A proposed command opens as its
  ordinary sentence, filled in from the tail, under a line saying so ("Proposed by the assistant from “…”: check it,
  then Enter to run"); nothing runs until Enter. Words that fit no one command go to `ask`. Without a model the bar is
  quiet, reads _Set up an AI assistant_, and `setup()` (the Assistant tab) closes the menu.
- **Assist inside a field** is in the core, not wired in the mocks: with a `CommandAssistant` on the manager and
  `assist` on the command, words a field cannot read as a value (no completion fits, or its format says `unreadable`)
  set `canAssist`, and Tab asks `CommandAssistant.propose`; it gives values (marked `proposed` until touched), an answer
  (shown as the result), or `escalate` (the words are more than this command).
- Footer: the keys for the mode (↑↓ move, ↵ run, ⌘↵ ask while the assistant is ready, → open group, ? help, esc close;
  in a sentence tab, ↵, ⌫, esc).
- ARIA: the input is a combobox and keeps focus; rows are `role="option"` in a `listbox`.

**State (mocks, `views/command-menu/command-menu.tsx`, Shell › Viewer mode › Command Menu Open):** all of the above. Not
yet: "open as form" and the dialog for params that are not flat; `opens` as a declared field; color as an inline
argument of its own; previews while browsing (no `PreviewSink` in the mocks); counts on demand (§10.1).

## 6. Keybindings

A `KeybindingService` on the plugin root element:

- **Scopes**: `global` (anywhere in the plugin), `viewport` (canvas focused or hovered), `view:<id>` (a focused view),
  `command menu`, `dialog`. The most specific scope wins; text-entry targets only receive ⌘/Ctrl chords (any-ui's rules
  for IME, repeat, and `e.code` fallback are kept).
- **Coexistence with canvas bindings.** Trackball keys (W/A/S/D/R/F/Q/E…) and behavior bindings (`I`, `O`, `T`, `G`,
  Ctrl+arrows) keep working through `InputObserver`. The service reads those bindings and reports conflicts in the Help
  view instead of overriding them. Later, behavior key bindings can be re-expressed as commands with `viewport` scope.
- **Platforms.** Chords are written once with `mod` (⌘ on macOS, Ctrl on Windows and Linux); `formatChord` in
  `commands/keys.ts` shows them per platform (`⌘K` / `Ctrl+K`, `⇧⌘P` / `Ctrl+Shift+P`).
- **Defaults** are part of command definitions; the UI spec can override them (`keybindings: [{ command, key }]`). User
  remapping is out of scope for the first version.
- Help is generated from the registry, behavior bindings, and trackball bindings, so it is always complete.

**State (mocks):** no `KeybindingService` and no scopes yet. `CommandManager.handleKey(event, mac)` runs a command's
`shortcut` from a keydown the host forwards while the command menu is closed: repeats, IME composition and handled
events are ignored; in a text-entry target only ⌘/Ctrl chords without Alt run (Alt types characters on macOS; Ctrl+Alt
is AltGr); the key is matched first, then its `code` when the key is not a plain letter or digit (⌥T on macOS is "†"),
so AZERTY's ⌘A does not run a ⌘Q command; among commands with one chord, one that can run wins, and a refused one is a
`blocked` outcome. `chordKeys` gives the keycaps (⌥⇧⌘ order on macOS, Ctrl+Alt+Shift elsewhere) and `formatChord` the
text. The mocks bind ⌥I (illumination) and ⌘S (snapshots); help is the command menu's `?` list (the prefixes and every
command with a shortcut), not yet a Help view with behavior and trackball bindings.

## 7. Menus and buttons are commands

A command can also just open UI: a panel, a tab, or a panel at a section (`Export animation…` opens the Export panel at
its `animation` section). Features that live in a panel are thus still found by the command menu, and extensions add
sections and commands rather than buttons.

Context menus (scene node, state-tree cell, viewport right-click), toolbar buttons, and view headers reference command
ids, optionally with fixed arguments:

```ts
{ command: 'structure.component.add', args: { structure: node.ref }, group: 'add' }
```

A button is enabled per the command's `when`, and its tooltip shows the title and keybinding. This keeps every action
discoverable in the command menu and avoids duplicate logic.

**State (mocks):** not yet. Toolbar items, the selection bar's actions (`SelectionModel.actions`) and the panels' own
buttons run closures that do what the matching mock command does; no menu entry references a command id. The one link so
far is the other way: the Structure tab's Preset and Styles buttons open the command menu inside the `preset` and
`styles` groups (`StructureTabProps.openCommands`).

## 8. Built-in command map (first pass)

The planned ids, and the ids the mock registers (`mockCommands`). They differ where the mock settled a different shape:
one command per preset or quick style in a group instead of one command with the preset as its argument, an `export`
group for the Export panel's sections, one command per panel (`panel.<tab>`), and a single `measure.open` /
`superpose.open` that opens the section holding the custom UI (§2.1). Which ids stay is for the spec review.

| Today                              | Commands                                                                                                                                                                                                                                           | In the mocks                                                                                                                                                                                                                                                                                      |
| ---------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Home tab actions                   | `data.open-files`, `data.download-structure` (presets: PDB, AlphaFold DB, PDB-IHM, URL), `data.download-volume`, generated state actions                                                                                                           | `data.download-structure` (ID, `from` PDB / AlphaFold DB / PDB-IHM; presets `alphafold`, `pdb-ihm`)                                                                                                                                                                                               |
| Source card                        | `structure.apply-hierarchy-preset`, `structure.set-model-index`, `structure.set-type`, `structure.select-in-scene`                                                                                                                                 | —                                                                                                                                                                                                                                                                                                 |
| Components card                    | `structure.apply-preset`, `structure.component.add`, `.add-representation`, `.set-color`, `.remove`, `.toggle-visibility`, `.rename`, `structure.component-options` (dialog)                                                                       | `preset.<id>` in the `preset` group, one per representation preset (`empty`, `auto`, `atomic-detail`, `polymer-cartoon`, `polymer-and-ligand`, …), right away                                                                                                                                     |
| Quick styles                       | `view.style.default`, `view.style.illustrative`, `structure.preset.cartoon` / `spacefill` / `surface`                                                                                                                                              | `styles.<id>` in the `styles` group: `repr-default`, `cartoon`, `spacefill`, `surface`, `style-default`, `illustrative`                                                                                                                                                                           |
| Selection toolbar                  | The `select` and `selection` groups, §10                                                                                                                                                                                                           | `select.residues`, `select.pymol` / `.vmd` / `.jmol` (`select.expression`), `select.ligands` / `.helices` / `.water` (`select.query`); `selection.color`, `selection.transparency`, `selection.clear`                                                                                             |
| Measurements                       | `measure.distance`, `.angle`, `.dihedral`, `.label`, `.orientation`, `.plane`, `measure.options`                                                                                                                                                   | `measure.open` (Selection tab at `measure`), `measure.distance` (blocked: "Select two atoms first")                                                                                                                                                                                               |
| Superposition                      | `superpose.chains`, `.tm-align`, `.ligands`, `.atoms`, `.db`                                                                                                                                                                                       | `superpose.open` (Selection tab at `superpose`; blocked below two structures)                                                                                                                                                                                                                     |
| Procedural animation               | `animation.wiggle.dynamics`, `.uncertainty`, `.clear`                                                                                                                                                                                              | `animation.wiggle` (by dynamics), `animation.open` (the Animation overlay)                                                                                                                                                                                                                        |
| Volume / particles                 | `volume.add-representation`, `volume.set-isovalue`, `volume.streaming.enable`, `particles.*`                                                                                                                                                       | —                                                                                                                                                                                                                                                                                                 |
| Viewport controls                  | `camera.reset`, `.reset-zoom`, `.orient-axes`, `.reset-axes`, `view.toggle-panels`, `view.expand`, `view.fullscreen`, `view.illumination.toggle`, `view.xr.toggle`, `view.settings`, `screenshot.copy`, `screenshot.download`, `screenshot.dialog` | `camera.reset`, `.reset-zoom`, `.orient-axes` (`camera` group), `view.expand`, `view.illumination` (⌥I), `view.xr` (blocked), `view.sequence`, `canvas.settings`, `canvas.background`; `export.image`, `export.animation`, `export.geometry` (`export` group, the Export overlay at that section) |
| Snapshots / animation / trajectory | `snapshot.add`, `.next`, `.previous`, `.play`, `.save` (molj/molx), `.open`, `animation.play`/`.stop`, `trajectory.next`/`.previous`/`.first`                                                                                                      | `snapshots.show` (⌘S, the Snapshots tab), `snapshot.next`, `snapshot.previous`                                                                                                                                                                                                                    |
| State tree                         | `state.undo`, `state.remove`, `state.toggle-visibility`, `state.set-current`, contextual actions                                                                                                                                                   | —                                                                                                                                                                                                                                                                                                 |
| Focus                              | `focus.target` (ligands etc. as completions), `focus.clear`, `focus.recenter`                                                                                                                                                                      | `focus.target` (the `ligand` format)                                                                                                                                                                                                                                                              |
| Shell                              | `view.show:<view-id>`, `view.theme`, `command-menu.open`, `help.keybindings`                                                                                                                                                                       | `panel.structure`, `panel.home`, `panel.log`, `panel.settings`                                                                                                                                                                                                                                    |

## 9. Extensibility

Extensions extend the command layer in four ways, all through their `UIContribution` ([shell.md](shell.md) §5) or the
plugin registry, never by patching built-in definitions:

| Mechanism                    | What an extension adds                                                                                                         | Example                                                                                                                         |
| ---------------------------- | ------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------- |
| `commands`                   | New commands, with `group` placing them anywhere in the tree, including built-in groups                                        | `geo-export.export` under `export`; `assembly-symmetry.enable` under `structure`                                                |
| `commandGroups`              | New groups, top level or nested under built-in ones                                                                            | `volseg` group with its segment commands; `select.validation` under `select.query`                                              |
| Registry entries             | Nothing UI-specific: registered selection queries, presets, themes, state actions, and animations appear as generated commands | rcsb's `hasClash` query under Select › By query › its category; dnatco presets under Apply preset; `DownloadTunnels` under Data |
| `formats`, `menus`, sections | Argument types, context-menu and toolbar entries that reference commands, and the sections custom UI lives in (§2.1)           | A `segment` format with completions; a toolbar button for the PAE plot; the MP4 export flow as an export section and `opens`    |

- **Registration follows activation.** A contribution's commands and groups are registered when it is activated and
  unregistered when it is deactivated (a run of theirs is cancelled, their open form closed), the same rule as views.
- **Ids are namespaced** by the contribution (`<contribution>.<name>`); a duplicate id is an error at activation.
  Placing a command in another contribution's group by id is allowed and is how built-in groups are extended.
- **Apps decide what is visible.** The UI spec can hide commands or whole groups and rebind keys
  (`ui.commands: { hide: ['superpose.*'], keybindings: […] }`), so an app shows a subset without forking definitions.
- **Script languages** follow the plugin: importing a transpiler module (`@molstar/model/script/transpilers/pymol`)
  makes its language available (`Script.getAvailableLanguages()`), and the Expression group lists one command per
  available language. A custom language registered with `registerTranspiler` appears the same way.

**State (mocks):** the manager has what this needs at its level: `register` and `registerGroups` each return the
function that takes them away again (for commands, cancelling a run of theirs and closing their form), groups nest by
`parent` id, a command joins any group by id, and formats are looked up by name through the `CommandFormats` port. Not
yet: the `UIContribution` fields, the duplicate-id check (the manager replaces a command registered again under the same
id), hiding by the UI spec, and languages from the transpiler registry (the mock lists PyMOL, VMD and Jmol by hand).

## 10. Selection commands

Today's selection toolbar (selection mode only) has four mode buttons (union, subtract, intersect, set) that each open
the same nested query menu, a residue-list helper, a Theme panel whose form mixes target, action, and representation,
plus create/remove component, undo, help, and close. It is already command-like but each piece is its own UI. In
prism-ui it becomes two command groups, used from the command menu, from context menus, and (`selection`) from the
selection bar's actions (§10.3).

### 10.1 `select`: changing what is selected

One group with a **shared `mode` argument** (`add` | `remove` | `intersect` | `set`, as Mol*'s selection manager names
them; default `add`) shown once as a segmented control in the group header, instead of four buttons each opening a copy
of the menu. Typed tails accept the words too ("select helix remove"). Modifier keys on Enter choose the mode for one
run: Enter adds, Shift+Enter removes, Alt+Enter intersects, ⌘/Ctrl+Enter replaces.

```
Select ›                                  shared: mode = add | remove | intersect | set
  All · Clear · Invert
  By query ›                              generated from the selection-query registry, grouped by `category`
    Structure ›            Polymer, Protein, Nucleic, Ligand, Water, Ion, Lipid, Carbohydrate, …
    Structure property ›   Trace, Backbone, Sidechain, Sidechain with trace, Helix, Beta strand/sheet, …
    Bond property ›        …
    <extension categories> e.g. Validation › Clashes; Confidence › pLDDT > 70; Membrane › Transmembrane
  In the loaded structures ›              computed per structure from the scene (dynamic children)
    Entities ›             polymer and carbohydrate entities (today's "Polymer/Carbohydrate Entities")
    Chains ›               (new)
    Ligands ›              ligands and non-standard residues
    Elements ›             element symbols present
  Residues…                               form: "A 10-20, B 5" with id type (auth/label), via compileIdListSelection
  Expression ›                            one command per available script language, text argument with syntax check
    PyMOL…  ·  VMD…  ·  Jmol…  ·  MolScript…
  Around the selection…                   (new) radius in Å, whole residues on/off
  Expand to ›                             (new) residues · chains · entities · structures
```

- **Preview while browsing.** Moving through rows highlights the row's selection in 3D through the preview sink
  (`lociHighlights`), so the user sees what "Sidechain" or "Chain B" means before choosing. Typing an expression
  previews it once it parses; parse errors show inline with the language's message.
- **Counts on demand.** The highlighted row shows how much it would select ("Helix: 412 residues"), computed for that
  row only, cancelled when the highlight moves.
- **Granularity** is not part of `select`; it is a picking setting (`selection.granularity`, below): a badge on the
  canvas mode switch's Select button, changed with `PickLevel` beside it ([shell.md](shell.md) §4.6).
- The `:` command menu prefix is a shortcut into Select › Residues or Expression: `:A 10-20` reads as a residue list,
  `:pymol resn HEM` as a PyMOL expression.

**State (mocks):** the `select` group with `select.residues` (the `residues` format), Expression › `select.pymol`,
`select.vmd`, `select.jmol` (each `expression.<language>`: examples as `fill` completions to adapt, a light syntax check
whose message shows inline), and By query › ligands, helices, water (right away, fixed). Groups have no shared arguments
yet, so `mode` is a field of each select command with arguments, after "and" ("and add to the selection"; options add,
remove, intersect, set, labelled as what they do). Not yet: the segmented control and Enter's modifiers, generated
queries, In the loaded structures, All and Invert (Clear is `selection.clear`), Around and Expand to, MolScript,
previews and counts.

### 10.2 `selection`: acting on the current selection

Commands whose target defaults to the current selection; `when` greys them out with "Nothing is selected".

```
Selection ›
  Focus  ·  Zoom to  ·  Clear  ·  Undo last change
  Color…                 color (named colors, recent colors, palette swatches, hex/RGB), representations (default all)
  Set transparency…      value, representations
  Set material…          material params, representations
  Clip…                  clipping params, representations
  Set emissive…          value, representations
  Reset color  ·  Reset material
  Create component…      label, representation (registered types)
  Remove from components
  Measure ›              Distance · Angle · Dihedral · Label · Orientation · Plane  (from the selection history)
  Granularity ›          Element · Residue · Chain · Entity · Model · Structure …  (radio; picking setting)
  Selection mode         toggle
```

- Today's Theme panel (`selection` + `action` mapped select + `representations`) splits into one command per action.
  Each has a small flat form that fits inline in the command menu, and the search finds "transparency" directly instead
  of through an Action dropdown.
- **Color** is a `color` argument format whose completions are swatches (recent colors first, then the palette shown
  today) and which accepts names and hex. Lighter and darker weights and the channels live in the color picker that
  opens under the field (`ColorPicker`, [styling.md](styling.md) §4; its `parseColor` reads typed colors). When the
  field is active, moving through swatches previews the color on the selection.
- The commands reuse `StructureComponentManager.applyTheme`, `add`, and `modifyByCurrentSelection` through
  `model/actions/selection.ts` (planned), so the inspector and context menus run the same code.

**State (mocks):** the `selection` group with `selection.color` (a text field in the `color` format, then "on" all
representations, cartoon, or ball & stick), `selection.transparency` ("to" a value from 0 to 1) and `selection.clear`;
they report what they would do and have no `when` yet. Not yet: the other actions, swatches and the color popover,
previews.

### 10.3 The selection bar

The mocks settled a different bar from the one first drawn here (a mode control and the `select` group as an inline
command menu, in selection mode only). As built (`views/canvas/selection-bar`, [shell.md](shell.md) §4.6):

```
12 residues   [ 🖌 ] [ 📏 ] [ ⌖ ] [ ✕ ]
 summary       Color, theme, component… · Measure… · Focus · Clear
```

- `SelectionBar` shows while anything is selected, at the bottom right of the canvas: the summary ("12 residues"), then
  `SelectionModel.actions` as icon buttons with tooltips, then Clear. In the mock, one button opens the Selection tab
  (color, theme, component), Measure… opens it at its `measure` section, and Focus is one click.
- Selecting (the mode, add/remove/intersect/set, queries and expressions) is the command menu's `select` group, not the
  bar; the picking mode and granularity are on the canvas mode switch (`ModeSwitch`, `PickLevel`).
- The actions are meant to be `selection` commands, so the command menu offers the same ones; in the mocks they are
  closures (§7). Help goes to the command menu's `?` and the Help view; the bar has no `?`.
