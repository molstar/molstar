# Commands, the command menu, keybindings

prism-ui's command layer (`src/commands/`, DOM-free) and its one view, the command menu
(`src/views/command-menu/command-menu.tsx`). What commands an app registers is the app's (the viewer's are over its
operations: the viewer `_spec/api.md` §9). The code is the truth for how.

## 1. The layer

- **Commands are user-facing**: string ids, titles, categories, arguments, enablement, keybindings; the command menu,
  keys, buttons and the assistant reach operations through them. Mol\*'s `PluginCommands` are not used for this.
- **Generic over its context** (`CommandManager<Ctx>`): no DOM, no plugin; runs in Node.

## 2. Definition

`CommandDefinition<Ctx>` (`commands/types.ts`): `id`, `title`, `category`, `icon?`, `group?` (by id; browsing only),
`alias?` (other names, and words that start it with arguments), `shortcut?` (`KeyChord`, `mod` = ⌘ or Ctrl),
`when?(ctx): true | string` (a reason greys it and says why), `params?(ctx, values): PD.Params`, `fields?` (inline
before-words and formats by param), `action?` (the run button's label), `reshapes?`, `presets?` (search-only ways in),
`defaults?`, `preview?`, `check?`, `run(args, ctx, signal)`, `assist?`, `view?` (where it lives in the UI).

- **Arguments are PD params.** `formFields` (`forms/flat-fields.ts`) reads the inline ones: text, string selects,
  booleans, numbers; the rest are form-only and run at their PD defaults (`argsFor`). `CommandValues` are the flat
  inline values; `CommandArgs` add the defaults.
- **Formats** are named argument types looked up through the `CommandFormats` port: `hint`, `check` (`{ ok, value }` or
  a message, `unreadable` for words), `complete` (a completion with `fill` fills the text for editing).
- **Results** (`CommandResult`: `message`, `detail`, `ops`, `undoable`, `data`) reach the host through `outcomes`
  (`done`, `failed`, `cancelled`, `blocked` for a refused key); the message goes to the status column.
- **Groups** (`CommandGroup`: `id`, `title`, `icon?`, `parent?`, `order?`) nest commands for browsing only: a search
  from the top still finds every command at any depth and shows its path.

### 2.1 How a command runs

- **Right away**: no `params`; it runs with its `defaults`.
- **As an inline sentence**: with `params`, `open` gives a `CommandInvocation`: the title, each inline field as a quiet
  before-word and a value, the run button. Form-only params show as "+ N more options at their defaults".
- **In its place in the UI**: a command with `view` offers it beside Run (⌥↵): the menu closes and the card or tab opens
  there, for what the sentence leaves out (a preview, the full form). What needs several steps or picking in 3D is a
  panel or section, not something inside the menu.

## 4. The manager

`CommandManager<Ctx>` (`commands/manager.ts`); options: `context` (a `ReadonlyState<Ctx>`), `formats`, and optionally
`preview` (a `PreviewSink`), `assistant` (`CommandAssistant`), `recents` (`RecentsStore`), `maxRecents`, `previewDeps`.

- **`register(defs)` returns its undo**, which cancels a run of theirs and closes their open form; an outcome that
  settles after it reports to no one. Registering an id again replaces the earlier command. `registerGroups` likewise;
  `groupPath`, `children`, `inGroup` serve browsing.
- **`query(text)`**: `recent` (empty text), `matches` (fuzzy subsequence over title, then aliases a little lower; equal
  scores keep categories together), `direct` (an alias followed by more text: that command with the rest as its tail). A
  preset is listed only when the text says more than the command alone.
- **`open(id, tail?, preset?)`** never throws: `form`, `ran`, or `blocked` with a reason (a run in flight, `when` says
  no). `close()`, `cancel()` (aborts the run's signal), `running`, `active`.
- **`CommandInvocation`** (`commands/invocation.ts`): states `editing`, `resolving` (the assistant reads words),
  `running`, `cancelling`, `done`, `closed`; a failed run goes back to `editing` with the error. Values start from
  `defaults`, then a preset, then the typed tail; a `reshapes` field rebuilds the params keeping what fits. The preview
  shows once the sentence is complete and `when` allows. Before a run `when`, the formats and `check` are asked again.
- **Reading the tail** (`commands/tail.ts`): words after a field's before-word; then words naming an option or a format
  value; then what is left, read by an open field's format.

## 5. The command menu

`CommandMenu` props: `manager`, `onClose`, `initial?`, `group?` (open inside it), `prefixes?`, `assistant?`, `direct?`
(opened for one command: Esc closes the menu).

- **Opened by** ⌘K / Ctrl+K and ⇧⌘P / Ctrl+Shift+P (the host listens), the rail's button, and buttons that open it
  inside a group (the viewer's Preset and Styles). Modal over the shell; the canvas stays live.
- **Browsing** (empty input): Recent, the top-level groups, then the ungrouped commands by category; inside a group, its
  subgroups then its commands. **Searching**: every command in scope, a direct alias match first, then matching groups,
  then commands. Rows render a screenful at a time (styling.md §5).
- **Rows**: icon, title with the matched letters, the reason it cannot run or its path or category, its shortcut
  (`formatChord`). Picking a blocked row says why above the list.
- **Prefixes** (`CommandMenuPrefix`: `char`, `title`, `placeholder?`, `items(text)`) turn the menu into another list;
  the host supplies them. `?` is built in: the prefixes, then every command with a shortcut.
- **The sentence**: ↑↓ through completions, Tab takes one or moves on (Shift+Tab back), Enter runs once no required
  field is empty (else goes there), Backspace on an empty input clears the field, then goes back a field, then to the
  list; Esc back to the list, Esc there closes. "50%" reads as 0.5. A run that is done closes the menu; a failure keeps
  the sentence open with the error.
- **The assistant** (`CommandMenuAssistant`: `ready`, `propose(text)`, `ask(text)`, `setup`), only when the app gives
  one: an Ask row takes what is typed (⌘↵, or ↵ when nothing matches); a proposed command opens as its sentence, filled
  in, and nothing runs until Enter; other words go to `ask`. Without a model the row offers the setup.
- **Assist inside a field** (core): with a `CommandAssistant` and `assist` on the command, words a field cannot read set
  `canAssist`, and Tab asks `propose` (values marked proposed, an answer, or `escalate`).
- **ARIA**: the input is a combobox and keeps focus; rows are `role="option"` in a `listbox`.

## 6. Keybindings

- **`CommandManager.handleKey(event, mac)`** runs a command's `shortcut` from a keydown the host forwards while the menu
  is closed. Ignored: repeats, IME composition, handled events. In a text field only ⌘/Ctrl chords without Alt run; a
  control that steps with arrows keeps plain keys. The key is matched first, then its `code` when the key is not a plain
  letter or digit (⌥T is "†" on macOS). Among commands with one chord, one that can run wins; a refused one is a
  `blocked` outcome, except a plain key, which is then not taken.
- **Keys go to one viewer**: the host checks `ownsKey(root, event)` (`dom/events.ts`, [library.md](library.md) §3)
  before forwarding.
- **Chords are written once with `mod`**; `chordKeys` and `formatChord` (`commands/keys.ts`) show them per platform.
- **Canvas keys coexist**: the plugin's trackball and behavior bindings keep working on the canvas; commands do not
  override them, and the Help card lists both (`views/shell/help-panel.tsx`).
