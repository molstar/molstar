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
`alias?` (other names, and words that start it with arguments), `shortcut?` (its default `KeyChord`, `mod` = ⌘ or Ctrl;
the keymap may replace it, §6), `scope?` (where the key works, §6), `when?(ctx): true | string` (a reason greys it and
says why), `params?(ctx, values): PD.Params`, `fields?` (inline before-words and formats by param), `action?` (the run
button's label), `reshapes?`, `presets?` (search-only ways in), `defaults?`, `preview?` (its form's, and its row's while
browsed, §5), `count?` (how much it would do, on its highlighted row, §5), `check?`, `run(args, ctx, signal)`,
`assist?`, `view?` (where it lives in the UI).

- **Arguments are PD params.** `formFields` (`forms/flat-fields.ts`) reads the inline ones: text, string selects,
  booleans, numbers; the rest are form-only and run at their PD defaults (`argsFor`). `CommandValues` are the flat
  inline values; `CommandArgs` add the defaults.
- **Formats** are named argument types looked up through the `CommandFormats` port: `hint`, `check` (`{ ok, value }` or
  a message, `unreadable` for words), `complete` (a completion with `fill` fills the text for editing).
- **Results** (`CommandResult`: `message`, `detail`, `ops`, `undoable`, `data`) reach the host through `outcomes`
  (`done`, `failed`, `cancelled`, `blocked` for a refused key); the message goes to the status column.
- **Groups** (`CommandGroup`: `id`, `title`, `icon?`, `parent?`, `order?`, `shared?`) nest commands for browsing only: a
  search from the top still finds every command at any depth and shows its path.
- **A shared argument** (`shared`: `name`, `title`, `options`, `default`; Select's mode) is set once for a group and
  goes to whichever command in it runs, at any depth: in the args of one without params, as the start of a form whose
  params have that name. So four mode buttons do not each open a copy of the same list. Put only commands that take it
  in such a group (an operation refuses an unknown value).

### 2.1 How a command runs

- **Right away**: no `params`; it runs with its `defaults`.
- **As an inline sentence**: with `params`, `open` gives a `CommandInvocation`: the title, each inline field as a quiet
  before-word and a value, the run button. Form-only params show as "+ N more options at their defaults".
- **In its place in the UI**: a command with `view` offers it beside Run (⌥↵): the menu closes and the card or tab opens
  there, for what the sentence leaves out (a preview, the full form). What needs several steps or picking in 3D is a
  panel or section, not something inside the menu.

## 4. The manager

`CommandManager<Ctx>` (`commands/manager.ts`); options: `context` (a `ReadonlyState<Ctx>`), `formats`, and optionally
`preview` (a `PreviewSink`), `assistant` (`CommandAssistant`), `recents` (`RecentsStore`), `maxRecents`, `previewDeps`,
`keymap` (the app's keys and the store of the person's, §6), `browseDelay` (§5).

- **`register(defs)` returns its undo**, which cancels a run of theirs and closes their open form; an outcome that
  settles after it reports to no one. An id already registered (or repeated in `defs`) throws, naming it, and nothing is
  added: an override is made before registering (the viewer's `overrides`). `registerGroups` likewise; `groupPath`,
  `children`, `inGroup` serve browsing.
- **`registerLazy({ make })`**: commands made when wanted, not at each change of what they come from (the Select group's
  per-structure commands). `stale()` marks them; the next `query`, `children`, `groupPath`, a `get` of an id they may
  hold (one they made, or one not registered) or a key the app's or the person's layer binds makes them, or a microtask
  while the menu holds `want()` (its mount), so a render never registers. An unbound key makes nothing (`handleKey` runs
  at every keystroke). `make` returning `undefined` keeps what was made (an open form of one survives); one that throws
  is logged, and the next change tries again.
- **`query(text)`**: `recent` (empty text), `matches` (fuzzy subsequence over title, then aliases a little lower; equal
  scores keep categories together), `direct` (an alias followed by more text: that command with the rest as its tail). A
  preset is listed only when the text says more than the command alone.
- **`open(id, tail?, preset?)`** never throws: `form`, `ran`, or `blocked` with a reason (a run in flight, `when` says
  no). `close()`, `cancel()` (aborts the run's signal), `running`, `active`.
- **`CommandInvocation`** (`commands/invocation.ts`): states `editing`, `resolving` (the assistant reads words),
  `running`, `cancelling`, `done`, `closed`; a failed run goes back to `editing` with the error. Values start from
  `defaults`, then a preset, then the typed tail; a `reshapes` field rebuilds the params keeping what fits. The preview
  shows once the sentence is complete and `when` allows. Before a run `when`, the formats and `check` are asked again.
- **Shared values**: `sharedValue(group)` (the person's pick, else the default), `setShared(group, value)`, `shared`.
- **Browsing** (`browse(id | undefined)`, `browsing`): the command menu names its highlighted command; one without
  params that can run is previewed (`preview` into the sink) and counted (`count`) after `browseDelay` (80 ms), so
  moving through the list computes only where it rests. A new highlight, none, or `open` aborts the count's signal and
  clears the preview.
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
- **Rows**: icon, title with the matched letters, the reason it cannot run or its path or category, its key after the
  keymap (`formatChord`). Picking a blocked row says why above the list.
- **Inside a group with a shared argument** (Select), its options show once under the input as toggle buttons (Tab
  reaches them; the input keeps the focus otherwise); the choice holds for the session.
- **The highlighted row** is previewed (the viewer lights its atoms in 3D) and its count ("1 residue · 20 atoms")
  replaces its path: that row only, through `browse` (§4). Closing the menu or opening a form clears it.
- **Prefixes** (`CommandMenuPrefix`: `char`, `title`, `placeholder?`, `items(text)`; DOM-free, in `commands/types.ts`)
  turn the menu into another list; the host supplies them. `char` is matched as the text's start, so a word with its
  space (`pymol `) is a prefix too. An item's `run` returns true to keep the menu open, or `{ text }` to put that text
  in the menu and keep it open: a completion edits the line rather than runs. `?` is built in: the prefixes, then every
  command with a shortcut.
- **The sentence**: ↑↓ through completions, Tab takes one or moves on (Shift+Tab back), Enter runs once no required
  field is empty (else goes there), Backspace on an empty input clears the field, then goes back a field, then to the
  list; Esc back to the list, Esc there closes. "50%" reads as 0.5. A run that is done closes the menu; a failure keeps
  the sentence open with the error.
- **Own forms** (§5.1) replace a command's option list with its own view.
- **The assistant** (`CommandMenuAssistant`: `ready`, `propose(text)`, `ask(text)`, `setup`), only when the app gives
  one: an Ask row takes what is typed (⌘↵, or ↵ when nothing matches); a proposed command opens as its sentence, filled
  in, and nothing runs until Enter; other words go to `ask`. Without a model the row offers the setup.
- **Assist inside a field** (core): with a `CommandAssistant` and `assist` on the command, words a field cannot read set
  `canAssist`, and Tab asks `propose` (values marked proposed, an answer, or `escalate`).
- **ARIA**: the input is a combobox and keeps focus; rows are `role="option"` in a `listbox`.

### 5.1 A command's own form

- **`CommandMenuProps.forms`**, by command id (`CommandFormView`: `fields?`, `render(props)`), from the app (the
  viewer's Show as… and Color by theme…): drawn under the sentence. While a field it names in `fields` is active, its
  options are not listed (the view draws them); other fields list theirs as before. "+ N more options" is not shown (the
  view shows what it offers).
- **The input keeps the focus and the keys**: `props.keys(handler)` sees each key first (true takes it: arrows through a
  grid, Enter on what is typed); the rest is the sentence's (Tab, Enter, ⌫, Esc). The view edits the invocation: `set`
  for inline fields, `setParam(name, value)` for a form-only or hidden param (in the run's args and the preview, kept in
  `state.params`, dropped when a reshape removes the param), `run()` as Enter.
- **Its own controls by Tab**: Tab on the sentence's last field (Shift+Tab on its first), with nothing typed to take,
  moves the focus into the view's controls that are Tab stops (a type's size, the uniform color's swatches; its choice
  tiles are not: the input's arrows move among them); Tab past the last control (Shift+Tab before the first) gives the
  focus back to the input, on the sentence's first (last) field. A view with no such control leaves Tab the sentence's.
  Why: what a view adds is offered nowhere else (hidden params, no "+ N more"), and a keyboard must reach it without the
  input giving up its keys.
- Why: a choice among many drawn things (types, themes) reads as what it draws, and keyboard flow stays the menu's. A
  hidden param is offered nowhere else, so `formOnly` leaves it out.

## 6. Keybindings

- **`CommandManager.handleKey(event, mac)`** runs a command's key from a keydown the host forwards while the menu is
  closed. Ignored: repeats, IME composition, handled events. In a text field only ⌘/Ctrl chords without Alt run; a
  control that steps with arrows keeps plain keys. The key is matched first, then its `code` when the key is not a plain
  letter or digit (⌥T is "†" on macOS). A refused command is a `blocked` outcome, except a plain key, which is then not
  taken. The commands are indexed by their key's letter once per list of commands and keymap, so a keydown checks only
  the commands on its key and allocates no list. Why: every keydown in the app comes here (17 µs at 437 commands
  before).
- **Scopes** (`CommandDefinition.scope`): `global` (the default) works anywhere in the app; a scope id (`canvas`,
  `panel:<tab id>`, any name a view gives) works only while the key's target is inside an element with that
  `data-msp-scope` (read from the event's composed path, across shadow roots). The shell marks the canvas (`canvas`) and
  each dock panel (`panel:<id>`). Among commands on one key, the nearest scope comes first, global last, then the first
  registered; the first that can run takes it, else the nearest is refused. Why: a key can mean one thing in a panel and
  another elsewhere without either losing it, and a scoped command that cannot run never swallows the global one.
- **The keymap** (`commands/keymap.ts`, `manager.keymap`), three layers, last wins: the command's `shortcut`, the app's
  (`keymap.app`: the viewer's `ui.keybindings`), the person's remaps (`remap(id, chord | null)`, `reset`, `resetAll`),
  saved through the host's `KeymapStore` (the viewer: its preferences). `null` removes a key. `manager.shortcut(c)` is
  the key now; the menu, Help, `?` and command buttons show it. What a store gives back is checked (chords and `null`
  only); a store that throws keeps nothing.
- **Conflicts** (`keymap.conflicts`): commands on one key in one scope; a scoped one over a global one is not a conflict
  (it wins inside its scope, by design). Help lists them; a remap that makes one says so.
- **Remapping** (Help, `views/shell/help-panel.tsx`): Change keys shows each row's Change (the next chord pressed is
  taken, `chordFromEvent`; Escape cancels; a modifier alone waits; Tab and the platform's other modifier are never
  keys), Clear and Reset, "Give a key to…" for a command without one, and Reset all. The recorder takes its keys before
  anything else (no shortcut runs while it records); the result is said in a polite live region.
- **Keys go to one viewer**: the host checks `ownsKey(root, event)` (`dom/events.ts`, [library.md](library.md) §3)
  before forwarding.
- **Chords are written once with `mod`**; `chordKeys` and `formatChord` (`commands/keys.ts`) show them per platform,
  `chordId` compares them (a symbol's Shift is part of it).
- **Canvas keys coexist**: the plugin's trackball and behavior bindings keep working on the canvas; commands do not
  override them, and the Help card lists both (`views/shell/help-panel.tsx`).

## 7. Menus and buttons by command id

- **`commandItem(manager, id, { label?, hint?, values? })`** (`commands/items.ts`): a menu entry (the context menu's,
  the canvas toolbar's) as the command is when the menu is built: its title or the given label, icon, key, `disabled`
  with its `when` reason, and `run` through `manager.open` (with `values`: the row it is for; a command with params
  opens its form). None when the command is not registered, so an override that removes a command takes its entries
  (`commandItems` filters). The viewer's camera Reset menu is made this way.
- **`CommandButton`** (`components/button/command-button.tsx`): the reactive button: icon, tooltip with the title and
  key, greyed with the reason while it cannot run, gone while the command is not registered.
- Why: a menu or button repeats no action's logic, every action stays findable in the command menu, and keys, overrides
  and enablement apply everywhere at once.
