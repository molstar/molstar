# Commands: the specification

`@molstar/ui`'s command layer (`src/commands/`, DOM-free, barrel `commands/index.ts`). The rules the code relies on and
why; code comments cite its sections (`commands §3.3`). The code is the truth for how. The command menu that draws it is
§5.

| File                   | What                                                                                        |
| ---------------------- | ------------------------------------------------------------------------------------------- |
| `types.ts`             | `CommandDefinition`, groups, presets, results and outcomes, the ports (formats, preview, …) |
| `registry.ts`          | `CommandRegistry`: commands, groups, owners, lazy commands (§4.1–4.2)                       |
| `query.ts`, `match.ts` | what typed text matches (§4.3), fuzzy scoring                                               |
| `manager.ts`           | `CommandManager`: open, run, outcomes, recents, shared arguments, keys (§4)                 |
| `browse.ts`            | the menu's highlighted row, previewed and counted (§4.5)                                    |
| `invocation-state.ts`  | `InvocationState`, `InvocationInit` (§3.1)                                                  |
| `invocation-model.ts`  | `InvocationModel`: values, derived views, `set`, reshaping, moving, the preview (§3.1–3.7)  |
| `invocation.ts`        | `CommandInvocation`: typing, choosing, Enter/Tab, the run, the assistant (§3.3, §3.5)       |
| `sentence.ts`          | pure: building the params, the sentence's fields, the choices, the options (§3.2–3.4)       |
| `tail.ts`              | reading the words typed after an alias (§3.6)                                               |
| `keys.ts`, `keymap.ts` | chords, the keymap's layers, the per-key index (§6)                                         |
| `items.ts`             | `commandItem`, `commandItems` (§7)                                                          |

## 1. The layer

- **Commands are user-facing**: string ids, titles, categories, arguments, enablement, keys. Menus, keys, buttons and an
  assistant reach an app's operations through them. Why: one definition, every way in agrees (enablement, keys,
  overrides).
- **Generic over its context** (`CommandManager<Ctx>`): no DOM, no app; runs in Node. A key event, a target and a scope
  are read through the few fields they have (`KeyEventLike`; `getAttribute`, `tagName`, `parentNode`).
- **States are core's** (`@molstar/core/reactive`): `State`, `derive`, `Subject` for outcomes; the manager and the
  invocation are `ReactiveModel`s (`retire` ends them). No rxjs.

## 2. Definition

`CommandDefinition<Ctx>`: `id`, `title`, `category`, `icon?` (an `IconName`), `group?`, `alias?`, `shortcut?`, `scope?`,
`when?(ctx): true | string` (a reason greys it and says why), `params?(ctx, values)`, `fields?`, `action?`, `reshapes?`,
`presets?`, `defaults?`, `preview?`, `count?`, `check?`, `run(args, ctx, signal)`, `assist?`, `view?`.

- **Arguments are PD params.** `formFields` (`forms/flat-fields.ts`) reads the sentence's: text, string selects,
  booleans, numbers. Every other param that is not hidden is an option (§3.4). `CommandValues` are the sentence's flat
  values; `CommandArgs` what a run gets (§3.5).
- **Results** (`CommandResult`: `message`, `detail`, `ops`, `undoable`, `data`) reach the host through `outcomes`
  (`done`, `failed`, `cancelled`, `blocked` for a refused key): events, never replayed.

### 2.1 How a command runs

- **Right away**: no `params`; it runs with its `defaults` and its groups' shared values.
- **As a sentence**: with `params`, `open` gives a `CommandInvocation` (§3).
- **In its place in the UI**: `view` (a card, a tab) for what needs a panel (picking in 3D, several steps).

### 2.2 Fields and formats

- `fields` by param name (`CommandFieldHint`): `before` (a quiet word: "as", "in"), `format` (a named input type),
  `presets` (a number's ready values, `n` or `[n, label]`).
- **Formats** are looked up through the `CommandFormats` port: `hint(ctx)`, `check(text, ctx)` (`{ ok, value }` or a
  message, `unreadable` for words), `complete(text, ctx)` (a completion with `fill` fills the text for editing). Why a
  port: what a "selection" or a "color" reads as is the app's.
- **Every call into a format is guarded**: `check` through `checkFormat` (`format.ts`), which reads a throw as text the
  format does not take, its message the error; `complete` and `hint` as nothing. Why: a format is the app's code, and
  every reader (typing, the tail, Enter, the run's last check) must leave the invocation editable; a throw out of the
  run's check once left it running for good.

### 2.3 Groups and shared arguments

- **Groups** (`CommandGroup`: `id`, `title`, `icon?`, `parent?`, `order?`, `shared?`) nest commands for browsing only: a
  search from the top finds every command at any depth.
- **A shared argument** (`shared`: `name`, `title`, `options`, `default`) is set once for a group (`setShared`,
  `sharedValue`) and goes to whichever command in it runs, at any depth: in the args of one without params, as the start
  of an invocation whose params have that name. Why: a few mode buttons, not a copy of the list per mode.

## 3. The invocation

### 3.1 State

`invocation.state` is one `State<InvocationState>` (compared shallowly): `status` (`editing`, `resolving` while the
assistant reads, `running`, `cancelling`, `done`, `closed`; a failed run goes back to `editing` with the error),
`params`, `fields` (the sentence's), `formOnly` (the options' names), `action`, `values`, `optionValues`, `field` (the
active one), `text`, `pick`, `proposed`, `missing`, `error` (`{ field?, message, hint? }`), `canAssist`, `note`,
`escalated`, `blocked`, `result`.

- **Values start** from `defaults`, then a preset's (or a shared argument's) values, then the typed tail (§3.6); those
  that do not fit their field are dropped, not thrown.
- **The views are derived** from narrow slices of the state (`sentence`, `choices`, `options`): typing tells the
  choices, not the sentence or the options. Why: the menu redraws only what changed (reactive §4).
- **Edits are refused unless `editing`**: while the assistant resolves or a run is in flight, `set`, `type`, `focus` do
  nothing.

### 3.2 The sentence

`invocation.sentence`: the fields in order (`SentenceField`: the `FormFieldView` with `value`, `display`,
`defaultValue`, `defaultDisplay`, `empty`, `active`, `proposed`, `differs`, `error`).

- **The active field** is the first empty one at the start (the last when all are filled). Taking a value moves on: the
  next empty required field after it (wrapping), then the next empty one.
- **Moving**: `focus(id)`, `step(±1)` (wrapping), `tab(±1)` (assist, else take a choice, else step), `back()` (an empty
  text: clear the value, else the previous field; false at the start: the menu goes back to its list).
- **Required** is the param's (not optional). An empty optional field stands for its default (`defaultDisplay`).
  `differs` compares with the param's default (`PD.isParamEqual`).

### 3.3 Choices

`invocation.choices` (`FieldChoices`, none without an active field): always the active field's, whatever it holds.

- A select: its options, an exact value or label first, then by fuzzy match on the label; a boolean: on and off; a
  format: its completions, what the text reads as first (`detail: 'as typed'`), and its `hint`; a number: its presets
  that fit and, for a slider (min and max), its `range`; a text without a format: nothing. The field's value is marked
  `current`.
- **What it follows**: the fields, the active field and its own value, the text, the context; a set of another field or
  of an option does not list the choices again. Typing is one change (the text and its `pick` in a `batch`), listing the
  choices once.
- **`invocation.range`** (`FieldRange`: the field, min, max, step): the active slider's bounds, the same object while
  they stand (compared shallowly). A view keyed on it (the slider under the pointer) is not built again as the value
  changes; keyed on `choices`, which are made anew with each value, a drag ended after its first step.
- **`pick`** is the highlighted item: typing highlights the first (not for a number: a typed number is taken as typed),
  `move(±1)` wraps, `point(i)` follows the pointer. The preview shows the highlighted item (§3.7).
- **`accept(i?)`** takes the highlighted item (a `fill` item goes into the text instead), else reads the text: a number
  ("60%" is 0.6), an option by value or label, a format's `check`, plain text. A refusal is the field's `error`.

### 3.4 Every param is editable

- **`set(name, value)` takes any param**: a field of the sentence (checked by its kind; moves on unless `stay`), or an
  option, a form-only or hidden param (kept in `optionValues`, given to the run, the preview and `check` over its
  default). `undefined` or `clear(name)` empties a field or puts an option back. Why: the menu's Options fold (and a
  command's own form) edit the same invocation; no param needs a panel just to be set.
- **`invocation.options`** (`InvocationOptions`): `params` (the form-only params, for a params form), `values` (set over
  the defaults), `changed` (`{ name, label, value }` in order: what differs from the default, the fold's summary).
- **Reshaping**: a param in `reshapes`, set or cleared (a field or an option), builds the params again from the
  sentence's values and the primitive options, keeping what still fits; a build that throws keeps the params and says
  why on the field.

### 3.5 Running

- `enter()`: take what is typed or highlighted (and run when that leaves no field empty), else `confirm()`.
- `confirm()`: refused with `blocked` or a missing field (the active field goes there); then `when`, the formats and
  `check(args)` are asked again; then `run(args, ctx, signal)`. Done closes it (through the manager) and clears the
  preview; a failure stays open with the error; `cancel()` aborts the signal, and a run that settles aborted reports
  `cancelled` (one that lands anyway is `done`: its effect is there).
- `args()` = the params' defaults, the sentence's values, the options set (`argsFor`). Optional fields left empty are
  absent.
- **The assistant** (`CommandAssistant`, with `assist` on the command): words a field cannot read set `canAssist`; Tab
  asks `propose` (values marked proposed, an answer that is `done`, or `escalate` with the words kept).

### 3.6 The tail

Words after a direct alias (§4.3), read in three passes (`tail.ts`): a field's before-word starts its value up to the
next before-word; then words naming an option or a listed format value; then what is left, read by an open field's
format as one run. What no field takes is typed into the first empty field (and taken when its format reads it).

### 3.7 The preview

Drawn into the `PreviewSink` while editing, once the required fields are filled (with the highlighted choice) and `when`
allows; cleared otherwise and on close. Drawn again only when the values, the options, the highlighted choice or the
context (`previewDeps`, item by item) changed. A run keeps it as it was (its own edit must not redraw it); done clears
it.

## 4. The manager

`CommandManager<Ctx>`; options: `context`, `formats`, and optionally `preview`, `assistant`, `recents` (`RecentsStore`),
`maxRecents`, `previewDeps`, `keymap`, `browseDelay`.

### 4.1 Registration

- **`register(defs)` is one change of `commands`** and returns its undo (one change too), which cancels a run of theirs
  (it then reports to no one) and closes their invocation. An id already registered (or repeated) throws, naming it, and
  nothing is added. `registerGroups` likewise. Why one change: a screen registers hundreds; each change re-derives every
  reader (and the key index, §6.1).
- `get`, `groupPath`, `children`, `inGroup` serve browsing; lists are indexed by id once per list.

### 4.2 Lazy commands

`registerLazy({ make })`: made when wanted, not at each change of what they come from. `stale()` marks them; the next
`query`, `children`, `groupPath`, a `get` of an id they may hold (one they made, or one not registered) or a key the
app's or the person's layer binds makes them, or a microtask while a view holds `want()`, so a render never registers.
`registered(id)` looks a command up without making anything: what a derived value or a render reads (`CommandButton`
holds `want()` while shown and derives with `registered`). A make's swap is one change (`batch`). `make` returning
`undefined` keeps what was made; one that throws is logged, and the next change tries again. An unbound key makes
nothing (`handleKey` runs at every keystroke).

### 4.3 Query

`query(text)`: `recent` (empty text), `matches` (fuzzy subsequence over the title, aliases a little lower; equal whole
scores keep categories together, then by title), `direct` (an alias followed by more text: that command with the rest as
its `tail`; the longest alias, one that can run first; none while the text is still the start of a title or alias). A
preset is listed only when the text says more than the command alone. A blocked direct match lists what the whole text
matches after it.

### 4.4 Open, close, cancel

`open(id, tail?, preset?)` never throws: `form` (an invocation, `active`), `ran` (no params), or `blocked` with a reason
(a run in flight, `when`, an unknown preset, params that throw). `close()` (false while it runs), `cancel()`, `running`,
`outcomes`, `recents`.

### 4.5 Browsing

`browse(id | undefined)`, `browsing`: the menu names its highlighted command; one without params that can run is
previewed and counted (`count`) after `browseDelay` (80 ms), so moving through the list computes only where it rests. A
new highlight, none, or `open` aborts the count's signal and clears the preview.

## 5. The command menu

`components/command-menu` (`CommandMenu`; files: `command-menu.tsx` the shell, `list.tsx` + `rows.tsx` the list,
`sentence-mode.tsx` + `sentence.tsx` the sentence, `choices.tsx`, `options-fold.tsx`, `footer.tsx` keys, `types.ts` the
props and ports). Props: `manager`, `onClose`, `initial?`, `group?`, `prefixes?`, `assistant?`, `direct?`, `forms?`,
`fieldViews?`. `CommandButton` is in `components/button/command-button.tsx` (§7's rules, as a button).

### 5.1 The shell and the list

- **A modal dialog over its positioned container** (a scrim fills it; what is behind stays visible for previews). Tab
  stays inside; closing gives the focus back to what opened it. A `done` outcome closes it. Open, it holds `want()`
  (§4.2). Why a container, not the viewport: a host places it over the part of the app it commands.
- **One input, a combobox that keeps the focus**: rows are `option`s in a `listbox`, captions `presentation`; the input
  names the highlighted one (`aria-activedescendant`). A polite status says what a search found once typing pauses.
- **Browsing** (empty): Recent, the top-level groups (Browse), then the ungrouped commands under their categories'
  captions; inside a group (a crumb before the input; ⌫ or ← on an empty input goes up, → on a group row opens it), its
  subgroups then its commands. Rows show the category only while searching (the captions say it while browsing).
  **Searching**: a direct alias first, then matching groups, then commands at any depth, the matched letters marked.
- **Rows** are 24px list rows (styles §4–5): the 12px icon in the slot at the inset (the head's 16px icon above it), the
  title (a UI label, never cut), the meta (a path, a reason, a count: variable, it fades), the key as caps (`kbd` per
  key in a `kbd`). Picking a blocked row says why above the list. Rows are made again when the context changes, so a
  row's reason it cannot run follows the app while the menu is open; rows that stay the same keep their objects (and
  their nodes), and the same rows are the same list.
- **The highlighted row** goes to `manager.browse` (§4.5): its count replaces its meta. A group's shared argument shows
  as toggles under the head.
- **A screenful at a time**: rows render 24 at first, more as the list scrolls or the highlight passes them; a spacer
  keeps the rest's height. Not `VirtualList`: rows and captions mix, and the highlighted option must exist for
  `aria-activedescendant` without a layout pass.
- **The list's filter is a frame state** (dom §3): keys faster than frames filter once.
- **Prefixes and `?`**: a prefix's items replace the rows (`{ text }` from an item edits the line); `?` lists the
  prefixes (picking one types it), then every command with a key.
- **The Ask row** (`CommandMenuAssistant`: `ready`, `propose`, `ask`, `setup`), only when the host gives one: ⌘↵, or ↵
  when nothing matches (it is then drawn as the default). A proposal opens its sentence filled in, under a line saying
  where it came from; nothing runs until Enter. Not ready: a quiet row to the setup.

### 5.2 Direct

Opened for one command (`direct`): Esc, and ⌫ past the sentence's start, close the menu rather than go back to the list;
the footer says which.

### 5.3 The sentence is the editor

- **The head is the sentence**: the command's icon in the slot, its title on the rows' text column, the fields as tokens
  a row tall, then ⌥↵ "Open in …" (with a `view`) and the run button (`action`, ↵). Why one row: a field is edited where
  it reads, not in an input away from it.
- **The caret is in the active field's token**: the one input moves into it (keeping the focus), sized to its text by a
  hidden copy. Not active, a token is a button (named `label: value`): its value; required and empty, its name dashed
  and italic, muted; optional and empty, its default, muted. Active, the input's placeholder is what the token reads as
  (a set value in the text's color: typing replaces it).
- **Moving**: a click on a token, Tab / Shift+Tab (§3.2's `tab`), ← / → at the text's edge (not past the first or last
  field), ⌫ on an empty input (`back`: clear, then the field before, then the list). Esc: the list (or close, direct).
- **The list is the active field's choices, always** (`invocation.choices`): the current one checked in the slot and in
  the accent; ↑↓ move, ↵ or a click takes (and moves on). One line above it: the field's error (an alert, with its
  hint), else the assistant's state or note, else the format's hint. A slider's range is a row above the list (dragging
  sets the value and keeps the field). A press in the list keeps the focus in the input.
- **The Options fold** holds the params the sentence does not (`state.formOnly`), closed by default; its head says how
  many, or what differs from the defaults (`options.changed`). Open, a `ParamsForm` whose edits go to `set` (so the run,
  the preview and `check` see them); made again only when the command reshapes, and told of values set elsewhere.
- **Outcomes**: done closes the menu; a failure keeps the sentence with its error.

### 5.4 Custom UIs

- **A command's own form** (`forms[id]`: `CommandFormView`: `fields?`, `options?: false`, `render(props)`), under the
  choices. While a field it names is active, no choices are listed (it draws them). `options: false` hides the fold (it
  draws those params itself).
- **A field editor by format** (`fieldViews[format]`: `CommandFieldView`), in place of the choices while a field of that
  format is active: `field`, `choices`, `set(value)` besides the invocation. Why by format: a picker (a color's
  swatches) serves every command with such a field, not one.
- **The input keeps the focus and the keys**: a view's `keys(handler)` sees each key first (true takes it); the rest are
  the sentence's. Edits go through the invocation (`set`, any param), `run()` is Enter.
- **Tab reaches what is below**: Tab on the last field (Shift+Tab on the first), nothing typed to take, moves the focus
  to the first (last) Tab stop below the choices (the own form's controls, the fold); past the last (before the first)
  it comes back to the input on the first (last) field. Controls with `tabIndex` -1 (tiles the arrows move among) are
  not stops. Nothing below: Tab stays the sentence's.

### 5.5 Footer

The keys for what is active (the list's, or the sentence's), as caps.

## 6. Keys

### 6.1 `handleKey(event, mac?)`

Runs a command's key from a keydown the host forwards while the menu is closed; true when it took the event.

- Ignored: repeats, IME composition, handled events. In a text field only ⌘/Ctrl chords without Alt run; a control that
  steps with arrows keeps plain keys. The key is matched first, then its `code` when the key is not a plain letter or
  digit (⌥T is "†" on macOS). A refused command is a `blocked` outcome, except on a plain key, which is then not taken.
- **The per-key index** (`KeyIndex`): commands by their key in lower case, made once per list of commands and keymap. A
  keydown looks up its key and code, and reads the target and its path only when a command is on that key; it allocates
  no list. Whether a lazy contribution is stale is a count the registry keeps, not a walk of them. Why: every keydown in
  the app comes here.

### 6.2 Scopes

`scope` on a command: `global` (the default) works anywhere; a scope id works only while the key's target is inside an
element with that `data-msui-scope` (`ScopeAttribute`; read from the composed path, across shadow roots). Among commands
on one key the nearest scope comes first, global last, then the first registered; the first that can run takes it, else
the nearest is refused. Why: a key can mean one thing in a panel and another elsewhere, and a scoped command that cannot
run never swallows the global one.

### 6.3 The keymap

`manager.keymap` (`Keymap`), three layers, last wins: the command's `shortcut`, the app's (`setApp`), the person's
(`remap(id, chord | null)`, `reset`, `resetAll`), saved through the host's `KeymapStore`. `null` removes a key.
`manager.shortcut(c)` is the key now. What a store gives back is checked (chords and `null` only); a store that throws
keeps nothing. `chordFromEvent` records a remap (none for a modifier alone, Tab, the platform's other modifier).

### 6.4 Conflicts

`keymap.conflicts(commands)`: commands on one key in one scope; a scoped one over a global one is not a conflict.
`sharing(id, chord, commands)`: what a remap would collide with.

### 6.5 Chords

Written once with `mod` (⌘ on macOS, Ctrl elsewhere); `chordKeys` (one cap per `kbd`) and `formatChord` show them per
platform, `chordId` compares them (a symbol's Shift is part of it). **The platform** is `KeyPlatform.mac`, detected once
from the browser by core's `isApplePlatform` (false outside a browser) and settable by a host (tests, a forced
platform); every default `mac` reads it. Why: one detection for the whole UI, in core with the other browser checks.

## 7. Menus and buttons by command id

`commandItem(manager, id, { label?, hint?, values?, id?, onBlocked?, mac? })`: an entry as the command is when the menu
is built (title or label, icon, key after the keymap, `disabled` with its `when` reason, `run` through `manager.open`:
with `values`, the row it is for). None when the command is not registered, so an override that removes a command takes
its entries; `commandItems` keeps those that are. Its shape fits `ContextMenuItem` as it is. Why: a menu repeats no
action's logic, and keys, overrides and enablement apply everywhere at once.
