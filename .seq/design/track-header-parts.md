# Track header parts: actions as registry entries, in the header's fixed slots

## Intent

The author, 2026-11-19 ("Fourth slice"): header actions are registry entries (icon, label, when, command), drawn in the
header's fixed slots, hit and announced. The built-ins (hide, color, filter) use the same entries. A host adds or
removes its own. A hidden-tracks chip in the section header shows a hidden track again. The headers' geometry is an API
that a host's own widget can use later. The story and its test, also the author's: a "color by this track" header
action on a fake 3D viewer. The click hands the track's color runs to the host, which builds its viewer's state. The
action is a host's registry entry, not built into Seq\* ([structure-syncing.md](structure-syncing.md)). On drawn against
DOM, the author: "can we support both?": the entries drawn on the canvas by default, or a host's DOM parts over the same
geometry; and on the secondary click, "configurable behavior again so it's the composing app's decision".

What is built today:
- the commands `track.hide`, `track.show`, `layer.hide`, `layer.show` and `layer.filter`;
- the "n filtered" chip, which is the header's only control. It is drawn in the label column and hit as
  `target: 'filter'`. A click on it clears the filter.

A hidden track is left out of the resolved view, so nothing on screen brings it back, and a hidden layer is the same.
[extensibility.md](extensibility.md) kept chrome parts "internal until an outside user". This design opens the track
header as the first of them. [ux-foundations.md](ux-foundations.md) section 9 already says what a header shows: "a track
header names its data and shows its state (hidden layers, filters, a pending load) as icons in fixed slots".

The design covers:
- what a track's header and a section's header hold, and where;
- the entry a host writes;
- the menu that lists every entry;
- how each is reached by pointer, keyboard and touch, and what is said;
- the geometry a host reads.

It does not cover:
- legends' library surface (a design of its own in this slice);
- the hover card's sections (internal, as built);
- what "color by this track" builds ([structure-syncing.md](structure-syncing.md)).

## Precedents

- **The current Seq\***:
  - a track menu (`TrackActions`, a DOM overlay on the hovered header line) built from a config list of built-ins
    (`filter`, `retry`, `hide`, `remove`, `focus`, `rows`, `colors`) and host objects (`id`, `label`, `when`, `run`,
    `checked`, `danger`);
  - `HiddenChips` and switches in a `SectionHeader` part;
  - a header hover that lights the track's shown items (EC-146);
  - a menu from the keyboard, deduplicated against the browser's (EC-156).

  Taken:
  - host entries with a `when`;
  - built-ins and a host's entries in one list;
  - hidden chips in the section header;
  - the keyboard's menu.

  Not taken:
  - a `run` function: an entry names a command here, so a key, an agent and a test run what the header runs;
  - the header lighting every item (Tensions);
  - `remove` as a built-in: hiding is the reversible form;
  - the menu as an overlay that follows the hovered line: the slots are drawn and hit with the picture (Decision).
- **IGV (desktop and igv.js)**: a track's name panel, a right-click or the gear opening the track menu (rename, color,
  data range, autoscale, display mode, remove), and custom menu items per track type (igv.js; details unverified).
  Taken:
  - a right-click anywhere on a track opens its menu;
  - color and the data range in the menu.

  Not taken:
  - display modes: the density steps adapt ([ux-foundations.md](ux-foundations.md) section 4);
  - rename and remove as built-ins.
- **JBrowse 2**: a label on each track with a drag handle, the name, a kebab menu and a close button. Closed tracks come
  back through the track selector, a drawer of checkboxes. A plugin's display model adds `trackMenuItems()`. A
  right-click on a feature opens the feature's menu. Taken:
  - the close as a hide;
  - menu items from extensions;
  - a feature's own menu (`on: 'element'`);
  - a list for showing tracks again: the chip's menu.

  Not taken:
  - the track selector as a panel: tracks are in the spec, and a hidden one is listed where it was;
  - drag to reorder now (planned).
- **Nightingale and ProtVista**: no controls in a track; the host draws its labels, category toggles and tooltips
  beside the components. Taken: a host's own widget beside a header, here from the geometry. Not taken: chrome left
  entirely to the host.
- **Mol\*'s structure component rows**: an eye toggle and a "⋯" menu of actions per row. Taken: the eye in the first
  slot and "⋯" for the rest.
- **Jalview**: hidden rows and columns marked where they were, "show all hidden" in the menu. Taken: the chip where the
  tracks were, its click showing all.
- **This codebase**:
  - the "n filtered" chip (core's Render, Filters) and the hidden-rows chip ("3 hidden", `rows.show`): one count chip
    whose click undoes the state;
  - `@molstar/ui`'s `CommandDefinition` (`title`, `icon`, `when` giving a reason, `params`) and its context menu;
  - the hover card's port on its behavior ([ports.md](ports.md));
  - the layout's `chrome` and `hitChrome` and the keyboard cursor's `part: 'header'` on rows;
  - `viewer.register(entry)` with its undo ([viewer-composition.md](viewer-composition.md)).

### Interactions in the precedents

| What a user does | Where | Here |
|---|---|---|
| hides a track from its header | JBrowse (×), Mol\* (eye), the current Seq\* (menu) | built: the `seq:hide` slot, `track.hide` |
| shows hidden tracks again | JBrowse (track selector), the current Seq\* (hidden chips), Jalview | built: the hidden-tracks chip (a click shows all; its menu shows one) |
| opens a track's menu (gear, kebab, right-click) | IGV, JBrowse, the current Seq\* | built: the "⋯" slot, a secondary click on the track (the `headers` behavior's `contextClick`), the menu key on its header |
| colors a track, a scheme from its menu | IGV, Jalview, the current Seq\* | built: `seq:color`, "Color ▸" with the fitting schemes, each previewed while hovered ([themes.md](themes.md)) |
| hides one layer of a track, and shows it again | Mol\* (a representation's eye) | built: "Layers ▸" in the menu (`layer.hide`, `layer.show`); the hidden-layers chip |
| sets a data range, autoscales | IGV | built as a command (`scale.set`); in the menu: planned with the legends design's lane (a numeric field) |
| clears a filter from the track | — | built: the "n filtered" chip, now the `seq:filtered` state entry |
| reorders tracks | JBrowse (drag handle), IGV | built: "Move up" and "Move down" in the menu (`track.move`); a drag from the header: planned, a later lane (a preview ending in one patch) |
| renames, removes, sets a display mode | IGV, JBrowse | no built-in: a label and the track list are the writer's, and density replaces display modes; a host's entry may add any of them over a patch |
| sees what a track is (About) | JBrowse | built: the pointer resting on the label shows the hover card (its label, its layers, its data and its states) |
| adds its own menu item | JBrowse plugins, igv.js, the current Seq\* | built: a registry entry; "color by this track" is one |
| a feature's own menu | JBrowse, IGV | built: `on: 'element'` entries in the menu opened on the element; the hover card shows none (Tensions) |
| header hover lights the track's items | the current Seq\* (EC-146) | refused by default: Tensions |
| collapses a group of tracks | ProtVista, the current Seq\* | no: the language has no track groups; a section is the group |

## Decision

### Two headers, one entry

- **A track's header** is its line in the label column: its label, its state chips and its action slots. A **section's
  header** is a line above the section's first track, in the same column. It holds the section's label, the
  hidden-tracks chip and the section's action slots. By default both are drawn by the layout's chrome and hit by its
  `hitChrome`, as the filter chip is today; a host may place DOM parts over the same boxes instead (below). Why drawn
  by default:
  - they scroll, clip, mirror and export with the picture;
  - one hit model serves the pointer, the keyboard cursor and the announcer;
  - no overlay has to follow the layout.

  The menu is DOM: a part behind a port (below).
- **A header entry** is a registry entry of kind `header-action` (kebab-case and singular as the registry's other kinds, `layer-type`; a `SeqExtension` lists them as `headerActions`), namespaced (`seq:hide`, `myapp:color-by`), listed in
  `spec.registry` or registered at run time with its undo ([viewer-composition.md](viewer-composition.md)). Its fields:

  ```ts
  interface HeaderEntry {
    readonly name: string;                         // namespaced; two entries under one name: a registration error
    readonly on: 'track' | 'section' | 'element';  // whose header, or an element's menu
    readonly kind?: 'action' | 'state';            // action (default): a slot and a menu item; state: a chip
    readonly icon?: IconName;                      // @molstar/ui's; an action without one: the menu only
    readonly label: string | ((f: HeaderFacts) => string);
    readonly when?: (f: HeaderFacts) => boolean;   // shown at all (absent: always)
    readonly command: string;                      // a command the viewer has; its `when` says enabled, or why not
    readonly args?: (f: HeaderFacts) => Record<string, unknown>;  // default { address: f.address }
    readonly text?: (f: HeaderFacts) => string | undefined;       // a state's chip text ("12 filtered")
    readonly order?: number;                       // slots by order, then name
    readonly slot?: boolean;                       // false: the menu only (default: true when it has an icon)
  }
  ```

  - **`HeaderFacts`** are plain, DOM-free data read from the resolved view and the last layout. They include the
    address and label of the track, section or element hit; the track's layers (type, address, hidden, filtered count,
    fitting schemes); its hidden tracks or layers; whether its data is pending; and its section's system. The functions
    are pure over the facts, run once a layout per shown header, never per frame or at registration.
  - **Shown or enabled**: the entry's `when` decides whether it shows at all. The command's `when` decides whether it
    runs now. A disabled entry is drawn muted, its reason in its card and in the menu, as the toolbar does. Why two:
    "this track has no color channel" hides the entry, while "nothing to undo yet" only greys it.
  - **The label and the icon default to the command's** (`title`, `icon`), so an entry over a known command is one line
    of data.
- **Built-ins on the same entries**, in `DefaultSeqViewerSpec().registry` (`BuiltInHeaderEntries`). None is drawn by
  code of its own:

  | Entry | On | Kind | Icon | Command | When |
  |---|---|---|---|---|---|
  | `seq:hide` | track | action, slot 1 | `eye-off` | `track.hide` | always |
  | `seq:color` | track | action, slot 2 | `palette` | `layer.channel` (its scheme param a submenu, below) | a layer with a color channel and more than one fitting scheme |
  | `seq:layers` | track | action, menu only | — | `layer.hide` and `layer.show` as checked items | more than one layer |
  | `seq:move-up`, `seq:move-down` | track | action, menu only | — | `track.move` | a sibling there |
  | `seq:filtered` | track | state | — (a text chip) | `layer.filter` (no value: cleared) | a layer with a user filter ("12 filtered", as today) |
  | `seq:hidden-layers` | track | state | — (a text chip) | `layer.show` (every hidden layer, one patch) | a hidden layer ("1 layer hidden") |
  | `seq:pending` | track | state, no command | `loader-circle` | — | data pending in view (no click; it says "loading") |
  | `seq:hidden-tracks` | section | state | — (a text chip) | `track.show` (every hidden track, one patch) | a hidden track ("2 hidden") |
  | `seq:reverse`, `seq:collapse` | section | action, menu only | — | `section.reverse`, `section.collapse` | a linear section; pieces to collapse |

- **A host adds, replaces or removes**:
  - adds by listing its entry (or `viewer.register`);
  - removes a built-in by leaving it out of the registry it spreads;
  - removes its own by the registration's undo.

  Registration is atomic and counted, and a conflict is one error before anything changes (two entries under one name;
  an entry naming a command the viewer lacks, `seq.viewer.header.command`). Why entries and not a config list: the
  registry is already how a host extends the viewer, and its rules (atomic, undoable, namespaced) come for free.
- **A command whose one param is a select is a submenu** of its options, each previewed while hovered when the command
  has a preview (`seq:color` through the theming behavior's `preview`). Why: no entry kind for submenus, and the
  palette, the menu and an agent see the same choices.

### The track header's geometry: fixed slots

- **The line**: the label column (`labelWidth`, `inset`), the track's label line (centred on its well, as built). At
  the column's right end the "⋯" slot; left of it:
  - the **action slots**: `header.slots` (2) columns of `header.slotPx` (20, an icon of 14). Their x is fixed by the
    column, not by the label, so every track's eye is in one column of the view, and the "⋯" is where a kebab menu is
    looked for (Mol\*'s order; the author, 2026-11-20: "Amend the design").
  - the **state chips**, ending where the action columns begin. On a track of two lines' height they move to a line of
    their own under the label (as the filter chip does today), and the label stays whole.
- **Actions show on hover** (`header.actions.show: hover`): while the pointer is over the track (its header or its
  data), or the keyboard cursor is on its header, its slots are drawn over a backdrop in the label column's
  background (the column has no well). The backdrop
  covers the label's tail and does not move it; where it cuts the label it makes room for an ellipsis and draws one,
  on the canvas and in the DOM parts (the novice pass, F13). Why:
  - a figure at rest shows names and states, not buttons;
  - drawing them is an emphasis pass over the label column, with no layout.

  States are always drawn: they are facts about the picture (ux-foundations section 9). Without the `headers` behavior a listed
  state entry's click runs its command over its layers; the minimal composition lists `seq:filtered` (its clear needs
  `layer.filter`, the legends behavior's, where a filter is set); other states are facts (the author, 2026-11-20:
  "Chip clears in minimal").
- **Slot assignment**: the action entries that apply, by `order` then name, take the slots. The rest, and every entry,
  slotted or not, are in the menu. "⋯" is drawn when the track has any entry. Why every entry in the menu: one place
  lists all, and a slot is a shortcut, never the only way.
- **Short tracks**: a track shorter than `header.minHeightPx` (16) draws no slots, and its states become marks (a 6 px
  accent squircle per state, its words in its card). Its entries are reached by the menu: a secondary click on it, or
  the menu key on its header. Why: an 8 px feature row cannot hold a 14 px icon, and a slot overlapping the next
  track's label is worse. **A narrow column** (one where the slots leave the chips no room, an 80 px column with actions
  holding one slot) draws its states as marks too (the author, 2026-11-20: "Marks when narrow"): all chips when every
  one fits, else all marks, never a chip beside a mark (a mix reads as two kinds of state).
- **Layouts**:
  - **tracks, text** (the first line's bands) **and grid** (its track's label) have a label column and take the slots;
  - **blocks** take them once a row, as their labels;
  - **rings** have no label column: their entries are reached through the menu (a secondary click on a ring's band,
    the menu key), their states said in the band's card;
  - **a rows track**: the label column is its tree and row headers, so its header line is its title row (shown when it
    has columns or any entry applies: `rows.titleRow: auto`).
- **Icons as paths**: the viewer reads the entries' `@molstar/ui` icons once (the Lucide markup's `path`, `line`,
  `circle`, `polyline`, `rect`) into path data for the draw lists, stroked in `msui-text-muted` (`msui-accent` while
  hovered). The core draws the paths it is given and imports nothing of `@molstar/ui` ([packages.md](packages.md)). An
  icon with another element is `seq.viewer.header.icon` and drawn as its label's first letter.

### Canvas slots or a host's DOM parts

- **One set of entries, two renderings** (decided: "can we support both?"). The `headers` behavior's `parts` param:
  `canvas` (the default) draws slots and chips as above; a `HeaderPartsPortFactory` hands them to a host's DOM. Why
  both: drawn parts scroll and export with the picture; DOM parts are focusable and styled by the host's design system.
  One entry list, one geometry and one command path, so the two never disagree.
- **With a parts port**, the chrome still lays out every box (`viewer.headers` unchanged) and draws nothing in the slot
  and chip boxes; the port gets `headerParts(viewer)` (DOM-free, published on change only): for each shown header, its
  geometry, whether its actions show now (`header.actions.show`), and each entry's item (name, kind, label, icon, chip
  text, enabled or why not, checked): the same items `headerMenu` builds. A press calls `parts.run(name, target)`, the
  slot click's path (`viewer.execute` with the entry's args). The port places parts at the boxes in the canvas's CSS
  px; it mounts actions only for headers whose actions show, so a 200-track view holds a few buttons, not 600.
- **`headerPartsDom()`** (`seqx-viewer`, on `@molstar/ui` buttons and icons): the reference parts, shipped (decided,
  the author, 2026-11-19: "Ship it"); a host passes its own. Why ship one: the Storybook and the twins need it, and `@molstar/ui` is already the menu's dependency.
- **Keyboard with DOM parts**: the cursor on a header moves DOM focus into its parts (one tab stop, roving: ← and →
  between them, Enter or Space runs, Esc back to the canvas with the cursor where it was); Tab is never captured
  (EC-155). Hits on the parts are the DOM's; the canvas sees no `header-slot` or `header-state` there.
- **Export is the canvas's in both**: chips drawn as states, slots never, so a figure does not depend on the rendering.

### The section header

- **Drawn when it has something to say** (`section.header: auto`): a hidden track, or a section entry with a slot that
  applies. Then it is a line of `header.sectionPx` (20) above the section's first track, holding the section's label
  (when it has one), the hidden-tracks chip and the section's slots, in the label column. `always` draws it for every
  section, and `none` never draws it (its chip's work then falls to the menu). Why `auto`: no line appears that says
  nothing, and records without hidden tracks do not change. ux-foundations section 9's "section headers are
  `@molstar/ui` headers (32px)" changes to this rule when agreed.
- **A section whose tracks are all hidden** keeps its header line (with "3 hidden"), so it never vanishes beyond reach.

### Chips: one rule for hidden and filtered

- **A count chip** (a squircle pill, ux-foundations section 5, in `msui-accent` text, as the filter chip) says a state
  of what is not shown: "12 filtered", "1 layer hidden", "2 hidden". Its click undoes the state in one patch:
  - the filter cleared;
  - the hidden layers shown;
  - every hidden track of the section shown, in order (`header.hiddenChip.click: all`).

  One undo step brings it back.
- **Its menu** (a secondary click on the chip, or the menu key) lists each hidden thing by its label as an item that
  shows that one (`track.show` with its address), then "Show all". Why click-all with a menu for one: the rows' chip and
  the filter chip already undo all at a click; one rule for every chip. `menu` makes the click open the list instead.
- **`track.show` and `layer.show` take several addresses**, one patch, as `layer.channel` takes several layers. A track
  hidden by a host's patch is counted like any other.

### The menu: a part behind the `headers` behavior's port

- **`headers`** (`HeadersProvider`, in the default spec) owns the slots' gestures, the menu and the geometry state. Its
  `menu` param is a `HeaderMenuPortFactory`, as the hover card's port is its behavior's ([ports.md](ports.md)). The
  default, `headerMenuOverlay()`, is `@molstar/ui`'s context menu in the canvas's element; a host passes its own. A
  viewer without the behavior draws headers with states only, with no slots and no menu.
- **The menu's data** (`headerMenu(viewer, target)`, DOM-free) has these groups:
  1. the element's entries (`on: 'element'`, when opened on one);
  2. the track's entries;
  3. the section's.

  Each item carries its label, its key when the command has one, enabled or why not, checked for the layers' items, and
  a submenu for a select param. The port shows it at a point, `show(menu, place)`, then `hide()`. A choice runs the
  command through `viewer.execute` with the entry's args; a hovered scheme previews it and Esc cancels.
- **Opened by**:
  - a click on "⋯" (drawn, or the DOM parts' button);
  - a secondary click, as the `headers` behavior's `contextClick` param says (decided: the composing app's): `everywhere`
    in the default spec, at what it hits (an element's entries first, then its track's and its section's); `header`,
    the label column only; `none`, the browser's everywhere;
  - the menu key or Shift+F10 at the keyboard cursor, once per ask (EC-156).

  **Shift+secondary click is always the browser's**, whatever the param. `everywhere` in the default spec is decided
  (the author, 2026-11-19: "Everywhere"). Why: IGV and
  JBrowse open a track's menu on a right-click over its data, and the canvas's browser menu ("save image") offers
  nothing a reader needs that Shift does not still reach.
- **Keyboard in the menu** is `@molstar/ui`'s (focus on the first item, arrows, Enter, Esc, a press outside). Closed, it
  gives the focus back to the canvas, where the cursor stays.

### Hit, cursor, announced

- **Hits**: a slot is `target: 'header-slot'`, a chip `target: 'header-state'` (the filter chip's `filter` becomes
  one), with the entry's `name` in the hit's `state`, the track's address in `layer`, the section's in `section`. A
  section header's line is `target: 'section-header'`, and a track's label line (the label, outside its slots and
  chips) `target: 'header-label'` (the author, 2026-11-20: "Add a label hit"): its hover shows the track's card and
  drives which header's actions show (no emphasis field of its own); on touch, a tap on it shows the slots, and a tap
  on a slot then runs it. A slot is hit only while its actions show. The bindings take them as any target ([ux-foundations.md](ux-foundations.md)
  section 8): click on a slot or a chip runs its entry, and the hover shows its card. They are never selection or
  navigation targets: a click there never clears the selection, and a drag from there neither selects nor pans.
- **The keyboard cursor reaches headers**, as it reaches a row's header:
  - `cursor.prev` at a section's first shown position (display order) goes onto the header of the cursor's track (the
    cursor gains `part: 'header'` and the track's address);
  - ↑ and ↓ step between the section's track headers and its header line, then on to the next section's;
  - → goes back to the data at the first shown position;
  - Enter or the menu key opens the menu.

  The track's actions are drawn while it is there. Why the menu and not a cursor through the slots: the menu is real
  DOM in the tab order of the moment, readable by a screen reader, with every entry in it. Tab is never captured
  (EC-155).
- **Said** (`describe.ts`, one description for the live region and the card):
  - a resting pointer on a slot says the entry's label and key ("Hide track");
  - a disabled slot also says why;
  - on a chip, its text and what a click does ("2 hidden: click shows them");
  - the cursor on a header says the track's label, its states, and "Enter for actions";
  - a run says its command's message ("Conservation hidden. 1 hidden in p53").

  The card for a slot or a chip is the hover card with a title alone. The same words are given to a host's card port.
- **After a hide by keyboard** the cursor moves to the next track's header (else the previous one's, else the section
  header's chip), so it never stands on what is gone.

### Touch

- With no hover, a tap on a track's label column shows its slots (as the pointer resting there would), and a tap on a
  slot runs it. A long press opens the menu (`header.menu.longPressMs`, 500, the platform's when it has one). Why: a tap
  on the label selects nothing today, so it is free, and a long press is the touch platforms' context menu.

### The geometry, for a host's parts and widgets

- **`viewer.headers`** is a reactive state: `HeaderGeometry[]`. Each one holds:
  - its target (the `section`, and the `track` address or none for a section header);
  - its line's box;
  - its label's box;
  - each slot's and chip's box, with its entry's name;
  - the well's box;
  - `shown` (on the page, not clipped).

  All are in the canvas's CSS px, as the hits are. It is published on change only, compared structurally
  (§relevant-rerenders): after a layout that moved a header, and during a text section's scroll, once a frame. Why
  publish it: the DOM parts (above) and a host placing a checkbox, a drag handle or a sparkline over a header read
  where, without the layout's arithmetic.
- **A widget's slot**, an entry of `kind: 'widget'` that reserves a box a host's DOM fills, is planned: the author's
  "later", when a host asks. The geometry already gives it a place: a widget is an entry whose slot draws nothing and
  whose box is published.

### Undo, history, export

- The built-ins are patches, one undo step each, as their commands' kinds say (`CommandChangeKinds`); the color's
  preview ends in one patch or none. A host's command declares its own kind (structure-syncing's "color by" changes no
  Seq\* state).
- **Export**: action slots and "⋯" are never exported (chrome for a pointer). State chips are exported as drawn
  (`header.states.export: include`, as the filter chip is today), and `none` leaves them out of a figure. The section
  header line is exported when it is drawn. With DOM parts the export still draws the chips on the canvas.

## The interaction matrix's cells

Rows this design fills or changes, with the columns of [ux-foundations.md](ux-foundations.md) section 8.

| | select | extend | add / remove | focus | back | zoom in / out | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|---|
| **a track or layer** | no: a track is not an element; its items by a legend entry or `address.select` (a host's entry may offer "Select items") | no: as select | no: as select | no: a track's extent is its section's; focus an item | Alt+Left (the view's) | as built: Cmd/Ctrl+wheel over the label column, about the first visible | as built | the `seq:hide` slot (`track.hide`); a layer from "Layers ▸" (`layer.hide`); shown again by the chips (`track.show`, `layer.show`) | "Color ▸" (`layer.channel`, previewed); "Move up/down" (`track.move`); the filter cleared by its chip (`layer.filter`); a host's entry ("color by this track") |
| **a track header** (its label line) | no: chrome; a click there selects nothing and clears nothing | — | — | the cursor onto it (← at the first position); Enter opens the menu | Alt+Left | as a track | as a track | shows its slots on hover, a tap, or the cursor | the menu: every entry; a slot runs one |
| **a slot or a chip** (drawn or a DOM part) | click (or Enter on a focused part) runs its entry; a chip's undoes its state (one patch) | — | — | — | Alt+Left; undo restores a chip's state | — | — | a hidden thing comes back by its chip | secondary click on a chip: the list, each shown alone |
| **a section header** (drawn when it has something to say) | no | — | — | — | Alt+Left | — | — | "n hidden" (`track.show`, all; its menu one) | its menu: `section.reverse`, `section.collapse`, a host's section entries |
| **an element** | as the grammar | as the grammar | as the grammar | as the grammar | as the grammar | — | — | no: hide its track | secondary click (under `contextClick: everywhere`) or the menu key: its `on: 'element'` entries, then its track's |
| **a ring's band** | as the grammar | as the grammar | as the grammar | as the grammar | as the grammar | as its layout | as its layout | the menu's `seq:hide` | the menu (no label column to slot) |

## Configuration

Items are `ViewerConfig` unless the table says otherwise.

| Item | Default | Why | Alternatives |
|---|---|---|---|
| registry `header-action` | `BuiltInHeaderEntries` in the default spec | the built-ins are entries like a host's | a host's list; any built-in left out |
| `headers` behavior | listed in the default spec | slots and the menu work out of the box | left out: states only, no slots, no menu |
| `headers.menu` (behavior param) | `headerMenuOverlay()` (`@molstar/ui`'s context menu) | as the hover card's port | a host's adapter; none: the menu's state moves and nothing is drawn |
| `headers.parts` (behavior param) | `canvas` | drawn parts scroll and export with the picture | a `HeaderPartsPortFactory`: `headerPartsDom()` or a host's |
| `headers.contextClick` (behavior param) | `everywhere` in the default spec | IGV's and JBrowse's right-click over a track; the composing app's call (decided) | `header` (the label column), `none` (the browser's) |
| `header.actions.show` | `hover` (the pointer over the track, the cursor on its header, a tap) | a figure at rest shows names and states, not buttons | `always`; `never` (the menu only) |
| `header.slots` | 2 | the eye and the color fit beside a label in 112 px | 0 to 4 |
| `header.slotPx` | 20 (an icon of 14) | a target of 20 px on a 24 px row | any |
| `header.minHeightPx` | 16 | below it an icon overlaps the next track | any |
| `header.menu.longPressMs` | 500 (the platform's when it says) | touch's context menu | any |
| `header.hiddenChip.click` | `all` (one patch) | one rule with the rows' and the filter chip | `menu` (the list) |
| `header.states.export` | `include` | as the filter chip is exported today | `none` |
| `section.header` | `auto` (when hidden tracks or a slotted section entry) | no empty line; records without hidden tracks unchanged | `always`, `none` |
| `header.sectionPx` | 20 | a line of the label font | any |
| `rows.titleRow` | `auto` (columns, or an entry applies) | a rows track's header line | `always`, `none` |
| Shift+secondary click | fixed: the browser's | the page's menu always reachable | — |
| a slot's and a chip's hit targets | fixed: `header-slot`, `header-state`, `section-header` | the bindings name them | — |
| entries' functions | fixed: pure over `HeaderFacts`, once a layout | no DOM, no per-frame work | — |

## Tensions

- **Drawn slots against DOM buttons.** Decided: both, drawn by default. Drawn slots are not focusable, so under
  `canvas` a screen reader reaches the entries only through the cursor and the menu; DOM parts are focusable but follow
  the layout from the geometry, a frame behind during a scroll. Standing: which a host picks is its trade.
- **Two renderings to keep alike.** A DOM part could show what the canvas would not (a stale `enabled`). Held by one
  item list (`headerParts` and `headerMenu` build the same items) and a test that both renderings show the same items
  for every example view. Standing.
- **Functions in entries against data.** A data-only `when` (layer types, a param present) could be printed, checked and
  sent to an agent. Functions cover what a host knows and Seq\* does not (a structure bound, a licence). Recommended:
  pure functions over `HeaderFacts`, never run at registration (decided). Standing.
- **Actions on hover against discoverability.** A novice may not find hover-only slots, and touch has no hover.
  Mitigations: "⋯" and the states are always drawn where something is stated, a tap shows the slots, the keys panel
  lists the menu key, and `always` is one item. Standing until a novice pass says otherwise.
- **A chip click shows all, against picking one.** A user who hid three tracks and wants one back clicks twice
  (show all, hide two) or opens the menu. All, for one rule across chips (decided); the menu shows one. Standing.
- **The secondary click is Seq\*'s in the default spec.** Decided: the `headers` behavior's param, Shift+secondary click
  always the browser's; `everywhere` in the default spec (decided). Standing: a page's own menu needs Shift there.
- **The section header line moves content.** Hiding a section's first hidden track removes its row and adds a 20 px
  line, so the section shrinks by less than the track. Accepted: `auto` keeps the line absent in every view that hides
  nothing. Standing.
- **Element actions in the menu, not the card.** The current Seq\* put item actions in the hover card (and
  structure-syncing's "show in 3D" is now an `on: 'element'` entry). The card is `aria-hidden` and goes when the pointer moves to reach it.
  Recommended: the menu, opened on the element. A host's card port may still draw actions from `headerMenu` data.
  Standing.
- **A header hover lighting the track's items** (EC-146). It reads well on a feature track, but it floods linked
  sections and a linked structure viewer with every item on a pass of the pointer down the label column. Refused by default.
  A host may bind a header's hover to `address.hover` of the track. Standing.
- **Template instances.** `track.hide` on an expanded instance has no node to patch. Recommended: its `when` says "from
  a template: hide the template" and the menu offers the template's hide (all instances). Hiding one instance (an
  exclusion in the template's `filter`) is planned with the language's next template change. Standing.

## Self-check

The claim: **this design covers what a user does with a track's header.**

1. **Assume it holds.** Every action on a track or a section is a registry entry naming a command. It is drawn in fixed
   slots on hover, listed in one menu, reached by pointer, keyboard and touch, and said. Every state of what is not
   shown is a count chip whose click undoes it. A host adds and removes entries, and reads the headers' boxes.
2. **What it rests on, negated.**
   - (a) *The header is there to click.* Negated: a track too short, a ring with no label column, a rows track whose
     column is its tree, an export, a touch screen.
   - (b) *Hiding is reached from the header.* Inverse: showing again, when the header is gone with its track.
   - (c) *One track at a time.* Plural: every track of a section hidden; several hidden, one wanted back; a template's
     many instances.
   - (d) *The viewer owns the entries.* Symmetric: the host adds, removes, replaces a built-in, or puts its own widget
     on the header.
   - (e) *The pointer.* Negated: the keyboard, a screen reader, an agent.
   - (f) *Actions change the spec.* Negated: a host's action changing another viewer's state ("color by"), an action
     with no command (pending).
   - (g) *The track is shown.* Negated: a link, a bookmark or a reveal names an item on a hidden track.
3. **The negated world.**
   - (a) gives the minimum height with marks and the menu, rings through the menu, the rows' title row, slots never
     exported, and the tap.
   - (b) gives the chips, and the section header line drawn when a track is hidden.
   - (c) gives the section header kept when all are hidden, click-all with a menu for one, `track.show` over several
     addresses, and the template's hide.
   - (d) gives the registry, built-ins as entries, the geometry state, and a host's DOM parts over it.
   - (e) gives the cursor's header part, the menu as DOM, focusable DOM parts when a host picks them, the words, a cursor kept off what was hidden, and entries
     naming commands (an agent runs `track.hide` as the slot does).
   - (f) gives commands of any change kind, and `state` entries with no command.
   - (g) the design did not cover. A focus or a reveal of an item on a hidden track draws nothing and says nothing.
4. **Held both.**
   - Integrated in the Decision: (a), (b), (c) except one instance of a template, (d) except the widget's slot, (e),
     (f).
   - Added by this check: (g). A jump whose target lies only on hidden tracks says so ("on a hidden track:
     Conservation; its section's chip shows it"), and the section's chip is drawn with the hover's outline for the
     announcer's rest. Why not show it: a hide is the writer's patch, and a reveal never edits the spec
     ([addressing.md](addressing.md)'s "into a hidden stretch, its break and the reason", the same rule).
   - Planned: the widget's slot (when a host asks); drag to reorder from the header (a later lane); one template
     instance hidden (the language's next template change); the data range in the menu (the legends design's lane).
   - Standing tensions: the header out of the tab order under `canvas`; two renderings kept alike; hover-only actions
     against discoverability; functions in entries.

## Consequences

The fourth slice builds it, beside the structure lanes. ST3's "color by this track" story needs HP2.

- **`seqx-core`**:
  - The layout's chrome draws the header line: slots as given (path data, the entry's name, enabled), chips (text,
    box), the section header line under `section.header`, and the short-track marks. Hits are `header-slot`,
    `header-state`, `section-header`. The filter chip becomes the `seq:filtered` state through the same path.
  - `TrackLayout.filter` is replaced by `TrackLayout.header` (slots, states, line). `SectionLayout.header` is the
    section's line.
  - Tests:
    - property tests: slots never overlap each other, a chip or the "⋯"; a slot's x is independent of the label; every
      hidden track is counted once in its section's chip; a chip's click then undo is the identity on the spec;
    - drawing against the data: for each view of the examples, every shown track has its header line; a hovered one has
      a slot for each slotted entry that applies, and each is hit as its entry; every state an entry's `text` gives is
      drawn and hit;
    - records: records with a filter chip change once (it moves before the action columns); nothing else changes under
      `auto`.
  - `_spec`: Render, Tracks layout (the header line), Filters (the chip through the state entry), Lanes and Hits (the
    targets).
- **`seqx-viewer`**:
  - the `header-action` registry kind with its checks;
  - `BuiltInHeaderEntries`;
  - `HeaderFacts`, made once a layout;
  - icons read into paths;
  - the `headers` behavior: hover and tap reveal, slot and chip clicks, the menu and its port (`headerMenu`,
    `headerMenuOverlay()`), the parts and their port (`headerParts`, `headerPartsDom()`), the long press, the secondary
    click under `contextClick` with Shift left to the browser;
  - the cursor's header part on every track (`cursor.prev` at the first position, ↑ ↓ across headers);
  - `describe.ts`'s header words;
  - `viewer.headers`;
  - `track.show` and `layer.show` over several addresses;
  - the config items above.
  - Tests:
    - each built-in through its slot, its menu item and its command alike (one patch, one undo);
    - a host's entry added and removed with its undo, and a conflict reported atomically;
    - a disabled entry says why;
    - the menu's data in Node, and the overlay in jsdom (focus in, Esc, focus back to the canvas);
    - the cursor onto a header and back, and off a hidden track;
    - a secondary click and the menu key once per ask (EC-156); each `contextClick` value; Shift+secondary click never
      taken;
    - the parts: `headerParts` and `headerMenu` give the same items for every example view; in jsdom, `headerPartsDom()`
      places a part at each box, a press runs the slot's command (one patch), roving focus and Esc back to the canvas,
      actions mounted only where they show; the export alike under both;
    - the geometry published only on change (a hover publishes nothing);
    - a quick representative: a pointer move over 200 tracks' label column allocates nothing and lays nothing out (the
      slots are an emphasis redraw).
  - `_spec`: Input (headers, the cursor's header part), Behaviors (`headers`), Ports (Header menu), Commands, Bindings
    (the targets), Configuration.
- **`seqx-examples`**: the p53 views with a hidden track and a hidden layer, and a filtered track (fixtures for the
  stories and the records).
- **Storybook**:
  - `Features/Track headers`: the built-ins on p53, drawn and as DOM parts (a story arg), a host's entry added and a built-in removed, the hidden-tracks chip
    with every track of a section hidden, short tracks, a ring through its menu, and a host's DOM widget placed from
    `viewer.headers`;
  - "Color by this track" with the fake structure viewer (shared with structure-syncing's ST3).

  Each has a headless twin (a fake clock, the menu's data, the patches, the geometry), and is looked at in light and
  dark. There is a `Seq* (next)/Design/Track headers` page.
- **Lanes** (named for the fourth slice's plan):
  - HP1 (strong): core's header line, hits and chips;
  - HP2 (strong, after HP1): the viewer's behavior, registry kind, menu port, parts port and `headerPartsDom()`, cursor
    and words;
  - HP3 (fast, after HP2 and ST2): the stories, twins and the Design page.
- **Designs that change when agreed**:
  - [ux-foundations.md](ux-foundations.md): section 9 (section headers under `section.header`, the slots, the menu);
    section 8's matrix (the rows above);
  - [extensibility.md](extensibility.md): track header actions open as a registry kind; context menu items as
    `on: 'element'` entries;
  - [ports.md](ports.md): the header menu's and the header parts' ports on the `headers` behavior;
  - [themes.md](themes.md): the "n filtered" chip as the `seq:filtered` entry; "the track menu's Color" as `seq:color`.

## Decided by the author

2026-11-19; each answer is folded into the sections named.

1. **Changed**: both renderings, drawn by default and a host's DOM parts over the geometry, one set of entries: Two
   headers, Canvas slots or a host's DOM parts, the geometry.
2. Actions on hover, states always. 3. The section header line under `section.header: auto`. 4. A chip's click shows
   all. 6. Pure functions over `HeaderFacts`. 7. The cursor onto headers. 8. The geometry now, a widget's slot later.
5. **Changed**: the secondary click a `headers` behavior param, Shift+secondary click the browser's: The menu.
9. The default spec's `headers.contextClick: everywhere` ("Everywhere"). 10. `headerPartsDom()` shipped ("Ship it").
