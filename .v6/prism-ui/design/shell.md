# Shell, layout, and the UI spec

The shell is the frame around the viewport: the rail and its tab panel, optional docks, the chrome over the canvas,
overlay cards, the command menu, and the status column. It is Prism's own layout (core `PluginContext.layout`,
`PluginLayout`, and the classic five-region `Layout` stay the classic UI's), and it is where contributions (built-in and
extension UI) are placed.

## 1. Principles

- **Canvas first.** The viewport fills the shell; persistent chrome is the exception.
- **Summoned, not resident.** Short tasks (export, canvas settings, animation, remote states) open as cards in the
  panel's place and close when done. Docked UI is for what the user keeps looking at while working.
- **One home for persistent UI.** A tabbed panel on the right holds Home, Structure, Selection, Snapshots, State tree,
  Log, Settings, and extension tabs. Left, top, and bottom docks exist but are optional; the top one holds the sequence,
  the others are empty by default; apps can place tabs there.
- **Appear when relevant.** Home shows when nothing is loaded; the camera and rendering tools only once something is;
  trajectory controls only for multi-model trajectories; snapshot controls only when there are snapshots; the selection
  bar only while something is selected; the Sequence button only while there is a sequence and its card is hidden; the
  log as one fading line.
- **Lean views.** Commands carry the operations ([commands.md](commands.md)); views show state and the common actions,
  not a form for everything.
- **Same spirit, not the same screens.** The workflow (load, style, select, measure, save) and the capabilities stay;
  the layout, density, and icons are redesigned.
- **Modular and embeddable.** Every tab, overlay, toolbar item, and section is a contribution; the same pieces work in
  the full app, a bare viewport, or single elements in a host page.

## 2. Structure

> **2026-10-15:** the pinned panel and the rail no longer shrink the canvas; the camera's framing keeps the scene in
> the visible area instead ([../spec/README.md](../spec/README.md) §3.5).

```
┌──┬────────────────────────────────────────────────────────┬───────────────────┬──┐
│  │ ┌ 1TQN ▾ Chain ▾ Cytochrome P450 3A4 ▾ Chain A ▾ ⌃ ┐   │                   │◆ │ head: Home (the Mol* mark)
│ l│ │ MALYGTHSHGLFKKLGIPGPTPLPFLGNILSYHKGFCMFDMECHKK…  │   │  panel: the active│──│
│ e│ └──────────────────────────────────────────────────┘   │  tab, or a card in│⬡ │ strip: Structure, Selection,
│ f│ ┌ Active site    ●○○○  ‹ › ▶ ┐                   [↖]   │  its place (full  │◌ │ Snapshots, State tree
│ t│ │ markdown description…      │ Picks ⬡ Residue ▾ [◎]   │  height; header:  │▤ │
│  │ └────────────────────────────┘                   ┈┈┈   │  title, pin, ✕)   │≡ │
│ (│ ┌ 1TQN  4 / 20  ━━●━━━━━━  ▶ ┐                   [⌖]   │                   │──│
│ o│ └────────────────────────────┘                   [✧]   │                   │▶ │ scene: Animation,
│ p│                                                  [⤢]   │                   │◐ │ Canvas settings
│ t│                               ┌ GLU 213      ┐         │                   │──│
│ .│             canvas       ⌖ ───┤ Chain A·1TQN │         │                   │⤓ │ export: Export,
│ )│                               └──────────────┘         │                   │☁ │ Remote states (when on)
│  │                                                        │                   │──│
│  │ ✓ Loaded 1TQN (3,624 atoms)                            │                   │▦ │ foot: History (+ log),
│  │ ⟳ Computing surface ━━━━── ✕                           │                   │⚙ │ (Assistant), Settings,
│  │ ✦ Show the active site · Working…  12 residues ✎ ⌏ ⌖ ✕ │                   │☀ │ theme switch,
│  │                                                        │                   │⌘K│ command menu
│  ├────────────────────────────────────────────────────────┤                   │  │
│  │ bottom dock (optional, canvas width)                   │                   │  │
└──┴────────────────────────────────────────────────────────┴───────────────────┴──┘
```

```
  canvas toolbar (the stack), a row in the bottom-right corner: mode (Interact ↖, Select ◎ with the picking levels),
                canvas (Reset camera ⌖, Illumination & XR ✧, Expand ⤢; each with a menu on hover, opening upward), view
                (Sequence: its menu, or only while its card is hidden)
  bottom: the status column on the left (message, tasks, the assistant's turn); on the right the corner: the selection
          bar to the left of the canvas toolbar, on its row
  also: the command menu (modal), overlay cards (in the panel's place), the busy overlay, the drop overlay
```

- **Right panel, floating.** By default the right dock floats over the canvas instead of taking space beside it
  (`DockState.floating`, on for `right` and `top`): the canvas gets the full width, the rail is docked to the right
  edge, and the active tab's content is a glass card beside it (more opaque than the small controls, for reading).
  Collapsed, only the rail stays. The panel runs the full height; its header is the tab's title (or its `heading`:
  Home's is the wordmark), the pin, and ✕ (collapses the dock). Its inner edge drags to resize (`setSize`, at least 180
  px; a double-click goes back to the default 340). `setFloating(dock, false)` docks it beside the canvas instead (the
  rail's buttons then sit in a frame at the canvas's corner); compact mode makes it a bottom sheet. Everything else on
  the right of the canvas (the sequence card's end) follows the panel's edge (the toolbar and the selection in the corner stay at the rail),
  without easing (`--msp-canvas-right-inset`). One margin frames the canvas (`--msp-canvas-edge`, 8px): every floating
  thing keeps it from the canvas edges, the docked rail, the open panel, the sequence card, and the cards stacked beside
  it.
- **Rail.** A full-height column at the right edge, no margin, buttons edge to edge, a hairline between sections, no
  labels (tooltips on the left name each button). Top to bottom: Home as the Mol* mark (`place: 'head'`, a section of
  its own); the tab strip (Structure, Selection, Snapshots, State tree, by `order`); the rail's toolbar sections, each
  named by its items' `group` and placed in the order of its first item: `scene` (Animation, Canvas settings) and
  `export` (Export, and Remote states when the spec turns them on); and at its foot, pushed down with a divider above,
  the tabs that are not working panels (`place: 'foot'`: Log, the Assistant with the agent module, Settings), then, a
  shade apart, the items without a group: the theme switch and the command menu (the main way in, its ⌘K under the icon:
  a toolbar item's `caption`, joined to its title in the tooltip and accessible name). Rail buttons are toolbar items
  (`CanvasToolbarItem`, `place: 'rail'` by default): a click runs `run` or toggles its `overlay` anchored to the button;
  it shows pressed while its card is open or `pressed` holds, shows only while `visible` holds, and may change its face
  (`face`: the theme switch shows the current theme's icon and name). Which buttons show is the app's choice
  (`ui.layout.buttons` in the spec, `ShellOptions.buttons`, `ShellModel.setButtons` at runtime; ids not listed show).
  Extensions do **not** add rail buttons for their features (the rail would grow with every extension): they contribute
  sections to a panel's slot instead (Tabs and sections, below).
- **Rail on a short shell.** The rail steps down with its own height (it is a size container, `container-type: size`;
  the queries are CSS, nothing is measured for them): **normal** from 620px (buttons 40px, icons 20px; the command menu
  with its caption 46px), **compact** from 530px (34px, icons 18px, the captioned 40px), **small** below (30px, icons
  16px, the captioned 36px; tooltips and accessible names are the same at every step). The steps are `--msp-rail-*`
  variables set on the rail's children in `dock.module.css`, which the buttons (`.tab`, `.tool`) read; the numbers are
  in `RailSteps` (`views/shell/rail.ts`) too, and `_test/rail.test.ts` keeps the CSS and the table equal. When even the
  small step is too tall, the lowest-priority buttons are left out and a **More** button (the ellipsis, below the
  toolbar sections, above the foot) opens a menu beside the rail, to its left, its bottom edge level with the button
  (outside the rail's clipping, as the canvas menus are): a row for each left-out button, in the rail's order, running
  what the button runs (a tab shows, an overlay opens anchored at More). Priority is a field on the registration
  (`ShellTab.priority`, `CanvasToolbarItem.priority`): the lower the number, the longer it stays, defaulting to the
  button's position in the rail from the top, so the foot goes first; the command menu says -1 (the last to go) and the
  theme switch 100 (the first). What fits is computed (`fitRail`) from the rail's height, the step's button sizes and
  the hairlines, only when the height or the set of buttons changes: one `ResizeObserver` on the rail (`RailLayout`),
  its callback sets a state only when the rounded height differs, and the left-out buttons are not rendered (their
  section's hairline goes with the last of them).
- **Compact rail.** An option (`ShellOptions.compactRail`, off by default; `ShellModel.setCompactRail`; in the stories the
  `compact` arg, `?args=compact:true`; not the narrow layout's compact `mode`): the rail column is hidden and the canvas
  runs to the viewport's edge. A button with the rail's head glyph (the Mol* triangle), a card of its own, sits at the
  bottom of the rail's column (centered on the rail's icons, at the viewport's right edge, level with the corner row)
  and shows and hides the column (`ShellModel.railOpen`, `toggleRail`; the command `shell.toggle-rail` does the same
  and can be given a key). The rail is an overlay: it slides in over the canvas's edge and takes no width, so opening
  it moves nothing (`--msp-canvas-end` stays 0; the corner toolbar keeps clear of the column whether it shows or not,
  to the left of the toggle), and it ends above the corner row as the floating panel does (the toggle and the toolbar
  stay visible and clickable); panels and cards open to its left as they do beside the normal rail, at whatever size
  step fits the shorter column. Hidden, the column is inert and paused (`activity`, §"Hidden UI"), slid out and faded
  (transform and opacity only). Hiding closes the overlay card and the right panel, floating or pinned (a pinned panel
  without its rail could not switch tabs; its pin is kept for when it opens again; a pinned panel beside the canvas
  does take its width). Only where the rail floats: the narrow layout and the docked one have none.
- **Theme switch.** At the rail's foot by Settings: a click opens the themes beside the rail (Light, Dark, Blue, High
  contrast, Auto (system): `UIThemes`; the rail's menu, `RailLayout.toggleMenu`, as More's), the current one checked;
  under More, a click goes to the next (`ThemeModel.next`). Its icon and tooltip are the current one's. `ThemeModel` (`shell/theme.ts`) keeps
  `data-msp-theme` on the element the UI is themed by. Settings › Interface offers the same choice. The rail, the panel
  and the cards follow the UI theme; the view's own chrome follows the canvas (Glass and canvas theme, below).
- **Command menu.** Opened by the rail's last button, ⌘K (Ctrl+K) or ⇧⌘P anywhere, Home's "⌘K Search every command"
  line, and the Structure tab's Preset and Styles, which open it inside their group (`CommandMenu`'s `group`, teaching
  the command menu as they are used). It is modal, over the whole shell (`ShellView`'s `commandMenu`); whether it is
  open is the host's state, not the shell model's. With the agent module, a bar under its list takes requests for the
  assistant ([agent.md](agent.md) §5). A command whose work has a place in the UI (`CommandDefinition.view`: Export ›
  Models, Selection › Measure, Snapshots) offers it in its form, beside Run (⌥↵): the menu closes and the place opens,
  for what the form leaves out. The rest is in [commands.md](commands.md).
- **Tabs and sections.** A tab (`ShellTab`: `title` for its tooltip, an optional `heading`, `icon`, `dock`, `order`,
  `place`, `compactBar`, `render`) goes into its dock by `order` among the tabs already there. A tab strip click
  activates a tab, or collapses the dock when it is the active one (`select`); `show(tab)` activates and expands. Panels
  and cards that hold several things are lists of **sections** (`ShellSection`: `slot`, `title`, `keywords`, `order`,
  `when`, `render`), shown by `SectionList` (`views/shell/section-list`) as blocks that fold: the Export, Canvas
  settings and Animation cards, the Selection tab, an inspector's `inspector:<kind>` slot. Extensions add sections, not
  buttons: geometry and model export are sections of the Export card (`export` slot), next to the built-in Image (with
  Animation as its choice, `export.picture`), Camera and Session & state; the card's sections start folded but Image
  (`SectionList`'s `collapsed`, `ShellSection.open`). Each section can be opened by a command at it: `show('selection', 'measure')` for a tab,
  `openOverlay('export', anchor, 'animation')` for a card; that section opens and scrolls into view and the others fold,
  so the command menu finds "mp4" (an alias of `export.animation`) and lands on it. The request is the tab's `tabFocus`
  until the tab is shown plainly or picked from the strip, which forget it; the section it opened is remembered as open.
  A section's `keywords` are declared for the command menu and a section filter; nothing reads them yet. Otherwise a
  section opens as the user last left it (`sectionOpen` / `setSectionOpen` by slot and section, kept while the shell
  lives, across tab switches and re-renders), else open by default, or folded where the list reads as a menu of what can
  be done (`SectionList`'s `collapsed`: the Selection tab).
- **Pages.** For what needs more room than a panel or a card (`ShellPage`: the agent's Inspector), a page over the
  canvas, beside the rail and the open panel (they stay; a panel opened meanwhile makes room): a card like the others,
  only full size (the app shows at its edge, its corners rounded), opaque (dense content would not read over glass), a
  header with its title, its own header content (tabs, facts) and ✕; Esc closes it and focus goes back. Opened by a
  command or a link at a focus (a tab of the page, an item); opening one closes an open card. The app is unchanged
  underneath.
- **Overlays.** Cards for short tasks (`ShellOverlay`: Export, Canvas settings, Animation, Remote states; §4.2), one at
  a time, opened by a rail button or a command. Every card takes the right panel's place, full height, whether the panel
  is open or not (one target for everything that opens on the right). Over the open panel it swaps in without a
  transition and the panel steps aside (no glass on glass, no tab shown selected); a tab click closes the card and shows
  that tab. It closes with Esc, its ✕, or a pointer down outside it and its button, but not on the canvas (the view is
  turned and zoomed with the card open: setting up an export); another card's button or a tab swaps
  it directly, without the panel underneath blinking for a frame. Cards do not block the canvas. Panels and cards share
  one header (title, pin, ✕).
- **Hidden panels.** A panel is mounted only while its tab is the dock's active one (`Show` in `DockView`: closed, or
  another tab picked, means gone, and what stays mounted by choice is paused). Whatever the shell hides while keeping it
  mounted is **inert and paused** ([styling.md](styling.md) §5, `activity` in `dom/control-flow`): the panel behind an
  overlay card (it only steps aside visually: `inert` on its card, its scope inactive, so the views under it follow no
  model until the card closes, and it returns as it was), the body of a folded section (`SectionList`: `inert`, scope
  inactive until it opens, so a list that starts folded binds nothing until opened), and a panel during its closing
  animation (`inert` from the first frame of the close, so it takes no pointer while it slides out). A dock's content
  is the card's `section[role=tabpanel]`; the `hidden` or class that does the hiding is bound outside the paused scope.
  Tooltips ignore elements in an inert part. A view that works on its own (a ResizeObserver, a rAF loop) pauses it with
  `onActive` ([rendering.md](rendering.md) §4.2).
- **Pinned panel.** The pin in every panel's and card's header (`setPinned('right', true)`, `DockState.pinned`) keeps
  the floating layout but moves the rail and the open panel beside the canvas: the canvas and its chrome end where they
  begin (`--msp-canvas-pinned`), the panel turns into a flat full-height column like the rail (hairlines on both sides),
  and cards take the same place and look. For workflows that keep a panel open the whole time and want the whole
  molecule in view beside it. Pinned, the panel and its cards are part of the workspace: Escape and clicks outside leave
  them; only their ✕, another card or a tab change them.
- **Escape** closes the innermost open thing, one per press: the assistant's turn card first, to its line (a capture
  handler; not when the focus is in a field, list, menu or card); then an open list or menu (a select's list, which
  handles it on itself; a canvas stack menu, the picking level); then the open card; then the open panel (`ShellView`
  listens on the window, after the document's handlers, and skips what they handled). The menus' and the card's handlers
  are bubble-phase document listeners that skip a handled event, so their order is the order they were added: the menus
  with the shell, the card when it opens, so a menu goes first. A search field empties first (`TextInput`'s `search`),
  and empty lets Escape through. The panel never closes from a field; a card does (its handler does not check the
  focus). Pinned cards and panels stay. The command menu handles its own Escape (from a sentence back to the list, then
  out).
- **Floating docks.** The top dock (sequence) floats too: a glass card along the top of the canvas, between the left
  edge and what floats on the right; hidden, it takes no room (the canvas stack's Sequence button brings it back). The
  canvas runs under everything, so the sequence reads as part of the visualization. Playback and the canvas stack sit
  below it (`--msp-canvas-top-inset`). The sequence card's header has its choices (Mol*'s sequence controls: structure,
  mode (Chain, Polymers, Everything), entity and chain in Chain mode) and its collapse.
  A top tab with `frame: 'none'` (the sequence in the viewer) is a **band** instead of a card (`FloatingDock`'s `band`
  class): it spans the canvas area (the left edge to where the canvas ends: the rail, or a pinned panel; a floating panel
  lies over it as over the canvas), over the canvas' top, in the canvas' own color (the ink) behind its lines, fading to
  transparent below them (`BandFadePx`, `BandDissolvePx`, a gradient, no filter); the faded part lets the pointer through.
  The canvas keeps its full height and ends at the same place (`--msp-canvas-end`); the chrome below it
  (`--msp-canvas-top-inset`) starts under the band's lines. The viewer asks the camera's default framing to keep the
  band clear.
- **Optional docks.** `left`, `top`, `bottom` take tabs placed by the app spec (mesoscale explorer puts its entity tree
  on the left) or moved there (`move`). The sequence is the only default occupant of `top`. Side docks run the full
  height; top and bottom docks span the canvas width only. A docked (not floating) top or bottom dock with one tab has
  no tab strip: its content carries its own header, and collapsed it is a slim bar with the tab's name that expands on
  click.
- **Canvas stack (the canvas toolbar).** What acts on the view sits on the canvas itself, in its **bottom-right corner**
  (anchored beside the rail, or beside a pinned panel's border: it does not move for a floating panel, which is transient
  and lies over the canvas and **stops above the corner**, leaving the toolbar's row and the edge gap, `--msp-corner-room`;
  an overlay card in the panel's place likewise), as **one glass row**: the mode (Interact, Select; §4.6), then the
  toolbar items placed on the canvas (`place: 'canvas'`) in their sections, separated by hairlines: `canvas` (Reset
  camera, Illumination & XR, Expand) and `view` (Sequence: while there is a sequence and its card is hidden, or, with a
  menu, always: the viewer's Sequence menu). Each button runs its action (to be its command, gated by the command's `when`
  and `PluginConfig.Viewport.*`; the mocks run plain functions); the camera and rendering tools show only once something
  is loaded (`visible`), and with nothing to show the row has no frame at all. A button can carry a small menu of
  variants (`menu`, `CanvasMenuItem`: a label, an icon, a hint under it, a check for a toggle, or `disabled` with its
  reason in place of the hint): it **opens upward** from the row, its right edge level with the button's, on hover or
  focus (a corner mark says so; no tooltip then, since both would open in the same spot), stays while the pointer travels
  to it, and closes on a choice, on Escape, or a moment after the pointer leaves; a click still runs the button's own
  action, or opens the menu when the button has none. Tooltips of the row's buttons open upward too. Reset camera (a focus
  glyph) resets on a click and offers Reset zoom, Orient axes and Reset axes. Illumination & XR (a lightbulb, pressed while
  illumination is on) has no action of its own: Illumination as a checked toggle, Enter XR. Canvas settings and Animation
  are on the rail: they set up the scene rather than act on the view. Where the row does not fit (a narrow shell) it
  scrolls across. (Until 2026-10-14 the stack was a column at the canvas's top right.)
- **Expand by context.** In an iframe, filling the page only fills the frame, so what the canvas stack's Expand button
  does depends on where Mol* runs (`FullscreenModel`, `shell/fullscreen.ts`; molstar#1663): on a page it is Expand (fill
  the page; fixed position on the plugin root, no `<body>` mutation), with Fullscreen in its menu (the whole screen, not
  just the page); in an iframe that allows fullscreen (`allow="fullscreen"`) the button is Fullscreen, pressed while on;
  in one that forbids it, Expand fills the frame and the menu's Fullscreen is greyed with why ("Not allowed by the page
  Mol* is embedded in"). Fullscreen is the Fullscreen API on the plugin root; leaving it by Escape or the browser's own
  controls is followed (`fullscreenchange`).
- **Playback** (top left of the canvas, §4.3): one column, two different things as two surfaces: the session's snapshots
  (a card), then a structure's trajectory (a bar). Each only when the data has it.
- **Bottom of the canvas**: the status column on the left; on the right **the corner** (`.corner`): the work in hand
  (the selection bar, `ShellView`'s `context`, then the status chips) to the **left of the canvas toolbar on the same row**,
  bottom-aligned (the toolbar stays pinned to the corner: the bar grows leftward and appears or goes without moving it;
  where the row does not fit, the chips wrap to a row above); its right edge is the toolbar's (`--msp-canvas-corner-inset`). Bottom-aligned: the status column grows upward,
  its last line level with the toolbar's row. When the bar is narrow (under 560px, a container query on the bar, not the
  window), the corner takes a row of its own above the status, so nothing overlaps. In the compact layout the status
  and the selection are docked under the canvas; the toolbar stays in the corner.
- **Status column** (`StatusModel`, `shell/status.ts`; `StatusLine`): one line per item, growing upward, frameless (ink
  in the canvas theme's text color with a soft halo of the background, so it reads over the molecule): the latest
  message (it fades after 4 s; an error stays until the next message or `dismiss`, the line has no ✕; a click opens the
  log, `onOpenLog`: History's lower half, or its Log tab when compact), then each running task (a spinner, its label, its progress, ✕ to cancel), then the `agent` slot
  (`ShellView`'s `agent`): the assistant's turn as its card or, closed, as a line like a task's (`StatusPill`;
  [agent.md](agent.md) §5). The model keeps the last 200 messages for the log. Chips (`StatusChip`: a label, a click
  to activate, ✕ to dismiss) follow the selection bar on the right.
- **Busy overlay.** A task in overlay mode (`StatusTask.overlay`: rendering an animation, a long export) holds the
  plugin: a scrim over everything and one card with what runs, its progress (or an indeterminate bar), and Cancel. Other
  tasks stay in the status column.
- **Empty canvas.** With nothing loaded (`ShellView`'s `empty`), the Mol* wordmark shows in the middle of the canvas,
  outlined in the canvas's subtle text color: only after the canvas has been empty a while (700 ms; clearing a scene to
  load the next never flashes it), fading in, then the logo's colors wipe in along the outline once and stay; under the
  pointer it comes alive (a vertical strip at the cursor: the logo's colors running across it, centred on the cursor,
  the line itself unchanged; fading out when it leaves); it fades out, faster (250 ms), when something loads. Clicks
  pass through.
- **Home's heading** is the Mol* wordmark, outlined in the text color: the logo's colors flash along it now and then,
  and under the pointer the same live strip as on the empty canvas runs across it, back to the plain outline on leave.
- **Dropping files.** Files dragged anywhere over the plugin bring up a drop overlay over everything (Mol*'s drag
  overlay: "Drop to open" and the formats it reads, `dropFormats`); only drags that carry files show it, and dropped
  files go to the plugin's drag-and-drop manager (`onDropFiles`). On unless the spec sets
  `components.disableDragOverlay` (the host then leaves out `onDropFiles`). Home's Open Files has a file field of its
  own that takes drops too.
- **Compact.** Below 720px of the shell's own width (`ShellOptions.compactWidth`), the panel becomes a bottom sheet
  (`min(45cqh, 420px)`: a share of the shell's height, not the window's) with its tabs along the bottom: the tabs marked
  `compactBar` (Home, Structure), then the rail's items without a group (the theme switch, with its static icon rather
  than its face, and the command menu), then More, which opens above it with the other tabs (Selection, Snapshots, State
  tree, Log, Settings) and, after a divider, the rail's actions (Animation, Canvas settings, Export, Remote states). The
  status column and the selection dock between the canvas and the sheet as plain full-width rows (progress first,
  selection below; the assistant's turn is not shown there yet); a card fills the canvas; the canvas stack keeps only
  the frame margin (no rail beside it) and scrolls on a short canvas; playback ends before the stack.
- **Glass and canvas theme.** Everything over the canvas uses flat glass (translucent fill, blur, a hairline edge drawn
  as a spread shadow, no drop shadows) and takes its theme from the canvas background automatically, not from the UI
  theme, so it stays visible (as stories-untold's stage ink does). `model/canvas-theme.ts` computes it from the Canvas3D
  props by relative luminance (threshold 0.179, where light and dark text have equal contrast): the background color, a
  gradient's colors averaged over it by opacity; images and skyboxes count as the background color for now; a
  transparent canvas uses the host's fallback (`canvasThemeOf`, and `canvasTheme` over a stream of the props). The glass
  is tinted from the canvas color too (`ShellView`'s `ink` state, a `CanvasInk` `{ theme, color }` from `canvasInk(props$)`, sets `data-msp-theme` and `--msp-color-bg` on the canvas layer, and every
  derived shade follows): cards read as frosted canvas, a shade off the background, rather than slabs of the UI's own
  (warm) panel color. Two layers sit over the canvas: the view's chrome (sequence card, playback, canvas stack, status
  column, selection bar, hover label) follows the canvas (the ink); the application's chrome (the rail, the panel,
  the cards) follows the UI theme, so switching the theme changes it.
- **Hover label** next to the cursor (§4.4).

**State (mocks):** the shell is built (`views/shell`: `ShellView`, `FloatingDock`, `DockView`, `CanvasStack`,
`OverlayHost`, `StatusLine`, `BusyOverlay`, `EmptyCanvas`, `FileDropOverlay`, `SectionList`, `LogView`) over
`ShellModel`, `StatusModel`, `FullscreenModel` and `ThemeModel`. The Shell/Shell stories (`views/shell/shell.stories`)
register the default UI by hand with mock models (`testing/`): the tabs, the cards and their sections, the rail and
canvas stack items, the mode switch, playback, the selection bar, the command menu and, opt-in, the agent module. Not
yet: Expand itself (the mock reports "Expanded"; Fullscreen works); a viewport model following
`events.canvas3d.settingsUpdated` (the stories pass the background); a UI to move tabs between docks (`move` exists);
default content for the left and bottom docks.

## 3. Shell model

`shell/shell.ts`, as built (abridged):

```ts
type DockId = 'left' | 'right' | 'top' | 'bottom';
type ShellMode = 'desktop' | 'compact' | 'embedded';

interface DockState {
  tabs: readonly string[];
  active: string | undefined;
  size: number;               // width (left, right) or height (top, bottom), px
  collapsed: boolean;
  floating: boolean;          // over the canvas (right and top by default); desktop only
  pinned: boolean;            // a floating side dock beside the canvas, the canvas making room
}

interface ShellLayout { docks: Record<DockId, DockState>; mode: 'auto' | ShellMode }

class ShellModel extends ReactiveModel {
  constructor(options?: ShellOptions);                      // mode, docks, compactWidth (720), buttons
  dock(id: DockId): DockModel;                              // state, tabs, active (none while collapsed)
  readonly mode: ReadonlyState<ShellMode>;                  // `auto` resolved by the shell's own width
  readonly layout: ShellLayout;                             // as it can be persisted
  readonly overlay: ReadonlyState<OpenOverlay | undefined>; // { overlay, anchor?, focus? }
  readonly toolbar: ReadonlyState<readonly CanvasToolbarItem[]>; // rail and canvas stack items, `buttons` applied
  readonly tabFocus: ReadonlyState<{ tab: string; section: string } | undefined>;
  sections(slot: string): ReadonlyState<readonly ShellSection[]>;
  readonly page: ReadonlyState<OpenPage | undefined>;    // { page, focus? }
  registerTabs(tabs): () => void; registerOverlays(overlays): () => void;   // each returns its undo
  registerPages(pages): () => void; openPage(id: string, focus?: string): void; closePage(): void;
  registerSections(sections): () => void; registerToolbar(items): () => void;
  setButtons(buttons: Record<string, boolean>): void;       // show or hide toolbar items by id
  show(tabId: string, section?: string): void;              // activate and expand; at a section
  select(tabId: string): void;                              // a tab strip click: activate, or collapse if active
  toggle(dock: DockId): void; move(tabId: string, dock: DockId, index?: number): void;
  setSize(dock: DockId, size: number): void; setFloating(dock: DockId, v: boolean): void;
  setPinned(dock: DockId, v: boolean): void; setMode(mode): void; setWidth(width: number): void;
  openOverlay(id: string, anchor?: Element, focus?: string): void;
  toggleOverlay(id: string, anchor?: Element): void; closeOverlay(): void;
  sectionOpen(slot: string, section: string): boolean | undefined;  // as the user left it
  setSectionOpen(slot: string, section: string, open: boolean): void;
}
```

Beside it: `StatusModel` (`shell/status.ts`: `current`, `history`, `tasks`, `chips`; `post`, `dismiss`, `setTasks`,
`setChips`), `FullscreenModel` (`shell/fullscreen.ts`: `mode`, `active`, `available`, `embedded`, `unavailableReason`,
`toggle`) and `ThemeModel` (`shell/theme.ts`: `theme`, `set`, `next`, `applyTo`).

- **Default layout**: the right dock floating, open, with Home active while nothing is loaded (the first tab registered)
  and Structure once something is (the host shows it); the top dock floating, with the sequence once a structure is
  loaded (its tab registered then); left and bottom docks empty. Default sizes: left 300, right 340, top 88, bottom 220
  px.
- **Modes** replace `controlsDisplay`. `desktop`: the right and top docks float over the canvas (or sit beside it,
  `setFloating(dock, false)`). `compact` (chosen by `auto` from the shell's own width via `ResizeObserver`, not window
  media queries; below `compactWidth`, 720px): the right panel becomes a bottom sheet, cards fill the canvas, the rail's
  actions go under More in the sheet's tab bar (§2 Compact). `embedded` (canvas, canvas stack, command menu and status
  column only) is a mode value with no layout of its own yet.
- **Persistence.** The layout (`ShellModel.layout`) is to be serialised to an injected store keyed by app id
  (localStorage in the viewer), not part of plugin state snapshots. Not yet: the store, restoring (`registerTabs` keeps
  a restored tab's place), and `reset()` back to the spec's initial layout.
- **Ownership.** The shell model is to be a field of `PrismUIContext` (`ctx.shell`, §6), so the layout belongs to the UI
  context, not to the plugin. Today the host creates `ShellModel` and `StatusModel` and passes them to `ShellView`.
- **Compatibility.** Prism neither reads nor writes core `PluginLayout` (`PluginCommands.Layout.Update`) or
  `behaviors.layout.leftPanelTabName`: they stay the classic UI's, unchanged and not deprecated in v6. Extension code
  that writes them (`leftPanelTabName.next('data')`) is ported to `ctx.shell` or shell commands, or does nothing in
  Prism ([../analysis/extensions.md](../analysis/extensions.md)).
- **Not yet:** hiding a tab (`hide`), dialogs (`openDialog`), Expand as layout state, and the canvas options of the
  first design (`canvas.statusLine`, `canvas.hoverLabel`).

## 4. Default UI

### 4.1 Tabs (right panel)

| Tab                  | Content                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | Replaces                                                                  |
| -------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------- |
| **Home**             | Head of the rail (the Mol* mark), headed by the wordmark; shown while nothing is loaded. The state actions applicable to the root (`model/home.ts`, from `plugin.state.data.actions` for the root: Download Structure with its sources, Open Files (a file field that also takes drops, format, default visuals), Download Density, Load Trajectory, Open Session, and whatever extensions register), each a block that folds open (title, a one-line description; the first one open) and ends in Mol*'s footer (`ActionFooter`): set defaults, then the action at the right. Then, when the spec turns them on, the state server's remote states as a quiet list (one line each, the date at the right; a click opens one; refresh, and the server address behind a gear; a long list scrolls on its own). At the bottom, always in view, "⌘K Search every command". No Recent section. Files dropped anywhere open too (the drop overlay). | Home tab, custom import controls                                          |
| **Structure** (main) | The scene tree over a fixed-height inspector; a header with Preset, Styles, Component, and a gear for the component options (below).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          | Source, Components, Volume, Particles, Measurements list, extension cards |
| **Selection**        | What to do with the selection, as sections that start folded (a list of what can be done): Color, transparency, material… (Mol*'s selection theme actions: color, reset color, transparency, emissive, material, reset material, clipping, wiggle, on all or some representations), Create component (label, representation), Measure (§4.1a′), Superpose (with two or more structures). The selection bar opens it, at Measure from its Measure button.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      | Selection toolbar panels, Measurements and Superposition controls         |
| **Snapshots**        | §4.1b.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        | Plugin State tab                                                          |
| **State tree**       | §4.1c.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        | Data tab                                                                  |
| **Log**              | At the rail's foot. Every message the status column showed, newest first, with kind and time (meant so: `LogView` reverses `StatusModel.history`, which is already newest first, so the mock lists oldest first); the status column's message opens it. With History shown, it is History's lower half (its Log tab when compact; one panel: what happened and what was said); without, a tab of its own ([agent.md](agent.md) §5).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | Bottom log region, log overlay                                            |
| **Assistant**        | With the agent module only, at the rail's foot (its tooltip says _Set up assistant_; headed Assistant): providers, keys, model ([agent.md](agent.md) §5).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     | (new)                                                                     |
| **Settings**         | At the rail's foot. §4.5.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     | Settings tab, viewport settings panel                                     |
| Extension tabs       | Segmentation (volseg), PAE plot (model-archive; also movable to the bottom dock), app-specific tabs.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          | Large extension cards                                                     |

Presets and quick styles are commands (the command menu's Preset and Styles groups); measuring and superposing are
sections of the Selection tab, procedural animation a section of the Animation card, the screenshot options the Export
card's Image section, the component options behind the Structure tab's gear; each is also opened by a command at its
section (`measure.open`, `superpose.open`, `export.image`, …).

The Structure tab is built on `SceneModel` (`model/scene.ts`): `roots`, `picked`, `pick`, and `inspect(node)`, which
returns the node's `ParamsFormModel`, its facts, its actions, and its section slot `inspector:<kind>`. Section
contributions for a node kind register into that slot, the same mechanism the Export card uses. Tree rows are keyed by
node id; hover highlights, a click picks, visibility and the row's ⋯ sit inline (the eye stays visible while a node is
hidden; Remove is the ⋯'s last entry). Only top-level rows (a structure, a volume, the measurements) carry an icon;
below, the indent says where a row belongs. The form shows the essential params first and the rest under Advanced;
mapped params (type, color theme) show their choice with a sliders button that opens the chosen case's params, laid out
by the same rules (v5's `MappedControl`). The tab fills the panel: the tree
scrolls on its own and the inspector is always at the bottom at a fixed height (scrolling inside; "Select an item in the
tree to edit its properties" while nothing is picked; "Nothing to set here" for a node without settings, facts, actions
or sections), so picking never moves the tree. The header has Preset and Styles, which open the command menu inside
their groups (teaching the command menu as they are used): Preset lists the representation presets as the registry
groups them (Empty, Automatic; Basic; Miscellaneous; extensions'), Styles the quick styles (representation Default,
Cartoon, Spacefill, Surface; style Default, Illustrative), each a command that a search anywhere also finds. Then
Component, and a gear for the component options (hydrogens, visual quality, lighting, material, clipping, non-covalent
interactions, interior, animation), which open under the header: a panel's options get a gear, while sliders stay for a
row's own params. Every structure starts with "Structure": its entities, then chains, computed when first expanded
(mesoscale-safe), each selectable, focusable, colorable, or made a component. A row's "…" opens its menu inline under it
(Mol*'s component menu: add representation, set coloring, modify by selection with include, subtract and intersect,
select this, edit label). Inspections carry facts (a structure's chains, residues, atoms; a map's grid, voxel size, σ),
shown compactly above the form. Volumes are top-level entries next to structures (an EMDB map beside its fitted model);
density streaming is an extension row under the structure (§4.1a) whose channels appear under it.

### 4.1a Structure: hierarchy presets and extension rows

A structure's inspector leads with **Hierarchy**, Mol*'s structure hierarchy presets (Default (Assembly), Unit cell,
Super cell, Crystal contacts), the structure transform's own params under Advanced; the structure row's ⋯ offers the
same presets (Hierarchy preset).

**Extension rows.** An extension can add a row under every structure it applies to (`structureRows` in its UI
contribution): its label and icon, when it applies, its options (an `enabled` switch among them), a description that
follows them, and what appears under it once on. The tree knows no extension by name. A row is off (hidden, its eye
crossed) until its eye turns it on; picked, the inspector shows its options. Volume streaming (X-ray and EM entries:
off, then 2Fo-Fc σ, Fo-Fc(+ve), Fo-Fc(−ve) or EM as volume rows with their iso values; region, radius, detail) and
assembly symmetry (RCSB: the symmetry, its axes as a row, coloring by it) are rows this way. Every row keeps its ⋯ slot,
so the eyes stay in one column.

### 4.1a′ Measurements

Measuring is done with the selection: a **Measure** section of the Selection tab (after Mol*'s measurement controls,
`views/structure/measurements-add.tsx`, `model/measurements.ts`), which the selection bar's Measure button opens
(`show('selection', 'measure')`: a tab can be shown at a section). The **points** are the latest selections in order
(Mol*'s selection additions history), the first four used ("+n older, not used" under them), each reorderable (a handle
on hover; arrows for the keyboard), removable, focused on click, highlighted on hover; then what can be **measured**
from them, each a row with how many points it takes and enabled once there are enough: Distance (2), Angle (3), Dihedral
(4), Label (1), Orientation and Plane (the selection); then the measurements **added**, adjustable here as in the
Structure tab (show or hide, their style inline, remove); then the **options** (distance unit, text color and size)
folded under Options. The Measurements row in the tree lists them (each picked shows its own style); its ⋯ has Add
measurement… (the Measure section) and Remove all, and its inspector a way to Measure.

**Superpose.** With two or more structures, the Selection tab also has a Superpose section (after Mol*'s superposition
controls, `views/structure/superpose.tsx`, `model/superposition.ts`): the method (Chains, TM-align, UniProt (SIFTS),
Ligands, Atoms) and its options (align sequences for chains; trace only), what it would superpose (one selection per
structure, the first the reference, a click focuses one; UniProt uses every structure; otherwise what to select), and
Superpose in the footer, enabled once there is enough.

### 4.1b Snapshots tab

A tab beside Structure (`images` icon, `snapshots.show`, ⌘S), in the Structure tab's shape: a header of buttons (**+
Add** adds the current state; **Save** opens two choices inline, Session `.molx` with its data or State `.molj`, the
description whose data loads from its sources; **Open** a file; a gear opens the save options, what a new snapshot
includes, Mol*'s `PluginState.SnapshotParams`), the snapshots as rows (position, thumbnail when saved with an image,
name, `#key`, the description's first line; a click applies one; hover shows the handle that drags the row to reorder,
Retake from the current view, and ⋯ for Move up, Move down and Remove), and under them at a fixed height the current
snapshot's editor: name, key, markdown description, how long it shows and the camera transition ("Apply a snapshot to
edit its name, key and description" until one is). Replaces the Plugin State tab's local snapshots.

### 4.1c State tree tab

The state tree for advanced use, in the Structure tab's shape and look: cells as rows (the transform's object, its
description, an error in the danger color; show/hide where there is a visual), each with a ⋯ menu inline: **Update
{transform}…**, the actions applicable to the cell (Mol*'s Apply Action), and Remove. A click selects a cell; the fixed
inspector below shows its updates (its own transform's, then those of the decorators folded into it: A → A transforms
such as custom structure properties, which Mol* shows as part of the cell) and its actions as Home-style action blocks;
a ⋯ entry selects the cell and opens its block.

### 4.2 Overlays

Every card opens in the panel's place (§2 Overlays); the ones with several parts are section lists that extensions add
to (`export`, `canvas-settings`, `animation` slots).

| Overlay             | Opened from                                                                                                | Content                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            | Replaces                                  |
| ------------------- | ---------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ----------------------------------------- |
| **Export**          | The rail's `export` section; commands at a section (`export.image`, `export.animation`, `export.geometry`) | Sections: **Image** (the preview at the output's shape with the crop, Auto-crop in its corner; the screenshot params: resolution, format and quality, transparency, axes, illumination; the footer: set defaults, the size, Copy, Download), **Session & state** (`.molx` with data or `.molj`, Download), and the extensions': **Animation** (video-export: the animation, resolution, transparency; Render holds the plugin through the busy overlay, with Cancel), **Geometry** (geo-export: GLB, STL, OBJ, USDZ), **Models** (model-export: mmCIF, BinaryCIF). | Screenshot panel, export extension panels |
| **Canvas settings** | The rail's `scene` section; `canvas.settings`, `canvas.background`                                         | Sections, Canvas3D's params by what people reach for: Background, Lighting & effects, Camera, Clipping & fog.                                                                                                                                                                                                                                                                                                                                                                                                                                                      | Viewport settings panel                   |
| **Animation**       | The rail's `scene` section; `animation.open`                                                               | Sections: **Play** (an animation: trajectory, camera spin and rock, state snapshots, assembly unwind; its options; Play, Stop while playing) and **Procedural animation** (wiggle by dynamics or model uncertainty; Apply, Clear).                                                                                                                                                                                                                                                                                                                                 | Animation viewport controls               |
| **Remote states**   | The rail's `export` section, next to Export, only when the spec turns remote states on                     | Upload the current state under a name (an action block, Upload in its footer); the server's states (a click opens one; ⋯ downloads or removes; refresh); the server address.                                                                                                                                                                                                                                                                                                                                                                                       | Remote states section                     |

Help is the command menu's `?`: its prefixes and every command with a shortcut ([commands.md](commands.md)). Not yet:
the old Help tab's mouse controls and selection help.

Dialogs (modal) are planned for flows that must not be lost by a stray click: command forms opened "as form", Zenodo
import. Not built yet. MP4 export, once planned as one, is the Export card's Animation section; its render holds the
plugin through the busy overlay instead.

### 4.3 Snapshots and timeline

Each part appears only when it applies.

- **Snapshot card** (top left, when snapshots exist; `SnapshotCard`, `views/canvas/playback.tsx`): the current
  snapshot's title, its position (a dot per snapshot up to 7, each jumping to its snapshot with its title as the
  tooltip; past that "3 / 12", which opens the list of titles, filterable when long), previous, next, and play for the
  slideshow. Keyboard: the existing Ctrl+←/→ bindings and the commands `snapshot.next` / `.previous`.
- **Snapshot description** (when the current snapshot has one): markdown under the card's head with the existing
  markdown extensions (`#key` links apply snapshots, command links, color swatches, images, audio;
  `SnapshotStepperModel.openLink`); it scrolls when long. Apps can replace it with their own renderer (a contribution),
  as mvs-stories and mesoscale-explorer do today.
- **Trajectory** (a bar under the snapshot card, when a structure has several models; `TrajectoryBar`): a different
  thing from snapshots. Snapshots step through saved states of the whole session; a trajectory steps through one
  structure's models (its model index, also in that structure's inspector). The bar names its structure (`1TQN 4 / 20`),
  a quiet scrubber, play; with several trajectories loaded it follows the focused structure.
- **Animation**: the rail's Animation card (play an animation, procedural wiggle); trajectory loop and duration live
  there.

Not yet: a progress bar during animated transitions, a `snapshot.play` command, collapsing the description, its text
size, and the replaceable renderer. The mocks render descriptions with `Markdown` (`components/markdown`), a small safe
stand-in (no HTML is parsed; headings, paragraphs, lists, emphasis, code and links, which go to `openLink`) until the
plugin's markdown extensions are ported.

### 4.4 Hover label

The label follows the cursor instead of sitting in a corner (`HoverLabel`, `views/canvas/hover-label.tsx`):

- Lines of text from the plugin's hover labels (`HoverLabelModel.lines`, the first one stronger: "GLU 213", then "Chain
  A · 1TQN"), positioned from the pointer over the canvas with an offset, flipped at the canvas edges,
  `pointer-events: none`. The sequence's hover feeds the same label.
- Positioned in `requestAnimationFrame` outside the reactive render path (one transform write per frame, no state or
  allocation per move); hidden while a button is held (rotating, dragging) and when the pointer leaves the canvas.
- Not yet: the setting `canvas.hoverLabel` (`cursor`, the default; `corner`, today's behavior; `off`), showing at the
  tap position on touch devices, and HTML or markdown labels from `behaviors.labels.highlight` (the mock's lines are
  text).

### 4.5 Settings tab

Today's settings panel lists every canvas parameter at one level with truncated labels ("Camera Rese…", "Transparent
…"), then the behaviors tree. The tab keeps the same parameters but organises them for use
(`views/settings/settings-tab.tsx`, `model/settings.ts`); it sits at the rail's foot:

- **Search** at the top filters across all settings by label and description (`ParamsForm`'s `filter`): only matching
  rows show, their groups open, a group whose name matches shows all its rows, a settings group with no match hides; "No
  setting matches" when none does. Escape empties it.
- **Groups by intent**, essentials first, each a block that folds (the first one open): Appearance (background,
  occlusion, shadows, outline, fog, marking), Camera (mode, field of view, clipping, reset), Interaction (pick padding,
  trackball, pointer), Performance (resolution, pixel and pick scale, transparency mode, multi-sampling, Hi-Z),
  Interface (theme, layout, pinned panel, sequence, tooltips), Behaviors (`plugin.state.behaviors` by category, each a
  switch with its own params behind the sliders: the on/off pattern), Advanced (renderer, handle, XR). Groups are to
  come from a mapping over the existing PD params, so no param is lost and new ones land in Advanced (the mock declares
  its own params).
- **Wider labels** than other forms (a full label as the row's tooltip); Behaviors give the label the row.
- Changes apply as they are made; no reset buttons.

### 4.6 Canvas mode and selection bar

What a click on the canvas does is a mode, and the difference matters, so the mode is always visible: the top block of
the canvas stack is two buttons, **Interact** and **Select** (`views/canvas/mode-switch.tsx`), once something is loaded.
Interact: drag rotates, a click focuses. Select: a click selects (Shift extends), drag still rotates. Select takes the
accent while on; each button's title explains its mode and switching posts a one-line hint to the status column. Select
carries the picking level as a small badge (its glyph, bottom right). **Hovering Select (or focusing it) shows the list of
picking levels directly, above it** (`PickLevel`, `ShellView`'s `canvasAside`, outside the row's frame so it is not
clipped; always on touch screens): Mol*'s granularities (atom, residue, chain, entity, model, operator, structure, then
the All instances levels: atom, residue, chain), the current one marked in the accent; one click from the toolbar picks it
and turns select mode on if it was off (Tab goes from Select into the list). Not yet: the mode as a command
(`selection-mode.toggle`, with a keybinding) and Escape leaving Select.

While anything is selected, the selection bar sits in the bottom-right corner, to the left of the
canvas toolbar on its row (moving with the panel edge; `SelectionBar`): what is selected ("12 residues"), then what to do with it as icon buttons: Color, theme, component…
(opens the Selection tab, beside Structure in the tab strip), Measure… (opens it at Measure), Focus, and Clear. One
button leads to the Selection tab rather than one per action; the tab's sections hold Mol*'s selection toolbar actions.
The command menu offers the same as commands: the `selection` group (`selection.color`, `selection.transparency`,
`selection.clear`), `measure.open`, `focus.target`. Selection operations (add, subtract, intersect, set) and queries
come from the command menu's `select` group ([commands.md](commands.md) §10).

### 4.7 Icons

One consistent modern line set (Lucide, ISC-licensed; the old UI's vendored Material icons are not carried over; only
the SVG path data that is used is copied into `components/icon/`, no runtime dependency), plus custom glyphs for
molecular concepts (cartoon, surface, ball-and-stick, ligand, volume, trajectory) that no generic set covers. Every icon
is chosen once per command, so buttons, menus, and command menu rows agree.

Implementation: `scripts/icons.mjs` extracts the inner markup of the used icons from `lucide-static` 1.48.0 (a
devDependency; `pnpm icons`, and `icons:check` to verify, not yet run in CI) into the generated
`components/icon/_generated/lucide.ts`; the custom molecular glyphs are hand-written in `components/icon/molecular.ts`
on the same 24×24 stroke grid. `Icon` parses each icon's markup once into a cached `<svg>` template and returns a clone
with size, class, and label (`aria-hidden` otherwise) set, so rendering an icon allocates one node tree and no parsing.

## 5. Contributions

> **Superseded (2026-10-15):** prism-ui has no contribution registry. Extensions place their UI in the viewer's named
> slots ([../spec/viewer-spec.md](../spec/viewer-spec.md) §3–4); `PrismUISpec` (§6) is replaced by `PrismViewerSpec`.

A contribution is a plain, declarative object, listed explicitly in the spec (like `registry` entries):

```ts
interface UIContribution {
  id: string;
  commands?: readonly UICommand[];
  commandGroups?: readonly CommandGroup[];        // new groups, or nested under built-in ones (commands.md §9)
  tabs?: readonly TabContribution[];               // panel tabs (right by default, any dock)
  sections?: readonly SectionContribution[];       // sections of a slot: export, selection, inspector:<kind>, …
  overlays?: readonly OverlayContribution[];       // cards in the panel's place
  canvasItems?: readonly CanvasItem[];             // rail and canvas stack buttons, playback parts, other chrome
  statusItems?: readonly StatusItem[];             // status chips
  menus?: readonly MenuItem[];                     // context menu entries (command refs)
  sceneEntries?: readonly SceneEntryContribution[];// scene-tree entries; replaces genericRepresentationControls
  paramEditors?: readonly ParamEditorContribution[]; // replaces customParamEditors (by transformer id or PD field)
  settings?: readonly SettingsGroupContribution[]; // groups or params added to the Settings tab
  structureRows?: readonly StructureRowContribution[]; // rows under applicable structures (volume streaming, symmetry)
  formats?: readonly [name: string, ArgumentFormat][];
  modules?: readonly PrismModule<any>[];           // its models (contexts.md §3), retired with it
}

interface TabContribution {                       // ShellTab, plus `when`
  id: string; title: string; icon: IconName;
  heading?: string | (() => Child);                // the panel's header when not the title (Home: the wordmark)
  dock?: DockId;                                   // default 'right'
  order?: number;
  place?: 'head' | 'strip' | 'foot';               // on a floating dock's rail (Home head, Settings foot)
  compactBar?: boolean;                            // in the compact tab bar rather than under More
  when?(ctx: CommandContext): boolean;            // e.g. only when a volseg entry exists (planned)
  render: ViewRenderer<ViewContext>;               // dom, React, Preact, or element (rendering.md §6)
}

interface SectionContribution {                   // ShellSection
  id: string; slot: string;                        // 'export', 'selection', 'inspector:structure', …
  title: string; icon?: IconName;
  keywords?: readonly string[];                    // found by these too ('mp4' for Export › Animation)
  order?: number;
  when?(ctx: ShellViewContext): boolean;          // ctx.subject: the inspected node for inspector:<kind>
  render: ViewRenderer<ShellViewContext>;
}

interface OverlayContribution {                   // ShellOverlay
  id: string; title: string;
  render: ViewRenderer<ShellViewContext>;          // opened by a command or a toolbar item; ctx.focus: a section
}
```

- **Placement of extension UI.** Most extension panels become rows under a structure (`structureRows`: volume streaming,
  assembly symmetry; §4.1a), inspector sections (kinemage), sections of a card (MP4, geometry and model export in
  Export; their commands open the card at them), or dialogs (Zenodo import). Only large, standalone views become their
  own tabs (segmentation, PAE plot). Extensions add no rail buttons.
- Contributions are activated after `plugin.init()`; each activation's registrations are undone on deactivation, as
  `plugin.register` does for registry entries.
- `when` replaces `isHidden` toggling inside cards; a tab, section, or item appears when its condition holds.
- An extension package exports its registry entry, its behavior, and (separately) its UI contribution, e.g.
  `@molstar/geo-export-extension/prism-ui`. Apps list what they want; nothing registers UI as a side effect of a
  behavior.

**State (mocks):** `UIContribution` and its registry are not built. The stories register straight into the shell
(`registerTabs`, `registerOverlays`, `registerSections`, `registerToolbar`, each returning its undo, which is what
deactivation will call); the extension rows are a scene mock's specs (`testing/structure-rows.ts`); a tab's `when` does
not exist yet (sections and toolbar items have theirs: `when`, `visible`).

## 6. The contexts and the UI spec

The contexts (`PrismContext`, `PrismUIContext`), modules, the mode, contributions, the command registry, history and the
app's lifecycle are specified in [contexts.md](contexts.md). In short: `PrismContext` holds what every mode shares (the
plugin behind sources and actions, the command registry, status, history, settings, the module registry);
`PrismUIContext` adds the shell, keybindings, contributions, theme and fullscreen; an app runs one mode, chosen before
it exists; a mode, an extension or the agent attaches its models as modules.

**The plugin is the engine, behind the context.** Prism keeps `PluginContext` for what it does well: the state tree and
transforms (sessions, snapshots, undo, the reproducibility story mode relies on), the canvas and its renderer, picking,
highlight and selection over loci, the camera, formats, presets, custom properties, the MVS runtime, the extensions. It
does not carry the plugin's UI-era surface:

- **One layer touches the plugin**: `sources/` and `actions/` folders, in any Prism package (core, a mode, an
  extension's prism entry); views, view models and commands see sources, actions and models only (contexts.md §7).
- **A slim plugin.** Prism's spec composes the plugin from v6's empty registries with only the behaviors and entries
  Prism uses.
- **Not used:** `PluginLayout` and `behaviors.layout.leftPanelTabName` (the shell is Prism's layout; Prism never sets
  `plugin.layout`, and an action calls `plugin.handleResize()` when the canvas inset changes), the `PluginCommands` bus
  for UI work (Prism's command registry), toasts (the status column), the UI-only maps (`customStructureControls`,
  `customImportControls`, `genericRepresentationControls`; contributions replace them), `customParamEditors`
  (`paramEditors` contributions), and behaviors that exist only to drive the classic UI. Each is a candidate for removal
  when the classic UI retires ([../../plans/breaking-v6-changes.md](../../plans/breaking-v6-changes.md)).
- **The engine is swappable.** As the UI sees sources only, a mode can run on another engine later (story mode could
  render MVS documents through a leaner path) without a UI change.

**State (mocks):** not built: the contexts, the specs, `DefaultUIContributions()` and `createPrismUI`. The mock app
(`PrismAppStory`, `ShellStory` with one mode) creates the models and passes them to `ShellView`; `ShellOptions.buttons`
is the spec's `layout.buttons`, and a host leaves out `ShellView`'s `onDropFiles` for `components.disableDragOverlay`.
`StatusModel` lives in prism-ui's `shell/` for now and moves to `@molstar/prism` with the context.

The UI spec (the spec's `ui`, besides the plugin's spec, the frame's and the mode's contributions; contexts.md §4):

```ts
interface PrismUISpec extends PrismSpec {
  ui?: {
    contributions?: readonly UIContribution[];     // the frame's, explicit and ordered; DefaultUIContributions()
    layout?: ShellOptions & {                       // the shell's own options (docks, mode, buttons)
      tabs?: Record<string, { dock?: DockId; order?: number; hidden?: boolean }>;
    };
    theme?: UITheme;                                // 'light' | 'dark' | 'blue' | 'auto'; the rail's theme button changes it
    commands?: { hide?: readonly string[] };        // command or group ids, globs allowed ('superpose.*')
    keybindings?: readonly { command: string; key: KeyChord | null }[];
  };
}
```

- Persistence of the layout, the mode and the theme is the spec's `persist: { store: PrismStore; key }` (contexts.md
  §8). The shell keeps where a tab was placed when it is unregistered, resolves a stale active tab to the first
  registered one, offers `layoutChanges` (debounced), `restore(layout)`, `reset()` and `hideTab`, and clamps sizes per
  axis.
- `DefaultUIContributions()` in `@molstar/prism-ui/default-spec` lists the built-in contributions; apps spread it and
  add extension contributions, mirroring `DefaultRegistry` and the v6 composition rule.
- Entry point: `createPrismUI({ target, spec, mode? })` returns `{ ctx, dispose }` (contexts.md §1); `<prism-viewer>`
  wraps it with a shadow root. `prism-ui/react` offers `<PrismUI>` for React hosts. Headless:
  `createPrism({ spec, mode })`.
- Today's spec fields map as follows: `components.controls.*: 'none'` → hide the tab or dock in `layout`;
  `components.viewport.view` → canvas parts and overlays; `components.structureTools` → contributions;
  `remoteState: 'none'` → omit the remote-states contribution (Home's list, the rail's button, the card);
  `customParamEditors` → `paramEditors`; `selectionTools.hide` → omit or reorder selection-bar items;
  `PluginConfig.Viewport.*` → toolbar commands' `when` (kept, so the viewer options keep working).

## 7. Embeddable elements

Top-level custom elements over the same models, for hosts that do not want the shell:

| Element                  | Content                                                                    | Used by                               |
| ------------------------ | -------------------------------------------------------------------------- | ------------------------------------- |
| `<prism-viewer>`         | The whole app (its mode, or the splash), in a shadow root (contexts.md §1) | `apps/prism`, hosts embedding Prism   |
| `<msp-viewport>`         | Canvas + overlays + toolbar (+ command menu)                               | mvs-stories, docking viewer, examples |
| `<msp-sequence>`         | Sequence view                                                              | Hosts with their own layout           |
| `<msp-params-form>`      | PD form bound to a `ParamsFormModel`                                       | alpha-orbitals example, host settings |
| `<msp-snapshot-stepper>` | Snapshot i/N, previous/next, description                                   | mvs-stories                           |
| `<msp-pae-plot>`         | PAE plot (model-archive contribution)                                      | alphafolddb-pae example               |

Each takes a `ui` or `model` property; several elements on a page can share one `PrismUIContext`. Not yet: the elements
themselves; the base class exists (`MspElement`, `defineElements` in `dom/element`), and the stories mount `ShellView`
and the views directly.
