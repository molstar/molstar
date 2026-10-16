# The shell

The frame around the canvas that prism-ui gives an app: the rail and its panel, the docks, the cards, pages, the chrome
over the canvas, the status column, the command menu's place. The model is `src/shell/shell.ts` (with
`managers/view-state.ts` and `managers/navigation.ts`), the views `src/views/shell/`. What goes in it (tabs, cards,
buttons) is the app's. The code is the truth for how.

## 1. Principles

- **Canvas first.** The canvas fills the shell; panels float over it; persistent chrome is the exception.
- **Summoned, not resident.** Short tasks open as cards in the panel's place and close when done; docked UI is for what
  stays in view while working.
- **One home for persistent UI**: a tabbed panel on the right. Left, top and bottom docks exist; the app places tabs
  there.
- **Appear when relevant.** Chrome shows only when it applies (camera tools once something is loaded, playback when
  there are snapshots, the selection bar while something is selected). The empty canvas shows the wordmark after 700 ms
  empty and idle (any work restarts the wait; a load never flashes it), with a splash behind it (`PrismSplash`:
  `prism-field.ts`, light through the star; `quantum-field.ts`, the same as a wave, `field: 'quantum'`); the word is
  glass (`glass.ts`, `GlassWordmark`) that assembles, leans under the pointer and shatters on a press; the light stays
  out of where the glass covers it (`glass-footprint.ts`). Which splash: random each time it shows, the wave 4 in 5
  (`EmptyCanvas` `field`). The wave: colored by wavenumber in the logo's colors; enters as a broad front, already
  flowing at the start, as strong as the glass is assembled. The canvas starts quiet (transparent, its axes too) and is
  quiet while the mark is up; whatever is in it brings it back gradually (`onCanvasShown`, `.canvas[data-quiet]`, 0.9 s,
  over the canvas's own color; same guards as the mark). `ShellView`'s `splash` turns it off.
- **Motion.** Open and close by `Show`'s `reveal`: grow (`data-reveal`) for lists and sections, `'side'` for panels,
  `'fade'` (opacity only) for what must not change size, like the sequence band drawn on a canvas (re-laid out every
  frame if it grew). One value replacing another (a right panel by another) can swap as a view transition scoped to the
  dock (`swap`, motion.ts `swapIn`), where the browser has scoped ones; never the document's (it would stop the 3D
  canvas while it runs).
- **Lean views.** Commands carry the operations ([commands.md](commands.md)); views show state and the common actions.

## 2. Structure

- **The right panel floats** over the canvas (`DockState.floating`, on for `right` and `top`): the canvas keeps its full
  width, the rail sits at the right edge, the active tab is a glass card beside it. Its header: title (or `heading`),
  pin, ✕. Its inner edge resizes it (`setSize`, at least 180px; double-click back to 340). `setFloating(dock, false)`
  docks it beside the canvas. One margin frames the canvas (`--msp-canvas-edge`, 8px).
- **Pinned** (`setPinned`): the rail and the open panel move beside the canvas, which ends where they begin; Escape and
  clicks outside leave a pinned panel or card.
- **The canvas never shrinks for the floating panel or the rail**: the app's camera framing keeps the scene in the
  visible area.
- **Rail.** Full height at the right edge, buttons edge to edge, no labels (tooltips name them): the `head` tab, the tab
  strip by `order`, the toolbar sections by `group` (in the order of their first item), then at the foot the `foot` tabs
  and the items without a group (the theme switch, the command menu with its ⌘K `caption`). Rail buttons are
  `CanvasToolbarItem`s (`run`, or toggle `overlay`; `pressed`, `visible`, `face`). Which show is the app's
  (`ShellOptions.buttons`, `setButtons`). Extensions add sections to panels, not rail buttons.
  - **Steps by its own height** (a size container; CSS, nothing measured): normal from 620px, compact from 530px, small
    below (`RailSteps` in `views/shell/rail.ts`, kept equal to `dock.module.css` by `_test/rail.test.ts`).
  - **More**: when even the small step is too tall, the lowest-priority buttons (`priority`; default their position, so
    the foot goes first; the command menu -1) go under a More menu. `fitRail` computes it on one `ResizeObserver`, only
    when the rounded height or the buttons change.
  - **Compact rail** (`ShellOptions.compactRail`, `toggleRail`): the column hides behind a button in the canvas's corner
    and slides in over the canvas, taking no width; hidden, it is inert and paused. Hiding closes the card and the
    panel.
- **Theme switch**: the themes beside the rail (`UIThemes`), `ThemeModel` (`shell/theme.ts`) keeps `data-msp-theme`. The
  rail, panel and cards follow the UI theme; the canvas chrome follows the canvas (Glass, below).
- **Command menu**: opened by the rail's last button and ⌘K / ⇧⌘P (the host listens); modal over the whole shell
  (`ShellView`'s `commandMenu`); whether it is open is the host's state. With an assistant, an Ask row under its list
  ([commands.md](commands.md) §5).
- **Tabs and sections.** A tab (`ShellTab`: `title`, `heading?`, `icon`, `dock`, `order`, `place`, `priority`,
  `compactBar`, `frame`, `render`) goes into its dock by `order`. A strip click activates it or collapses the active one
  (`select`); `show(tab, section?)` activates and expands. Panels and cards that hold several things are **sections**
  (`ShellSection`: `slot`, `title`, `heading?`, `keywords`, `order`, `when`, `open`, `render`) in a `SectionList`;
  `show(tab, section)` or `openOverlay(id, anchor, section)` opens that one and scrolls to it (`tabFocus`, the card's
  `focus`). A section otherwise opens as the person left it (`sectionOpen`/`setSectionOpen`).
- **Cards** (`ShellOverlay`: `title`, `render`, `persistent?`, `fill?`), one at a time, take the right panel's place,
  full height, whether the panel is open or not; over the open panel they swap in without a transition and the panel
  steps aside. A card closes on Esc, its ✕, or a pointer down outside it (not on the canvas, which stays usable; not a
  `persistent` one). Panels and cards share one header.
- **Pages** (`ShellPage`, `openPage`): a full-size opaque card over the canvas for what needs room (the agent's
  Inspector), with its own header; Esc closes it, focus goes back; opening one closes a card.
- **Hidden panels.** A panel is mounted only while it is its dock's active tab (`Show` in `DockView`). What the shell
  hides but keeps mounted is inert and paused ([styling.md](styling.md) §5, `activity`): the panel behind a card, a
  folded section's body (binds nothing until opened), a panel while it closes. The `hidden` or class that hides is bound
  outside the paused scope.
- **Escape** closes the innermost thing, one per press: a list or menu, then the card, then the panel (`ShellView` on
  the window, after the document's handlers, skipping handled events; not from a field). A search field empties first.
  Pinned panels and cards stay. The command menu handles its own.
- **Top dock** (floating): a glass card along the canvas's top, or with `frame: 'none'` a **band** over the canvas in
  its ink, fading out below its lines (`BandFadePx`, `BandDissolvePx`; the faded part lets the pointer through). The
  chrome below starts under it (`--msp-canvas-top-inset`).
- **Optional docks** `left`, `top`, `bottom` take the app's tabs (`move`); side docks run full height, top and bottom
  span the canvas.
- **Canvas toolbar**: one glass row in the canvas's bottom-right corner (`place: 'canvas'` items, sections by hairlines,
  after `ShellView`'s `canvasLead`); a floating panel stops above it (`--msp-corner-room`). A button's `menu`
  (`CanvasMenuItem`) opens upward on hover or focus; a click still runs its action. A narrow row scrolls across.
- **Expand by context** (`FullscreenModel`, `shell/fullscreen.ts`): on a page Expand fills it (fixed position, no
  `<body>` change) with Fullscreen in its menu; in an iframe that allows it the button is Fullscreen; one that forbids
  it greys Fullscreen with why.
- **Bottom of the canvas**: the status column on the left; the corner on the right: `context` (the selection bar) and
  the status chips left of the toolbar on its row. A narrow bar (container query, under 560px) puts the corner on a row
  of its own.
- **Status column** (`StatusModel`, `managers/notifications.ts`; `StatusLine`, `StatusPill`): one line per item, growing
  upward, frameless ink with a halo: the latest message (fades after 4 s, held while hovered; an error stays until its
  ✕; a click opens the log, `onOpenLog`), then each running task (spinner, progress, ✕ to cancel), then the `agent`
  slot. The model keeps 200 messages.
- **Busy overlay**: a task with `overlay` holds the app: a scrim and one card with its progress and Cancel.
- **Dropping files** anywhere shows the drop overlay (`dropFormats`) and hands the files to `onDropFiles`; without
  `onDropFiles`, no overlay.
- **Compact** below `compactWidth` (720px) of the shell's own width (`ResizeObserver`, not media queries): the panel is
  a bottom sheet (`min(45cqh, 420px)`) with its `compactBar` tabs, the items without a group, and More; the status and
  the selection dock between canvas and sheet; a card fills the canvas.
- **Glass and canvas ink.** Everything over the canvas is flat glass whose theme follows the canvas background
  (`model/canvas-theme.ts`: relative luminance, threshold 0.179; a gradient averaged by opacity; transparent: the host's
  fallback) and is tinted from the canvas color (`ShellView`'s `ink`, a `CanvasInk`, sets `data-msp-theme` and
  `--msp-color-bg` on the canvas layer). Two layers: the view's chrome follows the ink; the app's (rail, panel, cards)
  the UI theme.

## 3. The shell model

`ShellModel` (`shell/shell.ts`), DOM-free: `dock(id)` (`DockModel`: `state`, `tabs`, `active`), `mode`, `layout`,
`overlay`, `page`, `toolbar`, `tabFocus`, `sections(slot)`; `registerTabs`, `registerOverlays`, `registerPages`,
`registerSections`, `registerToolbar` (each returns its undo); `show`, `select`, `toggle`, `move`, `setSize`,
`setFloating`, `setPinned`, `setMode`, `setWidth`, `openOverlay`, `toggleOverlay`, `closeOverlay`, `openPage`,
`closePage`, `setButtons`, compact rail and `expanded`.

- **Two managers hold its state**, given or made: `viewState` (`ViewStateManager`: the layout, serializable: docks'
  tabs, sizes, float, pin, mode, compact rail, buttons) and `navigation` (`NavigationManager`: what is open: each dock's
  tab and collapse, the card, the page, a tab's section, `expanded`). An app keeps and restores the layout through
  `viewState`.
- **Modes**: `desktop`, `compact` (chosen by `auto` from the shell's width), `embedded` (a value, no layout of its own).
- **Default sizes** (`DefaultDockSizes`): left 300, right 340, top 88, bottom 220 px; minimums in `MinDockSizes`.
- **Prism's own layout**: the plugin's `PluginLayout` and `leftPanelTabName` are neither read nor written.
