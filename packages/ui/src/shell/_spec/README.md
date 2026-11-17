# Mol\* UI shell: the specification

`@molstar/ui/shell` (`packages/ui/src/shell/`): a generic app shell. The app gives it a main view and views; the shell
docks them around the main area, floats them as windows, or shows them in a layer over it. This spec is the behaviour
and its reasons; how a view's content is laid out is the foundations' (`styles §8`). Code comments cite its sections
(`shell §2`). The code is the truth for how.

| File                                  | What                                                                                                 |
| ------------------------------------- | ---------------------------------------------------------------------------------------------------- |
| `layout.ts`                           | `ShellLayout`, its invariants (`normalize`), its index, operations, validation, placement: pure (§2) |
| `shell-model.ts`                      | `ShellModel`: registrations, the layout, layers, navigation, modes, persistence (§2); labels (§6)    |
| `rail.ts`                             | `RailSteps`, `railEntries`, `fitRail`: the rail's order, steps and More (§3.2)                       |
| `view-menu.ts`                        | `viewMenu`, `paneMenu`, `railMenu`: the layout menus (§2, §2.1, §3)                                  |
| `pane-map.ts`                         | `PaneMap<T>`: what the app keeps a pane, app side (§2.1)                                             |
| `view-frame.tsx`                      | `ViewFrame`, `ViewHost`, `viewBody`, the view primitives (`ViewStack` … `ViewHeading`, `styles §8`)  |
| `section-list.tsx`                    | `SectionList`: a slot's sections as `Expandable`s (§2)                                               |
| `shell-view.tsx`                      | `ShellView`: the trees (the shell's, the main area's region), rail, windows, layers, F6 (§2.1, §3)   |
| `shell.module.css`, `view.module.css` | the shell's sizes and the view frame (`styles §8`)                                                   |

## 1. Principles

- **App-independent.** The shell knows views, groups, windows, cards and pages, not molecules. The main area is the
  app's (a canvas, a document, a chart); nothing in the shell reads a plugin, a canvas or a command manager. Why: one
  shell for every Mol\* app, and apps that are not about structures.
- **One docking model.** Every panel is a _view_ in a _group_; a group is docked in a split tree around the main area or
  floats as a _window_. Docks, pinned and floating panels, and dialogs are all cases of it (§2). Why: one set of
  operations (show, move, float, dock, close), one stored value.
- **Docked is beside the main area.** A docked group takes its width; the main area ends where it begins. Nothing floats
  over the main area except windows (and the layers: card, page, modal, busy). Inside the main area, the same model
  docks views beside the app's content (§2.1). Why: a panel over the main area hides what it is about. A panel floating
  over the main area, chrome tinted by the main area's content, and a splash are not part of this shell.
- **The Mol\* star** (`molstar` icon) marks Home (an app's view, by convention `place: 'head'`) and the rail's toggle in
  the main area's corner (§3.2).
- **Hidden UI does no work**: a view body no frame shows is paused (`dom §1`, activity), and so is a group that does not
  show (a folded rail group, the right column while the rail is hidden, a folded tab group's body); the rest of the
  shell under a modal is `inert`.
- **The app keeps what is the app's** (§5): the main area's toolbar, the status, the log, storage, commands.

## 2. The model

`ShellModel` (`shell-model.ts`), DOM-free; the layout is a value (`ShellLayout`, `layout.ts`) and every operation is a
pure function on it ending in `normalize`, so the model only ever holds a valid layout.

- **A view** (`ShellViewDef`): `id`, `title` (its name: the rail's tooltip, a tab's label, its header), `heading?`
  (shorter text or a view for its header), `icon`, `render` (a `ViewRenderer`: this runtime's, a custom element,
  React's), `actions?` (rendered in its header), `order?` (among its group's views, across contributors), `place?`
  (`head` | `strip` | `foot` on the rail), `visible?` (a reactive flag), `priority?` (More, §3.2), `compactBar?`,
  `size?` (as a window; default 360 × 320), `minSize?`, `home?`, `instances?`, `remember?`, `persistent?` (a card that
  stays on a press outside), `fill?` (the content fills the body and scrolls its own parts). `registerViews` returns its
  undo (an extension unloads); an unregistered view keeps its place in the layout (it lands there when registered again)
  and leaves any layer it was open in.
- **Sections** (`ShellSection`: `slot`, `title`, `heading?`, `order`, `when?`, `whenChanged?`, `open?`, `render`): what
  extensions add to another view's slot instead of a view of their own (`registerSections`). A section shows while its
  `when` holds, asked again as the view's subject changes and as its `whenChanged` does (a value fetched for the
  subject). `SectionList` shows a slot's sections as `Expandable`s by order, each open as the person left it this
  session (`sectionOpen`), else open, or folded under `collapsed` unless its own `open`. `show(id, { section })` opens
  that section and scrolls it into view (`ctx.focus`, a new value at each such `show`).
- **Rail items** (`RailItem`: `run?(anchor)`, `menu?()`, `group?`, `order?`, `priority?`, `caption?`, `pressed?`,
  `visible?`): app-wide tools on the rail that are no view (the theme, the command menu); `registerRailItems` returns
  its undo. An item with a `menu` and no `run` opens it on a press, beside the button (`openMenuBeside`: its top at the
  button's, its bottom where it would leave the window), read anew each time (a choice marks the chosen row with
  `checked`). Never on hover: the pointer travels along the rail to other buttons, the hover is the tooltip's, and touch
  has none.
- **The app's own settings are at the bottom**: the theme (and the like: what changes the app, not the work) is an item
  with no `group`, so it sits in the rail's foot, under the foot views. **The command menu is in the corner**: an item
  with `corner`, the rail's last button, under its Hide button (the shell's bottom-right corner; it stays when the rail
  runs short). Why: one place in every app built on the shell, found without looking; the top of the rail is for the
  work.
- **A group** (`GroupNode`): views shown one at a time, its `active` one, `collapsed` (only its strip shows), a `size`
  in px along its parent split. Its strip is the rail (§3.2) or a row of labelled tabs (§3.3).
- **The layout** (`ShellLayout`, versioned): `root`, a tree whose leaves are groups and the one main area and whose
  inner nodes are splits (`row` | `column`, `size` in px; the node holding the main area takes the rest); `windows` (a
  group, a `rect` relative to the main area, a stacking `z`); `hidden` (views out of view by choice); `rail` (whether
  the rail and the docked column show); `seq`. Default: `row(main, group('right'))`, the right group 340 px wide
  (`ShellOptions.size`).
- **Invariants** (`normalize`, after every operation): exactly one main area, in the tree (put back at the left when
  missing); no empty group (one holding only hidden views stays, to come back to); no split with fewer than two children
  (the only child takes its place and its size); no split directly in a split of its direction (its children join the
  parent); every child but the one taking the rest has a size of at least `MinGroupSize` (180; default
  `DefaultGroupSize`, 280); a group's `active` one of its views that is not hidden; each view in one place; window
  rectangles at least `MinWindow` (240 × 160), windows renumbered by stacking order; ids unique. Why: the view renders
  any tree it is given; every rule it relies on is checked once, here. A stored size is clamped here too, by its level
  (`validateLayout` reads it as stored): a pane's least is not a docked group's.
- **The index** (`layoutIndex`): what readers ask of a layout value (a node, its parent, a view's or a group's place,
  the rail group, the right column, the region, a group's column), built at the first ask and kept with the value
  (immutable). `findView`, `findGroup`, `inRegion`, `railGroup`, `columnOf` … read it. Why: a drag writes a layout a
  frame, and every tree node, group, splitter and port asked again by walking it.
- **Instances**: a view is `single` (default: one place; `show` brings it where it is) or `many`: each instance is a ref
  `id#key` (`viewRef`, `viewIdOf`), its own place, size, `subject` and section focus, opened by
  `show(id, { key, subject })` (key `default` when none) and gone on close. Instances are not placed at registration;
  the layout keeps their keys and the app maps them back to subjects, or drops them.
- **Homes** (`home`): a group id (default `right`): placed there at registration by `order` (unless the layout has it: a
  stored layout keeps where the person put it); a group that is not there is made beside the main area, on its right,
  340 px. `float`: a window, at `show`. `card`, `page`, `modal`: a layer (§3), not in the layout. `show`'s `target`
  places it elsewhere for this once. Placing is not opening: a group keeps showing what it showed, folded or not.
- **Operations** (`ShellModel`, named for what a person does):
  `show(id, { key, subject, section, anchor, anchorRect, target })` (a layer opens; a placed view becomes its group's
  active one, the group opens, a window raises, a hidden view comes back, a rail group's view shows the rail; an
  unplaced one goes home; a docked `show` closes the card); `select(ref)` (a press on a strip button: shows it, or folds
  its group when it is the open one, not in a window); `collapse(group)`; `move(ref, target)` (`group` at an index,
  `edge` of a group or the main area: a new group on that side, `float` with a rect); `float(ref, anchorRect?)`;
  `dock(ref)` (back to its home group); `merge(from, into)` (a window's views into another group); `splitBelow(ref)` (a
  new group of 240 px under its own, when its group has another view); `moveInColumn(ref, ±1)`; `close(ref)`;
  `resize(node, px)`; `setWindowRect`; `raise`; `hide`/`unhide`; `toggleRail`; `resetLayout`. A group left empty goes; a
  split left with one child collapses.
- **Close** (`close`): a layer's view leaves its layer; a docked single view is hidden (its place kept); a single view
  with a docked home closed from a window goes back home hidden (closed, not docked); a `float` view, an instance or a
  view in the main area's region leaves the layout, and its subject and focus go.
- **Every layout write** (`#setLayout`, one place) keeps what refers to the layout true in the same batch: the maximize
  and the active pane (§2.1); what was kept for a ref that left the layout and every layer goes (its port's states, its
  own title): a port an app still holds keeps working, a new one is made on demand.
- **Visible flags** (`visible`): while false the view is `absent`: it leaves its strip and the rail, keeps its place,
  and comes back there. `shows(ref)` is registered, not absent, not hidden.
- **Hidden views** (`hide`, `unhide`, `hiddenViews`; `layout.hidden`): out of view by the person's choice, without
  unregistering; its place kept; part of the layout, so a stored one keeps it hidden. `show` brings it back: asked for
  by name, it shows.
- **Window geometry** (`layout.ts`): a new window opens at its last rectangle this session (`remember` not false), else
  beside `anchorRect` (to its left, top level with it, 8 px apart; without room, below it), else cascaded from the main
  area's top-left (16 px in, 24 px a step); embedded, at 0.8 of its size; always clamped (`clampRect`): no larger than
  the main area and inside it. The main area's size (`setArea`) clamps every window again, and so does `setWindowRect`.
  Why: a window's header stays reachable however the shell shrinks.
- **Persistence**: the app stores the layout; the shell never touches storage. `save()` hands the value out (views with
  `remember: false` taken out, so they open at their home and default size each session); `restore(value)` takes one
  back: validated field by field (`validateLayout`: another version or a wrong shape is refused and the current layout
  stays), windows clamped, registered views it lacks placed at their homes, unknown refs kept (a view registering later
  lands where it was left). `ShellOptions.layout` starts from a stored one. Why: whether a session remembers the layout
  is the app's choice.
- **Reset** (`resetLayout`): the default layout, every registered view at its home, nothing hidden, no windows, the rail
  shown, the session's last window rectangles forgotten.
- **Modes** (`mode`): `desktop`; `compact` (`auto`, the default, below `compactWidth`, 720, of the shell's own width,
  `setWidth`); `embedded` (by option only): the rail starts hidden, windows open smaller.

### 2.1 The main area's region: the split stage

The main area is itself docked: a region of the same groups and splits around the app's content. A stage is a view; a
sequence, a plot, a diagram dock beside or above it. Why: one layout engine (normalize, sizes, save and restore, the
kept bodies), and a split stage, a strip above the content and a panel beside it are the same operation.

**The layout and the app meet at ports** (why: the shell is generic, the app's state stays the app's; neither reaches
into the other):

- **Pane port, layout to pane** (`PanePort`, `panePort(pane)`, one a pane; a view gets it as `ctx.pane`, the app's
  content as `main`'s `pane`): `id` (`content` or a view's ref); `status` (`PaneStatus`: `active`, `hovered`, `focused`,
  `maximized`, `several`: several shown, not one alone or maximized), one state, told only when a field changes;
  `insets` (`PaneInsets`: how much of its bottom-left and bottom-right the shell's `over` chrome covers, for the pane in
  that corner of the main area, both for a pane shown alone; none in a window); and what a person does to it:
  `toggleMaximized`, `float`, `close`. **Pane to layout**: only `setTitle` (an instance's own title, "Stage D": frames,
  tabs and windows show it, `paneTitle`). What a pane shows on its status or insets (its toolbar while active, a bottom
  bar kept clear) is the pane's. Hover and focus are the view's inputs (`setPaneHovered`, `setPaneFocused`), the `over`
  chrome's widths too (`setOverChrome`), fed by the region and by each window (and a window's pane maximized); headless
  they stay off. A splitter is the layout's: pressing, focusing or hovering it changes no pane (resizing is not
  choosing).
- **What a pane holds is the app's** (`PaneMap<T>`, `pane-map.ts`, app side, DOM-free): the app keeps its own entry a
  pane (a stage's model, a plugin), keyed by pane id; `active` is the active pane's entry, from the shell's `activePane`
  (any state of an id). A shared toolbar acts on `panes.active.value` and posts its result with the pane as the
  message's scope (status §1). The shell never stores a pane's content.
- **Layout**: `MainNode.root?`, the region's tree: groups and splits around the one `content` leaf (`ContentNode`, the
  app's `main` render). No tree: the content alone (every layout stored before the region reads so). `MainNode.top?`:
  the group overlaid on the content's top (`TopGroup`). `mainNode`, `inRegion(group)`, `regionHas(node)`.
- **Invariants** (`normalize`, the main area's rule one level down): exactly one `content` in a region with a tree, none
  outside it, put back at the left when lost; the node holding it takes the rest; a tree that is only the content is
  none; panes are at least `MinPaneSize` (64: a strip above the content is a few rows, not a panel).
- **Sizes within the main area** (`setArea`, `resize` of a pane): in each split, the sized panes leave the one taking
  the rest its least (its even share of the room, at most `MinMainSize`); over it, the pane being resized (and its grid
  peers) gives way first, then the rest shrink in proportion, each at least `MinPaneSize`. Why: a main area made smaller
  must not squeeze the content to nothing; a lined-up grid shrinks line by line and stays one.
- **Stage views** (`ShellViewDef.stage`): a canvas, not a panel: alone in its pane it has no header (the frame is named
  by its title); with tabs, its tabs. One canvas a stage view or viewports of one canvas is the app's choice (a page has
  about 16 WebGL contexts); the shell only gives each its body and size.
- **Homes**: `stage` (beside the content: `stageSide`, default right, `stageSize`) and `top` (the shell's `top` option,
  `TopMode`: `split`, the default, above the content, 96 px; `overlay`, the top group; `both`, the top group, "Above the
  main view" offered too). Placed at `show`, not at registration.
- **The top overlay**: the app places it (`main`'s `top`: a stage's band slot, `<Stage top={top()} …/>`); its view shows
  its body there, no header (tabs when several).
- **Operations**: `split(ref, side, of = 'content', size?)` (beside the content, a region group, or a view's pane; a
  placed view moves, from the right column, a window, another pane); `toTop(ref)`; `maximize(pane)`,
  `toggleMaximize(pane)`, `maximized`: one pane alone in the main area, Escape back; a region pane (`content` or a
  group's id: the others away and paused), or a window's view (its ref: out of its window over the main area,
  `MaximizedHost`, its body kept and back in its window after; the other windows inert and paused). The shell's chrome
  stays: maximize is the pane's, expand the shell's (§3.9). Panes show and hide at once, no slide (the slide is a
  panel's, folding beside the main area; panes sliding out in turn drew every layout between). It ends at any layout
  write that leaves its id no region node nor windowed view (closed, floated, docked, restored without it), when a view
  is shown or moved into the region outside it (it would stay hidden behind it; the top overlay on a maximized content
  shows), and on reset. `activePane`, `setActivePane` (the pane last pressed or focused, in the region or a window;
  `paneNode` maps a region pane to its node): back to `content` at the write its ref stops being a pane (closed, docked
  out of the main area). What maximizing a pane maximizes is `maximizeTarget` (its region node, its ref in a window,
  else `content`): the port and `view.maximize-pane` (enabled with several panes or the active one in a window).
  Commands: `view.maximize-pane`, `view.close-pane`, `view.grid-panes` (the active pane).
- **Close**: a view closed in the region leaves the layout; shown again, it goes to its home.
- **Windows** may hold a stage view (a detached 3D window): its body moves between the window and the region without a
  remount, its chrome its own; a press or the focus in it makes it the active pane.
- **Menus**: a pane's (`paneMenu`, its header's right-click): Maximize or Show every pane, Split right and below (a
  group of several), onto or off the top (by the `top` option), Arrange as a grid, Float, Dock right (a view whose home
  is the right group), Close. A view elsewhere adds Beside, Below, Above the main view, Over the main view's top.
- **What the base leaves to the app** (the shell and the stage are a base an app composes, not an app): the rail
  toggle's mark (`railToggleIcon`, default the Mol* star), the shell's active-pane mark (`paneMark`, default on), the
  stage's empty mark (`Stage`'s `emptyMark`, default the Mol* glass), what a pane holds (`PaneMap`), what a message is
  about (status scopes), and whether one button maximizes a pane and expands the shell (both, called by the app).
- **Toolbars are the app's to place** (not decided here): each stage keeps its own slots (a toolbar per pane, shown from
  its port's status: `Stage`'s `toolbarShown`; its bottom bar clear of the shared chrome: `Stage`'s `insets`), and
  `ShellView`'s `over` has a status and a corner once over the whole main area (a shared toolbar follows the app's
  `PaneMap`). The shell's own marks: a hairline between panes, and an accent line along the active pane's top edge while
  several show (the open tab's mark). Storybook shows per pane, shared, mixed and active only.
- **A grid** (`SplitNode.grid`): a split whose children are splits of the other direction with as many children each,
  lined up: resizing a cell resizes the cell at its place in every line, so a line across the grid moves as one. Kept
  only while the shape holds (normalize drops it when a pane is split or closed). A cell's splitter lights its line's
  others with it (hovered, focused, dragged: `data-linked`, `gridPeers`), so the line across reads as one.
  `arrangeGrid(columns?)` (the pane menu's Arrange as a grid, `view.grid-panes`): every pane, the content first, in
  columns filled top to bottom, equal sizes over the main area (asked before the main area is measured: again at its
  first size); a full grid is lined up, one with a short column is not.
- Not yet: dragging a view onto a pane's edge (the layout menu moves views).

## 3. The views

`ShellView` (`shell-view.tsx`): the root (a size container), the docked tree or the compact layout, then the layers. The
main area (`[data-msui-shell-main]`) is made once and moved between the desktop and compact layouts: the app's view
keeps its state (a canvas is not rebuilt).

### 3.1 Views and their bodies

- **One frame** (`ViewFrame`): a header (the title or a group's tabs, the view's `actions`, a ⋯ layout menu, a ✕) over
  one body (`viewBody`: a `ScrollArea`, or with `fill` a box the view fills). The same in every place (`placement`:
  `docked`, `window`, `card`, `page`, `modal`, `sheet`; the view reads it from its context). Why: a view looks the same
  wherever a person puts it.
- **The tree is rebuilt only when its structure changes** (`sameShape`: nodes, order, direction, the rail's group);
  which nodes show and which are folded, which takes the rest, the sizes and the splitters' sides are written in place
  (`applyState`, a size only when it changed; whether a node shows asked once a node). `applyState` runs only when what
  it reads changes (its tree's nodes, sizes, folds and views, what shows, the rail, its sources): a window's move and
  the other tree's drag skip it. What the rebuild makes anew is marked again after it (the active pane, `onBuild`). Why:
  rebuilding moves the main area (the app's canvas), which drops the focus in it and ends a splitter's drag; folding a
  split column or hiding the rail must not.
- **What is derived follows every source it reads and changes only when its result does** (the rail's buttons and fit, a
  window's rectangle, the bodies kept): a drag's frames write the layout, and only what the drag moves is written.
- **Bodies are kept across moves** (`ViewBodies`): a view's body is made once, in a scope of its own, and handed to the
  frame that shows it now; released by every frame it is paused and detached, not disposed. It is disposed when the view
  leaves the layout and every layer. Why: a scroll, a draft or a plot's zoom survives popping out, docking back and
  folding.
- Every frame has a ✕: a docked view's folds its group, a window's, card's, page's and modal's close it. The layout menu
  (Float, Dock, Split below, …) is a right-click or Shift+F10 on the header (and on a tab or a rail button); a ⋯ button
  for it only with `ShellView`'s `viewMenu`. Why: closing is what a header's button is for; moving a view is rarer.

### 3.2 The rail

- **The rail group** (`railGroup`) is the top group of the docked column at the shell's right edge (the last child of
  the root row, or its column's first). Its strip is the rail: a column of icon buttons, full height at the shell's
  right edge, 40 px wide; tooltips name them (`data-tip`).
- **Order** (`railEntries`): `head` views, the `strip` views, rail items by `group` (in the order of each group's first
  item, by `order` within), `foot` views, items without a group, the rail's Hide button, then `corner` items (the
  shell's bottom-right corner; kept when the rail runs short); sections a hairline apart, the foot at the bottom.
- **Steps by the shell's height** (CSS container queries, nothing measured): buttons 32 from 620 px, 28 from 520, 24
  below; a captioned item a caption (12) taller. `RailSteps` (`rail.ts`) holds the same numbers for `fitRail`, kept
  equal to `shell.module.css` by its test.
- **More**: when even the small step does not fit, `fitRail` takes out the entries with the highest `priority` (default
  their position, so the foot goes first; later first on a tie) until the rest and a More button fit; More opens them as
  a menu. Measured by one `ResizeObserver` on the rail, recomputed when the rounded height or the entries change.
- A press on a view's button selects it, or folds the column when it is the open one (`select`); a folded rail group's
  node does not show, the rail stays. The open view's button is its background only (`aria-selected`).
- A view button's context menu (or Shift+F10, the menu key) is `railMenu`: Hide, Show for each hidden view, Hide the
  rail, Reset layout (Show brings an instance back as itself, its key kept). A tab's (Shift+F10 too) is the menu of
  where it is: the view's docked or in a window, the pane's (`paneMenu`) in the main area and its top overlay.
- **The rail's toggle** (`railToggle` prop, default on): a Hide button (») ends the rail's foot; with `layout.rail`
  false the rail and every group of the right column are gone, and a star button in the main area's top-right corner
  (`RailToggle`, not in compact) shows them again. The same in embedded, where the rail starts hidden. Showing a view of
  the rail group shows the rail too.

### 3.3 Docked groups and splitters

- The tree renders as nested flex splits; a sized node takes its size (`--msui-shell-size`), the node holding the main
  area (or, in a column without it, the first that shows open) takes the rest. A node with nothing to show is
  `display: none`.
- A docked group is a surface with a hairline on the main area's side, and above it when it is under another in a
  column. The rail group's frame has the rail header (`styles §8`); a **lower group** (any group not at the rail) shows
  a 32 px header with a row of labelled tabs when it has several views (a `tablist`; a press on the open tab folds the
  group to its header), its title otherwise.
- The view menu (`viewMenu`) of a docked view: Float, Split below (its group has another view), Move up / Move down (in
  a column), Dock right (away from its home `right`), Move to each open window, Hide (single), Reset layout.
- **Splitters** sit on a sized node's edge toward the node that takes the rest (`role="separator"`, focusable): a drag
  (in the frame lane) or the arrows (16 px, Shift 64) resize it; the size is at least `MinGroupSize` and at most what
  leaves the main area `MinMainSize` (240) beside the rail; a narrower shell, or a stored size, gives way the same.

### 3.4 Windows

- A window (`WindowView`) floats over the main area, never over a docked group: a raised surface with an edge and a
  float shadow, a 32 px header (its title, or tabs when it holds several views), `role="dialog"`, `aria-modal="false"`.
- Its header drags it (not from a button, tab or field), its right edge, bottom edge and corner resize it, both in the
  frame lane, the release's place committed through it too; it stays clamped in the main area (§2). A press or focus in
  it raises it; the windows' nodes keep a stable order (by id), so a raise moves none (moving one would drop the focus
  in it): the stack is each window's `z-index`. A window none of whose views can show is hidden.
- ✕ or Escape (with the focus in it) closes its views (`close`, §2); a press outside does not. Why: a window is kept
  beside the work, not dismissed by it.
- Its menu: Dock, Move with the keyboard, Resize with the keyboard (arrows by 4 px, Shift 16; Enter or Escape ends;
  outlined meanwhile), Reset position (centered at its default size), Move to another window, Reset layout.

### 3.5 Tooltips

- The shell installs one tooltip layer over itself (`installTooltips`, the `tooltips` prop, default on), so the app's
  views and extensions' get tooltips without installing their own.
- Its docked groups, windows, card and sheet are tooltip edges (`data-tip-edge`): a row's card (`ViewRow`,
  `data-tip-row`) goes outside the view it is in, on the side with room (a docked group's: toward the main area), never
  over the view; with no room on either side, none. A row's card shows from a pointer over its text or the keyboard's
  focus on the row.

### 3.6 Layers

- **Card** (`home: 'card'`): one at a time, in the right column's place (its width, the rail beside it; with no column,
  340 px at the right edge), full height, over the docked groups. Closes on its ✕, Escape, a press outside (not a
  `persistent` one), opening a page, or showing a docked view. The groups it covers are inert and paused while it is
  open. Opened from the keyboard, its first control gets the focus; else the card itself (§4).
- **Page** (`home: 'page'`): opaque over the main area (`role="region"`), its own header; Escape or ✕ close it; opening
  it closes the card. What it covers (the app's main view, the windows) is inert, the windows paused.
- **Modal** (`home: 'modal'`): centered over the whole shell on a scrim, at most 480 px wide, `role="dialog"`,
  `aria-modal="true"`; Tab stays in it (`trapTab`), the rest of the shell is `inert` (its layout wrapper, whatever is in
  it as it changes); Escape or ✕ close it. For a confirmation or a short task.
- **Busy** (`busy` prop, a `ShellBusy`): a scrim and one card holding the app (`role="alertdialog"`): the label, a
  progress bar (indeterminate without `progress`), Cancel when it has `cancel`. The rest of the shell is inert; the
  focus goes to it (Cancel, else the card), stays in it, and goes back when it ends.
- **File drop** (`onDropFiles`): files dragged over the shell show one overlay (`dropLabel`, and `dropFormats` under it:
  what can be dropped); the cursor says copy; dropped, they go to `onDropFiles`, unless a view under the drop took it
  (prevented its default). Without it, no overlay and no handling.
- **No status slot**: the main area's chrome (status, chips, its toolbar) is the stage's (`stage/`), placed in the main
  area by the app.

### 3.7 Compact and embedded

- **Compact**: the tree is not rendered; the main area sits over a bottom sheet (`min(45cqh, 420px)`) and a bar of 44 px
  buttons (a touch target; icon and label). The sheet shows the card, else the view picked last (from the bar or More,
  by a `show`, a window opened or raised: a window is a card here); a bar button picked again closes the sheet; its ✕
  closes the card or the window, or the sheet. The bar holds the `compactBar` views; the rest are under More. No windows
  layer, no rail toggle.
- **Embedded**: the desktop layout with the rail starting hidden and windows at 0.8 of their size; it does not switch to
  compact.

### 3.8 Motion

On the motion vocabulary (styles §7.1), none at the first render, none under reduced motion (every duration 0):

- A view replacing another in a group, a window or the sheet fades in (`data-reveal="fade"`).
- A docked group that opens slides in from the rail's side (`side`); one that folds plays it back (`Motion.close`)
  before its node hides. The main area takes or gives the room when the slide ends.
- A card slides in (`side`) and back out; a page, a modal and the busy overlay fade; the rail slides in when shown.
- A window opened after the first render floats in, rising a unit (`msui-float-in`).

### 3.9 Expand and Fullscreen

- **Expand** (`ShellModel.expanded`, `expand`, `collapseExpanded`, `toggleExpanded`; `view.expand`): the shell's
  presentation, global: its own chrome steps away (the rail, the docked column, the windows: hidden, inert, each
  window's frame paused) and the main area takes the shell, showing what it shows (every pane, or the maximized one); a
  card, page or modal still opens over it. Escape steps back (after anything inside took the key): out of Expand, then
  out of a maximized pane; showing a docked view or a window leaves it too. What is open, not the layout: not saved. Why
  separate from maximize (§2.1): one pane alone is the layout's about a pane, the shell's chrome away is the shell's
  about itself; an app wanting both from one button calls both.
- **How far** (`ShellView`'s `fill`): `shell` (default: the main area takes the shell's box) or `page` (the shell over
  the page while expanded, `position: fixed`, nothing changed on the page: a viewer embedded in a page).
- **Fullscreen** (`fullscreen`, `fullscreenRequest`, `fullscreenAvailable`, `requestFullscreen`; `view.fullscreen`): the
  model asks, the view does it on the shell's root and writes back what the browser did (its own Esc too). Where it
  cannot be (no Fullscreen API, an iframe that does not allow it) `fullscreenAvailable` says why, the command is greyed
  with it and asking does nothing. Headless: unavailable.
- **Commands** (`shellCommands(shell, category?)`, `shell-commands.ts`): `view.expand`, `view.fullscreen`,
  `view.toggle-rail`, `view.reset-layout` (and §2.1's pane commands), titled and answered in `shell.labels` (§6); the
  app registers them with its command manager. Every person-facing operation of the shell is a model method; these are
  the ones a person asks for by name.

## 4. Keyboard and accessibility

- **Where the focus goes when a layer opens** (`focusInitial`, WAI-ARIA dialog): a modal focuses the first control of
  its content (a confirmation lists its least destructive action first), never its header's ✕ (the next Enter would
  close it); a card opened from the keyboard likewise, one opened by a press focuses itself (no ring); a page focuses
  itself. Closed, the focus goes back to what opened it.

- **Roles**: the root `application` (`aria-roledescription="workspace"`, labelled by `label`); the main area a `region`
  ("Main area"); a view a `section` labelled by its `h2` title (or by `aria-label` when tabs replace it); the rail a
  `nav` whose view sections are vertical `tablist`s and whose item sections are `group`s; a group's or a window's tabs a
  `tablist` of `tab`s; splitters `separator`s; a window a non-modal `dialog` (named by its pane's title, as it changes);
  the modal a modal `dialog`; busy an `alertdialog` with a `progressbar`; `ViewToolbar` a `toolbar`; `ViewError` an
  `alert`.
- **Roving focus**: the rail (vertical) and each tab row (horizontal) are one Tab stop, the arrows between their buttons
  (`rovingFocus`).
- **F6 / Shift+F6** (on the root) cycle the regions that show: the main area, each docked group, the rail, each window,
  the card, the compact bar; the focus goes to a region's first focusable (the main area itself). The shell marks its
  own (`data-msui-region`, under its root's `data-msui-shell`): not found by a label (translated) or a tag (an app's
  `nav` in a view), nor in a shell nested in a view.
- **Ids**: an element id the shell makes starts with its model's `uid` and encodes a ref one to one (`idPart`): two
  shells on a page, or refs `a#b` and `a_b`, never share one.
- **Escape** goes through `dismissable` (`components/popover/dismiss.ts`): one layer per press, the innermost (last
  opened) first: a menu or list, then a window's keyboard adjust, the window (only from inside it), the card, the page,
  the modal.
- **Menus from the keyboard**: a rail button's menu by Shift+F10 or the menu key; a view's layout menu by its ⋯ button
  (`aria-haspopup="menu"`).
- **Focus return**: closing a card, page or modal returns the focus to where it was (`rememberFocus`, else the anchor it
  was opened from); closing a window returns it to the main area.

## 5. What is the app's

The main area's content and its toolbar (its corner is kept clear only by the rail toggle), the status and the log (the
status slot; a log is a view), Home and every other view, the rail items, the commands that call the operations, where
the layout is stored, and whether a session restores it.

## 6. Labels

Every word the shell shows or says (its chrome, its menus, its commands, the fullscreen reasons, the rail's sections:
`railViews`, `railTools`, never their internal names) is in `ShellLabels`; `DefaultShellLabels` is English. Strings, or
functions for those with a title in them (`close(title)`, `moveToWindow(title)`). The app passes `ShellOptions.labels`
(a `Partial`), merged once in the model (`shell.labels`, frozen), so the view, `viewMenu` / `paneMenu` / `railMenu` and
`shellCommands` read the same words, headless too. Why: the shell is a base library; an app in another language replaces
its words without forking it. A view's own words are the view's. Search aliases of the commands stay English.
