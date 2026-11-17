## UI: DOM runtime
- `@molstar/ui/dom` (a new package, `@molstar/ui`): a JSX DOM runtime on `@molstar/core/reactive` (`"jsxImportSource":
  "@molstar/ui/dom"`): no virtual DOM; components run once, states update the nodes they bind. Scopes with pause and
  resume, `watch` (one render a microtask, parents first), bindings, `Show`/`Switch`/`Dynamic`/`For` (keyed,
  time-sliced), `activity`, the frame lane (`FrameState`, one animation frame a page), work queues and `latest`,
  scoped listeners, observers, timers and `drag`, `ViewElement`, focus and tree keyboard helpers, `DomProbe`. A
  `watch` or a binding on states is a sink of the state graph: one object whatever the number of states, nothing
  allocated to mark it. Nothing uses it yet.
- `@molstar/ui/styles/foundations.css`: the UI's foundations, design tokens (`--msui-*`: color by theme with
  `data-msui-theme`, type on a 4px baseline grid, spacing, sizes, shape, motion) and alignment primitives (`.msui-*`:
  one inset, icon slots, a level that is a slot, one label column, values where controls write, numbers by place).
  Storybook: a new UI group with an introduction and the foundations' stories.
- `@molstar/ui` components on the runtime and the foundations: icon, button, tooltip, scroll area, text input, select
  (opens over its field), popover, context menu, color picker, choices (a radio group when it shows a selection,
  labels never cut), virtual and recycled lists, drag to reorder, legend, markdown, action footer, logo (the Mol*
  wordmark), a generic toolbar (groups as blocks or on one surface with separators, menu buttons, submenus), a
  segmented picker (icons or text in a well), skeleton placeholders for what loads (crossfading into their content),
  an iso level control and a transfer function editor over a shared histogram plot (values in σ or absolute) and the
  params form with its model (no rxjs). Styles compile from CSS modules into registered modules (`pnpm -C packages/ui
  styles`, `styles:check` keeps spacing on the 4px grid); corners are squircles where supported, nested corners
  concentric; transitions are one vocabulary (floats rise from their anchor, folds grow and shrink, swaps crossfade in
  place); icons take their size from their container.
- Fix `PD.ObjectList` ignoring its `ctor` option: a list's "Add" now creates the item it gives, not the element's
  defaults.
- `@molstar/ui/react`: React bindings over core's reactive sources (`useReactive`, `useReactiveProp`,
  `useReactiveModel`), `reactView` for React or Preact views, and `UiElementHost`. A state is read through its `value`
  and followed by its changes only, so subscribing never notifies React's store; React is an optional peer dependency.
- `@molstar/ui/commands` and the command menu (`components/command-menu`): commands over PD params with keys, groups,
  presets and formats; in the menu, the command's sentence is its editor, the list under it always the active field's
  choices, and the rest of the params in an Options fold.
- `@molstar/ui` review fixes: overlays close by one rule (`popover/dismiss.ts`: the innermost layer takes Escape, Tab
  leaves menus and lists, the focus goes back to what opened them); a Select works uncontrolled; popovers follow their
  column's size; disabled reasons are read by screen readers and arrows reach disabled menu rows; toolbars and tabs
  follow `dir=rtl`. The params form keeps an open list item on the item it edits, gives intervals without bounds number
  fields, clamps to a single bound and shows what it stores (the step's decimals). Typed colors apply on Enter or blur;
  markdown loads only data: and same-origin images unless `remoteImages`. The runtime disposes and resumes every child
  scope, writes bindings parents first, publishes frame states in one batch, and creates `a`, `title`, `style` and
  `script` in SVG in the SVG namespace.
- `@molstar/ui` expandable panels (`components/expandable`): a header that opens a body under it, open a well holding
  both; the body renders on first open and stays, inert and paused, while closed; opening grows it in, closing plays
  back; open state kept by the panel or by the host (a `State`); `flat` for inside a surface.
- `@molstar/ui/shell`: a generic, app-independent shell. Views (what were tabs, cards, dialogs, pages) are registered
  once and placed by one docking model: groups in a split layout beside the app's main area (the right column, its
  top group's views on a rail at the edge, a lower group with tabs: Split below, Move up/down) or floating as windows
  (move, resize, raise, keyboard Move/Resize, several views as tabs, several instances of a view); cards in the right
  column's place, pages over the main area, modals over everything. Views keep their state while they move and pause
  while hidden; the rail steps down with the height and overflows into More, hides behind the Mol\* star; compact (a
  bottom sheet) and embedded modes; F6 between regions; the layout is a validated value the app saves and restores.
  View primitives (`ViewStack`, `ViewToolbar`, `ViewEmpty`, `ViewError`, `ViewFooter`, `ViewHeading`, `SectionList`)
  and the shell foundations (styles §8) keep contributed views aligned. Storybook: UI › Shell (Shell, View, Rail,
  Window).
- `@molstar/ui` glass wordmark (`components/logo/glass-wordmark.tsx`, `glass.ts`): the Mol\* wordmark cut into glass
  pieces that assemble each time it shows (a view's heading: the shell's Home), lean under the pointer and shatter on
  a press. It asks for frames only while something moves (none at rest, none while a resting pointer holds the pieces
  leaned, none while hidden) and a frame allocates nothing (guarded by `test:perf`); the cut and the pieces' atlas are
  shared by every glass. The logo's colors are foundation tokens (`--msui-logo-*`).
- `@molstar/ui` app parts, every one a DOM-free model that works headless (guarded by `src/_test/headless.test.ts`) with
  an optional view:
  - `log/`: `LogModel` (a bounded, batched log; a filter of kinds and a source, combined, and a text query) and
    `LogView` (All and a toggle a kind, the sources once there are any; recycled one-line rows, kept at the bottom).
  - `status/`: `StatusModel` (messages that fade, errors that stay, task pills with progress and Cancel, chips) and
    `StatusLine`, `StatusPill`, `TaskPills`, `StatusChips`.
  - `history/`: `HistoryModel` (one undo for a person's calls and an agent's turns: steps, groups, undo and redo to a
    step, checkpoints that switch off), `HistoryView`, and `historyCommands` (`edit.undo`, `edit.redo`).
  - `components/playback/`: stepper and timeline models (timed, counted, and durationless: running time or steps, no
    end) with clocks injected, and `PlaybackDeck`, `StepperRow`, `StepCard`, `TimelineRow`, `MuteButton`.
  - `stage/`: a model and a canvas wrapper for the main area (bands such as a sequence on top, a toolbar, the hover
    label, the status in its bottom bar, the empty stage with the glass splash) and placeholders for Storybook (a space-filling three-chain assembly, a sequence band).
  - `forms/settings-view.tsx`: groups of forms with one search (its words may span a group's title and a setting's
    path: "perf qua" finds Performance › Quality; the params form's search does the same over nested groups); `commands/help-view.tsx` and `KeysModel`: keys listed,
    conflicts first, and changed; the command menu's `?` opens the help (`help` prop).
- `@molstar/ui/shell`: Expand (the main area alone, the rest paused; Escape comes back; `fill: 'page'` covers the page
  when embedded) and Fullscreen (available or why not), and `shellCommands` (`view.expand`, `view.fullscreen`,
  `view.toggle-rail`, `view.reset-layout`). A rail item may open a menu on a press (`menu`; the context menu gains radio rows,
  `checked`, and `openMenuBeside`); a `corner` item is the rail's last button, under its Hide button (the command
  menu), and the app's settings (the theme) sit in the rail's foot. The shell's `status` slot is gone: the status lives on the stage. The file
  drop overlay names the accepted formats (`dropFormats`).
- `@molstar/ui/shell` split stage (shell §2.1): the main area is a region of the same docking model around the app's
  content. A stage is a view (`stage: true`, no header alone, `instances: 'many'`): stages side by side or in a quad,
  in a window; views dock beside, below or above the stage (`home: 'stage'`, `split`), or overlaid on its top
  (`home: 'top'`, `toTop`; the shell's `top` option: `split`, `overlay`, `both`; `main` gets `{ top }` for the
  stage's band slot). A pane maximizes (Escape back), the active pane is tracked (`activePane`) and ringed, `over` puts
  chrome once over the main area (a shared toolbar), each stage keeps its own; `paneMenu`, `view.maximize-pane`,
  `view.close-pane`. Stored layouts read as before. Storybook: UI › Stage › Split.
- `@molstar/ui` stage hover label: `hoverPlacement: 'auto'` (the new default, `hoverAt`): the corner on a roomy stage,
  beside the pointer on a narrow or short one (a small pane, a window); a corner label takes its own width (a folded
  toolbar's narrow corner no longer squeezes it).
- `@molstar/ui` split stage: panes have a hairline between them and the active one an accent line along its top. The
  layout gives each pane its status (`paneStatus`, `ctx.pane`, `main`'s `pane`: active, hovered, focused, maximized,
  several) and the pane reacts (a stage's `toolbarShown`: its toolbar only while active or hovered). A grid
  (`arrangeGrid`, `view.grid-panes`, Arrange as a grid): equal cells, lined up, a line across moving as one.
- `@molstar/ui/shell` ports (shell §2.1): each pane gets a `PanePort` (`ctx.pane`, `main`'s `pane`: its status, its
  insets under the shared chrome, expand, maximize, float, close; back to the layout only its title, `setTitle`). What a
  pane holds stays the app's (`PaneMap<T>`: the active pane's entry for a shared toolbar). Maximize is the pane's (a
  region pane alone, a window's pane out of its window over the main area; the shell's chrome stays) and Expand the
  shell's (global: its chrome away). The app chooses the rail toggle's mark (`railToggleIcon`), the shell's active-pane
  mark (`paneMark`) and the stage's empty mark (`emptyMark`). A grid's line lights across as one.
  Status messages take a plain `scope` (`post(…, { scope })`) and a `StatusLine` a `scope` (`{ id }`, `{ rest }`): a
  result shows in the pane it changed, one line speaks. A `Stage` keeps its bottom bar clear of `insets`.
- `@molstar/ui` review fixes (shell, stage, status, views):
  - Words: every view takes `labels` (`ShellLabels` through `ShellOptions.labels`, `StatusLabels`,
    `StageToolbarLabels`, `HistoryLabels` and `HistoryViewLabels`, `LogViewLabels`, `HelpLabels`, `KeysLabels`,
    `SettingsLabels`), English by default and domain-neutral; history refusals carry a `code`.
  - History limits (`capacity`, `checkpoints`, `maxRedo`, `maxRuns`, `setLimits`): what goes past them is dropped with
    its checkpoint; `setRecording(false)` and `recordIf` keep runs out of the history (earlier steps then unsafe);
    `clear`/`invalidate` during a run no longer record a stale step; a long turn's view is linear (latest runs, the
    rest behind Show all).
  - Shell: every layout write keeps the maximized pane and the active pane valid and forgets a gone pane's state; a
    per-layout index replaces the per-frame tree walks, and window moves no longer touch the docked trees; region sizes
    are clamped as the main area shrinks; stored sizes keep their level's minimum; window panes report hover and focus;
    tab menus are the pane menu in panes; F6 finds regions by data attributes; drops a view took are not passed on.
  - Stage: the folded toolbar's submenus open; any CSS color for the background (the engine resolves what the model
    cannot); a hidden toolbar shows when focused; no empty mark by default (`GlassEmptyMark` is opt-in, not bundled
    otherwise); the cursor label reads layout only on resize or scroll; the placeholders are not published.
  - Views: the settings search skips what is not shown and walks once per query; Help's recorder lets Tab through;
    the log forgets sources it no longer holds; Undo's keys follow the keymap and the platform.

