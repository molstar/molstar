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
