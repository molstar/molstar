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
  segmented picker (icons or text in a well), an iso level control and a transfer function editor over a shared
  histogram plot (values in σ or absolute) and the params form with its model (no rxjs). Styles compile from CSS
  modules into registered modules (`pnpm -C packages/ui styles`, `styles:check` keeps spacing on the 4px grid);
  corners are squircles where supported, nested corners concentric; icons take their size from their container.
- Fix `PD.ObjectList` ignoring its `ctor` option: a list's "Add" now creates the item it gives, not the element's
  defaults.
