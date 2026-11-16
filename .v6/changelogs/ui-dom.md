## UI: DOM runtime

- `@molstar/ui/dom` (a new package, `@molstar/ui`): a JSX DOM runtime on `@molstar/core/reactive` (`"jsxImportSource":
  "@molstar/ui/dom"`): no virtual DOM; components run once, states update the nodes they bind. Scopes with pause and
  resume, `watch` (one render a microtask, parents first), bindings, `Show`/`Switch`/`Dynamic`/`For` (keyed,
  time-sliced), `activity`, the frame lane (`FrameState`, one animation frame a page), work queues and `latest`, scoped
  listeners, observers, timers and `drag`, `ViewElement`, focus and tree keyboard helpers, `DomProbe`. A `watch` or a
  binding on states is a sink of the state graph: one object whatever the number of states, nothing allocated to mark
  it. Nothing uses it yet.
- `@molstar/ui/styles/foundations.css`: the UI's foundations, design tokens (`--msui-*`: color by theme with
  `data-msui-theme`, type on a 4px baseline grid, spacing, sizes, shape, motion) and alignment primitives (`.msui-*`:
  one inset, icon slots, a level that is a slot, one label column, values where controls write, numbers by place).
  Storybook: a new UI group with an introduction and the foundations' stories.
