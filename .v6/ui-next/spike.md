# Renderer check (plan.md Phase 2)

The state tree at 5,025 cells (25 structures × 4 models × 49 components) on the `dom` runtime, in
`_stories/performance.stories.tsx` (Internals › Performance › State tree). Wall time in Chrome until the result is all
there, and the longest task in that time (how long the page could not respond).

## Results

| Case                                                       | Total  | Longest task |
| ---------------------------------------------------------- | ------ | ------------ |
| State tree tab: first render (sliced)                      | 513 ms | none > 50 ms |
| State tree tab: select a cell (5,025 rows re-check)        | 33 ms  | none         |
| State tree tab: 20 relabels as immutable records           | 481 ms | 66 ms        |
| State tree tab: collapse and expand all 25 structures      | 328 ms | 311 ms       |
| Plain rows that follow their record by id: first render    | 200 ms | none         |
| Plain rows: 20 relabels                                    | 33 ms  | none         |
| Virtualized tree (`VirtualList`), per-node records: any op | 33 ms  | none         |

(33 ms is the measuring floor: two frames.)

## What it took, and what it showed

- **Measure outside Storybook's root.** The accessibility addon checks every change under `#storybook-root`; at 5,000
  rows it turned a 0.6 s render into 25 s. The check renders into a fixed, contained area on the body.
- **Time slicing in `For`.** A long list (or a tree: short lists, long whole) renders in 12 ms slices shared by every
  `For` in the task, continued with Mol\*'s `setImmediate` (a MessageChannel task), so the page never blocks; existing
  entries update at once. 12 ms measured best: 5 ms slices took ~40% longer in all; frame-aligned 12 ms was the same.
  Off in jsdom (`timeSlicing`), where nothing renders.
- **Rows skip off-screen work.** `content-visibility: auto` on tree rows took layout from 554 ms to 1 ms and frames from
  357 ms to 55 ms.
- **Records, not rebuilt trees.** `For` renders an entry again when its record changes, so an immutable tree rebuilds a
  root's whole subtree for one relabel (~24 ms each); rows that follow their own record by id do 20 relabels in a frame.
  Sources hand out per-node records (view-models.md §3).
- **Virtualize long lists.** Only the window is rendered: 5,025 cells cost what 40 do; expanding recomputes the flat
  list; scrolling swaps a few rows. The production scene tree, state tree, and long lists (remote states, History) use
  `VirtualList`.

## Decision

Keep the `dom` runtime (option A). None of the costs above are the renderer's to remove (Preact would need the same
slicing, records, and virtualization at this size), and the mocks proved its authoring. The Preact side of the gate was
not built.

## Next

A renderer bench (the js-framework-benchmark operations: create 1,000 / 10,000 rows, replace, update every 10th, select,
swap, remove, append, clear) measuring script time apart from layout and paint, medians of several runs, to make the
runtime itself as cheap as it can be and keep it there (a regression guard).
