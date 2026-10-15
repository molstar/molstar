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

## Renderer bench

`_stories/runtime-bench.stories.tsx` (Internals › Performance › Runtime): the js-framework-benchmark operations on the
`dom` runtime, rows written as the app writes them (a record per row, keyed `For`, the selected row by a keyed
selector). Script is the operation alone; medians of 5 after a warm-up, slicing off; Chrome, outside Storybook's root.
Total adds layout and two frames (~33 ms of it).

| Operation                        | Script, first | Script, now | Total, now |
| -------------------------------- | ------------- | ----------- | ---------- |
| create 1,000 rows                | 9.9 ms        | 6.3 ms      | 50 ms      |
| replace all 1,000 rows           | 14.8 ms       | 9.9 ms      | 50 ms      |
| update every 10th row            | 0.3 ms        | 0.3 ms      | 33 ms      |
| select a row                     | 0.4 ms        | 0 ms        | 33 ms      |
| swap rows 2 and 999              | 4.8 ms        | 0.6 ms      | 33 ms      |
| remove a row                     | 0.8 ms        | 0.4 ms      | 33 ms      |
| create 10,000 rows               | 116 ms        | 75–82 ms    | 430 ms     |
| append 1,000 rows (to 10k)       | 13 ms         | 9.8 ms      | 100 ms     |
| clear 10,000 rows                | 47.7 ms       | 30 ms       | 36 ms      |
| relabel 1,000 rows 10× (a burst) | 12.9 ms       | —           | 34 ms      |
| the same burst in `batch`        | —             | 2.4 ms      | 33 ms      |

What changed, by what it bought:

- **Scopes allocate on first use** (children set, cleanup list), a disposing parent is not edited by its children, and
  bindings keep their subscription instead of a closure (`Scope.own`): clear −29%, create −20%.
- **`For` makes an index state only when `children` takes it**; an emptied list goes in one range (measured 2.5× faster
  than emptying the parent with `textContent`, which was tried and reverted).
- **Static props skip bindings** (a plain value is written once; a plain `class` is the attribute), `for…in` over props,
  text alone is `textContent`: create 10,000 −19%.
- **Reactive children start as one Text node** (anchors only when nodes arrive).
- **`keyedSelector`** (`reactive/keyed-selector.ts`): the selected row of thousands notifies two rows, not all.
- **Delegated events** (`dom/render.ts`): one listener per type per render root instead of one per handler; each handler
  sees its own element as `currentTarget`, `stopPropagation` ends the walk, nested roots never run a handler twice, an
  event made not to bubble reaches only its target; focus, enter/leave, touch and wheel stay native. Clear 10,000 −30%
  (elements without listeners go faster), no listener per row.
- **`batch(fn)`** (`reactive/observable-state.ts`): states set inside emit once, with their final value, when the
  outermost batch ends (synchronous; reads inside see the new values). A burst of ten relabels per row: 12.9 → 2.4 ms. A
  microtask that defers DOM writes was not taken: code reads the DOM right after a change (focus, measuring, tests).
- **`For` places with the fewest moves**: entries whose old places form the longest increasing run stay; a swap is two
  moves, not a cascade (swap −85%); covered by a 200-round shuffle/insert/remove test.

The floor: the same rows with no bindings create 10,000 in ~50 ms; what is left is element creation itself (a compiled
template could clone rows instead, a build step the runtime does not have).

## Next

The bench as a regression guard: run it before changing `dom/` or `reactive/` (Internals › Performance › Runtime; a
visible page for totals, a hidden one gives script times only, about twice as slow). What is left is element creation
itself (~50 ms for 10,000 rows without bindings); a compiled template that clones rows would be the next step, a build
step the runtime does not have.
