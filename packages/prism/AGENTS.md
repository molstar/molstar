# Prism

The new Mol\* UI. The repository's rules ([../../AGENTS.md](../../AGENTS.md)) apply; these are Prism's own. Specs, with
the decisions code comments cite: each package's `_spec/` ([ui](ui/_spec/README.md), [viewer](viewer/_spec/README.md),
[assistant](assistant/_spec/README.md)). What is missing and what was deferred: [TODO.md](TODO.md) (an entry leaves when
it is done).

- Packages: `ui` (the common base: its own DOM JSX runtime, components, the shell, the command menu), `viewer` (Mol\*'s
  viewer on Prism), `stories` (the stories app's start; no modes: separate apps), `assistant` (the AI assistant's UI).
  They are opt-in: an app gets one only by depending on it, and `ui` never imports them
  (`ui/src/_test/opt-in-boundary.test.ts`).
- Not React: views are `.tsx` on the runtime in `ui/src/dom`; state lives in models. A component builds its DOM once and
  renders as a function of state: one `watch(sources, render)` writes what changed, batched in a microtask
  ([rendering.md §1–2](ui/_spec/rendering.md)); structure goes in regions (`Show`, `For`).
- **Reactive layering** (David, 2026-10-17; [rendering.md §2.2](ui/_spec/rendering.md)): models expose `State`s (the
  latest value plus listeners); rxjs only at the boundary (`toState` at a model's edge, `asObservable` out). Views never
  import rxjs or pipe a state (`ui/src/_test/views-without-rxjs.test.ts`): `derive` for computed values, `watch` for a
  component's own properties, `For`/`Show` for structure.
- **Narrowest source** (David, 2026-10-17; [view-models.md §3](ui/_spec/view-models.md)): a view listens only to what
  determines its output: one key of a shared value through `keyedSelector` or `keyedState`, one field through `derive`;
  a side effect that writes no DOM is an `effect()`. A change that writes nothing is waste, and is budgeted
  (`viewer/src/ui/_test/waste.test.tsx`, `waste-budget.json`).
- **Input in the frame lane** (David, 2026-10-17; [rendering.md §2.3](ui/_spec/rendering.md)): input-fed state (a scroll
  position, the hover label, a drag) is a `frameState`: `offer` at every event, published once a frame; a watch
  following one renders in the frame. Batch at the source, never per listener.
- **Subscriptions scale with what is shown** (David, 2026-10-17; [view-models.md §3](ui/_spec/view-models.md)): a
  collection exposes one change stream and keyed lookups (the scene's `changes`); a view of it listens once and updates
  only what changed; view state lives in the view, by key. Guarded per view (`subscriptions.test.tsx`, its budget
  shrink-only) and app-wide (`subscriptions-session.test.tsx`).
- **Work is wanted, not done per change** (David, 2026-10-17; [rendering.md §2.7](ui/_spec/rendering.md)): a
  recomputation many changes ask for (a tree's flattening) is a `workQueue` key, run once a burst in its lane; an answer
  that may come late (a search, a count, a fetch) goes through `latest`, so an old one never lands after a newer.
- **Scoped attachments** (David, 2026-10-17; [rendering.md §2.4](ui/_spec/rendering.md)): a component attaches
  listeners, observers, frames and timers only through `dom/scoped.ts` (`on`, `observeResize`, `observeIntersection`,
  `frame`, `frames`, `timeout`, `drag`): released with its scope, a listener detached while it is paused. A raw one in a
  component fails (`scoped-attach.test.ts`), as rxjs in a view does (`views-without-rxjs.test.ts`); the leak guard
  (`leaks.heap.tsx`) mounts and disposes views and counts what stays.
- **Perf conformance from day one**, as in Seq\* (David, 2026-10-17): base primitives and hot paths carry quick
  allocation guards (`test:perf`: bytes a step against `perf-baseline.json`, each with a positive control) and heap
  guards (`test:heap`). Opt in, by hand at checkpoints, not in CI. A new primitive or hot view lands with its guard.
- **Idle is idle; one frame tick** (David, 2026-10-17; [rendering.md §2.3](ui/_spec/rendering.md)): no Prism code keeps
  a loop running while nothing changes (a loaded viewer at rest asks for no frame and queues no microtask), and frame
  work joins the page's one frame tick (`requestFrameCallback`, or `frame()` in a component), never its own
  `requestAnimationFrame`: one frame request a browser frame, whatever the input. Guarded by
  `viewer/src/ui/_test/frames.quick.ts` (`test:perf`). The empty canvas's splash obeys it too: with no input, no glass
  animation and after its intro (12 s), its particles and waves fade out and it asks for no frame until the pointer
  wakes it (`ambient-pace.ts`).
- **Gestures attach once, and every end path releases** (David, 2026-10-17; [rendering.md §2.4](ui/_spec/rendering.md)):
  a drag, a slider, a splitter, a scrubber is `drag()` (`dom/scoped.ts`): its listeners attached once with the scope,
  the pointer captured, ended by a release, a cancel, the capture lost, a move with no button, the scope paused or
  ended. Never attach on the press and detach on the release (a lost release leaks them). Each gesture has a test that
  loses its release N times (`testing/gestures.ts`, `lostReleases`). No `MutationObserver` over the UI
  (`no-mutation-observers.test.ts`).
- **Bulk operations are one change** (David, 2026-10-17): an operation over many items (expand all, a select or a preset
  over many structures, an update of many cells, a load of many files) is one commit, one notification, one pass; never
  a copy of the whole collection, a scan of it (`find`, `includes`, `indexOf`) or a rebuild of an index per item: index
  once (a Map, a Set), batch (`dataTransaction`, one state build, a task's queue). Every bulk path has a scaling guard
  (`test:scaling`), time and allocation at n and 4n about linear.
- Styling: CSS modules compiled to `_generated/` by `styles`; every color and spacing from `--msp-*` tokens (the scale
  is `--msp-space-1..6`; `styles:check` rejects other values). The look is flat and quiet: few borders and shadows, the
  accent only for actions and the picked item, text that does not fit fades rather than ends in an ellipsis. Check light
  and dark, wide and narrow, in Storybook.
- **Snappy is a requirement** (David, 2026-10-14; the rules in [styling.md §5](ui/_spec/styling.md)): no layout reads
  after writes in the same frame, pointer-rate paths that only update what changed (once per frame, nothing rebuilt or
  allocated per event), animation of `transform` and `opacity` only, flat single-class selectors, contained and
  virtualized subtrees. Measure an interaction in the browser's performance panel before calling it done.
- **Hidden UI is inert and does no UI updates** (David, 2026-10-14; [styling.md §5](ui/_spec/styling.md)): what the
  shell hides but keeps mounted (a folded section, a panel behind an overlay card) is `inert` (no hover, tooltip, focus
  or screen reader) and its scope is paused, so its bindings hold no subscriptions until it shows again (`activity()`,
  `Scope.setActive`; the runtime does it once, views do nothing). A view that subscribes, observes or loops by itself
  pauses that with `onActive(active => …)`; content that is hidden by removing it needs none of this.

## Commands

`pnpm -C storybook storybook` or `mise run storybook` (port 6006; the top-level `storybook/` shows the stories of the
packages registered in `storybook/packages.ts`, every Prism package among them; stories go in `_stories/` folders),
`pnpm -C packages/prism/<pkg> styles` after editing a CSS module (`styles:check` verifies the generated files and the
spacing scale), `typecheck:storybook` for stories and tests.

### Tests

The repository's lanes hold Prism's unit tests (quick, medium, heavy: the root
[AGENTS.md](../../AGENTS.md#tests-three-lanes), Tests); the perf, heap and scaling tiers below are Prism's own, opt-in.

| command                                                       | what it catches                                                                                                                                                                                                                                                                                                                           | when                                                                       |
| ------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------- |
| `pnpm test:unit packages/prism/<pkg>`                         | the quick lane, jsdom, deterministic; with it the work-conformance checks: no rxjs in views (`views-without-rxjs.test.ts`), no raw attachments (`scoped-attach.test.ts`), the waste budget (`waste.test.tsx`), hidden UI writes nothing (`hidden-quiet.test.tsx`), a pick renders the rows it moves between (`scene-tree-keyed.test.tsx`) | every edit loop; CI                                                        |
| `pnpm test:medium packages/prism/<pkg>`                       | the medium lane: the axe sweeps past their representatives (`accessibility.test.tsx`), the hidden-UI setups past the first, the app's subscriptions at 24 and 48 structures (`subscriptions-session.test.tsx`)                                                                                                                            | when you change what they cover; CI                                        |
| `pnpm -C packages/prism/<ui\|viewer> test:perf`               | the quick allocation tier: bytes a step of each `*.quick.ts` case, read exactly (`v8.GCProfiler`), against `src/_test/perf-baseline.json` (`--update` for an intended move); the base primitives read 0                                                                                                                                   | when you touch a primitive or a hot path, in the edit loop; at checkpoints |
| `pnpm -C packages/prism/<ui\|viewer> test:heap`               | the heap guards (`*.heap.ts`, `*.heap.tsx`): growth over n and 4n steps (`heapScaling`; `retainedScaling` after collections for steps that allocate by design), each with a positive control; the leak guard (views mounted and disposed: heap and listeners)                                                                             | at checkpoints and before a release                                        |
| `pnpm -C packages/prism/<ui\|viewer\|assistant> test:scaling` | the bulk paths' scaling guards (`*.scaling.ts`, `*.scaling.tsx`; the harness `sequences/core/_test/model/scaling.ts`, shared with `packages/ai`): each at n and 4n, the fastest time and the allocation (`GCProfiler`) under 8× (a quadratic reads about 16×); a deliberately quadratic control caught                                    | when you touch a bulk path; at checkpoints                                 |

**Run what you touched.** In the edit loop, run the folders you changed
(`pnpm test:unit packages/prism/viewer/src/features/particles`, a test file), not the whole Prism suite, and the guards
only where their row says. The whole suite is for checkpoints. Its timing-sensitive checks (the waste audit, a drag's
lost releases) can fail when every file runs in parallel on a loaded machine: rerun the failing file alone, and only a
failure there is yours. Never raise a budget to quiet a run that passes alone.

The guards run on Seq\*'s harness (`packages/sequences/vitest.shared.ts`, `core/_test/model`; prism-ui's
`src/_test/perf.ts`); a `*.heap.ts` has its quick twin `*.quick.ts`, the two sharing their cases (`*-cases.ts`). CI runs
the quick and medium lanes, `pnpm test` the quick one; `test:perf`, `test:heap` and `test:scaling` are opt-in, never in
CI (as in Seq\*, packages/sequences/AGENTS.md).
