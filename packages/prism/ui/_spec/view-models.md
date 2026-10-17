# View models

The state layer views read: renderer-free models on prism-ui's reactive primitives (`src/reactive/`), and the rules an
app's sources and models follow (the viewer's are in its `sources/` and `model/`; its state contracts in the viewer
`_spec/state.md`). The code is the truth for how.

## 2. Primitives and conventions

- **`observableState(initial, { compare })`** (default `sameShallow`): `value`, `set`, `asReadonly` (a record state also
  `update`). A `set` made while subscribers are notified is queued, so every subscriber sees values in order.
  `ReadonlyState<T>` is `Observable<T> & { value }`. Also `mapState`, `combineBehaviors`, `keyedSelector`, `batch`,
  `createDisposableGroup`, the transaction coordinator (§5).
- **`ReactiveModel`**: `mount`/`unmount`/`retire` with `onMount`/`onUnmount`/`onRetire` hooks; mount-scoped and
  permanent subscriptions; one owner. `mountSignal` (protected) is an `AbortSignal` for this mount, made on first read
  and aborted when the mount ends (`unmount`, `retire`, a throwing `onMount`); a remount gets a fresh one; outside a
  mount it is already aborted.
- **`ParamsFormModel`** (`forms/params-form.ts`) is the form for PD params everywhere: `values`, `field(name)` (cached,
  distinct), `visible(name)`, `set`, `load(values)` (authoritative, not reported as a change), `reset`, `layout`
  (essentials, then advanced), `onChange`. Its owner decides what a change does (apply at once, or keep until a button);
  the same instance each time it is asked for, so a form keeps its state while views come and go.
- **Public fields are `ReadonlyState`**; writable states are private, changed through methods.
- **Lists take one of two shapes**: live nodes where rows are many and change in place (a scene node keeps its identity
  per ref, each field its own state, `children` a state of nodes); records where items are small and change as a whole
  (plain objects with an `id` in one `ReadonlyState<readonly T[]>`; a change replaces the record, never mutates it).
  Lists fixed for a model's life are plain arrays.
- **Optional parts are optional members**, present only when the app turns them on.
- **Expensive children are lazy** (`lazy`, `load()`: a structure's chains computed when first expanded).
- **State is a behavior, an event is a non-replaying `Subject`** (`CommandManager.outcomes`); a request a view acts on
  when it shows is state (`ShellModel.tabFocus`).
- **Errors are values** (a record's `status`/`error`), never stream errors.
- **Streams are made once per model** (fields, cached getters), never per render.
- **Dependencies come through the constructor**; upward calls through callbacks.
- **Combine, then debounce**: a view depending on several streams combines them, with a short debounce where bursts come
  (state-tree rebuilds, selection changes), so one change renders once.

## 3. Sources and records

Sources turn an engine's APIs (the plugin's managers and events) into `ReadonlyState`s; models and views read only them.

- **States at the boundary** ([rendering.md](rendering.md) §2.2): a state is its latest value plus listeners, no rxjs on
  its path. rxjs becomes a state once per value, where the engine's events come in (`toState(stream, initial)`, or an
  `observableState` the source sets), never per view binding. Only what consumes Mol\*'s and the plugin's events
  (`plugin.events`, cell updates, the viewer's `sources/`) keeps rxjs inside; models and sources expose states, and what
  they derive is `derive` (lazy and hot: it follows its sources only while something listens). Views receive states and
  derived values, and do no rxjs (rendering.md §2.6).

- **Records per node**: a tree's source hands out one record per node, keyed by id (the state ref), that changes on its
  own, not a tree rebuilt per change: a row follows its record by id, and `For` renders an entry again only when its
  record is replaced ([rendering.md](rendering.md) §2.5). A rebuilt tree costs a subtree per change.
- **Identity by id**: a rebuild that recreates the engine's objects keeps the same record for the same id; its fields
  emit only when what they show changed.
- **One change, one emission**: code that sets several states for one change wraps it in `batch`, so each emits once
  with its final value (views render once per task anyway: `watch`, [rendering.md](rendering.md) §2.3); long lists mark
  the current item with `keyedSelector` (a row watches `key(its key)`), so a pick renders the two rows concerned.
- **A view follows the narrowest source that decides its output** (David, 2026-10-17): one key of a shared value
  (`keyedSelector`, `keyedState`), one field of a record (`derive` the field, or split the state), never the whole for a
  part. A change that ends in no DOM write is waste. The waste audit (`viewer/src/ui/_test/waste.test.tsx`, the real
  shell in jsdom over a person's session: picks, hover, selections, a color, tabs, the command menu, the theme,
  trajectory frames, notifications) counts, per state change, the listener calls, the derived recomputes that came out
  the same and the writes and renders that changed no DOM, by source state and by view (where each was made;
  `RuntimeProbe.run`, `ReactiveProbe`); `waste-budget.json` holds today's numbers per source file, and a source over its
  budget fails. At the start the session made 16 writes and renders that changed nothing (7 the scene tree's
  scroll-into-view, a side effect now in `effect()`; 2 the dock's tab lists, now deduplicated; 7 the trajectory
  scrubber's `value`, a property the observer cannot see); 51 derived recomputes are the gates that keep the rest from
  waking (a navigation or trajectory field derived from its record).
- **A row never listens to a whole shared value it cares about one key of.** The current item of a list (the picked
  node, the selected cell, the highlighted option, the current model) is a `keyedSelector(state)`, made once per list; a
  row follows `key(its key)`, told only when it becomes or stops being current. Values per row (a tree's places and open
  states) are a `keyedState`; the flattening `replace`s the map and only the rows whose value changed are told. Measured
  on the scene tree with 48 rows mounted: a pick went from 97 listener calls and renders (82 on untouched rows) to 6
  (none), a reflow below the window from 101 (91) to 5 (none, no render)
  (`views/structure/_test/scene-tree-keyed.test.tsx`; rendering.md §1). Small fixed option groups (a theme's three
  choices) may read the shared value whole.
- **After the transaction**: sources update when the engine's transaction ends, not on intermediate trees.
- **Long lists stay cheap**: a list that can grow past a screenful of rows is virtualized (`components/virtual-list`):
  only the window is in the DOM. Rows are of one height, measured (compact mode's taller rows need nothing); what would
  make one taller goes elsewhere (a menu to a popover, a row's controls to the inspector), and a tree's extra lines
  (Show all, Loading…) are rows of their own. A tree is its shown nodes flattened (`flattenTree`, or the view's own walk
  when rows open by their own states), kept by key (`keepRows`: a row the same node at the same depth keeps its DOM);
  places (`aria-posinset`, `aria-setsize`) and open states come apart from the rows, so a sibling added remakes none of
  them. **A row is one `watch`** ([rendering.md](rendering.md) §2.3): its DOM made once, its render reading its node's
  states and the tree's (the picked node, the map of places and open states: `states.value.get(key)`), so a row holds
  one listener per source, not a chain per attribute (the scene tree's row: 16.1 → 4.0 KB). Its keyboard works on the
  flattened rows, not the DOM (`dom/tree-keys.ts`, `TreeNavigation`; it flushes before it reads a row). The viewer's
  scene and state trees are built so.
- **The scene tree recycles its node rows** (`RecycledList`, rendering.md §1): a fixed pool of slots re-pointed as the
  window moves, each placed by its own transform; a jump costs 1.7–2.4 ms against 4.7–4.9 with keyed rows, and half the
  allocation and GC. Other lists keep keyed rows (`VirtualList`). A slot shows every part a row may have (shown or
  hidden per item), follows its item's keys (`watch(...).retarget`), and on re-point closes its menu and moves its hover
  and keyboard registration.
- **An opened node comes in** (`reveal`, David, 2026-10-17): a change of a virtual list's rows that only inserts one run
  (`insertedRun`: a node opened, its children loaded) animates; any other change and every scroll are instant. Keyed
  rows grow in (`data-reveal`, the row height not measured meanwhile); recycled rows fade in at their place while the
  rows below glide down over them (`.glide`, `.enter`), the classes off after `--msp-motion-base`, or at once on a
  scroll (a slot re-pointed mid-glide would slide across the list).
- **A virtual list scrolls without work or blanks**: the browser scrolls before the scroll event reaches the list, so
  what it uncovers must already be there. The window stays while what shows (and a few rows) is in it, so most scroll
  events make nothing; past that it moves ahead of the scroll by a screenful; a jump (a scrollbar drag, a far wheel
  step) makes about a screenful and the window grows a screenful each way once the scrolling rests (making rows is what
  costs: a drag makes a window a frame). Its scroll handler reads the scroll and the viewport's height before it sets
  anything (setting them makes rows; a read after would lay the new rows out in the event: 1.7 of a jump's 4.4 ms). A
  row tells whether it has a menu without building it (`hasMenu`). Its area does not fade at its edges (`ScrollArea`
  `fade={false}`): the fade is a mask a scroll-driven animation repaints every frame.
- **Subscriptions scale with what is shown** (David, 2026-10-17): a view holds a constant and some listeners per row on
  screen, never per item of the model. So:
  - **A collection exposes one change stream and keyed lookups**: a view of the whole collection listens once. The
    scene's structure is `SceneModel.changes` (`viewer/src/model/scene-changes.ts`): the nodes whose children or
    expansion changed, by id, told once per Mol* transaction (the scene source's sync) or toggle, with a version; the
    roots are their own state, set only when they differ. A node's fields stay its own states: a row follows its node's
    (one row, its few listeners), so the fields need no keyed table of their own.
  - **Derived views update incrementally**: the scene tree flattens again only each changed node's subtree, in place,
    and keeps the rest (rows the same are the same objects); the state tree's toggle does the same. A full flattening is
    for what changes everything (a search, the roots, the extensions' rows).
  - **View state lives in the view**, by key: the state tree's collapsed cells are the view's set. The scene's expansion
    stays on its nodes, as the model's: the inspector opens a node's path (`reveal`), and the tree learns of it through
    the change stream.
  - Measured (jsdom, the scene tree with 1,000 components open, 4,000 rows): live listeners 8,457 → 454 (2,457 at 250
    components → the same 454); a mount 24–26 → 16 ms and 23.1 → 5.7 MB kept; a toggle's reflow 765 → 60 µs; expanding
    all 1,000 at once 10.3 → 3.7 ms, closing them 7.6 → 2.5. A thumb's drag tells no node state anything (only the
    list's own); the old listeners were idle during it, but held. The state tree held 269 to 271 listeners at 5,025 and
    20,100 cells already; its toggle, 9.0 and 11.5 ms, is 8.4 and 8.5 now (the rest is the rows' rendering).
  - **Guards**: `viewer/src/views/_test/subscriptions.test.tsx` (the trees' listeners the same for a model 4× larger,
    within `subscription-budget.json`, shrink-only; no node state told during a drag; a view listening to every item
    caught) and `viewer/src/ui/_test/subscriptions-session.test.tsx` (the real shell with 24 and 48 structures, each tab
    shown: no view's listeners grow; a planted view listening to every structure caught). The inventory
    (`waste-audit.ts`, `subscriptions()`: live listeners by the state and by the view that added them, through
    `ReactiveProbe.listened`; `SUBS_REPORT=<file>` writes it) found the scene tree's flattening the one view that grew
    with the model. It counts the reactive base's listeners, not rxjs subscriptions to Mol*'s own subjects.

### 4.3 Lifetimes

- **App models live as long as the app.**
- **View-only models** (a panel's list, an inspector form) are created when their view shows and unmounted when it
  closes, so hidden projections cost nothing.
- **Node models in keyed lists** are owned by their list and retired when their key goes.
- **One inspection**: the scene keeps the picked node's inspection only (another's is built again when picked), and
  builds it anew, its form kept, when what it was built from changes (a structure's data, a volume's representations);
  `inspection` is the one the inspector follows.

## 5. Results and state

- **Two channels.** A method that does work returns a promise of its result (never rejects; also posted to status): the
  outcome for whoever asked. State changes only through the reactive state fed by the sources from the engine's events,
  never from a promise's result: one source of truth whoever made the change (the person, the agent, a view).
- **Optimistic only when the future is known**: a change is projected at once when the request fully determines the
  target state and the source will report exactly that (a parameter edit): `createUiTransactionCoordinator` and
  `createOptimisticProjection` (`reactive/transaction.ts`); a contribution is absorbed when the source reports it and
  rolled back when the run fails. Synchronous changes (a visibility flag) and unpredictable ones (adding, loading) are
  not projected; the latter show progress.
- **Writing transform params**: forms emit the complete new params, each value of its PD type, and skip the update when
  nothing changed (`PD.areEqual`); `load` takes the authoritative params back after an update made elsewhere.
