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
- **A virtual list scrolls without work or blanks**: the browser scrolls before the scroll event reaches the list, so
  what it uncovers must already be there. The window stays while what shows (and a few rows) is in it, so most scroll
  events make nothing; past that it moves ahead of the scroll by a screenful; a jump (a scrollbar drag, a far wheel
  step) makes about a screenful and the window grows a screenful each way once the scrolling rests (making rows is what
  costs: a drag makes a window a frame). Its scroll handler reads the scroll and the viewport's height before it sets
  anything (setting them makes rows; a read after would lay the new rows out in the event: 1.7 of a jump's 4.4 ms). A
  row tells whether it has a menu without building it (`hasMenu`). Its area does not fade at its edges (`ScrollArea`
  `fade={false}`): the fade is a mask a scroll-driven animation repaints every frame.

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
