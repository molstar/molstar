# Reactive models; the view is a thin layer {#reactive-models}

- **Implementations:** [TypeScript](typescript.md)

Rules §frontend-state: state and logic live in nested reactive models; UI actions are methods on models, backed by commands (§command-reachable); views subscribe to model behaviors; the view library owns no application state, only transient local UI state; the main thread stays responsive.

## Problem

State and logic kept in components cannot be reached by agents, tests or other views, and leave with the view
library. A model with one reactive cell per field emits several times for one change. An optimistic write in a single
model ("set it, restore it on error") is seen by nobody else, and its restore can overwrite data that arrived
meanwhile. A binding that re-renders on subscribe or on unrelated changes wastes work, and long main-thread work freezes input.

## Shape

- **Model**: one aggregate state in a private `createState` (§api-surfaces), exposed as a `ReadonlyBehavior` `state`
  plus slices from `select` (`done`, `openCount`). State is plain data, so it has a declarative representation (§declarative-representations).
- **Actions** are methods with stable identity (bound once), taking only what the UI naturally provides
  (`toggle()`, `setText(text)`). Each submits a command as a transaction through `Send`.
- **`Send(label, command, input, touches)`**: submits one command as a transaction to the replica (§transactions);
  `touches` lists the entity ids it affects. It never raises for a rejection; the status arrives through `sync`.
  `sendVia(replica, origin)` builds one.
- **`sync(data, txStatus?)`**: the owner calls it with the *effective* data (confirmed ⊕ pending) and the latest
  transaction status for this entity; the model derives `pending` and `error` from it.
- **Nested models**: a list model reads a list port (`watchAll() -> Observable<list>` over the effective store, and the
  replica's transactions as a behavior) and keeps one item model per id (`item(id)`).
- **Lifecycle**: `mount()` and `unmount()` start and stop subscriptions and timers (a repeatable scope the view owns);
  `dispose()` is final and belongs to whoever created the model.
- **View binding**, one small module per view library:
  - reads the current value directly and subscribes to future changes only (`changes`, §api-surfaces);
  - selects a slice and compares it (shallow by default); an optional source with a fallback stops hidden UI updating.
- **Where state lives** (§thin-view-layer): hover, focus, open/closed and a text draft before commit stay in the view;
  anything saved, shared, replayed, driven by a command or inspected by an agent is in a model; a drag or draft that
  previews a document change is a model-level overlay (§transient-overlay).
- **Responsive main thread** (§responsive-ui): a task is a resumable computation with explicit safe points (a
  generator). `runChunked(task, options) -> async result`, with a budget (default 8 ms) and an optional signal and
  onSlice callback, runs it in time-boxed slices and yields to the host between slices; the same scheduler is the main-thread backend of a compute port, so moving work to
  a worker changes the composition, not callers (§location-agnostic-ui).

## Invariants

- One aggregate per model; independent concerns get nested models; no parallel cells duplicate parts of the state.
- The writable state is private; the outside sees only read-only behaviors.
- A model makes no optimistic write and no rollback: pending and failed state come from the transaction status, never
  from flags the model sets itself (§optimistic-from-store).
- A model never reaches around commands to a store for a mutation (§one-path-to-behavior).
- Item models keep their identity across updates (same id, same instance) and are dropped when their id leaves the data.
- An item that is new may already be pending; it is synced with its status when created.
- `mount` while mounted does nothing; `unmount` then `mount` resumes; `dispose` unmounts and drops the item models.
- Model files import no view library; models run headless on the in-memory composition (§headless-mode).
- A binding causes no render or store callback on subscribe, even with an `equals` that always says "different".
- A component re-renders only when its selection changes, and an equal selection keeps its previous reference.
- The source's identity is its lifetime: a new source per render starts over, so sources are memoized, and dev warns
  when a source changes on consecutive renders.
- Hidden UI stays mounted but holds no subscriptions or timers.
- `runChunked`: each slice stops at the first safe point after its budget is used; between slices control returns to
  the host; `onSlice` reports each slice's duration; on abort the task's cleanup runs and the call raises a
  cancellation that `isCancelled` recognizes (§error-codes).
- Status (pending, rejected, offline) is read from the replica and its link and projected by a model, never stored in one.

## Conformance cases

- List with item "a" added and toggled, server not yet answering → the list's `todos[0].done` is true, `openCount` is
  0, the item's `pending` is true.
- Server confirms both → the item's `pending` is false.
- Toggle rejected with a coded error while a remote edit of the text lands → `done` back to false, `text` is the remote
  edit, the item's `error.code` is the server's code.
- `item(id)` before and after a data change for that id → the same instance.
- `mount()` twice → one subscription to the list port.
- Binding with `equals` always "different" → subscribing causes no render.
- `runChunked` over a 200,000-item index build with `budgetMs: 8` → the longest slice is under 16 ms (guard).
- `runChunked` with an already aborted signal → raises a cancellation; the task's cleanup ran.

## Open

- The rejected-toggle case is described but its body is not written in the code.
