# API surfaces: sync / async / streaming {#api-surfaces}

- **Implementations:** [TypeScript](typescript.md)

Rule §api-boundaries: an API makes its surface visible in the type: a plain value, a promise, or an observable. Return an `Observable` rather than accepting a callback; never expose a transport (WebSocket, a provider's reactive query) directly; expose state as a read-only behavior.

## Problem

An API that takes a callback hides who owns the subscription and how to stop it; one that returns a socket or a
provider's reactive query ties every caller to that transport and makes it untestable without it. State delivered as a
plain event stream has no current value, so a late subscriber sees nothing until the next change. State kept as one
reactive cell per field emits several times for one change, and observers see half-applied states between them.

## Shape

Every operation uses one of three surfaces (§api-boundaries): a plain return value (cheap, pure, in-memory), an async
result (one request, one result), or an `Observable` (many values over time), returned rather than taken as a callback.
An example service: `describe(kind) -> schema` (sync), `submit(spec) -> async id`, `cancel(id) -> async`,
`progress(id) -> Observable<progress>`.

The minimal streaming contract the other patterns assume (any streaming library that fits it will do):

| Name | Shape |
|---|---|
| `Observable<T>` | `subscribe(observer) -> unsubscribe`; an observer is `next(value)`, optional `error(err)`, optional `complete()`, or just a `next` function |
| `observable(setup)` | builds one: `setup(emit, fail, done)` runs per subscription and may return a teardown |
| `createSubject()` | a multicast `Observable` with `next(value)`; no memory (events) |
| `ReadonlyBehavior<T>` | an `Observable` with a synchronously readable current `value` (state) |
| `createBehaviorSubject(initial, equals)` | a writable behavior; `equals` defaults to identity |
| `changes(behavior)` | the same behavior without the delivery of the current value on subscribe |
| `createState(initial, eq?)` | a behavior over one aggregate record, written only through `update(patch) -> changed?`; `eq` gives per-field comparison |
| `select(behavior, f, equals)` | a read-only behavior of a slice, emitting only distinct slices |
| `of(...values)`, `map`, `filter` | the small combinators the other patterns use |

**Events vs state.** A subject carries events (things that happen); a behavior carries state (things that are). Models
expose `ReadonlyBehavior` and keep the writable state private (§events-vs-state). Model state is one aggregate in one
behavior (§atomic-model-state); a standalone scalar may use a plain behavior.

**Transports stay in adapters.** A socket or provider query is opened inside an `observable(setup)` in an adapter;
unsubscribe, error and completion all close it. Callers see only the `Observable`, and tests hand them `of(...)` or a subject.

## Invariants

- After unsubscribe, error or complete, nothing more reaches the observer, and the teardown runs exactly once.
- A source that ends synchronously during `setup` still has its teardown run.
- A subject delivers to each subscription separately, even when the same callback subscribes twice.
- A subscriber added during an emission first receives the next emission; one removed during an emission is skipped.
- One throwing subscriber does not stop delivery to the others; the emitter then receives the first error.
- Subjects never error or complete; only sources that can end call `error` or `complete`.
- A behavior's `value` is always defined; a new subscriber receives the current value synchronously, then each change.
- A behavior does not re-emit a value equal (by its `equals`) to the current one.
- No subscriber of a behavior sees an older value after a newer one, including when a subscriber writes during delivery.
- `changes` delivers nothing at subscribe time; it forwards later values, errors and completion.
- `update(patch)` emits once if some field changed (per-field `eq`, default identity) and not at all otherwise, and
  returns which; it creates a new aggregate and never mutates a snapshot a reader holds.
- `update` costs O(size of patch) and allocates nothing when nothing changed.
- Comparison is shallow, per field: nested values compare by reference, so writers replace only what changed or
  give nested state its own model.
- `select` computes `f` on each read of `value`, so `f` must be cheap and not allocate; it forwards error and completion.
- `map` and `filter` forward error and completion.
- No API signature names a transport type or takes a callback for a stream.

## Conformance cases

- Behavior with initial `1`; subscribe → receives `1` before `subscribe` returns.
- Value `1`; `next(1)` → no emission; `next(2)` → subscribers receive `2`, `value` is `2`.
- Unsubscribe, then `next(3)` → nothing delivered to that observer.
- Two subscribers; the first calls `next(3)` when it receives `2`; then `next(2)` → no subscriber receives `2` after `3`.
- `changes(b)`; subscribe → nothing delivered; `b.next(5)` → `5`.
- `createState({a: 1, b: 2})`; `update({a: 1})` → returns false, no emission.
- Same state; `update({a: 3, b: 4})` → returns true, exactly one emission `{a: 3, b: 4}`; the earlier snapshot is still `{a: 1, b: 2}`.
- `select(state, s -> s.a)`; `update({b: 5})` → the slice emits nothing.
- Subject with two subscribers, the first throws on `x` → the second receives `x`; `next(x)` raises the first error.
- A subscriber added while `x` is being delivered → it does not receive `x`; it receives the next value.
- `observable` whose setup calls `done` synchronously → the observer completes and the teardown runs.
- `of(1, 2)` → `1`, `2`, complete; `map(of(1, 2), x -> 2x)` → `2`, `4`, complete; `filter(of(1, 2), odd)` → `1`, complete.
- `observable` that fails with `e`, then emits → the observer receives `error(e)` once and nothing after.

## Open

- §conformance lists an `observableSuite` for this pattern, but no suite is written here; the cases above come from the
  code's comments and notes, not from an existing suite.
