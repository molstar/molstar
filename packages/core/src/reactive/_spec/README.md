# Mol\* reactive: the specification

`@molstar/core/reactive` (`packages/core/src/reactive/`): streams in rxjs 7's shape and states, on one listener base,
built for few allocations. This spec is the rules the code relies on and why; code comments cite its sections
(`reactive §4`). The code is the truth for how.

| File              | What                                                                                  |
| ----------------- | ------------------------------------------------------------------------------------- |
| `types.ts`        | `Observer`, `Subscribable`, `OperatorFunction`, `ObservableInput`, `ReadonlyState`, … |
| `subscription.ts` | `Subscription`                                                                        |
| `observable.ts`   | `Observable`, `Subscriber`, `lift`, `from`, `isObservable`, the rxjs interop          |
| `listened.ts`     | `Listened`, the base of everything with listeners (§2)                                |
| `subject.ts`      | `Subject`, `BehaviorSubject`                                                          |
| `creation.ts`     | `of`, `EMPTY`, `NEVER`, `timer`, `merge`, `combineLatest`                             |
| `operators/`      | one file per operator                                                                 |
| `promise.ts`      | `firstValueFrom`, `lastValueFrom`                                                     |
| `state.ts`        | `State`, `RecordState`, `derive`, `batch`, `toState`, `isState`, `listen` (§4)        |
| `model.ts`        | `Scope`, `ReactiveModel` (§7)                                                         |
| `keyed.ts`        | `KeyedState`, `KeyedSelector` (§8)                                                    |
| `config.ts`       | `ReactiveConfig.onUnhandledError`, `reportUnhandledError` (§6)                        |

## 1. Two kinds

|                           | stream (`Observable`, `Subject`, `BehaviorSubject`) | state (`State`, `derive`)                   |
| ------------------------- | --------------------------------------------------- | ------------------------------------------- |
| current value             | `BehaviorSubject` only                              | always                                      |
| an equal value            | told again                                          | not told (`equals`, `Object.is` by default) |
| a re-entrant `next`/`set` | told at once, nested                                | queued, told after, in order                |
| ends                      | `complete`, `error`                                 | never                                       |
| derived                   | operators: cold, one chain per subscriber           | `derive`: shared, lazy, glitch-free         |

- **`ReadonlyState<T>`** = `Observable<T>` + `value` + `listen(fn) → unlisten`. A `State`, a derived value, a state's
  `asReadonly` view and a `BehaviorSubject` are each one; `derive` takes any.
- **Every state is an `Observable`**: `subscribe` replays the value, then each change.
- Why two kinds: streams keep rxjs's semantics, so code written against rxjs behaves the same; states give comparison,
  batching and ordered delivery, which derived values need to be consistent.

## 2. `Listened`

- An abstract `Observable` with a list of listeners: `listen`, `subscribe`, `observed`, protected `onFirst`/`onLast`
  (start and stop what feeds it), `emit`/`emitError`/`emitComplete` (at once, nested), `removeAll`.
- **Each listener is a `Subscriber` in a slot of a plain array**; the subscriber is the slot and the subscription at
  once (`_owner`, `_index`). Why: `subscribe(fn)` allocates one object, and removal needs no search.
- **Removal empties the slot**; the holes are dropped when they are more than half the list and no delivery runs
  (`_depth`). Why: O(1) amortized removal when many listeners go together, and a delivery's loop stays valid while
  listeners come and go.
- **During a delivery**, a listener added is not told the value being delivered; one removed and not yet told is not
  told it.
- The array is made on first use. A base class, not a list object beside it: one object per source.

## 3. Streams

### 3.1 rxjs shape

- Same names and shapes as rxjs 7 for what Mol\* uses: `Observable` (`subscribe(observerOrNext?)`, `pipe` typed to 9
  operators), `Subject` (`next`, `error`, `complete`, `observed`, `closed`, `isStopped`, `asObservable()`,
  `unsubscribe`), `BehaviorSubject` (`value`, `getValue`), `Subscription` (`closed`, `add`, `remove`, `unsubscribe`,
  `EMPTY`), `Subscriber`. rxjs 7's deprecated three-callback `subscribe` is not kept.
- **ours into rxjs**: `Observable` has `[Symbol.observable]()` (by type: the global `SymbolConstructor.observable`,
  declared as rxjs does; at run time under `Symbol.observable` when defined and always under `'@@observable'`), so
  rxjs's `from`, `combineLatest`, `switchMap` take ours, by type and at run time. A `Subject` is an rxjs `Observer`; a
  `Subscription` is an rxjs `SubscriptionLike`.
- **rxjs into ours**: `from` and the operators with inner sources take `ObservableInput` (ours, any `Subscribable`, a
  promise, an array, an iterable); an rxjs `Observable` is a `Subscribable`. `isObservable` is true for anything with
  `subscribe` and `pipe`.
- **Not compatible**: ours is not assignable to rxjs's `Observable` class (`source`, `operator`, `lift`, `forEach`,
  `toPromise`), and operators do not cross between the two `pipe`s: a mixed chain converts once with `from`.

### 3.2 Kept semantics

- `Observable` is cold: each `subscribe` runs the producer. A producer that throws is the stream's error.
- A `Subscriber` stops at the first `error` or `complete`, then tears down (its upstream with it).
- A `Subject` delivers at once, nested (§2). `error`/`complete` end it: later `next` is ignored, a later subscriber gets
  the end at once. `BehaviorSubject` replays its value nested; `getValue` throws the error after `error`; one that ended
  replays nothing.
- `combineLatest` emits a new array each time (a subscriber may keep it); `combineLatest([])` and `merge()` complete at
  once.
- Time operators (`debounceTime`, `auditTime`, `throttleTime` with `{ leading, trailing }`) time exactly as rxjs's, with
  `setTimeout` (no schedulers: fake timers drive them). `debounceTime` tells its pending value on complete; `auditTime`
  and a trailing `throttleTime` complete after the window's value.
- `switchMap` ends the previous inner before subscribing the next, and completes when the source and the inner have.
  `share()` is ref-counted and resets on end; there is no `shareReplay` (a shared, replayed value is `derive`, §4).
- Each case has a test with rxjs as the oracle (the same scenario through both, compared; `_test/streams.test.ts`).

## 4. States

- **`State`**: `set` compares with the latest value; `value` is always the latest set (inside a batch, while an earlier
  set is still being told). `asReadonly` is a view without `set`, made on first use.
- **`derive(source | sources, fn, equals?)`**: explicit sources (no dependency tracking). Why: a read costs a field
  load, not a registration; the graph is fixed when the derived value is made; the sources are visible where it is made.
  A conditional dependency derives from all its candidates.
- **Push dirty, pull values.** `set` stores the value, bumps the state's `_version` and a global version, queues the
  state if listened, and marks what derives from it (`_mark`: dirty, queued if listened, then its dependents; once per
  global version). Outside a batch it then flushes: each queued node in order is read and, if its value differs from
  what its listeners were last told (`_told`, by its `equals`), told. A derived value recomputes only when read, and
  only if a source's stamp (a node's version, a foreign source's value) moved; its own version moves only when the
  result differs. Why: a diamond computes and tells once with consistent values; a derived value nobody reads is never
  computed; a chain stops at the first unchanged value.
- **Hot and cold.** A node is hot while it has listeners or a hot derived value follows it (`_refs`); a hot derived
  value follows its sources (graph nodes by `_dependents`, a foreign source by a subscriber in its list, no replay).
  Cold, `value` computes from the sources and caches until the global version moves (or, with a foreign source, its
  value). Why: nothing is held or followed for a value nobody listens to.
- **Ordering.** A set during a flush joins it (queued after, so listeners get values in order and end on the latest).
  Within one flush, the sets to one node are told as its latest value. A subscriber set to a state while it gets the
  current value gets the new one after it (the replay runs in a batch).
- **`batch(fn)`** holds the flush to the outermost batch's end (also when `fn` throws); a node set back to what its
  listeners were told tells nothing.
- **A set back keeps the version.** A state remembers its value and version at its first set of a held change (a batch,
  or a flush running; the sets between two flushes are one change, an epoch); set back to that value, its version is
  restored. Versions come from one counter, so a version names one value of one state, never reused. Why: what compares
  versions (a derived value's stamps, `@molstar/ui/dom`'s bindings and watches) sees no change, and nothing recomputes,
  writes or renders for it.
- **A flush that does not settle** (a listener that sets its source forever) stops after 1,000,000 deliveries and
  reports.
- **`RecordState`**: a `State` of a record, `sameShallow` by default (a set of an equal copy tells nothing), with
  `update(partial)` (sets when any given field differs, building on the latest value).
- **Named sources**: `derive({ a, b }, ({ a, b }) => …)`; `fn` gets a new record of the values each compute (it may keep
  it).
- **`toState(stream, initial)`** is a state fed by one subscription, with `dispose`; a `BehaviorSubject` needs none.
- `_test/state.test.ts` has a property test: random graphs and set sequences, every listener's last value equal to a
  recompute from scratch.

### 4.1 Sinks

- The third kind of graph node (`Sink`, internal): sources and no value, for what renders or writes from states
  (`@molstar/ui/dom`'s `watch` and bindings). It is a dependent of each node it follows (`follow`, `unfollow`: a slot in
  the node's array, no object of its own; the node is hot while followed), marked with them, and run in the flush
  (`_flush`) after the change has spread, once per flush however many of its sources changed. Why: a reader of N states
  costs one object, its marking allocates nothing, and what it reads is consistent.
- `graphNode(source)` is the node behind a source (a state's view unwrapped), undefined for a foreign one; `versionOf`
  brings a node up to date and returns its version (a sink compares it to skip an unchanged value).
- **`observed`** on a node is true while anything follows it: a listener, a subscriber, or a dependent (a derived value,
  a sink). The flush tells a node's listeners only when it has some (`_hasListeners`).

## 5. Cost

- `subject.subscribe(fn)`: one object (§2). `listen(fn)`: that object and the closure it returns.
- `next`/`emit`: no allocation; one loop over the array.
- `state.set` with listeners or sinks: no allocation once the graph is hot. The flush queue is one reused array emptied
  by count (an array's `length = 0` frees its store, and the next push allocated it again: 150 bytes a set, found by
  `@molstar/ui`'s `test:perf`); edges exist.
- A pipe stage: one `Observable` at `pipe` time, one subscriber per subscription; an operator's state is its
  subscriber's fields, no closures per subscription. `pipe` takes a rest array: built at setup, not per value.
- `Subscription` finalizers and parents: `undefined`, one, or an array from two.

## 6. Errors

- What a listener or a subscriber's `next` throws, an `error` with no handler, and a derived value that throws while it
  is told go to `ReactiveConfig.onUnhandledError`, called at once; unset, the error is rethrown in a `setTimeout`. The
  producer and the other listeners go on. Why at once when set: tests and hosts see it where it happened.
- A derived value that throws on read throws to the reader and recomputes on the next read.
- An operator's callback that throws is that stream's error.

## 7. Model

- **`Scope`**: a `Subscription` whose teardowns run last added first (later resources may depend on earlier ones); every
  one is attempted and the first failure rethrown. `own(resource)` takes a function, any `Unsubscribable` or a
  `Disposable`. A `Subscription` added leaves the scope when it ends first, so a scope does not grow with subscriptions
  that end on their own.
- **`ReactiveModel`**: `mount` (its owner shows it), `unmount` (it may mount again), `retire` (for good). What `onMount`
  registers (`registerDispose`, `own`, `subscribe`, `listen`) ends with the mount, in reverse; what is registered
  permanently (`registerPermanentDispose`, `ownPermanently`, `subscribePermanently`, `listenPermanently`) ends with
  `retire`. The mount-scoped calls throw outside `onMount`; the permanent ones throw after `retire`.
- **`subscribe`** gives a state's or a `BehaviorSubject`'s current value first; **`listen`** gives only changes. On a
  `Listened` source each costs one object (the subscriber is the listener and the subscription).
- **Rules**: a failed `onMount` ends what it registered (best effort), aborts `mountSignal`, and leaves the model
  unmounted (a later mount may succeed); an `unmount` (or a `retire`) during `onMount` is carried out when it returns;
  `onUnmount` failing still ends the mount's resources, and `retire` still runs `onRetire` and the permanent ones (the
  first failure is rethrown); a second `mount` while mounted warns once and does nothing; after `retire` nothing mounts.
- **`mountSignal`**: an `AbortSignal` made on first read, aborted when the mount ends; outside a mount, an aborted one.
- **The lifecycle methods are the base class's**: a subclass that overrides them throws on construction (checked once
  per class, a `WeakSet`), and they are non-writable on the prototype. Why once per class: a model per row or panel
  costs its fields, not per-instance guards. State is one bit field and two scopes (the permanent one made on first
  use).

## 8. Keyed states

- **`KeyedState<K, V>`**: values by key (`get`, `has`, `values`, `set`, `delete`, `replace`); `key(k)` is a state of
  `k`'s value. A set tells only that key's followers; `replace` (takes the map as is) tells the followed keys whose
  value changed, came or went: its cost is the keys followed, not the map.
- **`KeyedSelector<K>`** over a source of the current key (a state, a `BehaviorSubject`, a stream): `key(k)` (or
  `is(k)`) is a state of whether `k` is current. A change of the current value marks the two keys it concerns, however
  many are followed. It follows the source only while a key is followed. With a source that has a `value`, a key reads
  it, so a key and anything else derived from the source agree within one change.
- **A key is a node of the state graph** (`derive` follows it as a node), registered while followed (one node, or an
  array when several `key(k)` of one key are followed) and dropped with its last follower. Unfollowed, it reads the
  value again on each read. Every keyed change moves the global version (`signalChange`), so a cold value derived from a
  key reads it again.
- Keyed changes batch with the states (`batch`).
