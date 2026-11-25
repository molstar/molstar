# API surfaces: sync / async / streaming: TypeScript

The concept, invariants and conformance cases are in [README.md](README.md); this file is the TypeScript implementation.

## Minimal Observable contract

Any library works (rxjs, Effect). This is the smallest shape the other patterns assume.

```ts
export type Unsubscribe = () => void;

export interface Observer<T> {
  next(value: T): void;
  error?(err: unknown): void;
  complete?(): void;
}

export interface Observable<T> {
  subscribe(observer: Observer<T> | ((v: T) => void)): Unsubscribe;
}

const toObserver = <T>(o: Observer<T> | ((v: T) => void)): Observer<T> => (typeof o === "function" ? { next: o } : o);

export function observable<T>(
  setup: (emit: (v: T) => void, fail: (e: unknown) => void, done: () => void) => Unsubscribe | void,
): Observable<T> {
  return {
    subscribe(o) {
      const obs = toObserver(o);
      let closed = false;
      let teardown: Unsubscribe | void;
      const end = () => { if (!closed) { closed = true; teardown?.(); } };   // error/complete also release resources
      teardown = setup(
        (v) => { if (!closed) obs.next(v); },
        (e) => { if (!closed) { end(); obs.error?.(e); } },
        () => { if (!closed) { end(); obs.complete?.(); } },
      );
      if (closed) teardown?.();                                              // ended synchronously during setup
      return end;
    },
  };
}

export interface Subject<T> extends Observable<T> { next(v: T): void }

export function createSubject<T>(): Subject<T> {
  const subs = new Set<{ next: (v: T) => void }>();                        // one entry per subscription, even for the same callback
  return {
    next(v) {                                                               // snapshot: subscribers added during emission start next time
      let failed = false, first: unknown;
      for (const s of [...subs]) {
        if (!subs.has(s)) continue;                                         // unsubscribed during this emission: skip
        try { s.next(v); } catch (e) { if (!failed) { failed = true; first = e; } }   // one throwing subscriber does not starve the rest
      }
      if (failed) throw first;
    },
    subscribe(o) {
      const entry = { next: toObserver(o).next.bind(o) as (v: T) => void };
      subs.add(entry);
      return () => { subs.delete(entry); };
    },
  };
}

export const of = <T>(...values: T[]): Observable<T> =>
  observable((emit, _fail, done) => { values.forEach(emit); done(); });
```

## Behavior subject: observe the current value

A `Subject` has no memory: late subscribers miss earlier values. State needs the opposite: always a current value, readable synchronously, and delivered to new subscribers. That is a **behavior**.

```ts
/** Read-only view: what models expose. Consumers read `value` and subscribe, but cannot write. */
export interface ReadonlyBehavior<T> extends Observable<T> { readonly value: T }

export interface BehaviorSubject<T> extends Subject<T>, ReadonlyBehavior<T> {}

export function createBehaviorSubject<T>(initial: T, equals: (a: T, b: T) => boolean = Object.is): BehaviorSubject<T> {
  let value = initial;
  const subject = createSubject<T>();
  return {
    get value() { return value; },
    next(v) { if (!equals(v, value)) { value = v; subject.next(v); } },     // distinct values only
    subscribe(o) {
      const f = toObserver(o).next.bind(o) as (v: T) => void;
      // register first, so a re-entrant next() is not missed; skip values a re-entrant next() already superseded,
      // so no subscriber sees an older value after a newer one
      const off = subject.subscribe((x) => { if (Object.is(x, value)) f(x); });
      f(value);                                                             // deliver the current value synchronously
      return off;
    },
  };
}

/** Future changes only: the same behavior without the synchronous delivery of the current value.
 *  For consumers that read `.value` themselves (view bindings, §reactive-models) and would otherwise redo work for a value they already have. */
export function changes<T>(src: ReadonlyBehavior<T>): Observable<T> {
  return observable((emit, fail, done) => {
    let initial = true;                                                     // the first delivery happens inside subscribe()
    const off = src.subscribe({ next: (v) => { if (!initial) emit(v); }, error: fail, complete: done });
    initial = false;
    return off;
  });
}
```

### Aggregated state with `update(partial)`

Model state is one aggregate in one behavior (§atomic-model-state), written with `update(partial)`; semantics: see [README.md](README.md#invariants).

```ts
export interface StateSubject<T extends object> extends BehaviorSubject<T> {
  /** Merge `patch`; emits once, only if some field actually changed. Returns whether it did. */
  update(patch: Partial<T>): boolean;
}

type FieldEq<T> = { [K in keyof T]?: (a: T[K], b: T[K]) => boolean };

export function createState<T extends object>(initial: T, eq?: FieldEq<T>): StateSubject<T> {
  const b = createBehaviorSubject<T>(initial);                             // identity: a new object means a change
  return Object.assign(b, {
    update(patch: Partial<T>): boolean {
      const cur = b.value;
      let next: T | undefined;
      for (const k in patch) {                                             // O(|patch|), not O(|state|); no allocation if unchanged
        if (!Object.hasOwn(patch, k)) continue;
        const v = patch[k] as T[typeof k];
        if (!(eq?.[k] ?? Object.is)(cur[k], v)) { next ??= { ...cur }; next[k] = v; }
      }
      if (!next) return false;
      b.next(next);                                                        // one new object, one emission
      return true;
    },
  });
}

/** Distinct slice of a behavior: consumers subscribe to what they use, not the whole state. */
export function select<T, S>(src: ReadonlyBehavior<T>, f: (t: T) => S, equals: (a: S, b: S) => boolean = Object.is): ReadonlyBehavior<S> {
  return {
    get value() { return f(src.value); },                                  // recomputes on read: keep f cheap and non-allocating
    subscribe(o) {
      const obs = toObserver(o);
      let first = true, last!: S;
      return src.subscribe({                                               // forward error and complete too, like map/filter
        next: (t) => { const s = f(t); if (first || !equals(last, s)) { first = false; last = s; obs.next(s); } },
        error: (e) => obs.error?.(e),
        complete: () => obs.complete?.(),
      });
    },
  };
}
```

Invariants (shallow per-key diffs, `update` as the only write path, events vs state): see [README.md](README.md#invariants). TypeScript notes:
- Supply a per-key `eq` where reference equality is wrong (small value objects).
- Stable snapshot references are also what `useSyncExternalStore` needs.
- `Unsubscribe | void` (and §representations' `boolean | void` visitor result) trips lint rules against `void` in unions; keep the type
  (a setup that returns nothing must typecheck) with a local suppression that says why.

## One API, three surfaces

```ts
export interface JobService {
  describe(kind: JobKind): JobSchema;                 // sync: pure, in-memory, no I/O
  submit(spec: JobSpec): Promise<JobId>;              // async: one request, one result
  cancel(id: JobId): Promise<void>;
  progress(id: JobId): Observable<JobProgress>;       // streaming: values over time, returned (not a callback parameter)
}

// BAD: a callback parameter, and the transport leaks out
// onProgress(id: JobId, cb: (p: JobProgress) => void): WebSocket;
```

## Wrapping a transport inside an adapter

```ts
function progressFromSocket(url: string, id: JobId): Observable<JobProgress> {
  return observable((emit, fail, done) => {
    const ws = new WebSocket(`${url}/jobs/${id}`);
    ws.onmessage = (m) => emit(JSON.parse(m.data) as JobProgress);
    ws.onerror = (e) => fail(e);
    ws.onclose = () => done();
    return () => ws.close();                          // unsubscribe, error and completion all close the socket
  });
}
```

The socket exists only in this function. Callers see `Observable<JobProgress>` and tests hand them `of(...)` or a `Subject`.

## Small combinators (enough for the other patterns)

```ts
export function map<A, B>(src: Observable<A>, f: (a: A) => B): Observable<B> {
  return observable((emit, fail, done) => src.subscribe({ next: (a) => emit(f(a)), error: fail, complete: done }));
}

export function filter<A>(src: Observable<A>, p: (a: A) => boolean): Observable<A> {
  return observable((emit, fail, done) => src.subscribe({ next: (a) => { if (p(a)) emit(a); }, error: fail, complete: done }));
}
```
