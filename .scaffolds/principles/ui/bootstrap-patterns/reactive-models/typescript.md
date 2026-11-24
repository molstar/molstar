# Reactive models; the view is a thin layer: TypeScript

The concept, invariants and conformance cases are in [README.md](README.md); this file is the TypeScript implementation.

## State is an aggregate, not a cell per field

A model keeps its state as **one aggregate object** in a `StateSubject` (see [§api-surfaces](../../../architecture/bootstrap-patterns/api-surfaces/typescript.md#aggregated-state-with-updatepartial)). Writes go through `update(partial)`; the outside sees a `ReadonlyBehavior`. Consumers subscribe to the whole state or to slices via `select`.

```ts
// BAD if propagated: one cell per field. A change spanning fields emits several times, observers can see
// inconsistent intermediate states, and every new field adds boilerplate and subscriptions.
class TodoModelBad { readonly text = createBehaviorSubject(""); readonly done = createBehaviorSubject(false); }

// GOOD: one aggregate, atomic updates, actions backed by commands submitted as transactions
interface Todo { readonly id: string; readonly text: string; readonly done: boolean }
interface TodoState extends Todo { readonly pending: boolean; readonly error?: ErrorInfo }

/** Submits one command as a transaction (§transactions). Never throws for a rejection: the status arrives through sync(). */
export type Send = (label: string, command: string, input: unknown, touches: readonly string[]) => void;
export const sendVia = (replica: Replica, origin = "ui"): Send => (label, command, input, touches) =>
  void replica.submit({ id: newId("tx"), label, origin, expect: [], steps: [{ command, input }], touches });

export class TodoModel {
  private readonly _state: StateSubject<TodoState>;
  readonly state: ReadonlyBehavior<TodoState>;
  readonly done: ReadonlyBehavior<boolean>;                            // slice for views that need only this

  constructor(todo: Todo, private readonly send: Send) {
    this._state = createState<TodoState>({ ...todo, pending: false }); // seeded from the data, not from defaults
    this.state = this._state;
    this.done = select(this._state, (s) => s.done);
  }

  /** Called by the owning list with the EFFECTIVE data (confirmed ⊕ pending, §transactions) and this item's transaction status. */
  sync(todo: Todo, status?: TxStatus) {
    this._state.update({ text: todo.text, done: todo.done, pending: status?.state === "pending",
                         error: status?.state === "failed" ? status.error : undefined });
  }

  // UI actions are bound closures with stable identity: views pass them directly (onChange={model.toggle})
  setText = (text: string) => this.send("Edit todo", "todo.setText", { id: this._state.value.id, text }, [this._state.value.id]);
  toggle = () => this.send("Toggle todo", "todo.setDone", { id: this._state.value.id, done: !this._state.value.done }, [this._state.value.id]);
}
```

Why there is no optimistic write or rollback in the model, and the other model invariants: see [README.md](README.md#invariants).

## Nested models, lifecycle

```ts
/** A list port's streaming surface over the EFFECTIVE store: the current list first, then each change (§ports-and-adapters, §transactions). */
export interface TodoListPort {
  watchAll(): Observable<readonly Todo[]>;
  txs: ReadonlyBehavior<ReadonlyMap<string, TxInfo>>;                 // the replica's pending and failed transactions (§transactions)
}

interface TodoListState { readonly todos: readonly Todo[]; readonly filter: "all" | "open" }

export class TodoListModel {
  private readonly _state = createState<TodoListState>({ todos: [], filter: "all" });
  readonly state: ReadonlyBehavior<TodoListState> = this._state;
  readonly openCount = select(this._state, (s) => s.todos.reduce((n, t) => n + (t.done ? 0 : 1), 0));
  private readonly items = new Map<string, TodoModel>();               // item models keep their identity across updates
  private off?: Unsubscribe;

  constructor(private readonly list: TodoListPort, private readonly send: Send) {}

  /** Item model by id, reconciled with the data: same id, same model instance. */
  item = (id: string): TodoModel | undefined => this.items.get(id);

  // repeatable mount scope (owned by the view, see below); final disposal belongs to the creator
  mount = () => {
    if (this.off) return;
    let todos: readonly Todo[] = [];
    const refresh = () => {
      const status = statusByEntity(this.list.txs.value);              // latest tx per touched id
      const seen = new Set<string>();
      for (const t of todos) {
        seen.add(t.id);
        let m = this.items.get(t.id);
        if (!m) this.items.set(t.id, (m = new TodoModel(t, this.send)));
        m.sync(t, status.get(t.id)?.status);                           // a new item may already be pending (add)
      }
      for (const id of [...this.items.keys()]) if (!seen.has(id)) this.items.delete(id);
      this._state.update({ todos });
    };
    const a = this.list.watchAll().subscribe((next) => { todos = next; refresh(); });
    const b = changes(this.list.txs).subscribe(refresh);                // §api-surfaces: later changes only; the first refresh ran above
    this.off = () => { a(); b(); };
  };
  unmount = () => { this.off?.(); this.off = undefined; };
  dispose = () => { this.unmount(); this.items.clear(); };

  setFilter = (filter: TodoListState["filter"]) => { this._state.update({ filter }); };
  add = (text: string) => { const id = newId("todo"); this.send("Add todo", "todo.add", { id, text }, [id]); };   // client id (§client-ids)
}
```

No React import appears in model files. Models are testable on the in-memory composition with no DOM.

## React as a view layer only

The binding is one small module; everything else in components is rendering. It must respect how behaviors work: a behavior delivers its current value on subscribe, and the hook already reads that value synchronously.

```tsx
// WASTEFUL: the subscription's immediate delivery fires the store callback on mount for nothing
// (redundant selector work; with useState/useEffect-based bindings it is an extra render)
export function useObservableNaive<T>(src: ReadonlyBehavior<T>): T {
  return useSyncExternalStore((cb) => src.subscribe(cb), () => src.value);
}
```

```tsx
// ui-react/use-observable.ts
type Sel<T, U> = { selector: (v: T) => U; equals: (a: U, b: U) => boolean };

export const sameShallow = <U,>(a: U, b: U): boolean => {
  if (Object.is(a, b)) return true;
  if (typeof a !== "object" || typeof b !== "object" || !a || !b) return false;
  const ka = Object.keys(a), kb = Object.keys(b);
  const ra = a as Record<string, unknown>, rb = b as Record<string, unknown>;
  return ka.length === kb.length && ka.every((k) => Object.hasOwn(rb, k) && Object.is(ra[k], rb[k]));
};

function createStore<T, U>(src: ReadonlyBehavior<T> | undefined, fallback: U, cb: { current: Sel<T, U> }) {
  let snapshot!: T, selection!: U, has = false;
  let selectedWith: ((v: T) => U) | undefined;

  const select = (v: T): boolean => {                                  // returns whether the *selection* changed
    const { selector, equals } = cb.current;
    selectedWith = selector;
    const next = selector(v);
    if (has && equals(selection, next)) return false;                  // an equal result keeps the old reference
    selection = next; has = true; return true;
  };

  return {
    // Read the current value directly; re-select only if the value or the selector changed, so the
    // returned reference is stable between changes (required by useSyncExternalStore).
    getSnapshot(): U {
      if (!src) return fallback;
      const v = src.value;
      if (!has || !Object.is(snapshot, v) || selectedWith !== cb.current.selector) { snapshot = v; select(v); }
      return selection;
    },
    // Future changes only (§api-surfaces's changes()): nothing fires on mount.
    subscribe(onStoreChange: () => void) {
      if (!src) return () => {};
      return changes(src).subscribe((v) => { snapshot = v; if (select(v)) onStoreChange(); });
    },
  };
}

export function useObservableProp<T, U>(src: ReadonlyBehavior<T>, selector: (v: T) => U, equals?: (a: U, b: U) => boolean): U;
export function useObservableProp<T, U>(src: ReadonlyBehavior<T> | undefined, selector: (v: T) => U, equals: ((a: U, b: U) => boolean) | undefined, fallback: U): U;
export function useObservableProp<T, U>(src: ReadonlyBehavior<T> | undefined, selector: (v: T) => U, equals: (a: U, b: U) => boolean = sameShallow, fallback?: U): U {
  const cb = useRef<Sel<T, U>>({ selector, equals });
  cb.current = { selector, equals };                                   // latest-ref pattern: selector must be pure
  const store = useMemo(() => createStore(src, fallback as U, cb), [src]);   // the source's identity is its lifetime
  return useSyncExternalStore(store.subscribe, store.getSnapshot, store.getSnapshot);
}

export function useObservable<T>(src: ReadonlyBehavior<T>, equals?: (a: T, b: T) => boolean): T;
export function useObservable<T>(src: ReadonlyBehavior<T> | undefined, equals: ((a: T, b: T) => boolean) | undefined, fallback: T): T;
export function useObservable<T>(src: ReadonlyBehavior<T> | undefined, equals?: (a: T, b: T) => boolean, fallback?: T): T {
  return useObservableProp(src, identity, equals, fallback as T);
}
const identity = <T,>(v: T) => v;
```

Binding invariants: see [README.md](README.md#invariants). In React:
- **Optional sources**: `useObservable(active ? model.state : undefined, undefined, fallback)` is how hidden UI stops updating.

```tsx
export function TodoRow({ model }: { model: TodoModel }) {
  const { text, done } = useObservable(model.state);                 // one subscription, one consistent snapshot
  return (
    <label>
      <input type="checkbox" checked={done} onChange={model.toggle} />
      {text}
    </label>
  );
}

const done = useObservableProp(model.state, (s) => s.done);          // or just the slice
```

### Where state lives

```tsx
// transient and local: stays in the view (nobody else needs it, §thin-view-layer)
const [hovered, setHovered] = useState(false);
const [draft, setDraft] = useState("");                               // text being typed; only the committed value goes to the model
<input value={draft} onChange={(e) => setDraft(e.target.value)} onBlur={() => model.setText(draft)} />
```

A drag or draft that previews a *document* change is not view-local: it is a model-level overlay (§transient-overlay, [§committed-transient-identity](../../../representations/bootstrap-patterns/committed-transient-identity/README.md)). Design-system components are presentational (props in, events out); only feature views connect to models ([§design-system-storybook](../design-system-storybook/typescript.md#two-kinds-of-components-presentational-and-connected)).

### Model lifecycle belongs to the host, not to React state

React owns a repeatable *mount scope* for a model (start/stop its subscriptions and timers); the model's creator owns final disposal (`dispose`).

```tsx
export function useReactiveModel<T extends { mount(): void; unmount(): void }>(model: T): T {
  useEffect(() => { model.mount(); return () => model.unmount(); }, [model]);   // repeatable (StrictMode, Suspense, tabs)
  return model;
}
```

Hidden UI does no work: the host passes an `active` flag (false for a hidden panel or folded section), and components skip subscriptions or timers while inactive but stay mounted.

Swapping React for another view library means rewriting this module and the components, not the models.

## Keeping the main thread responsive

Main-thread work is chunked and yields between slices (§responsive-ui). A task is a generator; the scheduler runs it in time-boxed slices. The same scheduler is the main-thread backend of a compute port; moving the work to a worker later changes composition, not callers (§location-agnostic-ui).

```ts
export type Task<T> = Generator<void, T, void>;                       // `yield` marks a safe point

const yieldToMain = (): Promise<void> => {
  const s = (globalThis as unknown as { scheduler?: { yield?: () => Promise<void> } }).scheduler;
  return s?.yield ? s.yield() : new Promise((r) => setTimeout(r, 0));
};

export async function runChunked<T>(
  task: Task<T>,
  opts: { budgetMs?: number; signal?: AbortSignal; onSlice?: (ms: number) => void } = {},
): Promise<T> {
  const budget = opts.budgetMs ?? 8;
  for (;;) {
    const start = performance.now();
    for (;;) {
      if (opts.signal?.aborted) {
        task.return(undefined as never);                               // run the generator's finally blocks
        throw opts.signal.reason ?? new DOMException("aborted", "AbortError");   // recognized by isCancelled (§error-codes)
      }
      const r = task.next();
      if (r.done) { opts.onSlice?.(performance.now() - start); return r.value; }
      if (performance.now() - start >= budget) break;                  // slice used up
    }
    opts.onSlice?.(performance.now() - start);
    await yieldToMain();                                               // let input and paint run
  }
}

function* buildIndexTask(items: readonly Item[], out: Index): Task<Index> {
  for (let i = 0; i < items.length; i++) {
    out.add(items[i]);
    if ((i & 255) === 255) yield;                                      // safe point every 256 items
  }
  return out;
}
// in a model: const index = await runChunked(buildIndexTask(items, new Index()), { signal });
```

```ts
// build-index.guard.test.ts (guard tier, see §test-tiers): the longest synchronous slice stays under one frame
it("never blocks longer than a frame", async () => {
  let worst = 0;
  await runChunked(buildIndexTask(makeItems(200_000), new Index()), { budgetMs: 8, onSlice: (ms) => { worst = Math.max(worst, ms); } });
  expect(worst).toBeLessThan(16);
});
```

## Tests need no browser

```ts
// todo.medium.test.ts: the model on the in-memory composition, through real commands and a real replica (§transactions)
test("a toggle is visible everywhere before the server answers, and settles when it does", async () => {
  const scope = new Scope();
  const app = createApp(local(scope), scope);                         // §ports-and-adapters
  const { replica, server } = replicaOver(app);                       // §transactions' FakeRemote: answers only on server.release()
  const list = new TodoListModel(appTodoPort(app, replica), sendVia(replica));
  list.mount();
  list.add("a");
  await server.received(1);                                           // applied locally, not yet confirmed
  const id = list.state.value.todos[0]!.id;
  list.item(id)!.toggle();
  await server.received(2);
  expect(list.state.value.todos[0]!.done).toBe(true);                 // the list sees it, not just the item
  expect(list.openCount.value).toBe(0);
  expect(list.item(id)!.state.value.pending).toBe(true);
  await server.release();                                             // confirm both
  expect(list.item(id)!.state.value.pending).toBe(false);
  list.dispose();                                                     // the test created the model, so it disposes it
  replica.dispose();
  await scope.dispose();
});

test("a rejected toggle disappears everywhere and leaves its error, keeping newer data", async () => {
  /* server.reject(next) with a coded error; a remote edit of the text lands while pending;
     expect done back to false, text = the remote edit, item.state.value.error.code = the server's code */
});
// helpers assumed: appTodoPort (watchAll over effective documents + replica.txs), replicaOver, statusByEntity (latest
// TxInfo per id in `touches`), newId (§committed-transient-identity), changes (§api-surfaces)
```

Binding tests need no component library: `react-dom/client` with React's `act` and a DOM shim is enough, so they stay in
the quick tier. One case every binding must pass: with an `equals` that always says "different", the immediate delivery
on subscribe still causes no render (that is what `changes()` is for, §api-surfaces). Status (pending, rejected, offline) is read
from the replica and its link and projected by a model, never stored in one (§optimistic-from-store).
