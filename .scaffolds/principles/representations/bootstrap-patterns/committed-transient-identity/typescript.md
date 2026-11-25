# Committed vs transient state, identity, hashing: TypeScript

The concept, invariants and conformance cases are in [README.md](README.md); this file is the TypeScript implementation.

## Stable ids and canonical hashing (`repr-core`, on top of §representations)

```ts
// helpers assumed: ulid (time-sortable, coordination-free); H, the one module hasher (sha-256 or xxhash from a library)
declare const H: { atom(text: string): string; list(childHashes: readonly string[]): string };   // domain-separates atoms from lists

// Entities carry an opaque, never-reused id attribute: (frame (id "f_8f3a") (name "sketch") ...)
export const newId = (prefix: string) => `${prefix}_${ulid()}`;
export const idOf = (n: Node): string | undefined => { const v = attrOf(n, "id"); return typeof v === "string" ? v : undefined; };

/** Memo for ONE immutable snapshot, collected with it. Core-owned; a query cache may extend it. */
export interface RootMemo { forRoot<T>(root: Node, key: string | symbol, build: () => T): T }

const ID_INDEX = Symbol("id-index");
export const indexById = (root: Node, cache: RootMemo): ReadonlyMap<string, Path> =>   // per immutable snapshot
  cache.forRoot(root, ID_INDEX, () => {
    const m = new Map<string, Path>();
    walk(root, [], (n, p) => { const id = idOf(n); if (id) m.set(id, p); });
    return m;
  });
```

Merkle hash over a canonical encoding, cached on the immutable node so unchanged subtrees cost O(1):

```ts
const hashes = new WeakMap<object, string>();        // keyed by node identity: structural sharing => cache hits

export function hashNode(n: Node): string {
  if (typeof n !== "object" || n === null || isSym(n)) return H.atom(print(n));   // print is canonical: strings quoted, symbols bare, lossless numbers
  let h = hashes.get(n);
  if (h === undefined) { h = H.list(n.map(hashNode)); hashes.set(n, h); }
  return h;
}

// = hash of the list node (parts...): typed like any node, so 1 and "1" differ (String(p) would collide them)
export const hashParts = (...parts: readonly (string | number)[]): string => hashNode(parts);
```

Two hashes (§merkle-hashing): `hashNode` is structural (sync, op checks, storage, mirrors); `valueHash` hashes the canonical spelling,
for equality, cache keys and dedup. Core does not know which spellings are equivalent; vocabularies supply a canonicalizer
(a shape vocabulary, for example), so core stays free of schemas:

```ts
/** Returns the canonical spelling of a subtree, or the same node when it has only one. Must be idempotent and pure. */
export type Canonicalize = (n: Node, path: Path) => Node;
export const identityCanon: Canonicalize = (n) => n;

const valueHashes = new WeakMap<Canonicalize, WeakMap<object, string>>();   // per canonicalizer, per immutable node

export function valueHash(n: Node, canon: Canonicalize = identityCanon, path: Path = []): string {
  const c = canon(n, path);
  if (typeof c !== "object" || c === null || isSym(c)) return hashNode(c);
  let memo = valueHashes.get(canon);
  if (!memo) valueHashes.set(canon, (memo = new WeakMap()));
  let h = memo.get(c);
  if (h === undefined) { h = H.list(c.map((x, i) => valueHash(x, canon, [...path, i]))); memo.set(c, h); }
  return h;
}
// with identityCanon, valueHash(n) === hashNode(n); a canonicalizer that returns new nodes should memoize them (per input
// node) so the WeakMap cache hits; path lets it find the slot (via the vocabulary) without a global registry
```

Canonical order: see [README.md](README.md#shape).

What hashes buy (§merkle-hashing): O(1) equality, change detection for reconcilers, cache keys, content-addressed storage, efficient sync:

```ts
/** Subtrees that differ, found by descending only into children whose hashes differ. */
export function diffByHash(a: Node, b: Node, path: Path = []): Path[] {
  if (hashNode(a) === hashNode(b)) return [];
  if (!Array.isArray(a) || !Array.isArray(b) || a.length !== b.length || tag(a) !== tag(b)) return [path];
  const y = b as readonly Node[];                                // Array.isArray narrows to any[]; keep the Node type
  return (a as readonly Node[]).flatMap((ac, i) => diffByHash(ac, y[i]!, [...path, i]));
}
```

## Id-addressed, versioned, idempotent operations

```ts
export type Op =
  | { op: "set-attr"; target: string; key: string; value?: Node }       // (key value); omit value = remove the attribute
  | { op: "insert";   parent: string; after?: string; node: Node }        // node carries an id; after = entity sibling id (omit = first)
  | { op: "remove";   target: string }
  | { op: "move";     target: string; parent: string; after?: string };

export interface OpEnvelope {
  id: string;                     // operation id: makes application idempotent
  base: string;                   // committed root hash it was made against
  result?: string;                // expected root hash after applying (convergence check)
  origin: string;                 // "ui" | "agent" | "executor:<name>" | "remote:<peer>"
  label?: string;                 // undo-step name
  tx?: string;                    // the transaction it belongs to (§transactions): acks and rebases match by it
  ops: readonly Op[];
}

export const opErrors = defineErrors("repr", {
  "repr.patch.target_missing": { description: "An operation targets an id that is not in the document.", retryable: false },
  "repr.patch.diverged":        { description: "Result hash differs from the expected one.",             retryable: true },
  "repr.patch.invalid_op":      { description: "Insert without an id or with a taken id, or a move into its own subtree.", retryable: false },
});

// Locates ids via an id index; throws repr.patch.target_missing / invalid_op. Entity children follow a node's attributes,
// so `after` (the previous entity sibling, omitted = first entity position) pins a position exactly. Inverses:
//   set-attr t k v -> set-attr t k <old value>, or with value omitted if k was absent
//   insert p a n   -> remove id(n)
//   remove t       -> insert <old parent> <old previous entity sibling, or omitted> <removed subtree>
//   move t p a     -> move t <old parent> <old previous entity sibling, or omitted>
// Pinned semantics: an entity is a child with an `id`; any other tagged child is a single-valued attribute. New
// attributes go before the first entity; insert without `after` goes to the first entity position (the end if none);
// `id` never changes (repr.patch.invalid_op); the root is never removed; every id inside an inserted subtree must be
// free. A no-op returns the same root reference and no inverse. Attribute inverses are exact for single values only.
declare function applyOp(doc: Node, op: Op): { next: Node; inverse: readonly Op[] };

/** The empty document: the root is an entity, so inverses of top-level removes have a parent id to insert into. */
export const emptyDocument: Node = [sym("doc"), [sym("id"), "root"]];

export function applyEnvelope(doc: Node, e: OpEnvelope): { next: Node; inverse: readonly Op[] } {
  let next = doc; const inverse: Op[] = [];
  for (const op of e.ops) { const r = applyOp(next, op); next = r.next; inverse.unshift(...r.inverse); }
  if (e.result && hashNode(next) !== e.result)
    throw opErrors.error("repr.patch.diverged", undefined, { hint: "resync from the server snapshot" });
  return { next, inverse };
}
```

`e.base !== hashNode(doc)` is not an error in itself: see [README.md](README.md#shape).

## Document with drafts: transient overlay, one commit

```ts
export interface Draft {
  update(ops: readonly Op[]): void;            // append to the overlay (applied incrementally)
  replace(ops: readonly Op[]): void;           // replace the overlay (drag: only the latest position matters)
  commit(): OpEnvelope | undefined;            // ONE envelope = one undo step, one sync message; undefined if empty or closed
  cancel(): void;                              // committed state untouched
}

export interface Doc {
  readonly committed: ReadonlyBehavior<Node>;  // confirmed ⊕ pending local ops: history and persistence
  readonly effective: ReadonlyBehavior<Node>;  // committed ⊕ open drafts: what executors render
  begin(label: string, origin: string): Draft;
  apply(e: OpEnvelope): void;                  // local commits, commands, remote ops
  resync(snapshot: Node, pending: readonly OpEnvelope[]): void;   // after divergence: adopt the snapshot, replay pending
}

export interface History { record(e: OpEnvelope, inverse: readonly Op[]): void }   // undo = apply `inverse` as a new envelope

export function createDoc(initial: Node, history: History, sync: DocSync): Doc {
  const committed = createBehaviorSubject(initial), effective = createBehaviorSubject(initial);
  type Open = { ops: Op[]; out: Node; closed: boolean };
  const drafts: Open[] = [];
  const seen = new Set<string>();              // applied op ids; grows until resync (trimming at acked log positions: sketch)
  const top = (): Node => { const d = drafts.at(-1); return d ? d.out : committed.value; };   // not ??: null is a Node
  const applyOps = (n: Node, ops: readonly Op[]) => ops.reduce((x, op) => applyOp(x, op).next, n);
  const close = (d: Open) => { if (d.closed) return false; d.closed = true; drafts.splice(drafts.indexOf(d), 1); return true; };

  // Replay drafts from index `from` (committed changed, or a draft was replaced). Appends to the top draft skip this.
  function rebuild(from = 0) {
    for (let i = from; i < drafts.length; i++) {
      const d = drafts[i]!, base = i === 0 ? committed.value : drafts[i - 1]!.out;
      try { d.out = applyOps(base, d.ops); }
      catch (err) {                            // its target vanished (e.g. a remote remove): auto-cancel that draft
        if (!(isCoded(err) && err.code === "repr.patch.target_missing")) throw err;
        close(d); i--;
      }
    }
    effective.next(top());
  }
  committed.subscribe(() => rebuild());

  const doc: Doc = {
    committed, effective,
    begin(label, origin) {
      const d: Open = { ops: [], out: top(), closed: false };
      drafts.push(d);
      return {
        update(ops) {
          if (d.closed) return;
          if (drafts.at(-1) === d) {                                    // O(new ops), not O(all ops)
            const out = applyOps(d.out, ops);                          // a bad op throws and leaves the draft as it was
            d.ops.push(...ops); d.out = out; effective.next(out);
          } else { d.ops.push(...ops); rebuild(drafts.indexOf(d)); }
        },
        replace(ops) { if (!d.closed) { d.ops = [...ops]; rebuild(drafts.indexOf(d)); } },
        cancel() { if (close(d)) rebuild(); },
        commit() {
          if (!close(d)) return undefined;
          if (d.ops.length === 0) { rebuild(); return undefined; }
          const env: OpEnvelope = { id: newId("op"), base: hashNode(committed.value), origin, label, ops: d.ops };
          const before = committed.value;
          try { doc.apply(env); } catch (err) { rebuild(); throw err; }
          if (committed.value === before) rebuild();   // no-op ops: committed.next did not emit, so drop the overlay here
          return env;
        },
      };
    },
    apply(e) {
      if (seen.has(e.id)) return;                                                  // idempotent: no history, no re-send
      const { next, inverse } = applyEnvelope(committed.value, e);
      seen.add(e.id);
      committed.next(next);
      if (e.origin === "ui" || e.origin === "agent") history.record(e, inverse);  // only local intents are undo steps
      if (!e.origin.startsWith("remote:")) sync.send(e);                          // remote ops are never echoed
    },
    resync(snapshot, pending) {
      seen.clear();                            // the adapter resumes `incoming` after the snapshot's log position
      committed.next(pending.reduce((n, e) => { seen.add(e.id); return applyEnvelope(n, e).next; }, snapshot));
    },
  };
  return doc;
}
```

Commands ([§commands](../../../architecture/bootstrap-patterns/commands/README.md)) go through `apply` or a one-shot draft; gestures use drafts. Only a commit reaches history, persistence and collaborators.

## Example: drawing on a whiteboard, reconciled back

```
pointer down → d = doc.begin("draw", "ui")              overlay opens; nothing committed
pointer move → d.replace([insert stroke with samples])   executors render the effective state immediately
pointer up   → data.put(samples); d.replace(... ref); d.commit()   ONE envelope: one undo step, one sync message
```

```ts
// Pointer input becomes INTENTS, never truth. Hit-testing is the renderer's screen→board mapping.
function onMove(ev: PointerEvent, ctx: DrawCtx) {
  ctx.draft ??= ctx.doc.begin("draw", "ui");
  ctx.samples.push(...ctx.pickPoint(ev), ev.pressure);
  // While transient, samples are inline (coalesce to one replace per frame in production); never committed inline.
  ctx.draft.replace([insertStroke(ctx, [sym("f32"), ...ctx.samples])]);
}

async function onUp(ctx: DrawCtx) {
  const d = ctx.draft; ctx.draft = undefined; if (!d) return;
  try {
    const ref = await ctx.data.put(encodeF32(ctx.samples), { type: "application/x-strokes" });   // bulk data, once (§bulk-data-not-nodes; the tree holds a reference)
    d.replace([insertStroke(ctx, refNode(ref))]);
    d.commit();
  } catch (e) { d.cancel(); throw e; }           // a failed upload must not leave the overlay open (no-op if committed)
}
const onCancel = (ctx: DrawCtx) => { ctx.draft?.cancel(); ctx.draft = undefined; };

const insertStroke = (ctx: DrawCtx, samples: Node): Op => ({ op: "insert", parent: ctx.frameId,
  node: [sym("stroke"), [sym("id"), ctx.strokeId], [sym("brush"), ctx.brush], [sym("samples"), samples]] });
// helpers assumed: DrawCtx (doc, data, draft, samples, frameId, strokeId, brush, pickPoint), encodeF32
```

Reconciling an executor with the effective state is pure and idempotent, and it skips unchanged subtrees by hash (§executors-not-truth):

```ts
// composition root: the executor port and bindExecutor are §representations'
const off = bindExecutor(doc.effective, canvasRenderer(canvas));        // or a flattened board derived from doc.effective (a derived form, never committed)
// inside reconcile(state, hint): changed subtrees = hint?.prev ? diffByHash(hint.prev, state) : [[]]
```

Echo suppression (`origin: "executor:canvas"`): see [README.md](README.md#shape).

## Collaboration: default strategy, behind a port

```ts
/** What a Doc needs: the op channel. Presence is a separate channel (§live-uncommitted-sharing), so a Doc never depends on it. */
export interface DocSync { send(e: OpEnvelope): void }          // committed operations only
export interface SyncPort extends DocSync {
  incoming: Observable<OpEnvelope>;                             // server-ordered; includes acks of our own ops
  snapshot(): Promise<{ doc: Node; hash: string }>;
}
export interface PresencePort {                                 // ephemeral, LWW (§live-uncommitted-sharing)
  set(at: readonly Locator[]): void;
  peers: Observable<ReadonlyMap<string, readonly Locator[]>>;
}
```

The adapter's contract with `Doc` (a sketch) and its open points: see [README.md](README.md#shape) and [Open](README.md#open).

## Conformance and tests

```ts
export const docSuite = suite<Doc>("Doc",
  kase("cancel leaves committed state and hash unchanged", "quick", (d) => { /* begin, replace, cancel */ }),
  kase("valueHash with identityCanon equals hashNode; with a canonicalizer, equivalent spellings agree", "quick", (d) => { /* ... */ }),
  kase("commit records exactly one history entry; empty or second commit records none", "quick", (d) => { /* ... */ }),
  kase("applying the same envelope twice is a no-op (no history, no send)", "quick", (d) => { /* ... */ }),
  kase("a remote remove auto-cancels a draft targeting it", "quick", (d) => { /* repr.patch.target_missing */ }),
  kase("id-addressed ops survive concurrent inserts", "quick", (d) => { /* reordered envelopes, same result hash */ }),
  kase("appending N ops to a draft costs O(N) applyOp calls", "guard", (d) => { /* count calls */ }),
);
```

Property tests (quick tier) and the cases as data: see [README.md](README.md#conformance-cases).
