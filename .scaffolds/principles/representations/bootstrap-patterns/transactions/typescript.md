# Transactions and app-wide optimistic state: TypeScript

The concept, invariants and conformance cases are in [README.md](README.md); this file is the TypeScript implementation.

## The form

```
(tx (id "tx_01J...") (label "Add shape") (touches "sh_42" "f_7")
  (expect (doc "board:b_1") (hash "h:9f2c..."))                            ; optional preconditions (a conditional write: compare-and-set on a hash)
  (step s (shape.add (frame "f_7") (id "sh_42") (kind "rect")))            ; `s` binds the step's result
  (step (notification.add (collaborator "u_3") (shape (use s id))))        ; (use <name> <field...>) reads it
  (step (board.count-shape (board "b_1") (shape (use s id)))))
```

```ts
export type Use = { readonly $use: string; readonly path: readonly string[] };   // (use m id) after parsing
export interface TxStep { readonly bind?: string; readonly command: string; readonly input: unknown }   // input may contain Use
export interface Tx {
  readonly id: string;                         // client-generated (§client-ids): acks match by id
  readonly label: string;                      // the undo step's name
  readonly origin: string;                     // "ui" | "agent" | ... (§committed-transient-identity)
  readonly expect: readonly { doc: string; hash: string }[];
  readonly steps: readonly TxStep[];
  readonly touches?: readonly string[];        // entity ids, so views can show pending/failed per entity
}

export const txErrors = defineErrors("tx", {
  "tx.precondition": { description: "A document changed since the transaction was made; re-read and retry.", retryable: true },
  "tx.unbound":      { description: "A step uses a name no earlier step bound.", retryable: false },
  "tx.not_document": { description: "A step wrote to a document this store does not own.", retryable: false },
  "tx.id_conflict":  { description: "You already used this transaction id for a different transaction.", retryable: false },
});
// printTx/parseTx convert between the form and Tx; (use ...) becomes a Use; step failures keep the step's own code
// and add details (step index), so callers branch on the real cause (§namespaced-codes)
```

## Commands write into the transaction, never commit

```ts
/** Given to commands in CommandContext.tx (§commands). Reads see this transaction's own earlier writes. Commands write every
 *  change through it, documents and outbox entries alike, also when ctx.dryRun is set; whether the collected ops are
 *  committed is runTx's `commit`, which commands never see. Nothing outside the documents happens inside a command:
 *  an effect on the world (mail, a provider call) is an outbox entry written here and performed after commit by a
 *  worker, so a preview, a local run or a refused transaction never sends anything (§all-or-nothing). */
export interface TxScope {
  read(doc: string): Node;                     // committed ⊕ this transaction's ops so far
  write(doc: string, ops: readonly Op[]): void;   // applied to the scope's view immediately; throws like applyOp (§committed-transient-identity)
}

// a mutating command, written once, used by UI, agents, tests and the server alike
export const countShape: CommandSpec<{ board: string; shape: string }, { board: string; count: number }> = {
  name: "board.count-shape", description: "Count a newly added shape in the board's shape count.", effect: "mutate",
  input: schema(), output: schema(),
  async run({ board, shape }, ctx) {
    if (!ctx.tx) throw commandError("command.invalid_input", "board.count-shape runs inside a transaction");
    const root = ctx.tx.read(`board:${board}`);                // sees this transaction's own earlier writes
    const count = Number(attrOf(entityById(root, board), "shape-count") ?? 0) + 1;
    ctx.tx.write(`board:${board}`, [
      { op: "set-attr", target: board, key: "shape-count", value: count },
      { op: "set-attr", target: board, key: "last-shape", value: shape },
    ]);
    return { board, count };
  },
};
// helpers assumed: entityById (the entity with this id in a tree), attrOf (an attribute's value; §committed-transient-identity)
```

## Running a transaction: all or nothing

`Actor` is the actor of §commands' `CommandContext` (§recorded-actors: who, of what kind, acting for whom).

```ts
export interface TxStore {
  doc(id: string): Doc | undefined;            // the documents this store owns (§committed-transient-identity); see "Documents that come into being"
  commands: TxCommands;                        // structural: a registry (§commands) satisfies it, tx does not import commands
  validators?: readonly TxValidator[];         // run on every changed document before commit, on server and replicas (§validators)
}
/** What a transaction needs from the command layer; a command registry plus its policy (§commands) fits it structurally. */
export interface TxCommands {
  invoke(command: string, input: unknown, ctx: CommandContext): Promise<unknown>;
}
/** Checks a document's next state against its previous one; throws a coded error to refuse the whole transaction. */
export type TxValidator = (doc: string, before: Node, after: Node, actor: CommandContext["actor"]) => void;
export interface TxResult { readonly results: Readonly<Record<string, unknown>>; readonly envelopes: ReadonlyMap<string, OpEnvelope> }
/** What a failed run had computed, per document: the state it started from, the overlay and the ops so far. Attached to
 *  the thrown error (`txPartialOf`), so a replica can rescue a locally refused change (§visible-tx-status). */
export interface TxPartial { readonly docs: ReadonlyMap<string, { readonly base: Node; readonly copy: Node; readonly ops: readonly Op[] }> }
export const txPartialOf = (e: unknown): TxPartial | undefined => (e as { txPartial?: TxPartial } | null)?.txPartial;

/** Runs the steps against private overlays, then (when `commit`) commits every document in one synchronous block.
 *  `commit: false` computes the envelopes only: the replica's local run and the server's run before its log append.
 *  `dryRun` is what commands see: this invocation is a preview. A preview never commits. */
export async function runTx(tx: Tx, store: TxStore,
                            ctx: { actor: CommandContext["actor"]; signal?: AbortSignal; dryRun?: boolean; commit?: boolean })
  : Promise<TxResult> {
  const commit = (ctx.commit ?? true) && !ctx.dryRun;
  const owned = (id: string): Doc => store.doc(id) ?? fail(txErrors.error("tx.not_document", `no document ${id}`));
  const checkExpect = () => { for (const e of tx.expect)
    if (hashNode(owned(e.doc).committed.value) !== e.hash) throw txErrors.error("tx.precondition", `${e.doc} moved`); };

  const ops = new Map<string, Op[]>();         // per document, in order
  const base = new Map<string, Node>();        // committed when the transaction first read or wrote the document
  const view = new Map<string, Node>();        // base ⊕ ops so far, for reads inside the transaction
  const scope: TxScope = {
    read: (d) => view.get(d) ?? base.get(d) ?? base.set(d, owned(d).committed.value).get(d)!,
    write: (d, more) => {
      let n = scope.read(d);
      for (const op of more) n = applyOp(n, op).next;          // throws before anything is recorded
      view.set(d, n);
      (ops.get(d) ?? ops.set(d, []).get(d)!).push(...more);
    },
  };
  const partial = (): TxPartial => ({ docs: new Map([...ops].map(([d, o]) => [d, { base: base.get(d)!, copy: view.get(d)!, ops: o }])) });
  const failWith = (e: CodedError): never => { throw Object.assign(e, { txPartial: partial() }); };

  try {
    checkExpect();                                             // fail fast; checked again at commit
    const results: Record<string, unknown> = {};
    for (const [i, s] of tx.steps.entries()) {
      const input = resolveUses(s.input, results);             // tx.unbound for an unknown name
      try {
        const out = await store.commands.invoke(s.command, input,     // commands see dryRun, never commit
          { actor: ctx.actor, signal: ctx.signal, dryRun: ctx.dryRun ?? false, tx: scope, idempotencyKey: `${tx.id}#${i}` });
        if (s.bind) results[s.bind] = out;
      } catch (e) {
        if (!isCoded(e)) throw e;                              // nothing was committed: overlays are dropped
        throw new CodedError({ ...e, message: e.message, details: [sym("tx-step"), i, e.details ?? null] });   // own code kept
      }
    }

    // Commit phase: synchronous, so nothing can interleave. Validate every document first, then apply them all.
    checkExpect();                                             // the steps awaited: committed may have moved since
    const envelopes = new Map<string, OpEnvelope>();
    for (const [d, list] of ops) {
      const doc = owned(d);
      const env: OpEnvelope = { id: `${tx.id}/${d}`, tx: tx.id, base: hashNode(doc.committed.value), origin: tx.origin, label: tx.label, ops: list };
      const { next } = applyEnvelope(doc.committed.value, env);  // pure; throws if committed moved under the overlay
      for (const v of store.validators ?? []) v(d, doc.committed.value, next, ctx.actor);   // a refusal commits nothing
      envelopes.set(d, env);
    }
    if (commit) {                                              // same input as validated: applying cannot fail, but a
      let failed = false, first: unknown;                      // subscriber may throw (§api-surfaces rethrows after notifying all):
      for (const [d, env] of envelopes)                        // still apply every document, then rethrow
        try { owned(d).apply(env); } catch (e) { if (!failed) { failed = true; first = e; } }
      if (failed) throw first;
    }
    return { results, envelopes };
  } catch (e) {
    return isCoded(e) ? failWith(e) : Promise.reject(e);      // coded failures carry what was computed
  }
}
// helpers assumed: resolveUses (Use -> value at path; tx.unbound), fail (throws its argument), schema (§commands)
```

Why it is atomic in-process: steps may await, and committed state may move meanwhile; the overlays are rebuilt only from
our own ops, and the commit phase re-validates every document against the *current* committed value before applying any.
Between validation and the last `apply` there is no `await`, so no other change can land, provided `committed` subscribers
never write to a document synchronously (write-backs are queued to a later task). A store backed by a database does the
same inside one database transaction (conditional writes: compare-and-set on the version), and a remote store runs `runTx` server-side (below).

## The replica: confirmed, pending, effective

```ts
export type TxStatus =
  | { readonly state: "pending" }
  | { readonly state: "confirmed" }
  | { readonly state: "failed"; readonly error: ErrorInfo; readonly rejected: boolean }    // rejected: the server said no;
                                                                  // otherwise the outcome is unknown (transport, local rebase)
  | { readonly state: "rescued"; readonly error: ErrorInfo };    // rejected, and the author's change is kept as a Rescue
/** What views show: pending / confirmed / rejected / rescued / failed (outcome unknown). */
export const statusLabel = (s: TxStatus) =>
  s.state === "failed" ? (s.rejected ? "rejected" : "failed") : s.state;

/** A rejected change the author may still want (§visible-tx-status), kept as data, never as a new document: per document, the state
 *  the change was computed on, that state with the change (what the author saw) and the ops. A three-way merge by id
 *  (base, the server's state now, copy; version history's, if the project keeps one) brings it back. Rescued
 *  automatically on a server rejection, and also when the local run refuses a change because a concurrent one
 *  already arrived (otherwise the draft is all that is left). */
export interface Rescue {
  readonly tx: Tx; readonly error: ErrorInfo;
  readonly docs: ReadonlyMap<string, { readonly base: Node; readonly copy: Node; readonly ops: readonly Op[] }>;
}
export interface TxInfo { readonly tx: Tx; readonly status: TxStatus }

/** The server side, behind a port (the project's wire carries it as one frame each way). */
export interface TxRemote {
  /** Runs the transaction authoritatively as the connection's verified actor; idempotent by (actor, tx.id): a
   *  resubmission returns the original outcome (§client-ids).
   *  Resolves only after its envelopes were delivered on `confirmed`: the reply frame is ordered after the broadcast
   *  frames on the same connection. */
  submit(tx: Tx): Promise<{ ok: true } | { ok: false; error: ErrorInfo }>;
  confirmed: Observable<{ doc: string; env: OpEnvelope; actor: Actor }>;   // every confirmed change, ours and others',
                                                               // with who ran it (§recorded-actors); a tx id is unique per actor only
  /** Optional: a document healed whole (a reconnect past the host's buffer, a restarted host, a document that became
   *  visible). The replica adopts it as confirmed and rebases what is still pending. A document that became invisible
   *  is **emptied**: resynced with `emptyDocument` (§committed-transient-identity). Known limit: a pending tx the host applied during the gap shows
   *  twice until its idempotent resubmission is answered. */
  resynced?: Observable<{ doc: string; root: Node }>;
}

export interface Replica {
  /** Optimistic everywhere at once; never throws for a server rejection: the status says so (§visible-tx-status). */
  submit(tx: Tx): Promise<TxInfo>;
  readonly txs: ReadonlyBehavior<ReadonlyMap<string, TxInfo>>;   // pending and failed, until dismissed or retried
  retry(id: string): Promise<TxInfo>;
  dismiss(id: string): void;
  undo(): Promise<TxInfo> | undefined;                         // the last confirmed local tx, as a new tx of inverses
  readonly rescues: ReadonlyBehavior<ReadonlyMap<string, Rescue>>;   // by tx id, until discarded
  discardRescue(id: string): void;
  dispose(): void;
}

export function createReplica(store: TxStore & { docs(): Iterable<string> }, remote: TxRemote,
                              actor: CommandContext["actor"],
                              autoRescue: (e: ErrorInfo) => boolean = (e) => e.retryable): Replica {
  const confirmed = new Map<string, Node>();                   // server state per document
  for (const d of store.docs()) confirmed.set(d, store.doc(d)!.committed.value);
  /** A document seen for the first time (a first-use document, 'Documents that come into being') starts empty. */
  const confirmedOf = (d: string): Node => confirmed.get(d) ?? confirmed.set(d, emptyDocument).get(d)!;
  const pending: { tx: Tx; envelopes: ReadonlyMap<string, OpEnvelope> }[] = [];   // in submission order
  const txs = createBehaviorSubject<ReadonlyMap<string, TxInfo>>(new Map());
  const rescues = createBehaviorSubject<ReadonlyMap<string, Rescue>>(new Map());
  const setStatus = (tx: Tx, status: TxStatus) => txs.next(new Map(txs.value).set(tx.id, { tx, status }));
  const forget = (id: string) => { const m = new Map(txs.value); m.delete(id); txs.next(m); };
  const rescue = (tx: Tx, error: ErrorInfo, docs: Rescue["docs"]) => {   // keep what the author made (§visible-tx-status)
    rescues.next(new Map(rescues.value).set(tx.id, { tx, error, docs }));
    setStatus(tx, { state: "rescued", error });
  };

  /** effective = confirmed ⊕ pending, per document; a pending tx that no longer applies fails here (§rebase). */
  function rebase(docs: Iterable<string>) {
    const todo = new Set(docs);                                // a Set visits entries added while iterating
    for (const d of todo) {
      let n: Node = confirmedOf(d);
      const keep: OpEnvelope[] = [];
      for (const p of [...pending]) {
        const env = p.envelopes.get(d); if (!env) continue;
        try { n = applyEnvelope(n, env).next; keep.push(env); }
        catch (e) {
          if (!isCoded(e)) throw e;
          // all or nothing in effective too: leave the queue and recompute its other documents. The server still decides:
          // a confirmation arrives on `confirmed` like anyone else's, and submit() then clears this status.
          pending.splice(pending.indexOf(p), 1);
          setStatus(p.tx, { state: "failed", error: e, rejected: false });   // the server may still confirm it
          for (const o of p.envelopes.keys()) if (o !== d) { todo.delete(o); todo.add(o); }   // (re)visit
        }
      }
      store.doc(d)!.resync(confirmedOf(d), keep);              // Doc.committed = confirmed ⊕ pending (§committed-transient-identity); drafts rebuild on top
    }
  }

  const lastConfirmed = new Map<string, { inverse: readonly Op[] }>();   // our envelopes' inverses, for undo
  const seen = new Set<string>();                              // our tx ids confirmed on the feed (bounded in practice)
  const off = remote.confirmed.subscribe(({ doc, env, actor: by }) => {
    const ours = env.tx !== undefined && by.id === actor.id;   // ids are per actor: another's equal id is not ours
    if (ours) seen.add(env.tx!);
    const { next, inverse } = applyEnvelope(confirmedOf(doc), env);   // a first-use document is first seen here
    confirmed.set(doc, next);
    if (ours) lastConfirmed.set(env.id, { inverse });         // bounded in practice: keep the undo depth
    const i = ours ? pending.findIndex((p) => p.tx.id === env.tx) : -1;
    if (i >= 0) pending[i] = dropDoc(pending[i]!, doc);
    rebase([doc]);                                             // others' changes move confirmed: rebase what is still ours
  });
  const offResynced = remote.resynced?.subscribe(({ doc, root }) => {   // healed whole, or emptied (root = emptyDocument)
    confirmed.set(doc, root);
    rebase([doc]);
  });

  // Undo (§optimistic-transactions: one transaction, one undo step): the inverse of a confirmed local tx, computed against confirmed state
  // when it was confirmed. Ops address ids, so later unrelated changes do not break it; a conflicting one fails it as a tx.
  const undoStack: Tx[] = [];
  const undoIds = new Set<string>();                           // an undo is not itself an undo step: a second undo()
  function remember(tx: Tx, docs: readonly string[]) {         // must undo the step before, not redo this one
    if (tx.origin !== "ui" && tx.origin !== "agent") return;   // only local intents are undo steps (§committed-transient-identity)
    if (undoIds.delete(tx.id)) return;
    const steps = docs.flatMap((d) => { const env = lastConfirmed.get(`${tx.id}/${d}`);
      return env ? [{ command: "ops.apply", input: { doc: d, ops: env.inverse } }] : []; });
    if (steps.length) undoStack.push({ id: newId("tx"), label: `Undo ${tx.label}`, origin: tx.origin, expect: [], steps, touches: tx.touches });
  }

  let turn: Promise<unknown> = Promise.resolve();             // local runs one at a time, in submission order
  async function submit(tx: Tx): Promise<TxInfo> {
    const queued = turn.then(async () => {
      // 1. compute the change locally (same commands as the server) on top of every earlier pending tx; commit nothing
      const local = await runTx(tx, store, { actor, commit: false });
      const seenBy = new Map([...local.envelopes.keys()].map((d) => [d, store.doc(d)!.committed.value]));   // what the author saw
      // 2. make it visible everywhere: pending, then rebase every touched document
      pending.push({ tx, envelopes: local.envelopes });
      setStatus(tx, { state: "pending" });
      rebase(local.envelopes.keys());
      return { local, seenBy };                                // captured now: by the time the reply arrives, `confirmed`
    });                                                        // has emptied the pending entry
    turn = queued.catch(() => {});
    let run: Awaited<typeof queued>;
    try { run = await queued; }
    catch (e) {                                                // refused locally (a concurrent change already arrived, ...)
      if (!isCoded(e)) throw e;
      const part = txPartialOf(e);
      if (autoRescue(e)) rescue(tx, e, part?.docs ?? new Map());   // the tx itself keeps the author's intent as data
      else setStatus(tx, { state: "failed", error: e, rejected: true });
      return txs.value.get(tx.id)!;
    }
    const docs = [...run.local.envelopes.keys()];
    // 3. one round trip; the server's envelopes arrive on `confirmed` before this resolves
    const r = await remote.submit(tx).then(
      (x) => (x.ok ? x : { ...x, rejected: true }),
      // transport: outcome unknown, unless our own broadcast already arrived (it is ordered before the reply)
      (e) => (seen.has(tx.id) ? { ok: true as const } : { ok: false as const, error: toErrorInfo(e), rejected: false }));
    const i = pending.findIndex((p) => p.tx.id === tx.id);
    if (i >= 0) pending.splice(i, 1);
    if (r.ok) { forget(tx.id); remember(tx, docs); }           // confirmed: nothing left to show; undoable
    else if (r.rejected && autoRescue(r.error))                // the server said no: keep what the author saw (§visible-tx-status)
      rescue(tx, r.error, new Map(docs.map((d) => { const env = run.local.envelopes.get(d)!, b = run.seenBy.get(d)!;
        return [d, { base: b, copy: applyEnvelope(b, env).next, ops: env.ops }]; })));
    else setStatus(tx, { state: "failed", error: r.error, rejected: r.rejected });
    rebase(docs);                                              // drop what is left of ours; server truth is in confirmed
    return r.ok ? { tx, status: { state: "confirmed" } } : txs.value.get(tx.id)!;
  }

  return {
    submit, txs, rescues,
    retry: (id) => { const t = txs.value.get(id);
                     if (t?.status.state !== "failed") return Promise.reject(new Error(`no failed tx ${id}`));
                     forget(id);
                     // unknown outcome: resend the SAME id, the server dedups it; a rejection: a new attempt, a new id
                     return submit(t.status.rejected ? { ...t.tx, id: newId("tx") } : t.tx); },
    dismiss: (id) => { const st = txs.value.get(id)?.status.state; if (st === "failed" || st === "rescued") forget(id); },
    discardRescue: (id) => { const m = new Map(rescues.value); m.delete(id); rescues.next(m); forget(id); },
    undo: () => { const u = undoStack.pop(); if (u) undoIds.add(u.id); return u && submit(u); },
    dispose: () => { off(); offResynced?.(); },
  };
}
// helpers assumed: dropDoc (removes one document's envelope from a pending entry: that part is now in confirmed),
// toErrorInfo (a transport failure becomes a retryable coded error; assumed), newId, emptyDocument (§committed-transient-identity). A rescue's
// `base` is what the author saw when the change was made (confirmed ⊕ earlier pending); merging it back is a
// three-way merge by id (base, the server's state now, copy).
```

What this gives (everyone sees it at once, no hand-written rollback, fresh data never undone, status as data): see [README.md](README.md#invariants).

## Server side and other stores

```ts
export type Outcome = { ok: true } | { ok: false; error: ErrorInfo };

/** The authoritative side. One transaction at a time per store; the log first, memory after; idempotent by tx id. */
export function txServer(store: TxStore, log: TreeStore, broadcast: (doc: string, env: OpEnvelope, actor: Actor) => void,
                         report: (e: unknown) => void) {
  let line: Promise<unknown> = Promise.resolve();
  return (tx: Tx, actor: CommandContext["actor"]): Promise<Outcome> => {
    const run = line.then(async (): Promise<Outcome> => {
      // ids are per actor (§client-ids): the log keys outcomes by (actor, id), so another actor's id never collides, blocks
      // or reveals anything; the same actor reusing an id for other content is tx.id_conflict
      const key = txKey(tx);                                     // hash of steps and preconditions
      const seen = await log.outcome(actor.id, tx.id);
      if (seen) return seen.key === key ? seen.outcome : { ok: false, error: toErrorInfo(txErrors.error("tx.id_conflict")) };
      let items: { doc: string; expected: string; env: OpEnvelope }[];
      try {
        const { envelopes } = await runTx(tx, store, { actor, commit: false });  // a real run, not a preview; nothing applied yet
        items = [...envelopes].map(([doc, env]) => ({ doc, expected: env.base, env: { ...env, id: `${actor.id}/${env.id}` } }));
        await log.appendTx(actor.id, tx.id, key, items);         // ONE database transaction: all envelopes + the outcome, each a conditional write on its expected base
      } catch (e) {
        if (!isCoded(e)) throw e;
        await log.recordOutcome(actor.id, tx.id, key, { ok: false, error: e });
        return { ok: false, error: e };
      }
      // committed: nothing below may change the outcome. A throwing subscriber is reported, never turned into {ok:false}.
      for (const { doc, env } of items) store.doc(doc)!.apply(env);   // memory follows the log; same input, cannot fail
      for (const { doc, env } of items)                          // then the reply is sent: ordered after these frames
        try { broadcast(doc, env, actor); } catch (e) { report(e); }
      return { ok: true };
    });
    line = run.catch(() => {});
    return run;
  };
}
// Server documents are created with a no-op History (§committed-transient-identity): undo is the client's transaction, not the server's.
// txKey: H(printed steps and preconditions); the log stores it with the outcome under (actor, id) (the store port's TreeStore.appendTx).
// Stored envelope ids carry the actor, so equal tx ids of two actors never meet in the log.
```

### Documents that come into being, and who may see them

Documents opened on first use: see [README.md](README.md#shape).

```ts
/** A TxStore's documents, opened on first use from §committed-transient-identity's emptyDocument. */
export function docStore(make: (id: string, root: Node) => Doc): Pick<TxStore, "doc"> & { docs(): Iterable<string> } {
  const docs = new Map<string, Doc>();
  return { doc: (id) => docs.get(id) ?? docs.set(id, make(id, emptyDocument)).get(id)!, docs: () => docs.keys() };
}
```

Who may see which documents (the hub's access hook): see [README.md](README.md#shape). The policy itself is the application's (§application-policies); the hub only holds the hook:

```ts
export interface TxAccess {
  visible(actor: Actor, doc: string): boolean | Promise<boolean>;
  reach?(doc: string, docs: readonly string[]): readonly string[];   // the documents a change to `doc` re-decides
}
```

The same `runTx` serves a workflow step (if the project has workflows: one transaction per command step), an agent's
write through a domain tool, agent start-up and remote batch commits (several commands in one transaction).

## Batches across ports, and save boundaries

Batches across ports and save boundaries: see [README.md](README.md#shape).

```ts
/** Buffer writes inside a loop; send them as ONE transaction at a save boundary (§buffered-writes). */
export interface UnitOfWork {
  add(step: TxStep): void;                     // nothing is sent
  readonly size: number;
  /** One transaction, one round trip; undefined if there is nothing to save. A failed save keeps the transaction, and
   *  the next save resends it under the SAME id, so the server applies it once (§client-ids). New steps wait for the next one. */
  save(label: string): Promise<TxResult | undefined>;
}

export function unitOfWork(submit: (tx: Tx) => Promise<TxResult>, origin: string): UnitOfWork {
  let steps: TxStep[] = [];
  let inFlight: Tx | undefined;
  return {
    add: (s) => { steps.push(s); },
    get size() { return steps.length + (inFlight ? inFlight.steps.length : 0); },
    async save(label) {
      if (!inFlight) {
        if (steps.length === 0) return undefined;
        inFlight = { id: newId("tx"), label, origin, expect: [], steps };
        steps = [];
      }
      const r = await submit(inFlight);        // throws: inFlight is kept for the resend
      inFlight = undefined;
      return r;
    },
  };
}
```

Read-your-writes inside the loop is the caller's in-memory view (a UI reads `effective`, which already includes pending transactions); where the boundaries go: see [README.md](README.md#shape).

## Conformance

```ts
export const txSuite = suite<TxStore>("Transactions",
  kase("a failing step commits nothing in any document", "quick", async (s) => { /* ... */ }),
  kase("a later step reads an earlier step's writes and its bound result", "quick", async (s) => { /* ... */ }),
  kase("committed moving under the overlay fails at commit, never half-applies", "quick", async (s) => { /* ... */ }),
  kase("dry-run returns the envelopes and changes nothing", "quick", async (s) => { /* hashes unchanged */ }),
  kase("a unit of work whose save failed resends the same tx id and is applied once", "quick", async (s) => { /* ... */ }),
  kase("expect with a stale hash fails with tx.precondition", "quick", async (s) => { /* ... */ }),
);

export const replicaSuite = suite<{ replica: Replica; remote: FakeRemote; docs: TxStore }>("Replica",
  kase("submit shows the change in every touched document before the server answers", "quick", async (r) => { /* ... */ }),
  kase("rejection removes it everywhere and records the code; no other change is lost", "quick", async (r) => { /* ... */ }),
  kase("a remote change during pending is kept; the pending tx is replayed on top", "quick", async (r) => { /* ... */ }),
  kase("a pending tx that no longer applies leaves every document at once and is failed", "quick", async (r) => { /* ... */ }),
  kase("after confirmation effective equals the server state, byte for byte", "quick", async (r) => { /* hashNode */ }),
  kase("order: pending txs apply in submission order across documents", "quick", async (r) => { /* ... */ }),
  kase("retry after a rejection submits a new id; after an unknown outcome it resends the same id", "quick", async (r) => { /* ... */ }),
  kase("a resubmitted id is applied once and returns the original outcome", "quick", async (r) => { /* txServer */ }),
  kase("undo submits the inverse of the last confirmed local tx as one tx", "quick", async (r) => { /* ... */ }),
  kase("dismiss clears a failed tx and never a pending one", "quick", async (r) => { /* ... */ }),
  kase("a tx confirmed on the feed before its reply is undoable (touched docs captured at submit)", "quick", async (r) => { /* ... */ }),
  kase("undo twice undoes two steps; an undo is never itself an undo step", "quick", async (r) => { /* ... */ }),
  kase("a lost reply after our own broadcast arrived counts as confirmed", "quick", async (r) => { /* ... */ }),
  kase("a throwing broadcast subscriber after commit leaves the outcome ok", "quick", async (r) => { /* txServer */ }),
  kase("the same actor reusing an id for other content gets tx.id_conflict", "quick", async (r) => { /* txServer */ }),
  kase("two actors using the same tx id are independent: neither blocks nor learns of the other", "quick", async (r) => { /* txServer */ }),
  kase("a rejected change is rescued with base, copy and ops; nothing typed is lost", "quick", async (r) => { /* ... */ }),
  kase("a change refused locally is rescued from the run's partial overlay, the tx kept as the intent", "quick", async (r) => { /* ... */ }),
  kase("a document first seen on the feed starts empty; a resync with the empty root empties it", "quick", async (r) => { /* ... */ }),
  kase("the server's run is not a preview: commands see dryRun false, and outbox entries commit with the tx", "quick", async (r) => { /* ... */ }),
);
// FakeRemote: an in-memory server running txServer, with controllable latency and rejection
```
