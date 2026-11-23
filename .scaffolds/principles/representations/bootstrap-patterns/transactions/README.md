# Transactions and app-wide optimistic state {#transactions}

- **Implementations:** [TypeScript](typescript.md)

Rules §optimistic-transactions, §safe-to-automate: a transaction is data; it runs all-or-nothing in one round trip; on a client it is applied
optimistically to the **shared store** (every view sees it), confirmed or rejected as a unit, and rebased when confirmed
state moves. One library (`lib/tx`) serves the UI, commands, agent tools, remote commits and background workers.

## Problem

When each component applies its own optimistic update and rolls it back in its own error handler, views disagree about
what happened, a rejection can overwrite newer data that arrived meanwhile, and a change spanning several documents
half-applies. When commands commit as they go, a preview, a failed step or a refused change leaves partial writes or
sends mail that should never have gone out. When every write is its own call, a loop costs a round trip per item.

## Shape

**The form.** A transaction (`Tx`) is data:

```
(tx (id "tx_01J...") (label "Add shape") (touches "sh_42" "f_7")
  (expect (doc "board:b_1") (hash "h:9f2c..."))                            ; optional preconditions (a conditional write: compare-and-set on a hash)
  (step s (shape.add (frame "f_7") (id "sh_42") (kind "rect")))            ; `s` binds the step's result
  (step (notification.add (collaborator "u_3") (shape (use s id))))        ; (use <name> <field...>) reads it
  (step (board.count-shape (board "b_1") (shape (use s id)))))
```

Fields: `id` (client-generated, §client-ids), `label` (the undo step's name), `origin` (as on an envelope,
§committed-transient-identity), `expect` (document, hash pairs), `steps` (an optional bind name, a command, an input
that may contain `use` references) and `touches` (entity ids, so views can show status per entity). Codes:
`tx.precondition` (retryable), `tx.unbound`, `tx.not_document`, `tx.id_conflict`. A failing step keeps its own code and
adds the step index to the details.

**Commands write into the transaction, never commit.** A command receives a `TxScope` with `read(doc) -> Node`
(committed ⊕ this transaction's ops so far) and `write(doc, ops)` (applied to the scope's view at once; fails like
`applyOp`). It writes every change through it, outbox entries included, also in a dry run. An effect on the world
(mail, a provider call) is an outbox entry performed after commit by a worker.

**`runTx(tx, store, {actor, signal, dryRun, commit}) -> { results, envelopes } | coded error`.** A `TxStore` gives the
documents it owns (`doc(id)`), the commands (`invoke(command, input, ctx)`) and optional validators
(`(doc, before, after, actor)`, throw to refuse, §validators). Steps run in order against private per-document
overlays; then, in one block with nothing interleaved, every `expect` is checked again, every document's envelope
(`id = <tx id>/<doc>`, `tx = <tx id>`) is computed against the current committed state and validated, and only then
(when `commit`) all are applied. This holds only while subscribers of committed state never write to a document
synchronously (write-backs are queued to a later task). `commit: false` computes envelopes only (the replica's local run, the server's run
before its log append). `dryRun` is what commands see; a preview never commits. A coded failure carries a `TxPartial`
(per document: base, the overlay, the ops so far). A store backed by a database does the same inside one database
transaction.

**The replica** (client). State per document: `confirmed` (server truth), `pending` (local transactions in submission
order, as envelopes per document) and the document's committed value = confirmed ⊕ pending (§committed-transient-identity).

| Member | Does |
|---|---|
| `submit(tx) -> TxInfo` | runs locally (commit false) on top of earlier pending, marks pending, rebases touched docs, sends one round trip; never throws for a rejection |
| `txs` | observable map of pending, failed and rescued transactions with status |
| `retry(id)` | rejected: a new attempt with a new id; outcome unknown: the same id again |
| `dismiss(id)` | clears a failed or rescued status, never a pending one |
| `undo()` | submits the inverse of the last confirmed local transaction as one new transaction |
| `rescues`, `discardRescue(id)` | rejected changes kept as data: per document base, copy (what the author saw) and ops |

`TxStatus` is `pending`, `confirmed`, `failed` (with the error and whether the server `rejected` it or the outcome is
unknown) or `rescued`. Rebase recomputes each touched document from confirmed plus what is still pending; a pending
transaction that no longer applies leaves the queue in every document at once and is marked failed (outcome unknown).
A transaction refused locally, or rejected by the server, is rescued when the error is retryable (the default policy).

**The remote** (`TxRemote`): `submit(tx) -> ok | error`, idempotent by (actor, tx id), resolving only after its
envelopes were delivered on `confirmed`; `confirmed` (every confirmed change with document, envelope and actor);
optional `resynced` (a document healed whole, or emptied with `emptyDocument` when it becomes invisible).

**The server** (`txServer`): one transaction at a time per store; looks up (actor, id) in the log and returns the stored
outcome if the content hash matches (`tx.id_conflict` if not); otherwise runs `runTx` (commit false), appends all
envelopes and the outcome in one database transaction (the store port's `TreeStore.appendTx`), records failures as
outcomes, then applies to memory and broadcasts. After the append nothing changes the outcome. Stored envelope ids
carry the actor. Server documents keep no undo history.

**Documents that come into being.** A store opens documents on first use from `emptyDocument`, so the first envelope
applies to the same root everywhere; a replica learns of one from the first confirmed envelope or its own first use,
and starts from a snapshot of the documents that exist when it connects. The hub decides per connection which
documents an actor receives through an access hook (`TxAccess.visible(actor, doc)`, optional `reach`), asked at
connect, at catch-up and on every change, and re-decides which documents a change reaches (a membership change takes
effect live); a document that becomes visible arrives whole, one that becomes invisible is emptied. Catch-up is no
existence oracle: a client naming a document it was never given gets nothing back.

**Unit of work** (§buffered-writes): `add(step)` buffers, `save(label)` sends all buffered steps as one transaction;
a failed save keeps it and the next save resends it under the same id. Where save boundaries go is a crash-loss
decision: everything an external effect depends on is saved before the effect, and everything after the last
boundary is redone under the same idempotency keys.

**One queue.** A draft commit in a replica is a one-step transaction, `(tx (step (ops.apply (doc ...) (ops ...))))`,
so gestures and commands share one pending queue, one status map and one rebase. There is one replica per store owner
(§one-queue-owner): a shared worker owns it and tabs subscribe, or each tab is a replica and the others are remote peers
on `confirmed`. Never a queue per component. Any port can contribute steps by registering its writes as commands
(§cross-port-batches).

## Invariants

- A transaction commits in every document it touched or in none; a failing step, a failing validator or a moved
  precondition commits nothing.
- A later step reads an earlier step's writes and its bound results; a `use` of an unbound name is `tx.unbound`.
- Commands never commit; a dry run changes no committed state.
- A precondition is checked before the steps and again at commit.
- Naming a document not visible to the actor returns nothing.
- Every view reads the same effective state, so a pending transaction appears, and a failed one disappears, everywhere
  on the same tick.
- Effective state is always confirmed ⊕ still-pending transactions in submission order; a remote change is never
  overwritten by a rejection.
- After confirmation, the effective state equals the server's, by hash.
- A transaction id is applied at most once per actor; a resubmission returns the original outcome; ids of different
  actors never interact.
- Once the server has appended a transaction, its outcome is ok, whatever later subscribers do.
- An undo is one transaction and is never itself an undo step.

## Conformance cases

- `runTx`: a failing step commits nothing in any document; a later step reads an earlier step's writes and bound
  result; committed moving under the overlay fails at commit, never half-applies; a dry run returns the envelopes and
  leaves hashes unchanged; `expect` with a stale hash fails with `tx.precondition`.
- A unit of work whose save failed resends the same tx id and is applied once.
- Replica: submit shows the change in every touched document before the server answers; a rejection removes it
  everywhere, records the code and loses no other change; a remote change during pending is kept and the pending
  transaction replayed on top; a pending transaction that no longer applies leaves every document at once and is failed.
- Pending transactions apply in submission order across documents.
- Retry after a rejection submits a new id; after an unknown outcome it resends the same id.
- Undo submits the inverse of the last confirmed local transaction as one transaction; undo twice undoes two steps.
- Dismiss clears a failed transaction and never a pending one.
- A transaction confirmed on the feed before its reply is undoable; a lost reply after our own broadcast arrived counts
  as confirmed.
- Server: a resubmitted id is applied once and returns the original outcome; the same actor reusing an id for other
  content gets `tx.id_conflict`; two actors using one id are independent; a throwing broadcast subscriber after commit
  leaves the outcome ok; the server's run is not a preview (commands see dry run false, outbox entries commit).
- Rescue: a rejected change is rescued with base, copy and ops; a change refused locally is rescued from the run's
  partial overlay, with the transaction kept as the intent.
- A document first seen on the feed starts empty; a resync with the empty root empties it.

## Open

- A pending transaction the host applied during a gap in the connection shows twice until its idempotent resubmission
  is answered (known limit of `resynced`).
- `TreeStore.appendTx` takes `(actor, tx id, key, items)` here; an older sketch of the store port may still show
  `appendTx(tx, items)` without actor or content key, and the keyed form is the one that holds.
- The undo depth and the set of confirmed ids the replica keeps are "bounded in practice"; no bound is stated.
