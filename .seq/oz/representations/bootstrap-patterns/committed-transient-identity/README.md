# Committed vs transient state, identity, hashing {#committed-transient-identity}

- **Implementations:** [TypeScript](typescript.md)

Rules §committed-vs-transient, §identity-and-hashing: the document is committed state, and an interaction is a draft overlay that ends in one commit or a cancel. Entities have stable ids, subtrees have cached canonical hashes, and operations are id-addressed, versioned and idempotent. Executors reconcile from state and are never the truth.

Where a piece of state lives (view-local, overlay, committed) is decided by §three-state-places; presence (§live-uncommitted-sharing) is a separate channel (see `SyncPort` below).

## Problem

If every pointer move writes to the document, a drag becomes hundreds of undo steps and sync messages, and a cancelled
gesture leaves traces. If edits address positions, two people editing at once corrupt each other's changes. Without
stable ids and cached hashes, every equality check, cache key and sync step walks whole trees, and a renderer that
writes back what it computed becomes a second source of truth.

## Shape

**Ids.** An entity is a list with an `id` attribute: `(frame (id "f_8f3a") (name "sketch") ...)`. Ids are opaque,
prefixed, unique, never reused and generated without coordination (`newId(prefix)`, time-sortable). `idOf(n)` reads
it; `indexById(root)` maps id to `Path`, built once per immutable snapshot (`RootMemo`).

**Two hashes** (§merkle-hashing):

| Name | What it hashes | Used for |
|---|---|---|
| `hashNode(n)` | structure: an atom hashes its printed form, a list hashes its children's hashes, domain-separated from atoms | sync, op checks, storage, mirrors |
| `valueHash(n, canon)` | the canonical spelling given by a `Canonicalize(n, path) -> n` function (pure, idempotent; `identityCanon` by default) | equality, cache keys, dedup |

Both are cached per immutable node, so an unchanged subtree costs O(1) and a one-leaf change rehashes only its path to
the root. Core knows no schema: vocabularies supply the canonicalizer (e.g. a shape vocabulary maps a box given by two corners to one given by origin and size). Where a vocabulary says
attribute order is not meaningful, nodes are normalized at construction or validation; the hash itself is
order-sensitive. `hashParts(...)` hashes a list of parts typed like any node, so `1` and `"1"` differ.
`diffByHash(a, b) -> paths` descends only into children whose hashes differ.

**Operations** (`Op`), addressed by id, never by position:

| Op | Fields | Inverse |
|---|---|---|
| `set-attr` | target, key, value (omitted = remove the attribute) | set-attr to the old value, or with value omitted if absent |
| `insert` | parent, after (an entity sibling id; omitted = first entity position), node (carries a free id) | remove the node's id |
| `remove` | target | insert into the old parent after the old previous entity sibling |
| `move` | target, parent, after | move back to the old parent and previous sibling |

Pinned semantics: an entity is a child with an `id`; any other tagged child is a single-valued attribute. New
attributes go before the first entity; insert without `after` goes to the first entity position (the end if none);
`id` never changes; the root is never removed; every id in an inserted subtree must be free; a move into its own
subtree is refused. A no-op returns the same root and no inverse. `applyOp(doc, op) -> { next, inverse } | coded
error`, with `repr.patch.target_missing` and `repr.patch.invalid_op`. `emptyDocument` is `(doc (id "root"))`.

**Envelope** (`OpEnvelope`): `id` (makes application idempotent), `base` (committed root hash it was made against),
`result` (expected root hash after, optional), `origin` (`ui`, `agent`, `executor:<name>`, `remote:<peer>`), `label`
(the undo step's name), `tx` (its transaction, §transactions) and `ops`. `applyEnvelope(doc, e) -> { next, inverse }`
applies ops in order, returns the inverse of the whole envelope, and fails with `repr.patch.diverged` when `result`
does not match. A `base` different from the current hash is not an error: id-addressed ops usually still apply; it
tells the sync layer to rebase pending local ops (§idempotent-ops).

**Doc and drafts.** A `Doc` has `committed` (confirmed ⊕ pending local ops: history and persistence) and `effective`
(committed ⊕ open drafts: what executors render), both observable. `begin(label, origin) -> Draft`;
`apply(envelope)` for local commits, commands and remote ops; `resync(snapshot, pending)` adopts a snapshot and replays
pending envelopes. A `Draft` has `update(ops)` (append), `replace(ops)` (only the latest matters, as in a drag),
`commit() -> envelope | nothing` (one envelope = one undo step, one sync message) and `cancel()`. A doc is built with a
`History` (`record(envelope, inverse)`; undo applies the inverse as a new envelope) and a `DocSync` (`send(envelope)`).

```
pointer down → d = begin("draw", "ui")           overlay opens; nothing committed
pointer move → d.replace([insert stroke])         executors render effective at once
pointer up   → put samples as a blob; d.replace([insert stroke with a ref]); d.commit()
```

**Executors reconcile** from `effective` and skip unchanged subtrees by hash (§executors-not-truth). An executor that
writes derived values back submits an envelope with origin `executor:<name>`, and skips subtrees whose hash equals what
it last rendered or wrote, so the loop settles in one round.

**Collaboration ports.** `DocSync` is all a doc needs. `SyncPort` adds `incoming` (server-ordered envelopes, including
acks of ours) and `snapshot()`. `PresencePort` (`set(locators)`, `peers`) is ephemeral and last-writer-wins, and a doc
never depends on it.

**The default adapter** (sketch: only its contract with `Doc` is fixed). The default for a small team is a
server-ordered log (§collaboration-adapters). The adapter keeps `confirmed` (server state up to the last entry it
received) and `pending` (local envelopes not yet acked, in order). An incoming envelope with the id of the head of
`pending` is an ack. Any other is applied to `confirmed` (checking its `result`); then, if `pending` is empty, it goes to
`doc.apply`, otherwise the adapter calls `doc.resync(confirmed, pending)`, because applying it after our pending ops
would reverse server order and non-commuting ops would diverge. On `repr.patch.diverged` it fetches a snapshot and
resyncs the same way. A CRDT adapter can implement the same port where concurrent sequence or text editing needs it.
The replica of §transactions is the default client side of this contract.

## Invariants

- Committed state changes only through `apply` (directly or by a draft's commit); drafts never reach history,
  persistence or collaborators.
- A commit produces exactly one envelope and one history entry; an empty, cancelled or second commit produces none.
- Applying an envelope whose id was already applied changes nothing: no history entry, no send.
- Only local intents (`ui`, `agent`) are recorded as undo steps; `remote:` envelopes are never sent back.
- `applyOp` either returns the next root and an exact inverse, or throws a coded error and changes nothing; untouched
  subtrees keep their identity.
- Applying an envelope's inverse restores the previous document exactly (attribute inverses are exact for single values).
- `hashNode(parse(print(x)))` equals `hashNode(x)`; distinct atoms of different kinds hash differently.
- `diffByHash(a, b)` is empty if and only if the hashes are equal.
- With `identityCanon`, `valueHash` equals `hashNode`; equivalent spellings under a canonicalizer share a `valueHash`
  and keep distinct `hashNode`s.
- A draft whose target vanished in a committed change is cancelled on its own; any other rebuild failure surfaces.
- Appending N ops to the top draft costs O(N) op applications, not a replay of all earlier ops.

## Conformance cases

Fixture `D` = `(doc (id "root") (title "T") (frame (id "f1") (name "sketch") (stroke (id "s1")) (stroke (id "s2"))) (frame (id "f2") (name "notes")))`.

- `set-attr f1 color "red"` on `D` gives `(frame (id "f1") (name "sketch") (color "red") (stroke (id "s1")) (stroke (id "s2")))`; replace, add and remove each invert back to `D`.
- `set-attr f1 name "sketch"` (no change) returns `D` itself and no inverse; removing an absent key returns `D`.
- `insert` into `f1` without `after` puts the node before `s1`; with `after s1` between `s1` and `s2`; into `f2`
  (no entities) at the end; the inverse removes it.
- `remove` of `s1`, `s2`, `f1`, `f2` each inverts back to `D` with positions kept; `move s1` to `f2` gives
  `(frame (id "f2") (name "notes") (stroke (id "s1")))`; `move s1` within `f1` after `s2` swaps them; both invert.
- Errors: remove `nope`, insert into `nope`, insert after `ghost` → `repr.patch.target_missing`; insert a node with
  the taken id `s1`, insert `(stroke)` without an id, move `f1` into `s1` or into itself, set `id` on `f1`, remove
  `root` → `repr.patch.invalid_op`.
- An envelope of insert, set-attr and move applies in order, and its inverse restores `D`; a wrong `result` gives
  `repr.patch.diverged`, the right one passes.
- Three envelopes (insert after `s1`, rename `f2`, insert into `f2`) applied in orders abc, cba and bac reach the same
  hash.
- Hashes: `1`, `1n`, `"1"`, the symbol `one`, `"one"`, `(1)`, `((1))`, `()`, `nil`, `true`, `"true"`, `0`, `-0`, NaN,
  `"nil"` give 15 distinct hashes; changing `(z 2)` to `(z 3)` in `(r (x 1) (y (z 2)) (w "k"))` keeps the hashes of
  `(x 1)` and `(w "k")` and `diffByHash` returns the single path `2 1 1`.
- Under a canonicalizer that sorts attributes, `(m (a 1) (b 2))` and `(m (b 2) (a 1))` share a `valueHash` and differ
  in `hashNode`; `(m (a 1) (b 3))` differs in both.
- 2000 ids from `newId("m")` are unique, sort in creation order and are prefixed `m_`.
- Drafts: update then replace shows the latest in `effective` while `committed` and its hash stay; cancel restores
  `effective` to `committed` and records and sends nothing.
- A commit of two ops records one entry and sends one envelope; a second commit and an empty draft's commit return
  nothing.
- A remote remove of `f1` cancels a draft renaming `f1` and keeps a draft renaming `f2`.
- A failing op in `update` on the top draft throws and leaves the draft as it was; a draft of no-op ops commits without changing
  `committed` and drops the overlay.
- `resync` adopts the snapshot, replays pending, rebuilds drafts on top, and treats a replayed pending id as applied.

## Open

- The pattern implements document-unique ids, the special case of §stable-entity-ids where the id alone is the
  address. Addresses as paths of sibling-unique ids, with operations rebased through moves, are not implemented here.
- The hasher is a module choice (sha-256 or xxhash) and the atom/list domain separation is not specified, so golden
  hashes are per implementation; the reference prototype uses sha-256 of `a:<printed atom>` and of
  `l:<child hashes joined by ,>`. Runtimes in different languages that sync must agree on one; the pattern does not say which.
- Adapter: `resync` throws if a pending op no longer applies (dropping or reporting it is the adapter's policy);
  history inverses recorded before a resync may be stale; the set of applied ids is only cleared on resync (trimming
  at acked log positions is a sketch).
- The O(N) draft-append guard is not yet tested in the reference prototype (it needs an instrumented `applyOp`).
- `update` on a draft below the top appends its ops before the rebuild, so a failing op stays in the draft (the rebuild
  throws), and a `target_missing` failure cancels that draft on its own instead; whether a nested `update` should behave
  like the top draft's is not said.
