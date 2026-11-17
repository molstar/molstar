# Representations: bootstrap patterns

Starting implementations, correct as written: take one in where it fits, and it becomes your project's stable pattern. What bootstrap, stable and sketch mean and how to use a pattern: [../../_core/system.md](../../_core/system.md#code-patterns).

| Pattern | Name |
|---|---|
| [§committed-transient-identity](committed-transient-identity/README.md) | Committed vs transient state, identity, hashing |
| [§representations](representations/README.md) | Representations |
| [§transactions](transactions/README.md) | Transactions and app-wide optimistic state |

## Shared names

The names these patterns own (one owner and one name per concept; see `_core/system.md`, "Shared names").

| Concept | Name | Owner |
|---|---|---|
| Tree | `Node`, `sym`, `print`/`parse`, `Path`, `nodeAt`, `attrOf`, `walk`, `TreeState`, `Executor`, positional `Patch` (wire only) | §representations |
| Edits and identity | `Op`/`OpEnvelope`, `hashNode` (structural), `valueHash`/`Canonicalize` (value), `Doc`/`Draft` (`effective` = committed ⊕ drafts) | §committed-transient-identity |
| Transactions | `Tx`, `TxScope`, `runTx`, `UnitOfWork`, `Replica`, `TxStatus`, `TxRemote`, `txServer` | §transactions |
