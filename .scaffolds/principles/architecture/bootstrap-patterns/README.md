# Architecture: bootstrap patterns

Starting implementations, correct as written: take one in where it fits, and it becomes your project's stable pattern. What bootstrap, stable and sketch mean and how to use a pattern: [../../_core/system.md](../../_core/system.md#code-patterns).

| Pattern | Name |
|---|---|
| [§api-surfaces](api-surfaces/README.md) | API surfaces: sync / async / streaming, behaviors |
| [§commands](commands/README.md) | Commands |
| [§error-codes](error-codes/README.md) | Errors with stable codes |
| [§hot-loop-allocations](hot-loop-allocations/README.md) | Allocations, arenas, data structures |
| [§ports-and-adapters](ports-and-adapters/README.md) | Ports and adapters, composition root |

## Shared names

The names these patterns own (one owner and one name per concept; see `_core/system.md`, "Shared names").

| Concept | Name | Owner |
|---|---|---|
| Streams and state | `Observable`, `createSubject`, `ReadonlyBehavior`, `createState`/`update`, `select`, `changes` | §api-surfaces |
| Wiring | `Scope`, `production`/`local`/`test`, `createApp` | §ports-and-adapters |
| Commands | `CommandRegistry`, `Invoker`, `bindInvoker`, `CommandContext`, `reprArg` | §commands |
| Errors | `CodedError`, `ErrorInfo`, `defineErrors`, `isCancelled` | §error-codes |
