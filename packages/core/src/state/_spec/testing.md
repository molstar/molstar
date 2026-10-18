# Testing

`pnpm test:unit packages/core/src/state` (Vitest, `_test/`). Test transformers are made in each file
(`StateTransformer.create` under a namespace of the file's own); a state runs tasks directly
(`runTask: (t) => t.run()`).

| File                     | What it checks                                                                         |
| ------------------------ | -------------------------------------------------------------------------------------- |
| `state.test.ts`          | `State.dispose` runs every transformer's `dispose`                                     |
| `dependencies.test.ts`   | links from the three sources, relinking, cycles, deferral                              |
| `cell-less.test.ts`      | cells below Null objects, waiters, deleted dependencies, stalls (dependencies.md §3-4) |
| `action-manager.test.ts` | the action registry                                                                    |

## How to test a change

- A change to the reconciler comes with a test that fails without it, through the public API: build with
  `state.build()`, commit with `updateTree`, then check the cells (`status`, `obj`, `errorText`, `dependencies`) and,
  where order matters, the events (reconcile.md §7).
- Whether an update throws is part of the result: a node's failure is its cell's `error`; an update that rejects is a
  failure of the whole pass (dependencies.md §2, §4).
- The two diffs (reconcile.md §2) are separate paths: a case that depends on what is found (cells below Null objects,
  cells left behind by a throw) needs a commit with one edit (the fast track) and one with more.
- The state logs every transformer error (`console.error`): a test that makes errors on purpose mutes it
  (`vi.spyOn(console, 'error')`).
