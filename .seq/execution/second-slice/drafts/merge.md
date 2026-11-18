# Merge in `@molstar/core/repr` (draft, for the author: lane K2's short agreement)

Split from `patches.md` at the author's word (2026-11-18, "let's go with 1-5"). It answers the second review's R2-B2 and
R2-B3. Once agreed, it moves to `.seq/design/merge.md`. Lane K2 builds it after lane K, and before lane N, whose rebase
test needs it.

## The function

```ts
function merge3(base: Node, ours: Node, theirs: Node, options?: MergeOptions): MergeResult;
interface MergeOptions {
  readonly language?: Language;
  readonly validators?: readonly PatchValidator[];
  readonly allowDangling?: boolean;
  readonly allowRebind?: boolean;
}
type MergeResult =
  | { readonly ok: true; readonly root: Node }
  | { readonly ok: false; readonly conflicts: readonly Conflict[]; readonly partial: Node };
interface Conflict {
  readonly kind: ConflictKind;
  readonly at: Address;
  readonly ours: readonly Op[];
  readonly theirs: readonly Op[];
  readonly issues?: readonly Issue[];      // the validators' issues, for `invalid`
}
type ConflictKind = 'set' | 'remove' | 'insert' | 'move' | 'replace' | 'order' | 'invalid';
```

- **Only `ok: true` has a `root`** (R2-B2): a caller cannot commit a merge and ignore its conflicts by accident.
- **`partial`** is `ours` with every change of theirs that is independent of every conflict. It is applied and
  validated as one patch, under the options. If that fails, `partial` is `ours`. A host shows `partial` with the
  conflicts, or rescues theirs as a patch to review.
- **Without a language** the merge is not validated, and it says so (an info).

## How

- `diff(base, ours)` and `diff(base, theirs)`, then each change is keyed by the node it touches (an address in `base`).
- **Changes that both apply**: to different nodes; to different keys of one node; equal changes, applied once.
- **Theirs is applied as one patch** on `ours`, under the options. Dangling refs, rebinding and the validators are
  judged on the merged root against `ours` (patches.md, Applying 6 and 7). This catches what no single op conflicts on:
  - a ref theirs wrote to a node only theirs inserted, where ours inserted another node with that id (R2-B2b);
  - two follows that make a cycle (R2-B2c).
- **Independence** is by node: a change of theirs depends on the changes of theirs to the nodes it names. That is its
  target, its parent, and the targets of the refs it writes. A change whose dependency conflicts is left out with it,
  and listed in the same conflict's `theirs`.

## Conflicts, complete for the six ops (R2-B3)

| Kind | Ours | Theirs |
|---|---|---|
| `set` | `set` of a key | `set` of the same key to another value, or `unset` of it |
| `remove` | `remove` of a node | any change at or under it: `set`, `unset`, `insert` under it, `move` of it or into it, `replace` of it or under it |
| `insert` | `insert` of an id | `insert` of the same id under the same parent with other content |
| `move` | `move` of a node | `move` of the same node to another parent or place, or another `rename` |
| `replace` | `replace` of a node | any change under it, or a different `replace` of it |
| `order` | moves of children of one parent | moves of the same child (not of different children) |
| `invalid` | — | the merged root its validators refuse, with their issues |

- Each row is symmetric: ours and theirs may swap.
- **Moves of different children** of one parent both apply. Ours are placed first. Theirs follow, each after its
  predecessor in theirs if that sibling still stands, else after its predecessor in ours. The result is deterministic.

## Laws (property tests)

- `merge3(b, x, b)` is `{ ok: true, root: x }`;
- `merge3(b, b, y)` is `{ ok: true, root: y }`, when `y` is valid;
- `merge3(b, x, x)` is `{ ok: true, root: x }`;
- with no conflicts and no moves of the same parent's children on both sides, the merge is symmetric;
- every `partial` is valid under the options.

## Rebase

A stale patch `p` made on `base` is rebased as `merge3(base, current, applyPatch(base, p).root)`. `base`'s root comes
from `SpecModel`'s recent roots, or from the host. It is the host's opt-in (patches.md, decision 4).

## For the author

1. **`partial` as ours plus theirs' independent changes** (recommended), or nothing but ours?
2. **Ours first on moves of different children** of one parent (recommended: deterministic, and ours is the live
   view)?
