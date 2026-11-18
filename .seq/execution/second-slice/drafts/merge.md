# Merge in `@molstar/core/repr` (draft, for the author: lane K2's short agreement)

Split from `patches.md` at the author's word (2026-11-18, "let's go with 1-5"). It answers the second review's R2-B2 and
R2-B3, and its re-review's N-5 to N-8 ("apply N-1 to N-14 as recommended"). Once agreed, it moves to `.seq/design/merge.md`. Lane K2 builds it after lane K, and before lane N, whose rebase
test needs it.

## The function

```ts
function merge3(base: Node, ours: Node, theirs: Node, options?: MergeOptions): MergeResult;
interface MergeOptions {
  readonly language?: Language;
  readonly validate?: ValidateOptions;      // the slots' languages, the units hook (patches.md, N-1)
  readonly validators?: readonly PatchValidator[];
  readonly allowDangling?: boolean;
  readonly allowRebind?: boolean;
}
type MergeResult =
  | { readonly ok: true; readonly root: Node; readonly changes: Changes; readonly issues: readonly Issue[] }   // changes against ours
  | { readonly ok: false; readonly conflicts: readonly Conflict[]; readonly partial: Node; readonly issues: readonly Issue[] };
interface Conflict {
  readonly kind: ConflictKind;
  readonly at: Address;
  readonly ours: readonly Op[];
  readonly theirs: readonly Op[];
  readonly issues?: readonly Issue[];      // the validators' issues, for `invalid`
}
type ConflictKind = 'set' | 'remove' | 'insert' | 'move' | 'replace' | 'order' | 'invalid';
// `move`: a node moved to two different parents, or renamed two ways; reorders within a parent are `order` (N-5).
// Reserved: merge3 works from diff, which infers no cross-parent move and no rename, so K2 never produces `move`
// and writes no branch or test for it; it is for a later merge of written patches.
```

- **Only `ok: true` has a `root`** (R2-B2): a caller cannot commit a merge and ignore its conflicts by accident.
- **`partial`** is `ours` with every change of theirs that is independent of every conflict. It is applied and
  validated as one patch, under the options. If that fails, `partial` is `ours`. A host shows `partial` with the
  conflicts, or rescues theirs as a patch to review.
- **Without a language** the merge is not validated, and it says so: an info in `issues` (N-7). Warnings kept under
  `allowDangling` are in `issues` too.
- **`changes`** on success is the net change from `ours` to `root` (patches.md, `Changes`), so `SpecModel` recompiles
  only what the merge changed.

## How

- `diff(base, ours)` and `diff(base, theirs)`, then each change is keyed by the node it touches (an address in `base`).
- **Changes that both apply**: to different nodes; to different keys of one node; equal changes, applied once.
- **Theirs is applied as one patch** on `ours`, under the options. Dangling refs, rebinding and the validators are
  judged on the merged root against `ours` (patches.md, Applying 6 and 7). This catches what no single op conflicts on:
  - a ref theirs wrote to a node only theirs inserted, where ours inserted another node with that id (R2-B2b);
  - two follows that make a cycle (R2-B2c).
- **Independence** is by node: a change of theirs depends on the changes of theirs to the nodes it names. That is its
  target, its parent, the targets of the refs it writes, and the sibling its `place` names when theirs made that
  sibling (N-6). A change whose dependency conflicts is left out with it, and listed in the same conflict's `theirs`.
- **Placement** (N-6): an insert or a move of theirs whose `place` names a sibling that no longer stands in `ours` goes
  after its nearest predecessor in theirs that still stands, else first. It never fails the merge.

## Conflicts, complete for the six ops (R2-B3)

| Kind | Ours | Theirs |
|---|---|---|
| `set` | `set` of a key | `set` of the same key to another value, or `unset` of it |
| `remove` | `remove` of a node | any change at or under it: `set`, `unset`, `insert` under it, `move` of it or into it, `replace` of it or under it |
| `insert` | `insert` of an id | `insert` of the same id under the same parent with other content |
| `move` | `move` of a node to another parent, or a `rename` | `move` of the same node to a different parent, or a different `rename` |
| `replace` | `replace` of a node | any change under it, or a different `replace` of it |
| `order` | a new order of a parent's children kept on both sides | a different new order of the same kept children |
| `invalid` | — | the merged root its validators refuse, with their issues |

- Each row is symmetric: ours and theirs may swap.
- **Orders are compared as orders, not as moves** (N-5), so the result never depends on which child `diff` chose to
  move. For each parent, take the children kept on both sides:
  - if only one side changed their order, it is taken;
  - if both changed it to the same order, it is taken once;
  - if both changed it differently, one `order` conflict at the parent, and `partial` keeps ours.
  
  Children only one side inserted are placed by the placement rule above.
- **`invalid`** (N-8): `at` is the address of the first refusing issue, or `[]` when it has none; `theirs` lists every
  change left out; `partial` is `ours`.

## Laws (property tests)

- `merge3(b, x, b)` is `{ ok: true, root: x }`;
- `merge3(b, b, y)` is `{ ok: true, root: y }`, when `y` is valid;
- `merge3(b, x, x)` is `{ ok: true, root: x }`;
- with no conflicts, the merge is symmetric up to the order of siblings that only one side inserted under the same
  parent (N-8);
- every `partial` is valid under the options;
- results are compared by hash; where nothing of theirs applied, `root` (or `partial`) is `ours`, the same object.

## Rebase

A stale patch `p` made on `base` is rebased as `merge3(base, current, applyPatch(base, p, options).root, options)`,
with the same options on both calls, so `p`'s moves rewrite refs as they would have (N-8). `base`'s root comes from
`SpecModel`'s recent roots, or from the host's `root(revision)`. It is the host's opt-in (patches.md, decision 4).

## For the author

1. **`partial` as ours plus theirs' independent changes** (recommended), or nothing but ours?
2. **Orders compared as orders** (recommended: one side's new order taken, two different ones a conflict), rather than
   merging moves child by child?
