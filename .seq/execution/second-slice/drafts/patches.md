# Patches and diff in `@molstar/core/repr` (draft, for the author's review)

A draft for Gate 0 of the second slice. Once approved, it is agreed into `.seq/design/patches.md`, and tree-core.md's
`repr/patch` line points at it. Lane K builds exactly the surface in [The API](#the-api); anything past it stops the
lane.

## Intent

Every change to a spec is a patch: the UI, commands, agents and a host alike (meso\*; changes.md). A patch is data. It
can be validated, applied, inverted, printed, logged, sent and diffed. An undo step is one patch batch, though the undo
stack itself is the host's (Seq\*'s history port).

## Precedents

- **meso\***: `set`, `unset`, `insert`, `remove` and `move`, each naming one node by address. A batch is all or nothing.
  Inverses are computed when a patch is applied.
- **The principles** (`.seq/oz/representations`):
  - operations address by id, never by position (§ops-address-by-id);
  - an operation carries the base hash it was made against, its origin and an id (§idempotent-ops);
  - a transaction is data (§transaction-is-data);
  - validators run on every commit, and only the issues a patch introduces refuse it (§validators).
- **This repository**: `repr`'s node, hash and address (`packages/core/src/repr/_spec`); addressing.md's refs, which
  follow a move; changes.md's change kinds and EC-81 (a patch made against another revision is refused).

## Decision

### A patch

```ts
interface Patch {
  readonly ops: readonly Op[];   // applied in order, all or nothing
  readonly base?: string;        // the revision it was made against; refused on any other (EC-81)
  readonly result?: string;      // the revision it gives; re-applying where the root already is `result` is a no-op
  readonly id?: string;          // the writer's own id (idempotence across a resend)
  readonly origin?: string;      // 'ui', 'agent', 'host', 'remote:<peer>': who wrote it
  readonly label?: string;       // for a history step
}
```

- **A revision** is the root's `hashNode` (32 hex digits), as the viewer's `SpecModel` already uses it. Why: the hash is
  already the spec's identity and is cached per node, so a revision costs nothing.
- **A patch is also a tree** in a small language, `PatchLanguage` (kind `patch`, one child per op), so it prints,
  encodes and decodes through the existing codecs: `(patch :base 3f… (set :at /views/main/t1 :key height :value 40))`.
  Why: §transaction-is-data, and no second encoding.

### The operations

Each op names spec nodes by **addresses of ids only**. A form is no target: what data makes (an item, a row, a dynamic
node) is not a spec node. A bookmark on a dynamic node is an `insert` of a ref whose `to` param holds that form's
address.

| Op | Fields | Does | Its inverse |
|---|---|---|---|
| `set` | `at`, `key`, `value`, `in?` (`params` default, or `custom`) | sets one key of a node's params or custom data | `set` with the old value, or `unset` if the key was absent |
| `unset` | `at`, `key`, `in?` | removes one key | `set` with the old value |
| `insert` | `parent`, `after` (a sibling's id, or `null` for first), `node` | inserts a node with an id unique among its new siblings | `remove` |
| `remove` | `at` | removes a node and its subtree | `insert` of the same subtree at the same place |
| `move` | `at`, `parent`, `after`, `id?` (a new id: a rename) | moves a node, with refs to it rewritten | `move` back, with refs rewritten back |
| `replace` | `at`, `node` | replaces a whole node, keeping its id and place | `replace` with the old node |

- **One key a `set`**: the inverse is exact, and two writers setting different keys commute. Why not a params merge: an
  inverse of a merge needs the old values of every key it touched anyway.
- **`replace` is for what has no ids**: children without ids (data forms in `args`, id-less children) cannot be
  addressed one by one, so they change with their nearest node that has an id.
- **`after`, not an index**: a position by sibling id survives other inserts (§ops-address-by-id).

### Applying

`applyPatch(root, patch, options)`:

1. **Base**: if `base` is given and is not the root's revision, the patch is refused with `repr.patch.stale` (EC-81).
   If `result` is given and is the root's revision, the patch was already applied: no change, no issue.
2. **Ops in order**, each on the result of the last, by structural sharing: only the path from the root to each touched
   node is new, and every other subtree is the same object with its cached hash. A missing target, an id taken, a
   malformed id, an `after` that is no sibling, or an op that moves a node under itself is refused with its code.
   Nothing is applied.
3. **Refs follow a move** when a language is given. Every `ref` param, and every address held where the language types
   an address, whose address has the moved node's address as a prefix is rewritten to the new prefix. Forms after the
   prefix are kept. The rewrites are part of the op, and so of its inverse. Without a language a move rewrites nothing
   and reports `repr.patch.refs_unknown` (info). A ref into a removed subtree is not touched: resolve reports it as
   dangling (addressing.md: reported, never silently re-bound).
4. **Validation** when a language is given: the touched nodes, their parents' rules (ids, parents, slots) and their
   subtrees are validated, normalized as on read. An error the patch introduces refuses it (§validators). Errors already
   present elsewhere do not.
5. **The result**: `{ ok: true, root, revision, inverse, issues }`, where `inverse` has `base` set to the new revision
   and `result` set to the old one. A refusal is `{ ok: false, issues }`, and the caller's root is untouched.

**Bounds**: at most `maxOps` ops a patch (default 10,000) and the security bounds on what a node holds (security.md).
Past either, `repr.patch.limit`.

### Diff

`diff(a, b)` gives a patch from `a` to `b`, with `base` = `a`'s revision and `result` = `b`'s:

- **Merkle descent**: equal hashes are skipped, so the cost follows what changed, not the spec's size.
- **Children are matched by id within a parent.** An id in both becomes the change of its subtree. Only in `b`, an
  `insert`; only in `a`, a `remove`. A reorder becomes `move`s, as few as the longest run kept in order allows.
- **Not inferred**: a node that changed parents is a `remove` and an `insert` (ids are local, so the same id under two
  parents may be two things). Nor are renames: a remove and an insert.
- **Id-less children and `args`** that differ become a `replace` of their owner.
- **Params and custom data**: a `set` or an `unset` a key that differs.
- Laws, held by property tests:
  - `apply(a, diff(a, b))` is `b` by hash;
  - `diff(a, a)` is empty;
  - `apply(apply(a, p).root, inverse)` is `a` by hash.

### Not in this design

- **Rebase and concurrent writers** (EC-93): a stale patch is refused, not rebased.
- **Transactions across documents**, presence, and composing patches into one. A history step groups patches; it does
  not merge them.

## The API

The agreed surface for lane K, entry `@molstar/core/repr/patch`:

```ts
type Op = SetOp | UnsetOp | InsertOp | RemoveOp | MoveOp | ReplaceOp;   // as the table above
interface Patch { … }                                                     // as above
type PatchResult =
  | { readonly ok: true; readonly root: Node; readonly revision: string; readonly inverse: Patch; readonly issues: readonly Issue[] }
  | { readonly ok: false; readonly issues: readonly Issue[] };

function revision(root: Node): string;
function applyPatch(root: Node, patch: Patch, options?: { language?: Language; maxOps?: number }): PatchResult;
function diff(a: Node, b: Node, options?: { language?: Language }): Patch;
const PatchLanguage: Language;                       // a patch as a tree, through the codecs
function patchToNode(patch: Patch): Node;
function patchFromNode(node: Node): { patch?: Patch; issues: readonly Issue[] };
const PatchCases: readonly PatchCase[];              // the conformance cases, as data
function checkPatches(apply: typeof applyPatch): readonly string[];   // failures, empty when all hold
```

Issue codes (`repr.patch.*`, cataloged): `stale`, `target`, `id`, `after`, `cycle`, `op`, `limit`, `invalid` (a
validation error the patch introduced, with the validator's issue beside it), `refs_unknown`.

## Tests

- The conformance cases as data (§suites-are-data): every op and its inverse, on ids at the root, in the middle and at
  the leaves; every refusal by code; a ref rewritten by a move and back by its inverse.
- Property tests on seeded trees: the three laws above.
- **Structural sharing**: an untouched subtree is the same object after a patch.
- **A refusal changes nothing**: the root's hash is the same, and no history step is recorded.
- **Two consumers**: a Seq\* spec patched by the UI and by an agent through one history (lane N's test). A stale
  agent's patch is refused, and the UI's step stays.
- **Performance representative**: apply a 10-op patch to the p53 spec and diff it back, under 1 ms. The edit target (a
  patch, resolve and redraw under 16 ms) is lane N's.

## For the author

1. **One key a `set`** (recommended) or a params merge?
2. **Refs rewritten on a move only with a language** (recommended: without one, a reader cannot tell a ref from a
   string), with an info issue otherwise?
3. **`diff` infers no cross-parent moves and no renames** (recommended: ids are local), so they come out as remove and
   insert?
4. **Validation of the touched nodes, their parents' rules and their subtrees** (recommended), not the whole spec?
   Refs across the spec stay resolve's pass.
5. **`replace` as the way to change id-less children** (recommended), rather than giving data forms ids?
