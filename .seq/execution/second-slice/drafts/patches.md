# Revisions, patches, diff and merge in `@molstar/core/repr` (draft 2, for the author's review)

A draft for Gate 0 of the second slice. Once approved, it is agreed into `.seq/design/patches.md`; tree-core.md and
format.md then point at it, and their "five operations" become six. Lane K builds exactly the surface in
[The API](#the-api); anything past it stops the lane.

- **Draft 1** was reviewed independently (Fable 5.2, 2026-11-18). Its six blocking findings (B1–B6) and six should-fix
  ones (S1–S6) are resolved here, each cited where it lands.
- **Draft 2** also takes the author's direction: "should we give this "git semantics" ... might be more future proof
  with concurrent editing too. and we would later bring similar concepts to meso as well."

## Decisions for the author first

1. **Git's model** ([The model](#the-model-snapshots-are-the-truth)): snapshots are the truth and patches are derived;
   a change is a commit (a root with its parent revisions). Lane K builds the pure primitives `applyPatch`, `diff` and
   `merge3`. Branches, a store of revisions and sync wait. Recommendation: yes.
2. **Seq's `layer` gets required ids** (its one kind with optional ids). Then hiding a layer is a one-key `set`, not a
   `replace` of its whole track, and every child in a Seq spec can be addressed. Recommendation: yes, in lane K's
   companion change to `seqx-core/lang` (a language change, so its version moves).
3. **A patch that leaves a ref dangling is refused by default** (link integrity, the built-in validator). A caller may
   allow it for one patch (`allowDangling`), and the result then lists what dangles. Recommendation: yes. Why: a delete
   that silently breaks a view is worse than one that asks first, and an undo then never passes through a broken state.
4. **A stale patch is refused, never rebased automatically** (EC-81). A host may rebase through `merge3` on purpose.
   Recommendation: yes (configurable over fixed: refusal is the default, rebase is a host's opt-in).

## The model: snapshots are the truth

What git gets right, and what `repr` already is:

- **A revision is a snapshot**: an immutable root, named by its `hashNode` (32 hex digits; what `SpecModel` uses
  today). Structural sharing keeps every unchanged subtree, so a snapshot costs only the path to what changed.
- **A commit** is `{ root, parents: revision[], origin, label, id? }`: data any holder keeps. A host's history keeps
  commits; Seq\* keeps none (the history port: its checkpoint is a root, its restore sets one).
- **A patch is derived**: what `diff(a, b)` says takes `a` to `b`. It is the wire and agent form of a change, and an
  agent or a user writes one to make a change. Applying it gives a new snapshot; the snapshot, not the patch, is what
  is recorded. Why: a recorded root is exact by construction. Ref rewrites, normalization and anything else applying
  did are in it, and its patch from the parent is `diff(parent, root)` whenever one is wanted (B1, B5).
- **Merge is three-way**: `merge3(base, ours, theirs)` by node identity, the path of ids, never position. Rebase is a
  merge whose result has one parent. Conflicts are data, never a silent pick (§transaction-is-data, the rescue of
  §visible-tx-status).
- **Not git**: there are no content-addressed blobs of text lines; the unit is a node and its params. Ids are local to
  a parent, so a node that changes parents is a remove and an insert, unless a patch said `move`.

meso\* gets this when it moves onto `@molstar/core/repr` (agreed at the first slice's Gate 0, migration later).

## A patch

```ts
interface Patch {
  readonly ops: readonly Op[];
  /** The revision it was made against: refused on any other (EC-81), unless it is `result` (already applied). */
  readonly base?: string;
  /** The revision it gives, when the writer knows it; a mismatch after applying refuses (`repr.patch.result`). */
  readonly result?: string;
  /** Who wrote it: 'user', 'agent', 'host', 'remote:<peer>'; the history port's origin kinds, one vocabulary. */
  readonly origin?: string;
  readonly label?: string;
  /** The writer's id for this patch; deduplicating a resend by (origin, id) is the holder's, not `applyPatch`'s. */
  readonly tag?: string;
}
```

- **Order of checks** (B3): if `result` is given and is the root's revision, the patch was already applied. It returns
  `{ ok: true, changed: false }` with the root as it was and no inverse, and the host records no step. Otherwise, if
  `base` is given and is not the root's revision, the patch is refused as `repr.patch.stale`. A patch with neither
  `base` nor `result` applies to whatever is there. Option `requireBase` refuses that, for agent and remote origins
  (default true when `origin` is `agent` or `remote:*`).
- **A patch is also a tree** in `PatchLanguage`: kind `patch`, one child per op, embedded nodes as `any` (their
  language varies), `tag` a param, not a node id (S6). It prints, encodes and decodes through the existing codecs:
  `(patch :base 3f… (set :at /views/main/t1 :key height :value 40))`.

## The operations

Written out (B6). A target is an `Address`, the path of ids (`readonly string[]`; `[]` is the root). Printed, it is an
address word of ids only. A form or `/*` is never a target or a parent: what data makes is not a spec node. A bookmark
on a dynamic node is an `insert` of a ref whose `to` holds that form's address.

```ts
type Where = 'params' | 'custom';
type After = { readonly after: string | null } | { readonly last: true };  // a sibling's id, null for first, or last

interface SetOp     { readonly op: 'set';     readonly at: Address; readonly key: string; readonly value: Value; readonly where?: Where }
interface UnsetOp   { readonly op: 'unset';   readonly at: Address; readonly key: string; readonly where?: Where }
interface InsertOp  { readonly op: 'insert';  readonly parent: Address; readonly place: After; readonly node: Node }
interface RemoveOp  { readonly op: 'remove';  readonly at: Address }
interface MoveOp    { readonly op: 'move';    readonly at: Address; readonly parent: Address; readonly place: After; readonly rename?: string }
interface ReplaceOp { readonly op: 'replace'; readonly at: Address; readonly node: Node }
type Op = SetOp | UnsetOp | InsertOp | RemoveOp | MoveOp | ReplaceOp;
```

| Op | Rules |
|---|---|
| `set` | One key, `where` defaulting to `params`. A record param (`for_each`) is set whole. An equal value (`valuesEqual`) is a no-op that keeps the node object. `key` may not be reserved (`kind`, `__proto__`, `$…`): `repr.patch.key`. |
| `unset` | One key; an absent key is a no-op. |
| `insert` | The node has an id unique among its new siblings. `place` is after a sibling, first, or last (the common append). |
| `remove` | The node and its subtree. |
| `move` | To a parent and place; `rename` gives it a new id there. Not under itself (`repr.patch.cycle`), not the root. |
| `replace` | A whole node at `at`, with the same id; its kind may change. On the root, the whole spec. |

- **On the root** (`[]`): `set`, `unset` and `replace`; never `remove` or `move`.
- **Siblings without ids** (B2): under a parent that has any child without an id, `insert`, `remove` and `move` are
  refused (`repr.patch.unaddressable`); that parent changes by `replace`. Why: a position after a sibling cannot be
  restored among siblings that have no names. With decision 2, no Seq parent has such children.
- **Ownership**: `applyPatch` owns the columns in the values and nodes it is given, as the tree does (their NaNs made
  canonical in place).
- **The inverse** is each op's inverse, taken against the state just before that op, in reverse order (the review
  checked that this composes).

## Applying

`applyPatch(root, patch, options)`:

1. **The checks above** (result, then base).
2. **The ops in order**, by structural sharing: only the path from the root to each touched node is new. A missing
   target, a taken id, a malformed id, a `place` naming no sibling, a cycle, a reserved key or an unaddressable parent
   refuses the whole patch with its code at `['ops', i, …]` (S4). Nothing is applied, and the caller's root is
   untouched.
3. **Refs follow a move or a rename** (S1), when a language is given:
   - **What is walked**: every value the language types `ref`, inside lists, records, maps, unions (by the member the
     validator would choose), tables and data forms' `args`, in params only.
   - **Not rewritten**: `any` params, `custom`, and addresses inside strings (text, links). A ref with `in` (a bare id)
     is rewritten on a rename within its scope.
   - **What matches**: only the leading run of id segments of a ref's address can be a prefix of the moved address.
     Forms after it are kept (`/views/main/(section :of /data/p53)` follows `/data/p53`'s move inside the form value
     too). Refs inside the moved subtree to itself move with it.
   - **Cost**: one index of the spec's refs, built once per apply, only when the patch has a `move`.
   - **Rebinding** (B1): an `insert`, `move` or `rename` that makes a dangling ref resolve again is refused
     (`repr.patch.rebound`), unless `allowRebind` is set. Why: addressing.md says a ref is never silently re-bound.
     The principles' "an id is never reused among its siblings" is carried as a tension: Seq's ids are authors' names,
     and an author may reuse one on purpose.
   - Without a language, a move rewrites nothing, and the result says `refs: 'not-rewritten'` once, not an issue per
     move.
4. **Validation** (B4), when a language is given:
   - Validate the new root, **skipping every subtree that is the same object as in the old root**: structural sharing
     marks exactly what the patch made.
   - Add the context rules of each op: a `set`/`unset` the param and its required-ness; an `insert`/`replace` the new
     subtree, its parent rule and its siblings' ids; a `move` the parent rule and sibling ids; a `remove` nothing in the
     schema (no kind requires children).
   - **Link integrity** (decision 3): every ref whose target the patch removed or moved without rewriting is
     `repr.patch.dangling`, an error unless `allowDangling`.
   - **The language's own validators** (`options.validators`: functions of the old and new root and the touched
     addresses, returning issues) run last. Seq's checks (`seq.lang.*`, follow cycles) join here.
   - Any error refuses the patch. Errors in untouched subtrees cannot appear, since those are skipped. This is what
     "introduced" means, with no issue identity needed.
   - `repr/schema` gains the entry this needs: `validateWithin(language, root, { skip: (node) => boolean })` (in the
     agreed surface, not a lane's surprise).
5. **The result** (B5):

```ts
type PatchResult =
  | { readonly ok: true; readonly changed: boolean; readonly root: Node; readonly revision: string;
      readonly inverse?: Patch;              // absent when nothing changed
      readonly touched: readonly Address[];  // every node made new (ops and ref rewrites), deepest first
      readonly dangling?: readonly Address[]; readonly refs?: 'not-rewritten'; readonly issues: readonly Issue[] }
  | { readonly ok: false; readonly issues: readonly Issue[] };
```

**Bounds** (S3): at most `maxOps` ops a patch (default 10,000). Depth is checked per op: depth(parent) + height(node)
is at most `MaxDepth` (512). security.md's tree depth of 256 is the Seq spec's tighter bound, its validator's, not the
core's. The spec's size after edits is the host's and Seq's bound (lane Q), not the core's.

## Diff

`diff(a, b)` gives a patch with `base` = `a`'s revision and `result` = `b`'s (S2):

- **Merkle descent**: equal hashes are skipped.
- **Children are matched by id within a parent.** A kept id with a changed subtree is descended into. A kind change
  under the same id, duplicate ids among siblings, or any child without an id makes the owner a `replace` (S2, B2).
- **Moves within a parent**: the fewest that keep the longest run in order, found by a longest increasing subsequence.
  Cost: the changed parents' fanout, O(n log n) each.
- **Emission order**, so every `place` is valid when applied in order:
  1. removes;
  2. `b`'s children left to right, each one not in the kept run placed after its `b` predecessor (inserted or moved);
  3. then the descents into kept children.
- **Not inferred**: changes of parent and renames (ids are local) come out as remove and insert, and the ref rewrites
  a `move` would have made come out as `set`s. Law 1 still holds.
- **Params and custom data**: one `set` or `unset` a key.
- **The laws**, held by property tests without a language (with one, they hold for valid, normalized trees):
  1. `applyPatch(a, diff(a, b)).root` is `b` by hash;
  2. `diff(a, a)` has no ops;
  3. `applyPatch(applyPatch(a, p).root, inverse)` is `a` by hash.

## Merge

`merge3(base, ours, theirs)` gives `{ root, conflicts }`:

- **By node identity** (the path of ids), from `diff(base, ours)` and `diff(base, theirs)`. Changes to different nodes,
  or to different keys of one node, both apply. Equal changes apply once.
- **Conflicts**, as data, each with its address and both sides:
  - the same key set to different values;
  - a node edited on one side and removed on the other;
  - the same id inserted with different content;
  - a node moved to different places;
  - sibling orders that interleave differently;
  - a merged tree its validators refuse.
- **On a conflict**, `root` is ours with every change that does not conflict applied. The conflicts say what was left
  out, so a host can show them or rescue them (keep the losing side as a patch to review).
- **Rebase** of a stale patch `p` made on `base`: `merge3(base, current, apply(base, p))`. It is a host's opt-in
  (decision 4), never automatic.
- **Laws** (property tests):
  - `merge3(b, x, b)` is `x`;
  - `merge3(b, b, y)` is `y`;
  - `merge3(b, x, x)` is `x`;
  - with no conflicts, the merge is symmetric up to sibling order (ours first on ties).

## The API

The agreed surface for lane K: entries `@molstar/core/repr/patch`, with one addition to `@molstar/core/repr/schema`.

```ts
// repr/patch
type Op = SetOp | UnsetOp | InsertOp | RemoveOp | MoveOp | ReplaceOp;    // as written out above
interface Patch { … }                                                      // as above
type PatchResult = …;                                                      // as above
interface ApplyOptions {
  readonly language?: Language;
  readonly validators?: readonly PatchValidator[];
  readonly maxOps?: number;            // 10,000
  readonly requireBase?: boolean;      // default: origin is 'agent' or 'remote:*'
  readonly allowDangling?: boolean;    // false
  readonly allowRebind?: boolean;      // false
}
type PatchValidator = (before: Node, after: Node, touched: readonly Address[]) => readonly Issue[];
function applyPatch(root: Node, patch: Patch, options?: ApplyOptions): PatchResult;
function diff(a: Node, b: Node): Patch;
function merge3(base: Node, ours: Node, theirs: Node, options?: { language?: Language }): MergeResult;
interface MergeResult { readonly root: Node; readonly conflicts: readonly Conflict[] }
interface Conflict { readonly kind: ConflictKind; readonly at: Address; readonly ours?: Op; readonly theirs?: Op }
interface Commit { readonly root: Node; readonly parents: readonly string[]; readonly origin?: string; readonly label?: string; readonly tag?: string }
const PatchLanguage: Language;
function patchToNode(patch: Patch): Node;
function patchFromNode(node: Node): { readonly patch?: Patch; readonly issues: readonly Issue[] };

// repr/patch/conformance (as repr/codec's suite is)
const PatchCases: readonly PatchCase[];
function checkPatches(apply: (root: Node, patch: Patch) => { ok: boolean; root?: Node; issues: readonly Issue[] }): readonly string[];

// repr/schema
function validateWithin(language: Language, root: Node, options: { skip: (node: Node) => boolean } & ValidateOptions): Issue[];
```

- **Dropped from draft 1**:
  - `revision()`: it was `hashNode` under a second name;
  - `MoveOp.id`, now `rename`;
  - `in`, now `where`;
  - `Patch.id`, now `tag`;
  - `diff`'s unused `language` option.
- **`Commit` is a type only**: who keeps commits (a host's history, later a store) is outside the core.

**Issue codes** (`repr.patch.*`, cataloged, each at `['ops', i, …]`; S4):
- `stale`, `result`;
- `target`, `id`, `place`, `cycle`, `key`, `unaddressable`, `op`, `limit`, `depth`;
- `dangling`, `rebound`.

A validator's own issue is carried as is, with the op that caused it in its `path` (one carrier, no duplicate).

## Tests

- **The conformance cases as data**: every op and its inverse, at the root, in the middle and at the leaves; every
  refusal by code; a ref rewritten by a move and back by its inverse; a rebinding refused; a dangling ref refused and
  allowed.
- **The laws** of diff and merge as seeded property tests.
- **Structural sharing**: an untouched subtree is the same object after a patch.
- **A refusal leaves the root's hash the same.** That no history step is recorded is lane N's test, over the
  composition.
- **Two consumers** (lane N): a Seq\* spec edited by the user and by an agent through one host history. A stale agent
  patch is refused, then rebased by the host's opt-in through `merge3`.
- **Performance representative**: a 10-op patch on the p53 spec, applied and diffed back. A quick-lane case, with its
  wall-clock bound (under 1 ms) in `itWithin` (the heavy lane).

## What else changes with it

- **format.md and tree-core.md**: six operations (meso\*'s five and `replace`), and the model above. A note for meso\*'s
  readers comes with its migration.
- **`seqx-core/lang`**: `layer` ids required (decision 2), the language version moved, and the examples given ids.
- **changes.md**: a spec patch's row says the history records the root and its parents, and the patch is derived.
  Concurrent writers (EC-93) are no longer "a design of their own": `merge3` is their primitive, and their policy
  (who rebases, when) is the host's.
- **Citations**: the lane's code and `_spec` cite their own sections, never `.seq/` or a principle's id. The `repr`
  `_spec`'s existing citations of `.seq/design/*` are fixed when the scaffold is retired.

## For the author, beyond the four decisions

From draft 1, with the review's agreement:
1. One key a `set`.
2. Refs rewritten only with a language.
3. No cross-parent moves or renames inferred by `diff`.
4. Validation by shared-subtree skip plus per-op context rules.
5. `replace` for siblings without ids.
