# Revisions, patches and diff in `@molstar/core/repr`

Agreed by the author, 2026-11-18: git's model, `layer` ids required, dangling refs refused by default, no automatic
rebase, and draft 1's five questions ("yes as recommended"). Revised after a second independent review (Fable 5.2) of
this design with the state semantics, as the author agreed ("let's go with 1-5"), then after its focused re-review of
the interfaces ("apply N-1 to N-14 as recommended"). Findings are cited where they land: R2-B1 to R2-B4, R2-S1 to
R2-S16, R2-C1 to R2-C10, then N-1 to N-14.

- **Merge** has its own short agreement (lane K2, `merge.md`); this design keeps only its place in the model.
- **tree-core.md and format.md** point here.
- **Lane K** of the second slice builds exactly the surface in [The API](#the-api); anything past it stops the lane.

## Decisions

1. **Git's model** ([The model](#the-model-revisions-are-the-truth)): revisions, immutable roots, are the truth, and
   patches are derived. Lane K builds `applyPatch` and `diff`; lane K2 builds `merge3`. Branches, a store of revisions
   and sync wait.
2. **Seq's `layer` ids are required**, and a Seq check requires an id on every child, extension kinds included
   (R2-S11). Hiding a layer is then a one-key `set`, and every child in a Seq spec can be addressed.
3. **A patch that leaves a ref dangling is refused by default**, with `allowDangling` per patch. Why: a delete that
   silently breaks a view is worse than one that asks first.
4. **A stale patch is refused, never rebased automatically** (EC-81). A host may rebase through `merge3` on purpose.

## The model: revisions are the truth

- **A revision** is an immutable root, named by its `hashNode` (32 hex digits; what `SpecModel` uses today).
  Structural sharing keeps every unchanged subtree, so a revision costs only the path to what changed.
  - The word used here is *revision*.
  - *Snapshot* means a story step everywhere in these designs (composition.md's `snapshots`, changes.md,
    [stories.md](stories.md)), never a revision (R2-S6).
- **A commit** is a revision's root with its parents (one, or two for a merge). The core has no commit type: the
  holder of a history defines its own (N-14). Seq\*'s is the history port's checkpoint (`seqx-viewer`, changes.md), and
  meso\* gets one in the core when it needs it.
  - An undo restores an earlier root. It is a move of the pointer, not a new commit (R2-S5).
- **A patch is derived**: what `diff(a, b)` says takes `a` to `b`.
  - It is the form a change travels in, from an agent, the wire, or a stored story step.
  - Applying one gives a new revision, and the revision is what is recorded. Why: a recorded root is exact by
    construction, with ref rewrites and normalization in it. Its patch from the parent is `diff(parent, root)` whenever
    one is wanted.
- **Merge is three-way**, by node identity (`merge.md`). Rebase is a merge whose result has one parent. To rebase, a
  holder needs the base revision's root: `SpecModel` keeps the recent roots by revision (`revisions.keepRecent`,
  default 64), and a host may resolve older ones through the history port's `root(revision)` face (R2-S5, N-14).
- **Not git**:
  - The unit is a node and its params, not lines of text.
  - Ids are local to a parent, so a node that changes parents is a remove and an insert, unless a patch said `move`.
  - The principles' "an operation made against an address a concurrent move changed is rebased through the move" is
    not held: a stale edit under a renamed node becomes an edit-against-remove conflict in a merge. Carried as a
    tension, beside the reuse of ids (addressing.md).
- meso\* gets this when it moves onto `@molstar/core/repr`.

## A patch

```ts
interface Patch {
  readonly ops: readonly Op[];
  /** The revision it was made against: refused on any other (EC-81), unless the root is already `result`. */
  readonly base?: string;
  /** The revision it gives, normalized, when the writer knows it: a mismatch after applying refuses. */
  readonly result?: string;
}
```

- **Metadata** (origin, label, a writer's tag) is not in the core's patch (R2-C6, R2-S4). It belongs to whoever holds
  the change: a history step, a commit, an agent's call.
- **Order of checks**:
  1. If `result` is given and is the root's revision, the patch was already applied. It returns `{ ok: true, changed:
     false }`, with the caller's own root object and no inverse.
  2. Otherwise, if `base` is given and is not the root's revision, it is refused: `repr.patch.stale`.
  3. If `requireBase` is set and the patch has no `base`, it is refused: `repr.patch.base`. The core's default is
     false; `SpecModel` sets it for agent and remote channels (R2-S4).

  The hash is MurmurHash3, which is not collision-resistant. It names revisions and is not a security boundary
  (R2-C1).
- **A patch is also a tree**. `patchLanguage(language)` gives the patch language for a language (R2-S6, N-4):
  - kind `patch`, one child per op, laid out by `patchToNode` the same way for every language; a place as `:after id`,
    `:first true` or `:last true` (lane K);
  - `at` and `parent` as plain address strings (not `ref`: the stored-patch rule below walks them, F-2), embedded nodes
    as slots of that language;
  - one object per language (memoized), its `id` the language's id with `.patch` after it (`seq.patch`), a valid
    language id, so slots can name it (F-1).

  Over no particular language, `AnyPatchLanguage` (id `patch`, version `1.0.0`) has `any` nodes and string addresses. Both print, encode and decode
  through the existing codecs: `(patch :base 3f… (set :at /views/main/t1 :key height :value 40))`.
- **A stored patch** (a Seq story step) has no `base`: it applies to whatever its spec's revision is. Its refs are
  walked by `repr/patch`'s own rule, not by the patch language's param types (N-12):
  - **what is a ref**: `at` and `parent`; the refs inside embedded nodes, by the target language; and a `set` value
    when the target language types that key as a ref in any kind that has it, so "set `layer.data` to `/data/other`"
    is a ref. Typed by the language alone, never by reading the spec (the author, 2026-11-18, after lane K3's review:
    typing by the dry run changed under a patch, cost a retype of every step at every move and weakened a refusal;
    drift `stored-refs.md`). A key that is a ref in one kind and not in another is taken as a ref;
  - **what they resolve against**: the spec with the step's earlier ops applied (a dry run), so a step may insert a node
    and then set a key on it. The dry run decides what a ref names, never whether it is one;
  - so a stored patch follows renames and moves, and counts in link integrity, like any other ref: a patch that leaves
    a step's ref dangling, or rebinds it, is refused as for any ref. An op of the step that no longer applies (its
    target gone, its sibling renamed, its own move now refused) is reported on the step and shown as not applying,
    never a refusal of the patch that caused it ([stories.md](stories.md));
  - **not rewritten**: `place.after`, a bare sibling id. A step whose sibling was renamed no longer applies, and is
    shown as not applying (changes.md). Accepted: a story step names its place by a sibling, as a writer did.

## The operations

A target is an `Address` (`readonly string[]`; `[]` is the root). Printed, it is an address word of ids only. A form or
`/*` is never a target or a parent. A bookmark on a dynamic node is an `insert` of a ref whose `to` holds that form's
address.

```ts
type Where = 'params' | 'custom';
type Place = { readonly after: string | null } | { readonly last: true };  // a sibling's id, null for first, or last

interface SetOp     { readonly op: 'set';     readonly at: Address; readonly key: string; readonly value: Value; readonly where?: Where }
interface UnsetOp   { readonly op: 'unset';   readonly at: Address; readonly key: string; readonly where?: Where }
interface InsertOp  { readonly op: 'insert';  readonly parent: Address; readonly place: Place; readonly node: Node }
interface RemoveOp  { readonly op: 'remove';  readonly at: Address }
interface MoveOp    { readonly op: 'move';    readonly at: Address; readonly parent: Address; readonly place: Place; readonly rename?: string }
interface ReplaceOp { readonly op: 'replace'; readonly at: Address; readonly node: Node }
type Op = SetOp | UnsetOp | InsertOp | RemoveOp | MoveOp | ReplaceOp;
```

| Op | Rules |
|---|---|
| `set` | One key, `where` defaulting to `params`. A record param (`for_each`) is set whole. A value equal to the old one (`valuesEqual`) is a no-op that keeps the node object. A reserved `key` (`kind`, `__proto__`, `$…`) is refused: `repr.patch.key`. |
| `unset` | One key. An absent key is a no-op. |
| `insert` | The node has an id unique among its new siblings. `place`: after a sibling, first, or last. |
| `remove` | The node and its subtree. Never the root. |
| `move` | To a parent and place; `rename` gives it a new id there. Never under itself (`repr.patch.cycle`), never the root. |
| `replace` | A whole node at `at`, with the same id, its kind free to change; on the root, the whole spec. A node equal by hash to the old one is a no-op that keeps the old object (R2-C3). |

- **On the root** (`[]`): `set`, `unset` and `replace` only.
- **Siblings without ids**: under a parent with any child that has no id, `insert`, `remove` and `move` are refused
  (`repr.patch.unaddressable`), and that parent changes by `replace`. Why: a position after a sibling cannot be
  restored among siblings that have no names. With decision 2, a Seq spec has no such parent.
- **Ownership**: `applyPatch` owns the columns in the values and nodes it is given, as the tree does (their NaNs made
  canonical in place).
- **The inverse** is each op's inverse, taken against the state just before that op, in reverse order.

## Applying

`applyPatch(root, patch, options)`. The steps, in order:

1. **The checks** above: result, then base.
2. **The ops**, in order, by structural sharing: only the path from the root to each changed node is new. Each op is
   refused, with the whole patch, for:
   - a missing target;
   - a taken or malformed id;
   - a `place` that names no sibling;
   - a cycle;
   - a reserved key;
   - an unaddressable parent;
   - a depth past the bound.

   Nothing is applied, and the caller's root is untouched.
3. **Normalization** (R2-B1, N-1), when a language is given. Each inserted or replaced subtree is normalized as on
   read (`normalize` of that subtree), and each `set` value against its kind's param (`normalizeValue`, a new entry of
   `repr/schema`): plain arrays in declared columns become typed columns. Slots resolve through `options.validate`'s
   languages, as validation's do. This costs only what the patch made. Without a language, nothing is normalized.
4. **Refs follow a move or a rename**, when a language is given:
   - **What is walked**: every value the language types `ref`, inside lists, records, maps, unions (by the member the
     validator would choose), tables, data forms' `args` and slots, in params only. The body of a stored patch's slot
     is skipped here: only the stored-patch rule ([A patch](#a-patch)) walks it, so nothing in it is rewritten twice,
     and an `at` naming a node the step itself inserts is not taken for a dangling ref (F-2).
   - **Not rewritten**: `any` params, `custom`, and addresses inside strings (text, links). Those go stale, and a
     language's validator may report them.
   - **What matches**: only the leading run of id segments of a ref's address can be a prefix of the moved address.
     Forms after it are kept, and addresses nested in form values are followed too. Refs inside the moved subtree to
     itself move with it.
   - **How refs are found**: through the revision's ref index, an internal seam written out now (N-2), so swapping its
     implementation changes nothing else:

     ```ts
     interface RefSite { readonly at?: Address; readonly path: readonly (string | number)[]; readonly value: string }
     interface RefIndex {
       under(prefix: Address): readonly RefSite[];         // refs whose id run names prefix or below
       dangling(): readonly RefSite[];                     // refs that name no node
       derive(after: Node): RefIndex;                      // the next revision's, by node identity and position
     }
     function refIndexOf(root: Node, language: Language, options?: ValidateOptions): RefIndex;  // cached per root and language
     ```

     It keeps its sites in chunks by node object, on plain `Map`s. A revision's index is derived lazily: `applyPatch`
     leaves a pending derive on the new root, pointing at the nearest built index, and the next apply that needs the
     index forces it. A node that is the same object at the same position keeps its chunk; a node whose sites read the
     root (a stored patch's typed `set`) is always walked again. The cache hits only for the same `validate` object,
     held weakly. Measured (lane R+K, 10,000 refs): a derive about 0.3 ms, a query 0.2 ms, so the persistent trie
     (immutable.md, In `repr`) is not needed. It is built only
     for ops that change an address: a move to another parent, a rename, an insert, a remove. A reorder within a parent
     does not.
   - **Not walked**: a ref with `in` (a bare id): its scope is the language's, and Seq has none.
   - Without a language, nothing is rewritten, and the result says `refs: 'not-rewritten'` once.
5. **Validation** (R2-B4), when a language is given. One rule, by how each new object came to be:

   | A new object that is | Is checked for |
   |---|---|
   | made by an `insert` or a `replace` | its whole subtree, nothing skipped inside it; old children kept under a node whose kind changed included |
   | a node whose own params or custom data changed | the changed keys, and that no key the patch unset is required by its kind (N-3) |
   | a node whose child list changed | its changed children's parent rule, and the uniqueness of its children's ids |
   | only a copy on the path to a change | nothing |

   Errors already present elsewhere, or on a path copy's own params, are never re-reported, so only what the patch
   introduced can refuse it. `repr/schema` gains the entry that checks exactly this: `validateChanges(language, root,
   changes, options)`, where `changes` is the result's `Changes`. `Changes` is declared in `repr/schema` and re-exported
   by `repr/patch`, so imports keep going `tree ← schema ← patch` (N-2).
6. **Link integrity** (R2-S1), judged on the old root against the new, never op by op, when a language is given. The
   core's own pass over refs (what a ref names: the node at its id run) is built by lane K.
   - **Dangling**: a ref in the new root that names no node, and either named one in the old root or was written by
     this patch. It is an error, `repr.patch.dangling`; with `allowDangling` it is a warning instead, kept in `issues`
     with its path and address (R2-S8).
     - **"Written by this patch"** (K's review, 2026-11-18): a ref the patch's ops wrote, a `set` value or a ref in an
       inserted or replaced subtree, **except** one that was already in the old root with the same value at the same
       place and already named no node there. Such a ref was not introduced: a `replace` that keeps it, or an inverse
       that restores removed content holding it, is not refused for it. So undoing a removal never fails on refs that
       dangled before the removal.
   - **Rebound**: a ref the patch did not change whose target changed identity. Either it named no node in the old root
     and names one in the new, or the node it named was removed or moved away by this patch and its address now names
     another node (remove `/b`, then rename `/a` to `b`: a ref to `/b/y` would now name old `a`'s `y`). It is an
     error, `repr.patch.rebound`, unless `allowRebind`. The hint offers both ways out: give the ref another target, or
     allow the rebinding (R2-C5; the identity case added after K's review, 2026-11-18: "never silently re-bound").
   - `diff`'s remove-then-insert of a node that changed parents passes, since only the end state is judged.
7. **The language's validators** (`options.validators`), last: `(before, after, changes) => issues`. Each costs in
   proportion to `changes`, and reports only what the changes introduced (R2-S9). Seq's checks join here, so `SpecModel`
   needs an incremental `checkSeqLang` (lane N).
8. **The result check**: if `result` was given and is not the new root's revision (normalized), it is refused,
   `repr.patch.result`. A writer that cannot normalize leaves `result` out.

**The result** (R2-S8). `Changes` is the **net** difference from the old root to the new (N-3):
- one entry a node, whatever sequence of ops reached it;
- `moved.from` in the old root, `moved.to` in the new;
- a node in `made` is not also in `params` or `children`;
- a node removed after other ops changed it is only in `removed`;
- a reorder within a parent appears in `children` only.



```ts
interface Changes {
  readonly made: readonly Address[];                                      // inserted or replaced, in the new root
  readonly params: readonly { readonly at: Address; readonly keys: readonly string[]; readonly where: Where }[];
  readonly children: readonly Address[];                                  // parents whose child list changed, in the new root
  readonly removed: readonly Address[];                                   // in the old root
  readonly moved: readonly { readonly from: Address; readonly to: Address }[];
}
type PatchResult =
  | { readonly ok: true; readonly changed: boolean; readonly root: Node; readonly revision: string;
      readonly inverse?: Patch;                    // absent when nothing changed
      readonly inverseNeeds?: { readonly allowDangling?: true; readonly allowRebind?: true };
      readonly changes: Changes;                   // ref rewrites included, as `params`
      readonly refs?: 'not-rewritten';
      readonly issues: readonly Issue[] }          // warnings and infos
  | { readonly ok: false; readonly issues: readonly Issue[] };
```

- **The inverse and the allowances** (R2-S2):
  - A patch applied with `allowDangling` has an inverse that makes those refs resolve again. `inverseNeeds` says that
    inverse needs `allowRebind`, and the same holds the other way.
  - A patch that removes content holding refs already dangling has an inverse that needs `allowDangling`
    (`inverseNeeds` says so): the inverse's insert makes that content again, and a patch carries nothing to tell it
    from new content. The allowance only ever lets through refs that dangled in the forward patch's old root, so law 3
    holds (lane K; confirmed by the author, 2026-11-18). A replace that keeps such a ref at the same
    place needs nothing: the ref is not "written" (above).
  - `inverseNeeds` is conservative: an inverse of a patch that rebound refs needs both allowances, and one that makes
    again content some ref of the old root named needs `allowRebind`. Precise needs would mean link integrity run on
    the inverse.
  - **A replace is an edit in place**: the node keeps its identity, and so do its children at ids the old node had (an
    inserted node has none). A ref into a node whose kind a diff changed is therefore not rebound.
  - **Every ref the inverse's moves would take along without the forward patch having rewritten it** gets an explicit
    restoring `set` in the inverse (K's review: not only refs that were rebound). So a move back never drags a ref the
    move forward did not touch: a ref under the move's destination address that dangled before and after, or one that
    named another node there.
  - Law 3 holds when the inverse is applied with `inverseNeeds`. A host's undo restores roots and needs none of this.
- **Issue paths** (R2-S9): an issue about the spec has the spec's path and address. One about an op has the patch
  tree's path (`['children', i, …]`, as `patchToNode` lays it out).

**Bounds**:
- at most `maxOps` ops a written patch (default 10,000). A patch from `diff` is applied with no op bound, so law 1
  holds on any diff (R2-C2);
- depth(parent) + height(node) at most `MaxDepth` (512), and a value's own nesting counted too;
- security.md's tree depth of 256 is the Seq spec's tighter bound, its validator's;
- the spec's size after edits is the host's and Seq's (lane Q).

## Diff

`diff(a, b)` gives a patch with `base` = `a`'s revision and `result` = `b`'s:

- **Merkle descent**: equal hashes are skipped.
- **Children are matched by id within a parent.** A kept id with a changed subtree is descended into. The fallbacks
  (R2-S3):
  - a kind change under the same id: a `replace` of that node;
  - duplicate ids among siblings, or any child without an id: a `replace` of the nearest ancestor that has an address,
    up to the root;
  - a root that differs in its own id or kind: a `replace` of the root.
- **Moves within a parent**: the fewest that keep the longest run in order, by a longest increasing subsequence. Cost:
  each changed parent's fanout, O(n log n).
- **Emission order**, so every `place` is valid when applied in order:
  1. removes;
  2. `b`'s children left to right, each one not in the kept run placed after its `b` predecessor (inserted or moved);
  3. the descents into kept children.
- **Not inferred**: changes of parent and renames come out as remove and insert, and the ref rewrites a `move` would
  have made come out as `set`s.
- **Params and custom data**: one `set` or `unset` a key.
- **The laws**, held by seeded property tests without a language, and with one on valid, normalized trees applied with
  `allowDangling` and `allowRebind` (two trees may differ in which refs dangle; N-4):
  1. `applyPatch(a, diff(a, b)).root` is `b` by hash;
  2. `diff(a, a)` has no ops;
  3. applying `p`, then its inverse with `inverseNeeds`, gives `a` by hash.

## Previews

A preview is a patch over the current root (R2-S7):

- **What is shown** is `applyPatch(root, preview)`. When the root changes (an agent's patch lands mid-drag), the
  preview is applied again on the new root. If it no longer applies, it is dropped with an issue, and the gesture ends.
- **At interactive rates** the preview applies with the language's normalization and validation but not the language's
  own validators (`options.validators`). Committing it applies it once more, in full.
- `viewer-models.md` holds the rest.

## The API

The agreed surface for lane K: entry `@molstar/core/repr/patch`, with one addition to `@molstar/core/repr/schema`.

```ts
// repr/patch
type Op = SetOp | UnsetOp | InsertOp | RemoveOp | MoveOp | ReplaceOp;    // as written out above
interface Patch { … }  interface Changes { … }  type PatchResult = …;      // as above
interface ApplyOptions {
  readonly language?: Language;
  readonly validate?: ValidateOptions;   // the slots' languages, the units hook: as validation and normalization take them (N-1)
  readonly validators?: readonly PatchValidator[];
  readonly maxOps?: number;            // 10,000; Infinity for a diff
  readonly requireBase?: boolean;      // false
  readonly allowDangling?: boolean;    // false
  readonly allowRebind?: boolean;      // false
}
type PatchValidator = (before: Node, after: Node, changes: Changes) => readonly Issue[];
function applyPatch(root: Node, patch: Patch, options?: ApplyOptions): PatchResult;
function diff(a: Node, b: Node): Patch;
function patchLanguage(language: Language): Language;   // memoized; id `<language id>.patch`
const AnyPatchLanguage: Language;                      // over no particular language
function patchToNode(patch: Patch): Node;
function patchFromNode(node: Node): { readonly patch?: Patch; readonly issues: readonly Issue[] };

// repr/patch/conformance (as repr/codec's suite is)
interface PatchCase {
  readonly name: string; readonly root: Node; readonly patch: Patch; readonly options?: ApplyOptions;
  readonly expect: { readonly ok: boolean; readonly root?: Node; readonly codes?: readonly string[] };   // codes: the exact set
}
const PatchCases: readonly PatchCase[];
function checkPatches(apply: (root: Node, patch: Patch, options?: ApplyOptions) => PatchResult,
                      cases?: readonly PatchCase[]): readonly string[];   // default PatchCases; Seq passes its own (N-4)

// repr/schema
interface Changes { … }                                                  // as above; re-exported by repr/patch
function validateChanges(language: Language, root: Node, changes: Changes, options?: ValidateOptions): Issue[];
function normalizeValue(type: ParamType, value: Value, options?: ValidateOptions): Value;
```

**Issue codes** (`repr.patch.*`, cataloged):
- the checks: `stale`, `base`, `result`;
- the ops: `target`, `id`, `place`, `cycle`, `key`, `unaddressable`, `op`, `limit`, `depth`;
- link integrity: `dangling`, `rebound`.

A validator's own issue is carried as is: one carrier, no duplicate.

## Tests

- **The conformance cases as data**, each with its options; `checkPatches` also runs Seq's own cases over
  `SeqLanguage` (lane N):
  - every op and its inverse, at the root, in the middle and at the leaves;
  - every refusal by code;
  - normalization of a plain array set into a column;
  - a ref rewritten by a move and back by its inverse;
  - a rebinding refused and allowed;
  - a dangling ref refused and allowed, with `inverseNeeds`;
  - a remove-then-insert from `diff` passing link integrity;
  - `Changes` net after insert-then-move, set-then-remove and a swap by three renames;
  - a stored patch that inserts a node and sets a key on it; a stored `set` of a ref value rewritten by a rename.
- **The laws** of diff as seeded property tests.
- **Structural sharing**: an untouched subtree is the same object after a patch; an equal `set` or `replace` keeps the
  node.
- **A refusal leaves the root's hash the same.**
- **Validation**: an error present before an edit is not re-reported by an edit below it; a kind-changing `replace`
  validates its kept children.
- **Lane N's** (over the composition): no history step for a refusal; two consumers (the user and an agent) on one
  host history; a stale agent patch refused, then rebased by the host's opt-in through `merge3` with the base root from
  `SpecModel`'s recent roots; a preview surviving an agent's patch mid-drag.
- **Performance representative**: a 10-op patch on the p53 spec, applied and diffed back. A quick-lane case, with its
  wall-clock bound (under 1 ms) in `itWithin` (the heavy lane).

## What else changes with it

- **format.md and tree-core.md**: six operations (meso\*'s five and `replace`), and the model above.
- **`seqx-core/lang`**: `layer` ids required, a check that every child of the spec has an id (`seq.lang.id.required`;
  not the bodies of slots, whose nodes, a stored patch's ops, have none), the examples and the builder given ids. The language version moves with no migration: the language is pre-release, `0.x`
  (R2-S11).
- **changes.md**: revisions and commits, story steps as stored patches, edits while a step is shown, previews as
  patches.
- **viewer-models.md**: `SpecModel`'s checkpoints are commits, it keeps the recent roots, applies previews as patches,
  and composes a shown story step over them (changes.md).
- **Citations**: the lane's code and `_spec` cite their own sections, never `.scaffolds/seq/` or a principle's id.

## In the specs

| Decision | Where |
|---|---|
| Git's model: revisions the truth, patches derived; no commit type or undo stack in the core | repr › A revision |
| Branches, a store of revisions, sync | design-only: not built (waiting) |
| `layer` ids required; every child of the spec has an id | core › Every child of the spec has an id |
| Dangling refs refused by default, `allowDangling` | repr › Link integrity |
| A stale patch refused, never rebased automatically | repr › A patch |
| The host's rebase through `merge3`, opt-in | viewer › A stale patch |
| Commits are `SpecModel`'s; an undo moves the pointer | viewer › SpecModel |
| The recent roots (`revisions.keepRecent`) and `root(revision)` | viewer › Recent roots |
| *Snapshot* means a story step only | core › Story steps |
| Not git: a stale edit under a renamed node is a conflict, not rebased through the move | design-only: a tension (with the reuse of ids) |
| The patch's shape, no metadata, the order of checks; the hash no security boundary | repr › A patch |
| `requireBase` for agent and remote channels | viewer › SpecModel |
| The patch as a tree, `patchLanguage`, `AnyPatchLanguage` | repr › The patch as a tree |
| A stored patch: what is a ref (by the language), what it resolves against (its dry run), `place.after` not rewritten | repr › A stored patch |
| The six ops, on the root, siblings without ids | repr › Six ops |
| `applyPatch` owns the columns it is given | repr › The tree owns a column it is given |
| The inverse | repr › The inverse |
| All or nothing; the bounds (`maxOps`, `MaxDepth`) | repr › All or nothing |
| The Seq spec's tighter depth (256) | core › Language › Bounds |
| Normalization of what the patch made | repr › Normalization |
| Refs follow a move or a rename; what is walked, what is not | repr › Refs follow a move or a rename |
| The ref index, derived lazily, measured | repr › The ref index |
| Validation of what the patch introduced; `validateChanges` | repr › Validation |
| Link integrity: dangling, "written by this patch", rebound, identity by lineage | repr › Link integrity |
| The language's validators; Seq's incremental checks | core › Checks of a change |
| The result check | repr › The result check |
| `Changes` net | repr › Changes |
| The inverse's allowances (`inverseNeeds`) | repr › The allowances turned around |
| The inverse's restoring sets | repr › The inverse's restoring sets |
| Issue paths and codes | repr › Issue paths |
| `diff`: Merkle descent, matching by id, fallbacks, moves by a longest increasing subsequence, emission order | repr › diff(a, b) |
| The three laws as seeded property tests | repr › Tests |
| Previews: re-applied on a new root, without the language's own validators, committed in full | viewer › Previews are patches |
| The API and the conformance suite | repr › Patches › The conformance suite |
| Seq's conformance cases and the 10-op p53 bound | core › Patches |
| Lane N's tests over the composition | viewer › SpecModel |
| What else changes (format, changes, viewer-models, citations) | design-only: pointers to other designs |
