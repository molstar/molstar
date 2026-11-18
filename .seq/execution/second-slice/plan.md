# Plan: the second slice

Work id `second-slice`. Version 1's first milestone, **core editing and addressing** ([version-1](../../design/version-1.md),
Milestones after the slice), plus the selection the first slice left out and what the first slice carried over
(`packages/seq/TODO.md`, Next). Designs: [changes](../../design/changes.md), [addressing](../../design/addressing.md),
[tree-core](../../design/tree-core.md), [viewer-models](../../design/viewer-models.md) (Hover and selection),
[security](../../design/security.md), [performance](../../design/performance.md), and those each stage names.

Starts once the first slice is closed (`.seq/records/first-slice/`), on `packages-seq`.

## Done when

- A spec is edited the way changes.md's table says:
  - patches validated and applied;
  - undo and redo through `@molstar/ui/history`, one step a batch;
  - a preview that ends in one patch or nothing;
  - snapshots entered and left;
  - refs and bookmarks, on static and dynamic nodes;
  - a change made against another revision refused (EC-81).

  Each one shown in a story and in tests.
- Everything in a view has an address. Addresses resolve both ways (to a node, to a Location), and links and actions
  act through commands on them ([addressing](../../design/addressing.md)).
- Selection, highlights and the cursor are built:
  - the models, the click and drag-select bindings in the default preset, the keyboard cursor;
  - the selection and the hover projected onto linked sections and drawn there;
  - the selection tint's place decided.
- Every bound in security.md's table is enforced with its coded issue, and the asset manager change has landed if the
  author agrees.
- The MessagePack codec passes the codec suite, and a spec round-trips in both encodings.
- Measured against performance.md:
  - an edit and its undo;
  - a selection step;
  - MessagePack decode;
  - **browser frame times and first paint** (a visible pane, Chrome);
  - the layout's allocation bound held.

  Where an edit path has no target, performance.md gets one.
- The author has used the stories (§used-before-done). Designs, `_spec`s, `status.md`, the features rows and the edge
  cases' outcomes say what was built, and the record is written.

## How this slice runs (the first slice's fold-back, applied)

- **The lane's review includes the look**: before review, the orchestrator opens the lane's stories in a **visible**
  pane, in light and dark, and hands the reviewer the screenshots. A hidden pane runs no frames, so it is never
  evidence.
- **A fix lands with a test shown failing on its parent**: the report names the command, and the reviewer reruns it.
- **The merge message names the review's refuted or unfixed findings**, and why.
- **Each lane keeps its package `_spec`'s Configuration table** (item, default, why, alternatives). gaps.md lists only
  the defaults the lane was unsure of.
- **Shared things are tested with two consumers live** (the viewer and an export, two views on one store, an agent and
  the user on one history).
- **A path's performance target gets a quick representative** in the stage that builds it, not at the end.

## Gate 0: before any lane

- [x] GATE 0 (2026-11-18): "rest as recommended", with one amendment: the asset manager. Recorded in `state.md`,
      Decided. The text below is what was asked.
  - **The first slice's open decisions** (its record, Open decisions; recommendation: confirm as built, except 10):
    1. undeclared typed columns tagged;
    2. a cyclic range with start > end wraps, an error elsewhere;
    3. the query decisions (`chain` as Nodes, quoted `.n` labels, the boundary rule, structural Location equality);
    4. the protocol and store details (malformed not retried, `list`, row identity, per-consumer want sets,
       cancellation, `AssetReader`);
    5. the `seq.lang.*` codes and alphabets;
    6. every lane's configurable defaults (render, viewer, canvas, the asset policy, summaries, nesting, templates):
       confirmed as a whole, any one moved later through its config item;
    7. the theme read from the UI root (`theme.source`);
    8. the selection tint drawn over the data (lane P builds on it);
    9. letters' density in CSS px, the other steps in device px;
    10. `values.coarser`'s default: recommendation **1** (the direct path makes level 1 cheap, and 1 keeps the picture
        exact), against the 2 decided at stage 7;
    11. the asset manager change: below, with the base library;
    12. the TP53 exon coordinates in the test geometry: checked against RefSeq in stage 0 (S0), or marked unverified.
  - **Public API changes in the base library** ([base-library](../../design/base-library.md)):
    - **Decided** (2026-11-18): `@molstar/core/repr`'s revisions, patches and `diff`, git's model: `design/patches.md`,
      agreed with its four decisions, revised after three independent reviews (the last: sound to build on). Its API
      section is the agreed surface for lane K. Nothing to confirm here.
    - **`@molstar/core/data/immutable`** (the author, 2026-11-18: "implement core/immutable and use it in repr"; moving
      Mol\*'s state tree and model off immutable.js is a later, separate PR): the draft `drafts/immutable.md` (draft 2,
      after the second review), its API section the agreed surface for lane R, with its four decisions (`repr` uses it
      behind `RefIndex` once measured; CHAMP; our own API close to the immutable.js subset Mol\* uses, checked by a
      `StateTree` spike; configurable per collection for hash and equality). Recommendation: as drafted.
    - **`merge3`** (lane K2's short agreement, split from the patch design): the draft `drafts/merge.md`, with its two
      questions. Recommendation: as drafted.
    - **The asset manager**: `credentials: 'omit'`, `redirect: 'error'`, the signal passed through, a streamed byte
      limit. Recommendation: yes, in this slice, with the policy then claiming all four.
    - **Range reads in the asset manager**: recommendation: wait for bigWig (milestone 3).
    - **The MessagePack codec**: a `repr` codec in io on `common/msgpack`'s `encodeMsgPack` and `decodeMsgPack` and
      their extension hooks, with no change to io's MessagePack itself (tree-core.md, base-library.md). Typed columns
      decode without a copy where the hooks allow it, else with one copy, measured. Recommendation: yes (io depends on
      core; `repr` stays dependency-free).
  - **Decided** (2026-11-18, "yes, go ahead"): history is a Seq\* port, and the host owns the undo stack. The designs
    say so (changes.md, tree-core.md, base-library.md, viewer-models.md). Lane N builds it:
    - the port in `seqx-viewer/ports`, with a participant face (`checkpoint`, `restore`, `changed`) and a run face (an
      edit as a host step with its label and change kind);
    - an edit made outside a run reported to the host;
    - adapters in the composition roots: `@molstar/ui/history`, later the Mol\* plugin's, none headless.

    `SpecModel` stops importing `@molstar/ui/history`.
  - **Selection history**: none by default; a host may opt in (changes.md). Recommendation: as written.

## Stages and lanes

Lanes run in worktrees, each on its own branch `seq/<lane>`, merged by the orchestrator (§commit-often,
§lanes-report-not-write). At most two heavy lanes and two light ones at once (§parallelism). Lanes continue the first
slice's letters from K; stage 0's carry-over lane is S0.

### Stage 0: carried from the first slice (orchestrator and one fast lane)

- [ ] S0b (fast) `values.coarser`'s default to 1 (Gate 0, item 10): `RenderConfig`'s default, the core `_spec`'s
      Configuration row, the test that names the default; the 3M-bars heavy bounds still hold (rerun).
- [x] S0 (fast) carry-overs, each small (merged 2026-11-18, `360043baf`; the TP53 exons were right, checked against
      NM_000546.6):
  - the stale comment in `canvas/backend.ts` (two canvases: three);
  - `vitest` in the Storybook package's dev dependencies, and the Vite plugin's test under `src`;
  - the TP53 exon coordinates checked against RefSeq NM_000546. Fixtures only: no network in tests. If unverifiable
    offline, mark it in the fixture and the TODO.
- [ ] The Gene by region story rechecked in a visible pane (light and dark), with screenshots.
- [ ] GATE: green; screenshots in the report.

### Stage 1: the base (three lanes in parallel, after Gate 0)

- [ ] R (strong) `@molstar/core/data/immutable`, as `drafts/immutable.md` (draft 2) agreed at Gate 0 (moved to
      `design/immutable.md` before the lane starts):
  - `HashMap`, `HashSet` (CHAMP, flat below a measured threshold), `OrderedHashMap`, `OrderedHashSet`, `Vector`
    (32-way with a tail), transients that read like the persistent collections, cursors, per-collection hash and
    equality (SameValueZero keys by default), `defaultHash` with a finalizer;
  - a type-check-only spike of `StateTree` on the API, before it is frozen;
  - model-based property tests against JS collections, persistence and transient tests (the colliding-keys
    regression), collisions, equality; allocation cases (the opt-in tier); ratios against immutable.js with 22-char
    keys, and the bundle size (heavy lane);
  - its `_spec`, base-library.md's row, a v6 changelog entry.
- [ ] K (strong) `@molstar/core/repr/patch`, as `design/patches.md` (agreed 2026-11-18, revised after the second
      review), beside R:
  - the six ops with their inverses and `inverseNeeds`, `applyPatch` (result then base, all or nothing, normalization
    of what it made, refs following a move, validation by how each object came to be, link integrity old root against
    new, the language's validators, the result check), `diff` with its fallbacks and emission order, `patchLanguage`
    and `AnyPatchLanguage` (no core commit type: `SpecCommit` is lane N's, in `seqx-viewer`);
  - the core's pass over refs (what a ref names), and `RefIndex` (its signature as written) on plain `Map`s;
  - `repr/schema`'s `normalizeValue` and `Changes`; the child index (a native `Map` per children array) in
    `resolveAddress`;
  - `repr/schema`'s `validateChanges`;
  - the conformance cases as data with their options (§suites-are-data), the laws of diff as property tests;
  - the companion change in `seqx-core/lang`: `layer` ids required, `seq.lang.id.required` for every child, the
    builder and examples given ids, the language version moved (no migration: `0.x`);
  - format.md and tree-core.md as they now read; the `repr` `_spec` and a v6 changelog entry.
- [ ] K2 (strong) `merge3`, as `drafts/merge.md` agreed at Gate 0 (moved to `design/merge.md` before the lane starts),
      after K: the result type, the conflict table, independence by node, the laws, rebase from a base root.
- [ ] R+K (fast, after R and K; gated): `repr`'s ref index on `data/immutable`'s persistent trie, behind `RefIndex`,
      only if the measured case (a spec of 10,000 bookmark refs, edited by moves) shows the per-apply rebuild matters;
      else recorded as not needed, with the numbers.
- [ ] L (fast) the MessagePack codec:
  - `codecs/msgpack` in io on `common/msgpack`, passing `checkCodec`;
  - typed columns decoded without a copy where aligned;
  - the p53 specs round-tripping JSON ↔ MessagePack;
  - decode throughput as a heavy-lane bound (performance.md, at least 200 MB/s).
- [ ] M (fast) the asset manager's interface, as amended at Gate 0 (security.md):
  - `AssetManager.ResolveOptions` gains `credentials`, `redirect`, `signal`, `maxBytes` (enforced while the body
    streams; `fetchData` gains the same); `DefaultAssetManager` honors them, its defaults unchanged for every Mol\* app;
  - `checkAssetManager` (data and a function, on a fake fetch, never the network): the four options honored or the
    read refused, references and release, the shared download, decompression, errors not cached;
    `DefaultAssetManager` runs it;
  - Seq\*'s asset policy passes its choices on every read (browser: `credentials: 'omit'` unless the composition
    says otherwise, `redirect: 'error'`, the signal, the size bound), and claims all four in its documentation;
  - the assets `_spec` (`util/_spec/assets.md`), a v6 changelog entry; a public API addition to `@molstar/core`, no
    break.
- [ ] GATE: merged, reviewed, the checks green; `pnpm check:workspace` after `pnpm build` (core and io changed).

### Stage 2: editing and addressing (two strong lanes, after stage 1)

- [ ] N (strong) editing in `seqx-viewer` (and `seqx-core` where the model compiles increments), after K2:
  - `SpecModel.applyPatch` with `SeqLanguage`, its validate options and an **incremental `checkSeqLang`** as the patch
    validator (cost in proportion to `Changes`); then resolve and the view recompiling only what `Changes` reaches
    (composition.md); Seq's own conformance cases through `checkPatches(apply, cases)`;
  - **commits and recent roots**: `SpecCommit` as the port's checkpoint, `revisions.keepRecent` (64), the port's
    optional `root(revision)`;
  - history through Seq\*'s history port: the port in `seqx-viewer/ports`, `SpecModel` free of `@molstar/ui/history`,
    and the `@molstar/ui/history` adapter in the default composition. One step a batch, an agent's turn as one step,
    an outside change reported to the host. Tested with two participants on one host stack (a Seq\* spec and a
    stand-in structure state, undone in order), and with no history at all;
  - **the rebase test** (why N waits for K2): a stale agent patch refused, then rebased by the host's opt-in through
    `merge3` with the base root from the recent roots;
  - **previews as patches** re-applied when the root changes, dropped with an issue when they no longer apply (a
    preview surviving an agent's patch mid-drag), ending in one patch or nothing (Esc cancels);
  - **story steps as stored patches** (`patchLanguage(SeqLanguage)`): entered and left (navigation, not history),
    shown over the base and the preview, edits while shown patching the base, masked edits reported, a step that no
    longer applies shown as such;
  - every command declaring its change kind, with the test changes.md asks for: every command's kind matches what it
    changes, undo touches only spec patches, an arrival never changes the spec's hash.

  Target: an edit (patch, resolve, redraw) under 16 ms, its undo the same; a quick representative in the lane.
- [ ] O (strong) addressing in `seqx-core` (and the viewer's commands):
  - ids and hierarchical addresses `path/(query :key value)` resolving to a node, a Location and back;
  - user-defined refs and bookmarks, on dynamic nodes too: a ref waits while its node is not held;
  - actions and links as commands on addresses, with no markup from data;
  - addresses of the session state (selection and focus).

  The address's two consumers live in one test: the viewer and an export resolving the same address.
- [ ] GATE: merged, reviewed with the look (an editing story, a bookmarks story), the checks green.

### Stage 3: selection and security (one strong lane, one fast; after stage 2)

- [ ] P (strong) selection, highlights and the cursor:
  - `SelectionModel` (element keys and a Location, the anchor, the cursor), `HighlightModel` (channels; the pass exists
    and is empty);
  - the default preset's click, Shift-click and drag-select bindings and the keyboard cursor;
  - the grammar's toggle;
  - `Emphasis.projected`, so the hover and the selection projected onto linked sections are drawn (the open item from
    the first slice);
  - the selection tint's place as Gate 0 decided, as a config item;
  - the selection holding its nodes;
  - selection as addresses (from O).

  Target: a drag or cursor step under 4 ms; hover and selection allocation-free warm (the opt-in tier).
- [ ] Q (fast) the security bounds:
  - every row of security.md's table enforced before allocating, each a config item with its coded issue (EC-58, 76,
    86, 91, 179);
  - SVG output escaping tested with hostile labels;
  - the headless policy's root and hosts.
- [ ] GATE: merged, reviewed with the look (selection on both p53 views, light and dark), the checks green.

### Stage 4: stories and use (orchestrator, light)

- [ ] S1 stories under `Seq* (next)/Examples/p53/`:
  - **Edit and undo**: add and hide a track, set a theme, undo, redo, a preview by dragging a track's height;
  - **Bookmarks and links**: a ref on a transcript that arrives later, a link that moves the view;
  - **Selection**: drag-select on the gene, projected onto the protein.

  Each has its headless twin composing it as the story does.
- [ ] S2 the Features and Examples pages list what is built.
- [ ] GATE (§used-before-done): the author uses the stories in a visible browser. What they find is fixed or listed
      before the close.

### Stage 5: measure (orchestrator, light)

- [ ] M2 the heavy lane and the opt-in tiers:
  - an edit and its undo;
  - a selection step;
  - MessagePack decode;
  - the layout's allocation bound and draw and hover at zero (regressions from the first slice);
  - **browser** frame times (draw included) and first paint, in Chrome at DPR 2, on the p53 stories and the 3M
    fixtures;
  - bundle sizes again (editing and addressing add to the viewer).
- [ ] GATE: the author confirms or moves the targets the slice touched.

### Finalize

- [ ] FINALIZE: reconcile:
  - the designs where the code taught otherwise (drift notes resolved);
  - the `_spec`s with what is built;
  - `status.md`, the features rows, the edge cases' outcomes, `TODO.md`;
  - the fold-back (what this slice taught, for the author).

  Then write the record in `.seq/records/second-slice/` and end `state.md` with its closing line.

## Order and parallelism

Gate 0 → stage 0 (beside Gate 0) → R, K, L, M in parallel (two heavy, two light at once: R and K first, L and M
light) → K2 (needs K) and O (needs K for refs as patches) → N (needs K2: the rebase test) → R+K (needs R and K, gated
by its measurement)
→ P (needs O's addresses for the selection) beside Q → stories and use → measure → finalize.

## Hot files (the orchestrator's only)

The first slice's list:
- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/`;
- `packages/seq/*/package.json` and `tsconfig*.json`;
- `.seq/`.

New in this slice: `packages/core/package.json` and `packages/io/package.json` (lanes R, K, L and M ask for export
changes in their reports).

## Out of scope

Moving Mol\*'s state tree and `@molstar/model` off immutable.js (a later, separate PR); the figures milestone (the text layout, alignments, rings, dot plots, PNG, the render CLI), structures, services,
bigWig and range reads, the other datasets, concurrent writers (EC-93), agent tools.
