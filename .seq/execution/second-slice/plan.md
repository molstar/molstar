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

- [ ] GATE 0: the author confirms, each with its recommendation:
  - **the first slice's open decisions** (its record, Open decisions). The orchestrator re-presents them with the
    first slice's evidence:
    - undeclared columns tagged;
    - cyclic ranges wrap;
    - the query decisions;
    - the protocol details;
    - the `seq.lang.*` codes;
    - the selection tint over the data;
    - letters' density in CSS px.

    Recommendation: confirm all as built.
  - **public API changes in the base library** ([base-library](../../design/base-library.md)):
    - `@molstar/core/repr`'s `patch` and `diff` (patch kinds, inverses, the revision a patch names). Recommendation:
      as changes.md and tree-core.md say.
    - **the asset manager**: `credentials: 'omit'`, `redirect: 'error'`, the signal passed through, a streamed byte
      limit. Recommendation: yes, in this slice, with the policy then claiming all four.
    - range reads in the asset manager. Recommendation: wait for bigWig (milestone 3).
  - **the MessagePack codec** builds on io's existing `common/msgpack`, extended for typed columns without a copy where
    aligned. Recommendation: yes, a `repr` codec in io (io may depend on core; `repr` stays dependency-free).
  - **history as a Seq\* port** (the author, 2026-11-18: the viewer lives inside a bigger app, so the app owns the
    undo stack). Seq\* never owns one. `seqx-viewer/ports` defines a history port in Seq\*'s own terms:
    - a participant face the host drives: `checkpoint()` returns the spec's root, `restore(c)` sets it back and drops a
      preview, `changed` tells of edits;
    - a run face Seq\* asks of the host: an edit runs as a step with its label and change kind, so a host's stack
      interleaves Seq\* edits with its own (a structure's, Prism's);
    - an edit made outside a run is reported to the host, not acted on: whether earlier steps are invalidated is the
      host's rule.

    Adapters, each in the composition root that wires it:
    - `@molstar/ui/history` for the standalone viewer and Prism;
    - the Mol\* plugin's state history, later;
    - none (no undo) headless.

    `SpecModel` stops importing `@molstar/ui/history`. changes.md's "through `HistoryPort`" becomes "through Seq\*'s
    history port; `@molstar/ui/history` is one adapter". Recommendation: yes.
  - **selection history**: none by default; a host may opt in (changes.md). Recommendation: as written.

## Stages and lanes

Lanes run in worktrees, each on its own branch `seq/<lane>`, merged by the orchestrator (§commit-often,
§lanes-report-not-write). At most two heavy lanes and two light ones at once (§parallelism). Lanes continue the first
slice's letters from K.

### Stage 0: carried from the first slice (orchestrator and one fast lane)

- [ ] K0 (fast) carry-overs, each small:
  - the stale comment in `canvas/backend.ts` (two canvases: three);
  - `vitest` in the Storybook package's dev dependencies, and the Vite plugin's test under `src`;
  - the TP53 exon coordinates checked against RefSeq NM_000546. Fixtures only: no network in tests. If unverifiable
    offline, mark it in the fixture and the TODO.
- [ ] The Gene by region story rechecked in a visible pane (light and dark), with screenshots.
- [ ] GATE: green; screenshots in the report.

### Stage 1: the base (three lanes in parallel, after Gate 0)

- [ ] K (strong) `@molstar/core/repr` patches:
  - `repr/patch`: the patch kinds (set, insert, remove, move by address), validation before apply, inverses;
  - `diff` of two roots by hash, structural sharing kept;
  - a patch naming its base revision, refused on another (EC-81);
  - streams of patches as one batch;
  - the conformance cases as data (§suites-are-data);
  - the `repr` `_spec` and a v6 changelog entry.
- [ ] L (fast) the MessagePack codec:
  - `codecs/msgpack` in io on `common/msgpack`, passing `checkCodec`;
  - typed columns decoded without a copy where aligned;
  - the p53 specs round-tripping JSON ↔ MessagePack;
  - decode throughput as a heavy-lane bound (performance.md, at least 200 MB/s).
- [ ] M (fast, if Gate 0 agrees) the asset manager's fetch: credentials omitted, redirects refused, the signal
      passed, a streamed byte limit. Then the asset policy claims them, and security.md drops "what waits". Tests run on
      a fake fetch, never the network.
- [ ] GATE: merged, reviewed, the checks green; `pnpm check:workspace` after `pnpm build` (core and io changed).

### Stage 2: editing and addressing (two strong lanes, after stage 1)

- [ ] N (strong) editing in `seqx-viewer` (and `seqx-core` where the model compiles increments):
  - `SpecModel.applyPatch` validated, then resolve and the view recompiling only what the changed addresses reach
    (composition.md);
  - history through Seq\*'s history port (Gate 0): the port in `seqx-viewer/ports`, `SpecModel` free of
    `@molstar/ui/history`, and the `@molstar/ui/history` adapter in the default composition. One step a batch, an
    agent's turn as one step, an outside change reported to the host. Tested with two participants on one host stack
    (a Seq\* spec and a stand-in structure state, undone in order), and with no history at all;
  - previews as a transient overlay, ending in one patch or nothing (Esc cancels);
  - snapshots entered and left (navigation, not history);
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

Gate 0 → stage 0 (beside Gate 0) → K, L, M in parallel → N and O in parallel (N needs K; O needs K for refs as patches)
→ P (needs O's addresses for the selection) beside Q → stories and use → measure → finalize.

## Hot files (the orchestrator's only)

The first slice's list:
- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/`;
- `packages/seq/*/package.json` and `tsconfig*.json`;
- `.seq/`.

New in this slice: `packages/core/package.json` and `packages/io/package.json` (lanes K, L and M ask for export changes
in their reports).

## Out of scope

The figures milestone (the text layout, alignments, rings, dot plots, PNG, the render CLI), structures, services,
bigWig and range reads, the other datasets, concurrent writers (EC-93), agent tools.
