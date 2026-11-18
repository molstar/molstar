# State: second-slice   (updated 2026-11-18, by the orchestrator)

Branch: `packages-seq`. The loop started 2026-11-18 (`/loop run .seq/execution/second-slice/loop.md`).

## Done

- Lane S0b merged: `values.coarser` defaults to 1; 3M bars within the heavy bounds (pan 3.0 / 5.8 ms median / max).

- Lane L merged: the MessagePack codec of the repr tree; reviewed, nine findings fixed in one round (a one-line
  BOM fix in io's `common/msgpack`, allowed); the merged tree 1108 tests green.

- Lane M merged (`4410e7738`): the asset manager's options and `checkAssetManager`; reviewed, a credential leak and
  nine more fixed in one round; the merged tree 863 tests green.

- The interfaces reviewed independently three times (Fable 5.2, 2026-11-18) and revised each time; the last check's
  verdict: `design/patches.md` (lane K), `drafts/merge.md` (K2), `drafts/immutable.md` (R) and the state rules (N) sound
  to build on, with its last fixes in `31e9b8780`.

- Plan, loop and state written; the loop rewritten for this slice.
- Before the loop (the author, 2026-11-18): the history designs corrected; the patch and diff design drafted
  (agreed the same day as `design/patches.md`, git's model); Gate 0 aligned with the first slice's record;
  stage 0's lane renamed S0; the first slice's worktrees and branches removed.

## In progress

- Stage 0: S0 merged (`360043baf`); S0b merged. The Gene by region recheck waits for a visible
  pane (the author asked to show it).
- Lanes M and L reported; checks green on the orchestrator's run (M: 789 tests; L: 305 and the heavy bound, lint green
  once io's `./repr/_test/*` export null is added); both reviewed (opus): M, no regression for existing callers but a credential leak through the asset-id path, among
  ten findings; L, round-trip breaks for a leading U+FEFF and lone surrogates, and decode writing into the caller's
  buffer, among nine. Both in their fix round.
- Lane K reported: checks green on the orchestrator's run (589 quick, 35 heavy; lint green once core's
  `./repr/patch/_test/*` null is added); under independent review (opus). Its layer keys are now addresses, so the p53
  stories need the look in a visible pane before the merge (the pane is hidden).
- Lane R reported: checks green on the orchestrator's run (320 quick in 0.7 s, 10 heavy, 12 allocation; lint green once
  core's `./data/immutable/_test/*` null is added); under independent review (opus).
- Stage 1, started from `f629c1374` (the briefs): R (`seq/r-immutable`, opus), K (`seq/k-patches`, opus), L
  (`seq/l-msgpack`, sonnet), M (`seq/m-assets`, sonnet). S0b waits for a light slot.

## Next

- On each lane's report: checks, the look where something is drawn, an independent review, one fix round, merge.

## Waiting for the author

None.

## Decided

- **Three decisions, as recommended** (the author, 2026-11-18: "Ok go as recommended"):
  - **Persistent collections against immutable.js** (lane R): an int32 key is hashed by its value; a string-hash memo
    is opt-in, inside a config's `hash`, documented; the target reads "no slower than immutable.js on collections of
    1,000 or more, and within 2× below". immutable.md says so. Lane R builds the first two in its fix round.
  - **MessagePack throughput** (lane L): the 200 MB/s target applies to documents of columns; node-heavy documents
    keep the 20 MB/s regression guard. performance.md says so. No trusted decode.
  - **`layer` gains `hidden: boolean`**, default false (gaps 5): language-description.md says so; lane K adds it to
    SeqLanguage in its fix round and its cases hide a layer with it.

- **Gate 0** (the author, 2026-11-18): "rest as recommended", with one amendment.
  - The first slice's open decisions 1–9 confirmed as built; 10 `values.coarser` default 1 (lane S0b); 11 below; 12 the
    TP53 exons checked by S0 (right, NM_000546.6).
  - `@molstar/core/data/immutable` as drafted (lane R); `merge3` as drafted, both questions yes (lane K2); range reads
    wait for bigWig; MessagePack in io on its hooks (lane L); selection history none by default.
  - **The asset manager, amended**: "we will need to be able to support credentials in the asset manager, but given
    how it's implemented, each app can get its own derived implementation. so what i am thinking: we make sure the
    interface has all we need, and then app composition still provides that. we should still have basic conformance
    tests for the asset manager that any implementation needs to satisfy." So `ResolveOptions` carries credentials,
    redirects, the signal and a streamed byte limit; each app composes its manager; `checkAssetManager` is the suite
    every implementation passes; Seq\*'s policy passes its choices per read (security.md). Lane M.

- History as a Seq\* port, the host owning the undo stack; `@molstar/ui/history` one adapter (the author,
  2026-11-18: "yes, go ahead"). Gate 0's history item is settled; the designs say so (changes.md, tree-core.md, base-library.md,
  viewer-models.md).

## Decisions taken (awaiting the author)

None yet.

## Verify

As loop.md, Checks.
