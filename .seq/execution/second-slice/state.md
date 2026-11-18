# State: second-slice   (updated 2026-11-18, by the orchestrator)

Branch: `packages-seq`. The loop started 2026-11-18 (`/loop run .seq/execution/second-slice/loop.md`).

## Done

- Lane K merged (`b3ea63511`): `@molstar/core/repr/patch` and `layer.hidden`; reviewed (opus), both law-3 failures and
  the rest fixed in one round (the fix's 4 tests fail on its parent, checked); the merged tree 14,380 quick tests
  green, the heavy lane 3 runs of 3 after doubling the single-sample bounds. Its drift is in patches.md, for the
  author to confirm. The look done after the merge (below).

- Lane R merged (`47b754e25`): `@molstar/core/data/immutable`; reviewed (opus), a lost-collision-bucket delete and the
  rest fixed in one round (the fix fails on its parent, checked); 14,333 quick, 30 allocation, 10 heavy green.

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

- Stage 2 started (2026-11-18, from `57b1d7c4e`; the lane worktrees were cut from master, each lane reset onto it):
  - K2 merged (`5927c9e4b`), R+K merged (`989292665`), each after one review round.
  - O reported (`3f4d0f128`, checks green on the orchestrator's run); looked at in a visible pane from its worktree
    (Storybook on 6008): the links story works in light and dark; in the bookmarks story a bookmark's focus says
    "Focused" but the view does not move, and "Bookmark this view" records the last focused address. Under review.
  - Law 3 with stored patches: case 572 fixed by the orchestrator (`ff2ca5476`); seed 4 case 43 (pre-existing) is
    lane K3's (`seq/k3-stored-law3`, opus).
  - N's brief written (`briefs/n-editing.md`); it starts from the merged tree once O lands.

## Next

- On each lane's report: checks, the look where something is drawn, an independent review, one fix round, merge.

## Waiting for the author

None.

## Decided

- **Integer keys and lane K's drift, as recommended** (the author, 2026-11-18: "ok as recommended"):
  - integer keys keep `fmix32`; their miss against immutable.js at 1,000 entries (up to 1.9×) is accepted, each case
    with its own heavy bound; an identity `hash` is a config's choice (immutable.md);
  - patches.md's "The inverse and the allowances" confirmed: undoing a removal of already-dangling refs needs
    `allowDangling`, `inverseNeeds` is conservative, a replace is an edit in place.

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
