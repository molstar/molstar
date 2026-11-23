# second-slice: core editing and addressing in Seq\* (next) (2026-11-18)

## Asked

- Run the second slice's loop (the work id `second-slice`): core editing and addressing,
  selection, security bounds, MessagePack, built in `packages/seq` and `@molstar/core`/`@molstar/io` through lanes.
- Standing rules carried in: the Mol\* header with the author only; no mention of the previous implementation or of the
  scaffold in code; the Storybook group "Seq\* (next)"; "if there is a way to make something configurable, then do that
  and choose your recommendation as the default"; the fold-back is approved before it is committed; Seq\* owns no host
  things (history a port, the animation loop the application's); persistent collections in
  `@molstar/core/data/immutable`; the asset manager carries credentials through its interface.
- At the gates, the author used the stories in a visible browser and reported what they found (stage 4).
- Close the slice: measure, finalize with the fold-back as a step ("make the fact that foldback is the step part of
  the principles themselves"), then plan the third slice.

## Decided

- Gate 0's three decisions and the stage gates "as recommended"; the integer-key miss in the immutable collections
  and lane K's drift accepted (design/immutable.md).
- Stored patches type their refs by key across the language's kinds, not by a dry run (drift note option 1 with (a));
  a step that no longer applies is shown so ([design/patches.md](../../design/patches.md)).
- Stories: a story is a sequence of full states; a step is stored as a full state or as a patch, keyframes and deltas
  (Braid-style, `story.maxChain` 8); transitions planned from two states' diff and sampled purely; the animation loop is
  the application's, behind an `AnimationPort` ([design/stories.md](../../design/stories.md); "ok, let's do it like
  this"). Built later: the library and the port.
- The Complex dataset is CDK2–cyclin A with an AlphaFill ligand ([design/examples.md](../../design/examples.md)).
- The selection's look, shaped at the author's look: one band per section across tracks and gaps (`column`, `tracks`
  configurable), a soft tint with one outline, the hover outside its mark, the cursor after a key; range squircles
  consistent with the highlight.
- From the author's use (stage 4's gate):
  - several ranges in one selection (Cmd/Ctrl+drag adds, Alt+drag removes, Shift extends the anchor's range): "you
    can build it, and the gestures are configurable anyway";
  - following both ways: "could we have following in both ways? .. i think we should" (`follow.followerMoves: lead`);
  - zooming out (a view history, back and forward, the zoom keys), and every drawn mark hit;
  - to the third slice: a reversed follower, the protein and its gene on one axis (letters through a codon map),
    collapsed introns ("we can do these in the next slice").
- The fold-back is a step of every close, always written; merging it is the author's decision ("not every foldback
  should be "merged" .. some exploration work for example should not merge its lessons"); checked by
  `_checks/closed-work.ts`.
- Stage 5: the change validators their own entry (`@molstar/seqx-core/lang/changes`), the read bundle's target 160 kB
  ("Yes, as recommended").
- The nine gaps as recommended; the fold-back merged ("Yes, all of it"), with B1 rewritten after the author asked to
  apply the self-check to designs.

## Delivered

Sixteen lanes, each through a brief, an independent review, one fix round and (where it draws) a look in light and
dark before the merge; 238 commits from `fc2bba6ce`.

- `@molstar/core`: `data/immutable` (R), `repr/patch` with apply, diff, `merge3`, the ref index derived lazily, stored
  refs typed by key, conformance (K, K2, R+K, K3); the asset manager's options (M). `@molstar/io`: the MessagePack codec
  of the repr tree (L).
- Seq\*: editing (`SpecModel`, the history port and its `ui-history` adapters, previews; N), addressing (refs, links,
  bookmarks, `followMoves`; O), selection, highlights and the cursor (P), the security bounds and the headless asset
  policy (Q), the p53 example stories and pages (S), several ranges and followers that lead and the far-zoom hover (T),
  hit boxes for merged runs, zooming out, the focus of an item through a path (U); `layer.hidden`, `track.hidden`, the
  language at 0.3.0 (refs, snapshots).
- Fixes found by looking: WRAP53's second exon never drawn (a tile-edge dedupe by id alone, `b73d3199b`), the Editing
  story's host stand-in, a domain focus zooming to one cell (U).
- The test harness: `bestOf`, and `itWithin`'s heavy bound the best of three runs of the whole test.
- Verified: the merged tree's quick lane (800 tests), medium and heavy lanes, the allocation tier, typecheck, lint,
  format, `check:workspace`, the build, the five `.seq` checks; every fix with a test shown failing on its parent.

## Measurements

Every time and allocation
target passes in Node and in Chrome at dpr 2 (an edit 0.44 ms median, its undo 0.3 ms at p95; selection steps under
0.2 ms; browser pan and zoom about 1 ms median, 5.5 ms worst, no frame dropped; node-heavy MessagePack 60 MB/s; draw and
hover 0 bytes). The read bundle missed (176 kB against 150): after the entry split 150.6 kB against 160; headless
rendering 267.9 kB, the viewer 504.9 kB. Not measured: first paint whole in the browser, the 3M fixtures in the
browser.

## Learned

- The author's use found missing design, not bugs: six of eight findings, five of them the negation of a design's own
  premise (focus and back, one-way follow, one range, one direction, separate sections). The self-check now runs on
  every design when it is agreed, interaction is a matrix, and a novice pass precedes the author's use.
- Canvas, SVG and hits can agree and all be wrong; drawing is checked against the data.
- A restriction with no reason ("letters on their own system only") blocks a figure; it names its reason now.
- An index that re-exports an editor's module makes every reader pay (26 kB); sizes are measured at every gate.
- Single-sample bounds flake under load; the best of three runs, and heavy runs never overlap.
- The look before the merge paid in four of five drawing lanes: the tension carried from the first slice is released.

## Open decisions

None left with the author: the gaps and the fold-back are decided.

## Left open

- In `packages/seq/TODO.md`: the stories library and the animation port, navigators, the reversed follower, letters
  through a map, collapsed introns, `track.add`, a link before load kept as waiting, the history's toolbar, exporting
  `projectPoint`, the codon-split follow edges, drawing checked against the data in every example view, and the rest
  of the followups.
- Defaults no story shows (store sizes, bounds) are never confirmed by use: the part of "configurable over fixed" still
  carried.
- Whether an early use gate pays when a slice adds figures, not interaction: counted at the third slice's close.

## Self-check

The fold-back's self-check (in the work's archive, found by its work id) on the claim "these
proposals make the principles truer to what was built". Delta: merged the self-check on designs, the interaction
matrix, drawing against data, a restriction's reason, the best of three, gaps at every gate, sizes at every gate;
released "look before merge"; settled "configurable over fixed" except for defaults no story shows; carried "does an
early use gate pay for figures".
