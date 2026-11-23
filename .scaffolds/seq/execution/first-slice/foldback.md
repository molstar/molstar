# Fold-back: what the first slice taught (draft, for the author)

The first slice built eleven lanes (A to J, C2) through briefs, an independent review and one fix round each, with gates
after each stage. Each proposal below says what happened (the evidence is the commit or the note), what it would change,
and where. **Nothing here is applied until the author approves it** (§authors-words-decide). Once approved:

- rules for how the code works go into the packages' `_spec/` (they outlive `.seq/`);
- rules for how work runs go into `.seq/principles/` as adaptations, and to the oz package's home as a fold-back;
- design changes go into `.seq/design/`.

Where to: **P** principles (the oz package's home; an adaptation here until it is merged there), **D** a design, **S** a
package `_spec`, **C** a check.

## Process

1. **Look at the picture before the merge, not only at the gate.** (P: §stories-are-tests, §used-before-done; the
   loop file)
   - *What happened*: every lane passed its own tests and an independent review. Even so, five bugs of what is drawn
     were found only by opening a story after the merge:
     - values drawn on the right third only (`26c37f6bc`, `c960f158a`);
     - the dark theme ignored (`c80e00a1d`);
     - the base canvas not redrawn after a re-hit (`4a790a573`);
     - nested domains drawn late, because exports detached their source (`edf068a8e`, `72bde5669`);
     - an SVG of height 0 (`cdf667452`).

     None was a design question; each had a one-line cause.
   - *Proposal*: the lane's review includes the look. The reviewer, or the orchestrator before the review, opens the
     lane's stories in light and dark and reports what it saw; the gate's screenshots stay as the record. The rule
     §stories-are-tests already says "looked at once before it is called done"; the loop file put the look after
     the merge.
2. **A fix's test is shown failing without the fix.** (P: testing; C)
   - *What happened*: the first test for the height-0 SVG passed with and without the fix, because the test's
     `settle()` waited for the frames that hid the bug. It was caught only because the run was repeated on the parent.
   - *Proposal*: a bug fix lands with a test, and the report names the command that shows it failing on the fix's
     parent (the reviewer reruns it).
   - *Check*: a script that reverts the fix's non-test files on a scratch worktree and expects the new test to fail.
     Optional, at gates.
3. **A visual or timing observation states its conditions.** (P: §used-before-done)
   - *What happened*: "the domains appear late" was partly an artifact. The browser pane was hidden, and a hidden page
     runs no animation frames. Time went into chasing it before a headless test on a fake clock separated the real bug
     from the artifact.
   - *Proposal*: an observation that a fix rests on names where it was made (a visible or hidden page, the browser,
     the DPR). It is reproduced headless on a fake clock before the fix starts. A hidden pane is inconclusive for
     anything that moves.
4. **Configurable over fixed works, but the defaults pile up unread.** (P: §gaps-queue; carried as a tension, below)
   - *What happened*: the author's rule produced about sixty config items and defaults over gaps 6 and 9–14. Each one
     is documented and tested. None has yet been confirmed, and gaps.md lists them in prose, hard to decide on in one
     pass.
   - *Proposal*:
     - each package's `_spec` keeps one Configuration table (item, default, why, alternatives);
     - gaps.md points at the tables and lists only the defaults a lane was unsure of;
     - the gate asks for the unsure ones, not all of them.
5. **Hot files with the orchestrator, lanes asking in their reports.** (P: §orchestrator-owns-lockfile; a verdict, no
   change)
   - *What happened*: eleven lanes ran with no conflict on a lockfile, the inventory or tsconfig. Every hot-file change
     came as a request in a report and was applied after the merge.
   - *Proposal*: none to the rule. Record the verdict and the count, so a later reading can rely on it.
6. **Refuted review findings are named in the merge.** (P: the standing tension "proportional review" in the oz
   README)
   - *What happened*: this slice did not record them either. So, like the earlier evidence, it cannot settle the
     tension.
   - *Proposal*: adopt the proposal already carried there. The merge message names the findings refuted or left
     unfixed and why. The loop file says so.
7. **Performance targets name their path and their test when they are written, and get a first measurement early.**
   (D: performance.md; P: §ui-performance)
   - *What happened*:
     - the targets were first measured at stage 7, after four stages built on them;
     - "a frame allocates nothing" was read two ways, the draw loop or the layout too, and needed the author;
     - the 3M-bars miss came from the gap between a source's levels, which no design spoke of.
   - *Proposal*:
     - each target names the code path it covers and the test that holds it;
     - the stage that builds a path runs a quick representative of its target (one scenario, one bound);
     - the design states the relation between a source's levels and the frame's cost, its density and the level
       choice's `k`.

   - *Lane J's results*:
     - both misses closed by the code, not by moving the targets: 3M bars at 2 ms max, draw and hover at 0 bytes;
     - the review then showed the decision's premise had moved: with the direct path, `coarser` 2 is no longer needed
       for time (gaps 15). A target decided on a measurement is re-read when the measurement moves.

## Code (into the `_spec`s)

8. **Reused buffers are compared against fresh ones.** (S: core Render; D: layout-and-draw-lists)
   - *What happened*: reusing buffers (the performance rule) brought in a correctness bug: `grow` did not copy.
     Unit tests on small inputs never grew a buffer mid-layout.
   - *Proposal*: every buffer reused across frames has a test that lays the same input out twice and compares the
     lists: once into fresh buffers, once into buffers warmed by a different, smaller input. Purity by contract is then
     a test, not a hope.
9. **What two consumers share is held per holder, and tested with both live.** (S: core Data; D: data-store)
   - *What happened*: the viewer and the story's SVG exports share one store. An export ending detached the viewer's
     nested sources (`edf068a8e`), and a detached source left tiles pending (`72bde5669`). The single-consumer tests
     could not see either.
   - *Proposal*: anything a store or a registry shares is reference-counted by holder. Its tests include two consumers
     interleaved, one ending while the other reads.
10. **Derived signals publish on change only.** (S: viewer; P: §relevant-rerenders)
    - *What happened*: `layoutIssues` was set every frame with an equal value. The story re-exported on every set, and
      the churn masked the height-0 bug.
    - *Proposal*: a state derived each frame is set only when it differs (structural equality for lists). The viewer's
      `_spec` says which signals compare how.
11. **"At rest" means the picture, not the input.** (S: viewer, `rested`; D: viewer-composition)
    - *What happened*: `rested` followed only the view's moves, so a picture changed by data arriving while the view was
      still never "rested" again (`cdf667452`).
    - *Proposal*: a signal named for a state of the picture (rested, drawn, ready) fires on every cause that changes the
      picture. Its test drives each cause: a move, data arriving, a resize.
12. **Each story has a headless twin that composes it as the story does.** (P: §stories-are-tests,
    §whole-composition-tested; S: storybook)
    - *What happened*: `gene-by-region.test.ts` composes the viewer, the delayed source and the SVG exports as
      `example-root.ts` does. That one test found three bugs the view's own tests could not see.
    - *Proposal*: a story that composes more than one consumer (a canvas and an export, two views) has a headless test
      of that composition on a fake clock, next to the story.
13. **Every asynchronous source has a synthetic twin with an injectable clock.** (D: testing, sources)
    - *What happened*: the delayed `seq` source made every async bug of stage 6 reproducible in milliseconds, with no
      real waiting.
    - *Proposal*: a source runtime lands with its synthetic, delay-injectable twin, used by its stories and its tests.
14. **An edge case cited in code exists in the catalog.** (C: a drift check; D: edge-cases)
    - *What happened*: lane H's followups recorded EC-143 as missing from the catalog. It was not missing: the
      catalog is the source's (`sources/seqstar/edge-cases.md`), and the design's `edge-cases.md` only refers to it.
      The finalize reconcile found the error by hand.
    - *Proposal*: the drift check resolves every `EC-<n>` in `packages/seq` against the catalog, like its other
      registry references (§mechanical-drift-checks). The ids are numbers only, so code citing them stays free of
      `.seq/` paths. When `.seq/` goes, the catalog's kept items need a home in a `_spec` (or the citations go).

15. **The rules found a bug while they were being written.** (evidence for 11, no new rule)
    - *What happened*: writing 11 into the viewer's `_spec` ("a signal named for the picture fires on every cause,
      each driven in its test") needed a test per cause. The test found `rested` firing twice a move: once at the
      settle timer, again after the settled frame's layout (`e16fe448c`).
    - *Proposal*: none beyond 11. Recorded as evidence that 11 catches what the story's test alone did not.

## Self-check (§self-check), on the claim "these proposals make the principles truer to what was built"

1. **Assume it holds.**
   - 1–3 move checks earlier: before the merge, and before the fix.
   - 8–14 turn bugs into rules, so the next slice cannot repeat them.
   - Every proposal has a commit behind it.
2. **What it rests on, negated.** It rests on the bugs being of kinds that recur. Negated: they were one-offs of a
   first slice (a new renderer, a new store, a new viewer), and rules made from them only add process.
3. **The negated world.** It explains why five visual bugs came in two stages (4 and 5) and none in stage 6's own
   code. Stage 6's bugs were all about sharing between consumers (9–11). By that reading, 1 and 8 are first-slice
   costs, and 9–12 are the recurring kind: every later slice adds consumers.
4. **Both held.**
   - *Integrates*:
     - 9–12 and 14 stand under either reading (sharing and composition recur; a check is cheap);
     - 2 is cheap, and its one miss was real;
     - 5 becomes a verdict with its count.
   - *Carried* (standing tensions, into the oz README beside "proportional review"):
     - **"configurable over fixed" against a config surface the author can confirm.** The rule grew about sixty
       defaults in one slice, and the gaps queue cannot be decided at that size. 4 eases it but does not settle it:
       does a default need confirming at all once it is documented, tested and easy to change?
     - **"look before merge" (1) against the cost of a look per lane.** If the next slices' visual bugs stay at zero
       before merge, 1 was a first-slice cost.

## Not proposed

- Rules from the `check:workspace` failure (the `@molstar/ui` export): it is the host's, recorded in gaps.md.
- Anything from the hidden-pane artifact beyond 3: it is the tool's behavior, not the code's.
