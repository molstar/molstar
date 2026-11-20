# Fold-back: what the third slice taught (draft, for the author)

The third slice built version 1's second milestone, figures, in over sixty lanes (C0 to LX; the list in `state.md`,
about 640 commits from C0's merge), each through a brief, an independent review and a fix round, a look before every
merge that draws, and gates after each stage; CE is parked unmerged for the fourth slice. It also ran the second
slice's new steps for the first time: a self-check on every design, the interaction matrix, a novice pass, the author's
use as lanes merged, test-only lanes for the edge cases (EO, ET, AT) and a finalize prepared by read-only lanes (SR,
EO). This file follows the principles' rule "The fold-back is a step of every close": what the slice taught, the
earlier fold-backs' rules checked against this slice's evidence, the count of what each gate found, the self-check, and
the author's decision on its own line. **Nothing here is applied outside this folder before that line.**

Where to: **P** principles (the instance's adaptations, and the oz package's home), **D** a design, **S** a package
`_spec`, **C** a check, **L** the next loop file.

## A. What each gate and use found, by kind

From `state.md` ("Findings by kind"), `novice-findings.md` and the use counts in `state.md`. The table in `state.md`
stops at stage 4; the reviews of stages 5 and 6 (FS-A1 to LX) were recorded in its prose, not counted, so their
findings are missing here (a miss of the orchestrator, B6).

| Gate or use | Bug | Missing design | Look |
|---|---|---|---|
| Lane reviews and looks, stages 0–4 (C0, AX, SL, TH, TX, TA, TV, RB, PG, NV, GD, RF, CH) | 82 | 3 | 17 |
| Stage gates (stage 2's, stages 3–4's) | 2 | 0 | 2 |
| Novice pass (stage 6, `novice-findings.md`; 21 findings, 2 good) | 5 | 8 | 8 |
| The author's use (stages 3–6, as lanes merged) | 6 | 7, and 1 after the count (LX: leaders crossing) | 6 |
| Edge-case lanes (EO's outcomes, ET's and AT's tests) | 3 (EC-135 not wired; the SVG's pairs arcs unnamed; a folded clade unnamed) | 2 (EC-63 slices against glyphs; EC-138 against saturated fills) | 0 |
| Finalize's reconcile (SR's list, this close) | 0 | 6 (configuration the designs name, the code lacks) | 0 |

Read: lane reviews find bugs (82 of 98 bugs); the novice pass and the author find missing design and looks (16 of 27
missing designs, 14 of 33 looks), and nearly all of the author's were **defaults**: how a figure looks before anyone
sets anything (B1).

## B. What this slice taught

1. **Defaults are design, decided with a picture before the build.** (P: §design-first; D: each layer's and layout's
   design; L: the design gate)
   - *What happened*: of the author's seven missing designs and six looks, eleven were a default nobody had drawn: the
     CDR tracks "feel empty" (heights by content, DB), "really good defaults for everything" (the defaults audit DA,
     twelve findings, after the build), the plasmid's labels inside and beside inconsistent (one label rule), the clade
     visuals "quite ugly" (RX), "2 columns hidden" misleading and misplaced, the page moving when a legend row appears
     (`channels: reserve`), the ruler blank over the constant domain (`unnumbered`), a leader through another label
     (LX). Four lanes (DB, DS, LX, RX) and an audit lane followed. The designs were self-checked (the second slice's
     rule), and the self-check's negation found inverses, not defaults.
   - *Proposal*: a design that adds a layer, a layout or a part has a **Defaults** section: what a reader sees with no
     setting (heights, labels, colors, legends, empty and dense cases), each with a sketch or a reference picture, and
     the design gate reviews it as the defaults audit did. The audit's questions (does it read at the fit view, at the
     closest zoom, empty, crowded, in dark) become the section's checklist.
2. **Configuration tables are checked, not kept by hand.** (C; S)
   - *What happened*: the first slice's rule (one Configuration table a `_spec`) held as a rule but not as a fact: SR
     found 17 missing or stale rows (the core's `rows.*`, `imageBytes`, heights for rows and logo, the viewer's
     `ui.canvas`, `ui.history`, `bounds.themeDepth`, no table for `renderSvg`'s options), and six items the designs name
     that the code never had (`rows.clades.cutLine`, `reveal.opensClades`, `navigator.idle`, `rows.reduce`, the
     `figure` preset, a `seqx-examples` table).
   - *Proposal*: a check in `.seq/_checks/` reads `DefaultRenderConfig`'s and `DefaultViewerConfig`'s keys (flattened)
     and fails on a key with no row in its `_spec` table, or a row with no key; a design's Configuration row names a
     key the code has, or says "not built".
3. **A drift note is a gap: it is asked at the next gate, never parked for the close.** (P: amends "Gaps are asked at
   every gate"; L)
   - *What happened*: SR's and EO's design questions (`drift/finalize-config.md`, `drift/ec-contradictions.md`) were
     written "for the next gate" and reached this close unasked; `gaps.md`'s "Open" heading stayed after GD had built
     every item under it.
   - *Proposal*: each gate's message to the author lists the drift notes' open items beside the gaps, each a yes or no
     with the recommendation; a note's item leaves the note when asked or resolved.
4. **An edge case a lane builds is tested by name in that lane.** (P: extends §whole-composition-tested; C)
   - *What happened*: EO reconstructed 56 outcome rows by reading tests for their behavior (12 kept with a stated gap,
     9 not met); test-only lanes ET and AT then added seven tests and found one edge case not wired (EC-135) and two
     accessibility faults. Tests that cite their EC (`(EC-43)`) were found by grep in seconds.
   - *Proposal*: a lane's test of a catalog item names it in its title; `.seq/_checks/ec-ids.ts` also lists the items a
     built feature's row covers (`features.md`, Cases) that no test names, as the close's worklist.
5. **The Features page moves with the merges.** (P: amends Gate 0's "The Storybook's Design pages follow the designs";
   L)
   - *What happened*: at this close the Design page "Features" still listed under Planned almost everything the slice
     built, and `features.md` had 41 stale rows; the Design pages of designs were kept in step lane by lane.
   - *Proposal*: a lane that builds a feature moves its row (`features.md`) and its line on the Features page to built,
     with its story's link, in its own change; the story-link test already checks the link.
6. **The findings table is kept for every gate and lane, not only the early stages.** (L)
   - *What happened*: `state.md`'s table stops at stage 4; stages 5 and 6 (about forty lanes) recorded their review
     findings in prose ("review findings fixed"), so this close cannot count them.
   - *Proposal*: the loop file's merge step adds the lane's row (bug, missing design, look) before the merge note.
7. **The quick lane's time is the lane's to keep.** (L; the repository's "A test slow by nature moves")
   - *What happened*: three lanes (QT, QT2, QT3) moved tests after gates found packages over the 2 s guide; each moved
     what earlier lanes had added.
   - *Proposal*: a lane's checks include its packages' quick-lane time against the guide; a lane over it moves its own
     slow tests (a quick representative kept) before its review.
8. **Parallel lanes need their own worktrees, reviewers too, and the lane count is the author's.** (L)
   - *What happened*: PG's reviewer switched NV's worktree (restored); the author set the lane count five times
     ("1heavy+2", "1+1", one at a time, a pause for usage, "1+1-2"); `resume.md` carried a hand-over to a new session
     without loss.
   - *Proposal*: every reviewer works in a worktree of its own; the loop file holds the author's current lane limit as
     one line, and `resume.md`'s shape is the template for any pause.
9. **A seeded generator is checked before fixtures depend on it.** (S: base-library; C)
   - *What happened*: CS found core's `PCG` biased (G 33%, C 20% over 1 Mb); the fix (PCG32) changed Mol\*'s distinct
     colors and a stored snapshot in `@molstar/model`, and broke the rows pyramid's stated tolerance (decided: small
     blocks exact).
   - *Proposal*: core's seeded generators carry a reference-vector test against their published algorithm; a
     generated dataset states the statistics its tests assert.

## C. The earlier fold-backs' rules, checked against this slice

From the first slice:

1. **The look is part of the lane's review, before the merge** — *held.* Looks found 17 faults across the stage 0–4
   lanes that draw (PG 5, TH 3, NV 2, RB 2, TV 2, C0, AX, SL 1 each) and two at a gate; RB's look also found a test changed to fit the
   code. No change.
2. **A fix's test is shown failing without the fix** — *held; amend for tests changed.* RB's review found a test
   bent to the code (tracks not wrapping on a cyclic axis) and restored it to the design; PG's own look-fix regressed
   and was reverted. **Proposal (P)**: a lane that changes an existing test's expectation says why in its report, and
   the reviewer checks it against the design.
3. **An observation a fix rests on states its conditions** — *held.* NV's blank pictures were at dpr 2 only; RW's look
   was said to be light only; the heavy bounds' failures were tied to load and reproduced alone. No change.
4. **Defaults are confirmed from tables** (amended by the second slice) — *amended again* by B1 and B2: defaults are
   designed with a picture, and the tables checked mechanically.
5. **A merge names the review findings refuted or left unfixed** — *held where recorded* (PG's A7 to GD, DC's public
   API accepted, RX's lows); the prose records of stages 5–6 seldom say "none refuted" (B6). No new evidence for the
   package's "proportional review" tension: the reviews were single and independent, and no merge names a refuted
   finding.
6. **A performance target names its path and its test, and is measured early** — *held.* The allocation tier ran with
   each lane; TA bisected a 42.6 → 72.4 KB pan regression to TH within the stage. Stage 7 still found the viewer's
   emphasis frame at 800 B (lane EA): its quick representative is in the core only. **Proposal (S)**: the viewer's warm
   emphasis frame gets a quick representative in the Storybook's perf tier at the stage that builds it.
7. **A state derived each frame publishes on change only** — *held, extended.* DS made the open ring group an input
   from the hover, not a state set by a draw (`32108c1e5`). No change.
8. **A story composing more than one consumer has a headless twin** — *held.* Twins exist for every dataset's composed
   stories; the novice still found the SVG twin a step behind on two stories (a story's wiring, not a twin's). No change.

From the second slice:

1. **The self-check runs on every design when it is agreed** — *held for interaction, insufficient for defaults.* The
   designs were self-checked as a set; the author's use found few missing inverses this time, and many missing
   defaults. Amended by B1 (a Defaults section the gate reviews).
2. **Interaction is a matrix** — *held; amend.* The plan named its cells; the cursor's rows were marked built in the
   design while column-only (NP found it; RK built it). **Proposal (L)**: a gate checks a cell built by the test that
   names its gesture, not by the design's word.
3. **Precedents cover interactions; a novice pass precedes the author's use; an early use gate** — *held for the
   novice pass, dropped for the early use gate.* The novice pass found 21 (bug 5, missing design 8, look 8) before the
   author; the early use gate was deferred twice and replaced by the author's use as lanes merged, which found what it
   would have. **Proposal (P)**: the use gate is continuous use counted by kind, with a look at what changed before
   the close; no fixed early gate.
4. **What is drawn is checked against the data** — *held.* TE's tile-edge item and the author's broken rail labels at
   insertions were both data-against-drawing faults, found and tested. No change.
5. **A restriction in a `_spec` states its reason** — *held.* RX found the rows rail idling on `whole` alone against
   the design; the `_spec` corrected with its reason. No change.
6. **Bounds are the best of three, and heavy runs never overlap** — *held; amend.* Bounds still failed in a full heavy
   run under load and passed alone (rows-time, pairs-grids-time); the author capped the heavy lane at two workers.
   **Proposal (P)**: the heavy lane runs alone and at two workers (the cap recorded with its measurement).
7. **Gaps are asked at every gate; a default is confirmed by use** — *amended* by B3 (drift notes are gaps). Defaults
   were confirmed by use throughout (the author decided eleven while using the views), as the rule says.
8. **Sizes are measured at every gate** — *held as a measurement, not as a bound.* Every gate measured; the bundles grew
   from 150.6 / 267.9 / 507.8 kB to 175.4 / 712.0 / 1142.4 kB (the minimal viewer 853.1; `measurements.md`), the author calling them guides; SF's `--files` attributes growth. No
   change; the targets are reset after RR (plan, stage 7).
9. **A fix of an unreachable state says so** — *no evidence.* No change.

From the third slice's Gate 0:

- **The Storybook's Design pages follow the designs** — *held for the designs' pages, missed for Features* (B5).

From the second slice's close:

- **The fold-back is a step of every close** — *held* (this file).

## D. Not proposed

- Rules from CE's parking: a lane whose design moved under it was parked with its review in `followups.md`; the loop
  handled it as written.
- Rules from the lockfile and hot files: OP's hot files were applied by the orchestrator after the merge, as before.

## Self-check (§self-check), on the claim "these proposals make the principles truer to what was built"

1. **Assume it holds.** B1 moves the defaults' decisions to the design gate, where the author's use found them; B2 and
   B4 turn two hand-kept lists into checks; B3, B5 and B6 close three places where state went stale between gates; C
   amends four rules by this slice's evidence.
2. **What it rests on, negated.** It rests on the defaults findings being a kind that recurs. Negated: this was the
   slice of figures, where defaults are the product; the fourth slice (integration: molecules, structure syncing,
   ports) adds few pictures, so a Defaults section is ceremony there, and its findings will be about protocols and
   APIs.
3. **The negated world.** It explains why the second slice's use found gestures and this one found looks: each slice's
   use finds what that slice made visible first. It predicts the fourth slice's use findings will be about what a host
   must wire (ports, adapters, the fake 3D viewer's syncing), which no picture shows and a novice pass on a host's
   task (not a reader's) would find.
4. **Both held.**
   - *Integrates*: B2, B4, B5 and B6 stand under either reading (they keep lists true, whatever the slice builds). B1
     is scoped to designs that add a layer, a layout or a part, so it costs nothing where none is added.
   - *Carried*: the novice pass's task script follows the slice: for the fourth, a host's integration task ("sync a
     selection with the 3D view") instead of a reader's; counted at its close.

## Questions for the author left by this close

The two drift notes kept in `drift/`: `finalize-config.md` (items 2–4: drop or keep planned the configuration the
designs name and the code lacks) and `ec-contradictions.md` (EC-63 slices against glyphs; EC-138's outline against
saturated fills), each with its recommendation.

Decided by the author:
