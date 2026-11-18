# Fold-back: what the second slice taught (draft, for the author)

The second slice built sixteen lanes (S0, M, L, S0b, R, K, K2, R+K, O, K3, N, Q, P, S, T, U) through briefs, an
independent review and one fix round each, with a look in the browser before every merge that draws, and gates after
each stage; 238 commits from `fc2bba6ce`. This file follows the principles' rule "The fold-back is a step of every
close": what the slice taught, the first slice's fold-back rules checked against this slice's evidence, the self-check,
and then the author's decision on its own line. **Nothing here is applied outside this folder before that line.**

Where to: **P** principles (the instance's adaptations, and the oz package's home), **D** a design, **S** a package
`_spec`, **C** a check, **L** the next loop file.

## A. The first slice's rules, checked against this slice

1. **The look is part of the lane's review, before the merge** — *held, and the tension it carried can be released.*
   - Evidence: the look before the merge found bugs that tests and the review had passed, in four of the five lanes
     that draw: O (a bookmark's focus did not move the view), P (a drag became a one-cell click), S (five small layout
     faults), U (a double-click on a domain zoomed to 17 nucleotides, from before the lane). It also found a renderer
     bug from the first slice: WRAP53's second exon never drawn (`b73d3199b`).
   - The carried tension "look before merge, against the cost of a look per lane" said: if later slices find no visual
     bug before the merge, the look was a first-slice cost. This slice found them in most lanes. **Proposal (P)**:
     release the tension with that count.
2. **A fix's test is shown failing without the fix** — *held; amend for defensive bounds.* Every fix this slice named
   the command that fails on its parent (WRAP53: the layout test; T: five tests; U: nine). One fix could have none: a
   depth bound on a recursion a spec cannot reach (T, finding 5). **Proposal (P)**: a fix that guards an unreachable
   state says so in its report and its commit, instead of a failing test.
3. **An observation a fix rests on states its conditions** — *held.* The far-zoom blink and the missing exon were each
   reproduced headless before the fix (T's `far-hover.test.ts`, the diagnosis of WRAP53 with `delayMs: 0`); every
   look was taken in a visible pane. No change.
4. **Defaults are confirmed from tables, the unsure ones first** — *not followed through; amend.* The gaps list stayed
   short (nine items), but seven of them were never put to the author at a gate: lane M's and L's since stage 1. The
   author decided gestures and looks when they saw them ("the gestures are configurable anyway"), and never asked for
   a table. **Proposal (P, L)**: each gate's message to the author ends with the gaps raised since the last gate, each
   a yes/no with the recommendation; the finalize step may not carry an unasked gap. See also C below on the carried
   tension.
5. **A merge names the review findings refuted or left unfixed, and why** — *held.* T's and U's merges name every
   finding and that none was refuted; S's merge message did not (a miss of the orchestrator). Evidence for the oz
   tension "proportional review": in stage 4's three reviews (S, T, U: 3 + 5 + 8 findings), one reviewer each, every
   finding cited a line and every fixed one was confirmed by a test failing on the parent; none was refuted.
   **Proposal (P, the package's home)**: add this as later evidence to that tension: one reviewer whose findings each
   come with a failing test needed no skeptic.
6. **A performance target names its path and its test, and is measured early** — *held for time, missed for size.*
   Every lane ran its time representative (`itWithin`). The bundle sizes were measured only in stage 5: the read
   bundle grew 26 kB at stage 2 (R+K: the immutable collections reached through `lang`'s index) and was found three
   stages later. **Proposal (P, L)**: the size targets are part of every gate's checks (`pnpm -C packages/seq/core
   size`, about 10 s), not only the measure stage.
7. **A state derived each frame publishes on change only** — *no evidence either way* (no incident this slice). No
   change.
8. **A story that composes more than one consumer has a headless twin** — *held; amend.* S's twins and T's and U's
   composed the stories as the stories do. But no twin caught the two bugs the looks found in drawn data (WRAP53's
   second exon, the merged runs with no hit box): the twins asserted canvas equals SVG, and both were equally wrong.
   See B2.

## B. What this slice taught

1. **The author's use finds missing design, not bugs.** (P: §used-before-done; D)
   - *What happened*: stage 4's gate (the author using the stories) produced eight findings. Two were bugs (the
     far-zoom blink, boxes with no hit). Six were design the slice had not made: several ranges in one selection, a
     follower that leads, zooming out, a reversed follower, the protein and its gene on one axis, collapsed introns.
     The first three cost two extra lanes (T, U); the rest went to the third slice. Two of the six were in the
     designs' later columns; four were in no design at all (`letters` on another system's axis was refused in a
     `_spec` with no reason, a decision made in passing).
   - *Read again*: five of the six were a missing **inverse or symmetry** of something built: focus without its way
     back (zoom out, back), a follow in one direction only (a follower that leads), a selection that replaces with no
     add or remove (several ranges), a direction one way only (a reversed follower), a protein's data on its gene but
     not its letters (one axis). The sixth (collapsed introns) was planned but not built. The author: "some things are
     often under-designed. how can we better identify and surface these gaps?"
   - *Proposal* (the author asked for ways, then: "the .seq/oz/_core/recurrence.md asks for a negation (.. which sort
     of is inverse). should we apply the recurrence to the designs too? the "4 phase cycle" is generally quite good at
     identifying residuals"):
     1. **The self-check runs on every design when it is agreed, not only at a fold-back** (P: §design-first,
        §self-check; the package's home: recurrence.md "Checking the loop against itself" and the
        `principles-self-check` skill name the design gate). The claim is the design's own ("this design covers what a
        user does with X"); step 2 negates what it rests on; step 4 builds what integrates or plans it with a slice,
        and keeps the rest in the design as a standing tension (§surface-tensions), visible. Run on this slice's
        designs it surfaces all five: focus rests on "the user wants to go to X" (negated: to leave it: zoom out,
        back); follow on "the user moves the leader" (the follower leads); selection on "one range" (several); the
        codon map on "genomic order is how a gene is read" (along its protein: reversed); protein on its gene on
        "separate sections" (one axis). For an action the negation is usually its inverse or its symmetric case, so
        no separate inverse list.
     2. **An interaction matrix in the UX design** (D: `ux-foundations.md`): rows the things a user acts on (a
        section, a follower, an item, a range, the selection, the view), columns the verbs (select, extend, add,
        remove, focus, back, zoom in and out, pan, hide, edit), each cell the gesture and command or "no, because".
        Where the self-checks' results for interaction land; the slice plan names the cells it fills, the gate checks
        them.
     3. **A precedent capability check and a novice pass before the author's gate** (P: §precedents, §used-before-done;
        L). §precedents is applied to interactions too: what a user can do in the precedents (the current Seq\*, IGV,
        Jalview, the genome browsers) that the slice's views cannot is listed, each built, planned or refused. Then a
        light agent plays a user in a visible browser with a short task script per view ("select two regions",
        "zoom out after focusing", "move from the gene") and reports friction before the author sees it.
     The use gate also moves earlier: a slice that builds interaction has a short one after its first drawn stage.
     The finalize step counts what each gate found, by kind (bug, missing design, look), so a fold-back can tell
     which of the three catches what.
2. **Check what is drawn against the data, not only one consumer against the other.** (P: §whole-composition-tested;
   S; C)
   - *What happened*: WRAP53's second exon was dropped by the renderer since the first slice (a tile-edge dedupe by id
     alone); merged runs at the fit view had no hit box (U). Canvas, SVG and the hit index agreed with each other, so
     every twin passed. The look found them by comparing the picture with what the data says (a transcript of two
     exons showed one).
   - *Proposal*: the examples' tests carry a data-to-drawing check: for each ranges layer at a few viewports (the fit
     view, a zoom across a tile edge), every item of the data in view has its marks drawn and each drawn mark is hit.
     U added the hit half as a property (`hits.test.ts`); the exon half is one assertion per view.
3. **A restriction in a `_spec` states its reason, and lanes question it.** (P: "Question the specs"; S)
   - *What happened*: "letters are drawn on their own system only" stood since the first slice with no reason; the
     author found the feature it blocked. The repository rule says to question specs with no reason recorded; no lane
     did.
   - *Proposal*: a `_spec` rule that refuses or restricts (an issue code, "only", "never") names its reason; the
     reviewer's checklist asks for it on the lines a lane adds. Applied already for `letters` (`4e026805b`).
4. **Timing bounds in the heavy lane need the machine quiet or a best-of.** (P: §performance-tested; L)
   - *What happened*: single-sample bounds failed under load four times this slice (load average 6.8 during a lane's
     heavy run; a 10-op patch at 7.2 ms against 5); each passed on rerun. `bestOf` (`test-harness`) was added and used
     in two suites.
   - *Proposal*: `itWithin` bounds take the best of three by default; a lane does not run the heavy lane while another
     heavy run is going (the loop's parallelism already says two heavy lanes at most: say "never two heavy runs").
5. **An entry's index re-exports only what its readers need.** (S; C)
   - *What happened*: `@molstar/seqx-core/lang`'s index re-exported the editing validator, and with no
     `"sideEffects": false` every reader of the language paid 26 kB for editing. Split as its own entry (`c49294ff7`).
   - *Proposal*: in the core `_spec`, the entries rule: an entry a reader imports holds no module only an editor or a
     viewer needs; the size check at each gate (A6) is its guard.
6. **The orchestrator's own edits need the same care as a lane's.** (L)
   - *What happened*: an over-wide edit of `state.md` removed three of its sections (restored from git); an unquoted
     heredoc ran a backtick span as a command and dropped a config name from a brief (`4c95e4ba9`).
   - *Proposal*: the loop file says the orchestrator writes briefs and work files with the editor (or quoted heredocs),
     and diffs a work file after a scripted edit. No principle: it is tool hygiene.
7. **The fold-back is a step of every close.** Already applied by the author's rule this slice (principles, "The
   fold-back is a step of every close"; `_checks/closed-work.ts`). Listed for the record; nothing more to apply.

## C. The carried tensions

- **"Configurable over fixed" against a config surface the author can confirm** (carried from the first slice).
  This slice added about twenty config items (`selection.extendScope`, `follow.followerMoves`, `view.historyLimit`,
  `view.backOnFocused`, `view.jumpMs`, `view.motion`, `view.keys`, `navigate.doubleClickZoomStep`, `story.*`, the
  emphasis items, the bounds). The author confirmed defaults only by using them (the selection's look, the gestures,
  the follow), never from a table, and said so ("the gestures are configurable anyway"). **Proposal (P)**: settle it
  by evidence: a default is confirmed by use at a use gate, or stands as documented and tested; the gaps list carries
  only the defaults a lane was unsure of **and** that change a picture or a gesture. Carry what does not integrate:
  defaults no story shows (store sizes, bounds) are never confirmed by use.
- **"Look before merge" against its cost**: release (A1).

## Self-check (§self-check), on the claim "these proposals make the principles truer to what was built"

1. **Assume it holds.** B1 moves the author's use earlier; B2 and B3 turn the looks' and the author's finds into
   checks; A4 and A6 turn two rules that were followed only in part into gate steps; C closes a tension by the
   author's own behavior.
2. **What it rests on, negated.** It rests on the second slice's finds being of kinds that recur. Negated: this slice
   was the first with interaction (selection, editing, navigation), so missing gestures were bound to surface once;
   the third slice (figures, layouts) will find layout bugs, not missing gestures, and an early use gate only adds a
   stop.
3. **The negated world.** It explains why six of eight use findings were gestures and navigation, and why they came
   all at once at stage 4: the first stories anyone could play with were built there. It also predicts that the third
   slice's use findings will be about how a figure reads (labels, wrapping, legends), which a look can catch and a
   lane's twin cannot.
4. **Both held.**
   - *Integrates*: B2 (data against drawing) and B3 (a restriction's reason) stand under either reading: both found
     faults from the first slice that no interaction exposed. A6 (size at every gate) is cheap and its one miss was
     real. A1's release is a verdict with its count (four of five drawing lanes).
   - *Carried*: B1 under the negated reading: whether an early use gate pays when a slice adds figures rather than
     interaction. Kept as a proposal for the third slice's plan (a use gate after its first layout), and counted at its
     close: what the early gate found, by kind.

## Gaps for the author (§gaps-queue, from gaps.md; asked now, per A4)

**Decided** (the author, 2026-11-18): "Nine gaps: ok as recommended".

1. Lane M: `SizeLimitExceeded` exported from `@molstar/core/util/data-source` (a typed error). Recommendation: keep.
2. Lane M: Seq\*'s asset policy `redirect` (now an option, `f6f012846`). Recommendation: keep as built.
3. Lane L: MessagePack extension codes `0x10` (special numbers), `0x19` (category), `0x1a` (bigint) beyond the agreed
   `0x11`–`0x18`. Recommendation: confirm.
4. Lane L: a one-line bug fix in io's `common/msgpack` (a leading U+FEFF kept), against "no change to io's MessagePack".
   Recommendation: keep (a bug fix, no API change, BinaryCIF shares it).
5. Lane K: `ChildIndexThreshold` 8, the patch language's id `patch` 1.0.0, a place written `:after id`/`:first`/`:last`.
   Recommendation: keep.
6. Lane R: `FlatThreshold` 16 (measured; the design said about 8); `Vector.indexOf` by `Object.is`. Recommendation:
   keep.
7. Lane N: `track.hidden`, the language staying 0.3.0 (unreleased). Recommendation: keep.
8. Lane P: the selection's look (`column`, a soft tint, one outline; the cursor after a key; the hover outside its
   mark), seen and shaped by the author. Recommendation: keep.
9. Lane U: a step zoom (-, =, Shift+double-click) keeps the steps forward, as a move by hand does. Recommendation: keep.

## Not proposed

- Rules from the `check:workspace` glass-wordmark failure: the host's, fixed by removing stale outputs.
- A rule from the frame-time flake: recorded as a followup; B4's best-of covers its kind.

## The author's answers so far (2026-11-18)

- The nine gaps: as recommended.
- B1: "yeah, this is real .. some things are often under-designed. how can we better identify and surface these
  gaps?" — the three ways above, for the author's yes.
- B2 to B6: "ok for other new lessons".
- A1 to A8 and C: not yet answered.

## Decision

(the author's line goes here: `Decided by the author, <date>: <merged | merged in part | not merged>: "<their words>"`)
