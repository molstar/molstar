# Fold-back: what the fourth slice taught (draft, for the author)

The fourth slice built version 1's third milestone, integration, in 35 merged lanes (AR to HT; the list and each
lane's findings in `state.md`, "Findings by lane"), each through a brief, an independent review and a fix round, the
orchestrator's look before a merge that draws, and gates after stages 1, 2 and 4. It was the first slice to run the
third slice's rules: Defaults sections, checked Configuration tables, drift notes as gaps, edge cases named by tests,
the Features page moving with the merges, the quick lane's time kept by the lane, §module-size, sizes as bounds and
continuous use. It ran overnight under "Running while the author is away" and beside the fifth slice's stage 0. This
file follows the principles' rule "The fold-back is a step of every close". **Nothing here is applied before the
author's line at the end.** Each proposed change is one line the author decides: "Proposed: … Why: … Where: …".

Where to: **P** principles (the instance's adaptations), **D** a design, **S** a package `_spec`, **C** a check, **L**
the next loop file.

## A. What each gate, lane and use found, by kind

From `state.md` (the findings table, 35 rows), `novice-findings.md`, `fixes.md` and this close's reconcile.

| Gate, lane or use | Bug | Missing design | Look | Other |
|---|---|---|---|---|
| Lane reviews (35 lanes; a strong reviewer in its own worktree) | 126 | 35 (with the lanes' own drift; decided by the author or the designs) | 2 | — |
| The orchestrator's check and look at merge (HP2 2, SV 2, FX5 1, PC 1, HT 3) | 9 | 0 | 16 | — |
| The author's use, this slice (2026-11-21: FX5, FX6, HT) | 8 (a submenu closing on the way, FX6's five, a submenu replacing its menu, the lollipops' missing hover) | 1 (hiding a track moves the view: HT) | 1 (the PAE's colours) | — |
| The author's check of the third slice, carried (fixes.md, FX1–FX3) | 3 | 4 (the ellipsis rule, a long card, the overview's double-click, submenus from a select) | 2 | — |
| Stage gates (1, 2, 4) | 0 | 0 | 0 | 3 quick-lane folders over 2 s (stage 1); a ratio bound missed once at load 4.9 (stage 4); 2 allocation cases red on the base (before FX5) |
| Novice pass (29 findings) | 5 (2 not bugs: F5 as designed, F7 not reproduced) | 11 | 13 | — |
| Edge-case outcomes (`ec-outcomes-draft.md`) | — | — | — | not met 5, changed 2, 8 kept without their number in a test |
| Finalize's reconcile (this close) | 0 | 0 | 0 | 28 stale rows in `features.md`, 6 groups on the Features page, 7 stale passages in track-header-parts, 4 precedent rows in other designs |

Read: reviews found the bugs (126 of 151, in 33 of 35 lanes; the two pure splits, SZ, none); the orchestrator's
looks found the lanes' looks (16 of 18); **seven of the author's eight bugs were interactions in a menu, the eighth a
hover**, all after lanes and the orchestrator had looked at screenshots (B1).

## B. What this slice taught

1. **A look is a gesture, not a screenshot.** (L; D: ux-foundations' matrix)
   - *What happened*: the header menu passed HP2's, HP3's and FX5's reviews, their jsdom twins and the screenshots of
     every look; the author then found in minutes a submenu closing on the way to it (FX5), the menu's items not
     applied, its click and hover reaching the view, `Color ▸` closed by its own preview and `⋯` opening nothing (FX6),
     and a submenu replacing its menu (HT). FX6's cause: the menu mounted in the view's container, a press starting a
     pan's pointer capture, which no synthetic event in jsdom reproduces. The novice pass opened menus (F17, F18) and
     missed them too: its browser tool cannot send a pointer's path.
   - Proposed: the look at a merge that adds or changes an interactive part drives each gesture of its matrix row in
     the visible pane with real pointer events (hover, press and release, a path onto a submenu, a press outside, Esc),
     and records what each did; screenshots show the result, never stand for the gesture. Why: 8 of 8 author bugs were
     gestures the screenshots could not show. Where: loop.md step 5 ("look"), the brief's matrix cells.
   - Proposed: a part that mounts DOM over the canvas has a jsdom test that the view receives none of its pointer
     events (`isViewEvent`, FX6). Why: the cause of five of the author's bugs. Where: S (viewer `_spec`, Input), the
     brief template.

2. **The orchestrator looks itself, every merge.** (L)
   - *What happened*: FX6 was merged on the lane's own look; the author's next use found HT's first bug, and the
     author's correction: "the merge-time look by the orchestrator, every time" (HT's row). The loop already said the
     orchestrator looks; the step was skipped under time.
   - Proposed: the findings row names who looked (orchestrator, with the stories and themes seen) and a merge with no
     orchestrator's look is not merged. Why: the rule existed and was skipped once, and that merge is where the author
     found a bug. Where: loop.md step 5 ("record"), the findings table's columns.

3. **The lane limit counts every running agent.** (L)
   - *What happened*: on 2026-11-21 SV, FX5 and PC ran together (merged 09:12, 09:17, 09:23) against "Two lanes at
     once", with the fifth slice's reviews beside them; overnight the fourth and fifth slices shared the limit by
     hand. The heavier load showed in the timed checks (B6).
   - Proposed: before starting a lane the orchestrator lists what runs (lanes, reviewers, fix rounds, another slice's
     draft lanes) in `state.md`'s In progress with start times, and starts nothing past the limit. Why: the limit was
     passed without anyone deciding it. Where: loop.md, "Lanes at once" and step 2.

4. **A pane's screenshot can be one step stale.** (L)
   - *What happened*: HP3's look and the novice pass both found the Storybook pane repainting only after a mouse move,
     so a screenshot right after an action showed the state before it; HP3 checked by the canvas's pixels and the DOM
     instead.
   - Proposed: after an action, move the pointer off the target and take a second screenshot; a look's evidence is two
     agreeing screenshots or one plus the DOM or canvas pixels. Why: a stale picture is a false pass or a false bug.
     Where: loop.md step 5 ("look"), the novice brief.

5. **The timed checks need a quiet machine, and the loop must make one.** (P: amends "The heavy lane runs alone";
   L)
   - *What happened*: stage 2's heavy run at load 3.7 and stage 4's at 2.6 and 4.9 (the chromosome view's pan ratio
     1.97× against 1.9× once); AC's core folder 2.28 s at load 8.9; a services story test failed once under load (FX6)
     and Prism's `commands` and `collapse` tests once each in a combined run (HT), all green alone. Two allocation
     cases (`headers.quick`, core's grid hit) were red on the base before FX5, found by FX5's run, not by a gate.
   - Proposed: a gate's timed runs (the quick lane's time, heavy, the allocation tiers) start only when no lane runs
     and `uptime` is under 3, recorded beside the numbers; a test that fails once under load is rerun alone and listed
     in the TODO's watch list, a second failure is a fix lane's. Why: every timing miss this slice was load. Where:
     loop.md, Gates.
   - Proposed: the orchestrator's check runs the allocation tier for a lane that touches a frame's path (its brief
     says so), so a red case is bisected to the lane that made it, not found on the base. Why: two cases went red
     between gates unattributed. Where: loop.md step 5 ("check").

6. **Defaults sections were a rule, not a practice.** (P: §defaults-designed; C)
   - *What happened*: no design in `.seq/design/` has a Defaults section (a heading search at this close), the stage 0
     designs and the amendments since included (track-header-parts' parts, the axis's connector, extension-openings'
     fallback). Two of the author's findings were defaults: the hidden-tracks chip moving the view (HT, a new design)
     and the PAE's colours (FX5).
   - Proposed: a check that a design whose Consequences name a layer, a layout or a part has a `Defaults` heading,
     with the existing designs listed in a baseline that only shrinks. Why: the third slice's rule held nowhere without
     a check, as its Configuration tables did not. Where: C (`.seq/_checks/`), §defaults-designed.

7. **`features.md` did not move with the merges; the Design page did.** (C; L)
   - *What happened*: lanes moved the Storybook Features page (the structure link, host data, track headers, legends,
     extensions) but this close still found 28 rows of `features.md` stale (Archives, Extension descriptors, Host data,
     A second consumer "planned"; Map kinds "offset and codon"), the same miss as the third slice's B5.
   - Proposed: a check fails a `features.md` row whose status says "planned" or "to come" for a lane `plan.md` marks
     merged (the lane's code in the row's status or the plan's line naming the row). Why: the rule failed twice by
     hand. Where: C (`ec-ids.ts` or a new check), loop.md step 5 ("features").

8. **An edge case built is not always named.** (P: §ec-tested-by-name)
   - *What happened*: EX built EC-235 (the fallback drawn and said) and its tests never cite it; seven other items are
     held by tests without their number (EC-36, 51, 110, 132, 175, 207, 238–244). The worklist cannot see them.
   - Proposed: the review checks each `EC-<n>` the brief names against the lane's test titles, and a missing one is a
     review finding. Why: the brief named EC-235, the review did not ask. Where: loop.md step 5 ("review").

9. **Independent reviews found real bugs in nearly every lane.** (P: held; L)
   - *What happened*: 126 review bugs in 33 of 35 lanes, from crafted inputs exhausting memory (MC, BW) to tests that
     could not fail (EX, HP3); the two lanes without one were pure splits (SZ), whose reviews found none.
   - Proposed: keep the strong reviewer for every lane but a pure split, which gets the fast model with the "no
     behavior change" check. Why: SZ's two reviews cost a strong run each and found nothing. Where: loop.md step 5.

10. **The novice pass came before the last lanes that draw.** (P: §use-continuous)
    - *What happened*: it ran at stage 4's gate; SV (services, skipped), FX5, PC, FX6 and HT merged after it, and the
      author's eight bugs were in what changed after. Its 29 findings were real (3 of 5 bugs; F1 and F22 the hover
      card again).
    - Proposed: the novice pass runs after the last lane that draws, its script including a pointer path onto each
      menu and submenu. Why: what changed after it went unseen by it. Where: the plan template, the novice brief.

11. **The slice reached finalize with two plan items unticked.** (L)
    - *What happened*: stage 3's gate and M4 are not ticked in `plan.md`; the author's "run slice 4 loop all the way to
      the close" moved the loop on after the use gates passed.
    - Proposed: Finalize starts with the unticked plan items listed, each run, cut by the author's word, or carried
      to the next plan, in `state.md`. Why: "When every item is ticked" was not true and nothing said so. Where:
      loop.md, Finalize.

## C. The earlier fold-backs' rules, checked against this slice

From the third slice:

1. **Defaults are designed with a picture** — *not applied* (B6). Amend by a check.
2. **Configuration tables are checked** — *held.* Both config-table tests green at this close; every new item
   (`hover.cardItems`, `labels.minChars`, `headers.submenuRestMs`, `navigator.rowsRail`, `ranges.connector`, …) had
   its row when merged. No change.
3. **A drift note is a gap** — *held.* AC's, ST1's, MX's and HD's drifts were asked at the next gate and decided
   (`state.md`, Decided); `drift/` is empty at the close. No change.
4. **An edge case is tested by name** — *held; amend* (B8).
5. **The Features page moves with the merges** — *held for the page, missed for `features.md`* (B7).
6. **A seeded generator is checked** — *no evidence* (no new generator). No change.
7. **A changed test expectation is explained** — *no evidence recorded.* No change.
8. **A matrix cell is built when a test names its gesture** — *held for the reviews, insufficient for the look* (B1):
   the menu's cells had tests and still failed under real pointers.
9. **Use is continuous** — *held.* The author's use passed the use gates and found 8 bugs, 1 missing design, 1 look as
   lanes merged; the novice script was fitted (a task a story); amended by B10 for its timing.
10. **The heavy lane runs alone, at two workers** — *held at two workers; its quiet machine not held* (B5).
11. **A hot loop's quick representative in every package** — *held.* FX5's review found the warm hover frame
    allocating through `items()` by the Storybook tier; amended by B5 for when the tier runs.
12. **§module-size** — *held.* SZ split five files (2871 → 865 the largest), the baseline is empty; ML2's and MK's
    `model.ts` (1410, 1483) and FX3's `resolve/model.ts` (1433) were raised by reviews and kept under 1500. At this close 18 files are over 800 lines, `render/navigator.ts` at
    1473 (CN's pictures) the nearest the bound: a split before the fifth slice's navigators grow it. No change to the
    rule.
13. **Sizes are bounds** — *held.* Four growths came to the author with the feature behind them (QE, MR, ST2, EX) and
    each bound moved by decision (230 / 830 / 1330 / 940 kB). No change.
14. **The quick lane's time is the lane's** — *held in part.* FX3 still fixed three folders after stage 1's gate; the
    viewer folder stays about 2.8 s by the author's budget decision. No change.
15. **Reviewers in their own worktrees; the lane limit is the author's** — *held for worktrees; the limit passed*
    (B3).

From the first and second slices: the look before a merge (*held in text, skipped once*, B2); a fix's test failing
without the fix (*held*: reviews checked it, e.g. FX5's and HT's); an observation states its conditions (*held*: the
load beside each timing miss); a merge names refuted findings (*held*: the merge commits and rows name them); a
performance target names its path and test (*held*: the link's hover bound of 2 kB a lit unit set by decision);
interaction is a matrix (*held*, B1 amends the look); the fold-back is a step of every close (*held*: this file).

## D. Not proposed

- The overnight runs' gaps taken as recommended (`gaps.md`): every one was accepted by the author the next morning; the
  loop's "Running while the author is away" held as written.
- Public API changes under the author's standing permission: listed in `gaps.md` per lane and accepted; no rule.
- The fifth slice's drafts beside this slice's lanes: covered by B3's count.

## Self-check (§self-check), on the claim "these proposals make the loop truer to what was built"

1. **Assume it holds.** B1, B2 and B4 make the look a test of gestures done by the orchestrator on fresh pictures,
   where the author found every bug; B3 and B5 make the limit and the quiet machine facts the loop checks; B6 and B7
   turn two hand-kept rules into checks, as the third slice did for Configuration tables; B8, B10 and B11 close gaps
   in the review, the novice pass's timing and the finalize's entry.
2. **What it rests on, negated.** It rests on the menu bugs being a kind that a scripted look finds. Negated: they
   were found because the author used the view as a person does, with intent and mistakes; a scripted gesture list
   finds only what the matrix names, and the menu's failures (a pointer's diagonal path, a press that starts a pan)
   were not in any matrix row.
3. **The negated world.** It explains why the novice pass, which did act with intent, still missed them (its tool
   could not send a path), and predicts that the fifth slice's author use will again find interaction faults off the
   matrix, in whatever new part it adds (a genome browser's track menus, a ring's drag).
4. **Both held.**
   - *Integrates*: B1's gesture list covers the matrix's named cells and makes the look cheaper to repeat; B10 moves
     the novice pass after the last drawing lane with pointer paths, which is the closest stand-in for use with intent.
     B2–B8 and B11 stand under either reading (they keep the loop's own rules true).
   - *Carried*: the third slice predicted this slice's use would find what a host must wire; it found menu
     interactions instead (the slice still added visible parts). Use's findings follow what a slice makes visible,
     not what it is about: the fifth slice's plan names its interactive parts and gives each a gesture look.

## Questions for the author left by this close

1. Stage 3's gate and M4 are unticked in `plan.md` (B11): run them in the finalize's gate, or cut M4's browser
   measurements to the fifth slice?
2. The edge-case outcomes (`ec-outcomes-draft.md`): paste as drafted into `design/edge-cases.md`?

Decided by the author, 2026-11-21: merged ("All 11"). B1–B5 and B7–B11 applied in the fifth slice's `loop.md`; B1's
overlay rule in the viewer `_spec` (Input, "Overlays own their events", FX6); B6's and B7's checks built by the fifth
slice's first fix lane (its plan). Questions: stage 3's gate and M4 were run before the close (M4's bigWig pan, real-vsync
frames and the split patch carried to the fifth slice); the edge-case outcomes pasted into `design/edge-cases.md`.
