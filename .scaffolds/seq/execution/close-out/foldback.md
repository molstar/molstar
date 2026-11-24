# Fold-back: what the close-out taught (draft, for the author)

The close-out finished version 1 of Seq\* (next):
- 26 lanes, in five waves and the fix lanes FX25 to FX36: search, track groups, values marks, streamed patches, the
  drag rule, the novice questions' answers, `height: 'content'`, the legend's layers menu, the measurements' overs;
- three gates, the author's use gate, a novice pass in three groups, and the measurements with a re-measure;
- the language `0.10.0` to `0.11.0`.

It was the first close to run three lanes at once throughout, to split a novice pass by groups, and to measure with
microtasks drained.

**Nothing here is applied before the author's line at the end.** Each proposal is one line: "Proposed: … Why: … Where:
…".

| Code | Where |
|---|---|
| **P** | principles (the instance's adaptations) |
| **D** | a design |
| **C** | a check |
| **L** | the next loop file (and `_method/`) |

## A. What each gate, lane and use found, by kind

| Gate, lane or use | Bug | Missing design | Look | Other |
|---|---|---|---|---|
| Lane reviews and the orchestrator's looks (26 lanes) | 67 | 19 | 9 | size raises, each with its cause |
| The author's use, during the lanes | 5 (a canvas bridge gap, a legend row hiding a track, the index pages' font, a Manhattan link that threw, the legend menu's styling) | 2 (`height: 'content'`, Show and Hide all) | 1 | |
| Gates (1, 2, last) | 2 (a 16 B warm frame; a leak control mis-scaled) | | | 2 flakes fixed (FX33) |
| Novice pass (A, B, C) | 9 (2 not reproduced, 1 by design) | 12 (4 decided, D1 to D4) | 22 | |
| Measurements | 4 overs (3 fixed by FX36, 1 carried) | | | the method amended |

Read: the reviews and looks found most of the bugs, again with every lane's own checks green. The author's own use
found the bugs no test or review saw: each was on a canvas or a page in a browser.

## B. What this close taught

1. **A headless test of a gesture or a look can pass while the browser disagrees.** (L, P)
   - Proposed: a lane that changes a gesture, a menu or a DOM part tests it through the story's own composition with
     the browser's event order (pointerdown, pointerup, click, then dblclick), and asserts where the element mounts
     (its container, the UI root), not only its classes. The orchestrator's look in the pane stays the merge's
     condition for such a lane.
   - Why: FX35's double-click passed its tests and opened nothing on a label in the browser; FX32's class-equality
     test passed on its parent (the menu sat outside `.msui-root`); FX34's fold-then-undo twin case passed under the
     fake clock.
   - Where: `loop.md`, Each step (the look); the instance's testing adaptation.
2. **The measurements' method missed work.** (L, D)
   - Proposed: the measurement method drains microtasks after each flush and stubs pointer capture for drags; a
     report gives "old method | with microtasks" until the old numbers are retired. A story's SVG twin never exports
     inside an edit's frame (now FX36's code; the rule in `performance.md`).
   - Why: the fifth slice's drags likely measured no drag (synthetic `setPointerCapture` throws); the split's ~27 ms
     tail and search's 12 to 25 ms were the story's twin exporting in the edit's microtasks.
   - Where: `performance.md`, the measurements' method; the next plan's measurement brief.
3. **A hidden pane stops every agent that uses it.** (L)
   - Proposed: before starting an agent in the Browser pane, the orchestrator reads `tabs_context`; if hidden, it asks
     the author to bring the window forward. Pane agents run one at a time, each in a foreground tab of its own.
   - Why: the novice pass's first start (two groups) and the re-measure stopped on a hidden pane; parallel tabs
     throttle the background one.
   - Where: `loop.md`, the novice pass and the measurements.
4. **A test on a fake clock that steps a fixed time over real I/O is a flake.** (P)
   - Proposed: a story test on a fake clock waits on the viewer's `settle()` (a real turn a step), never a fixed step,
     when a dynamic import or a read is on its path.
   - Why: FX33: the Manhattan test read the domain before the lazy reference-lines chunk landed (failed three times in
     full runs); the ring's read-by-region case turned the clock with no real tick.
   - Where: the instance's testing adaptation.
5. **A brief's shape is checked when it is written, not at the merge.** (L)
   - Proposed: the orchestrator runs `features-planned --brief` on a brief before starting its lane.
   - Why: FX34's brief had no "**Items**" list; the check failed only at its record.
   - Where: `loop.md`, Each step (the brief).
6. **Bundle budgets ran out at the close.** (L)
   - Proposed: the next plan opens with a trim pass on the default and minimal viewers (by module, against the close's
     growth) before any lane that adds to them.
   - Why: raises at most feature merges (from the fifth slice's bounds: spec 265 to 268, headless 980 to 1,000,
     default 1,545 to 1,592, minimal 1,090 to 1,124; `performance.md` holds each cause); every bundle ends within 1.3
     to 4.8 kB of its bound, the default +51.0 kB and the minimal +34.0 kB over the close.
   - Where: the next plan.

## C. The earlier fold-backs' rules, checked against this work

- The heavy lane runs alone, at two workers: held (GATE 2 and the last gate, 187 and 188 green alone).
- Lane models: held (Sonnet for scoped and fix lanes, Opus for the perf lane and the measurements).
- Configuration tables: held (the checks green at every merge).
- Edge cases by name: held as a check, but `--worklist` held 56 built rows' items untested at the close (the
  reconcile's).
- The Features list moves with the merges: amended by the reconcile (rows the lanes did not move).
- The Storybook's Design pages follow the designs: held (`designs-in-specs` green; FX29, FX30, FX32, FX35 moved their
  pages).
- Interaction matrix, defaults designed: held.
- Mocks: not exercised (the 14 questions were decided in text by the re-plan; no new view).
- A fold-back at every close: this file.
- The novice pass is standing: held, its bugs confirmed in the pane before placing; amended by B3.
- The sixth fold-back's B5 (never `cd` into a lane's worktree): held.
- The fifth fold-back's B10 (the use gate the author's; interactions mocked at design): held for the use gate; the
  mocks not exercised.

## D. Not proposed

- The 200-copy split and merge over the edit bound: a design question for the next plan (lays shared between
  identical copies), not a rule.
- The old Seq\* viewer's label churn over its baseline: its owner's, outside Seq\* (next).
- The privacy plugin's refusals: the author's configuration, resolved by turning it off; no rule.
- The TODO's new lines (the novice pass's looks, the lesser overs, the square's whole genome after a cut): the next
  plan's.

## Self-check (§self-check), on the claim "these proposals make the loop truer to what was built"

1. **Assume it holds.** B1 and B3 close the browser gap (a gesture passes in the pane, the pane is there to look); B2
   makes the numbers mean the work; B4 removes a class of flake; B5 and B6 move two late failures earlier.
2. **What it rests on, negated.** The browser gaps were chance, the headless suite already catches the class, and the
   old measurements were right.
3. **The negated world.** It explains FX29's and FX30's clean headless runs; it does not explain FX35's label (no test
   hit outside the bands), FX32's parent pass, or the fifth slice's 0.3 ms release that changed nothing.
4. **Both held.** Headless tests stay the lanes' proof and the browser the merge's; carried as a standing tension:
   the look costs the orchestrator's time on every lane that draws.

## Questions for the author left by this close

- B6: open the next plan with the trim pass, or keep raising under the standing leave?

## Decided by the author
