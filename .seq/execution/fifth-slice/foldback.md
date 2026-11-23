# Fold-back: what the fifth slice taught (draft, for the author)

The fifth slice built version 1's genome milestone in 56 merged lanes. Their list and each lane's findings are in
`state.md`, "Findings by lane". The milestone:
- coordinates over assemblies, and sections over a complex;
- BGZF and tabix, and the genome formats;
- region reads and their pairs;
- bridges between assemblies (synteny, the dot plot, liftover);
- the Genomes dataset;
- the genome ring as an extension, with its square and kinds of links.

It was the first slice to run the fourth slice's rules (B1–B11): a look is a gesture, the orchestrator looks itself, a
quiet machine for timing, Defaults checked, `features.md` checked, edge cases named, and the novice pass after the last
drawing lane. It also ran under the author's mid-slice decisions:
- Sonnet for scoped lanes;
- classes fixed by a rule and a sweep;
- the engine, not the game.

This file follows the principles' rule "The fold-back is a step of every close". **Nothing here is applied before the
author's line at the end.** Each proposed change is one line the author decides: "Proposed: … Why: … Where: …".

Where to:

| Code | Where |
|---|---|
| **P** | principles (the instance's adaptations) |
| **D** | a design |
| **S** | a package `_spec` |
| **C** | a check |
| **L** | the next loop file |

## A. What each gate, lane and use found, by kind

From `state.md` (the findings table), `novice-findings.md`, `measurements.md`, `fixes.md` and this close's reconcile.

| Gate, lane or use | Bug | Missing design | Look | Other |
|---|---|---|---|---|
| Lane reviews and the orchestrator's looks (56 lanes: 35 strong-model, 21 Sonnet fix or showcase lanes) | 427 | 69 | 42 | — |
| The author's use during the slice (fixes.md: the square, arcs and sub-bands, "+n", the ring's top labels, the select box on chr1, letters off-centre, …) | about 20 | 6 (the square as a real window, scattered sets, kinds of links, the engine and the game, centring letters on the whole cell, the Prism swap's data port) | 4 | — |
| Stage gates (1, 2–3, 4) | 0 | 7 decided (stage 1's five; stage 2–3's ClinVar and genotype tokens) | 0 | timing misses at load 4 (stage 4: 2, rerun quiet: pass); one quick-lane flake (fixed by FX24) |
| Novice passes (two, 36 findings) | 5 (2 fixed by FX24; 3 not bugs: a tool's click, a wheel by design, the square following a pair) | 14 | 15 | 2 unverified |
| M5 | 0 | 0 | 0 | 2 over target (the split patch 16.2 ms median; one Signal zoom frame of 25.7 ms); 2 heavy bounds that fail under the lane's contention; the minimal viewer +150.7 kB |
| Finalize's reconcile (this close) | see the reconcile commit | | | |

Read:
- **Reviews again found the bugs**: in every lane but three (NS, VP, FX19).
- **The author's use found what no review, test or look did**, and almost all of it was interaction on the ring and
  its square: an undraggable box, a two-click switch, arcs that couldn't be hovered, sub-bands lighting everything,
  "+n" not opening, labels overlapping.
- **The novice passes found mostly missing design** (14 of 36): the page gave no way or no hint.

## B. What this slice taught

1. **A lane's model follows its scope.** (L; P)
   - Proposed: a lane with a closed item list and named stories runs on a Sonnet subagent. That covers a fix lane over
     `fixes.md`, a storybook-only lane, a pure split and a fast review. A lane that builds a design's core (the
     language, resolve, layouts, runtimes on untrusted input) and its independent strong review keep the strong model.
     The brief names the model, and the findings row records it. The orchestrator's look and checks gate both alike.
   - Why: the author, 2026-11-21: "maybe use this more in general for well-scoped subtasks?", after NS and FX2. This
     slice ran 21 Sonnet lanes (FX2–FX24, NS, RL, both novice passes). Their reviews found 2 to 8 bugs each, as the
     strong lanes' did per item, and no Sonnet lane's merge was reverted.
   - Where: L (the loop's Lanes section), P (the instance's adaptations: who runs a lane).

2. **A finding that repeats a fixed one is a class: fix it by a rule and a sweep, not by the instance.** (L; C; P)
   - Proposed: when a finding matches one already fixed (the same miss on another layer, layout or story), the
     orchestrator writes the rule the instances break into its design or `_spec`. It adds a sweep test over every story
     and layer the rule covers, folds the open instances into that one item, and briefs it to a strong lane.
   - Why: the author, 2026-11-21: "we should not be fixing this one-by-one btw. since we already fixed a similar issue,
     we need a general rule for this and global fix". The hover contract and its sweep (FX9) closed four open
     instances at once. FX24's "a card over blocks lists the block under the pointer's items" is the same shape: one
     rule, and 17 stories' cards corrected together.
   - Where: L (the fix-point step: classify before logging), C (a sweep per contract), P.

3. **A brief's items are checked as landed, not as reported.** (L; C)
   - Proposed: at a merge, the orchestrator ticks each numbered item of the brief against a test or a story, by name,
     in the findings row. An item without one is "not landed" and goes back to the lane or into `fixes.md`. The
     design's "Built" note lists only ticked items.
   - Why: GC2 was merged with `section.center` in its brief, and the RR design's `PairsReader` was marked as planned.
     Neither existed: RX found both missing on 2026-11-23 and had to take on the language bump. The reports said done,
     and nothing compared them with the brief.
   - Where: L (the merge step), C (features-planned already checks the plan's names; extend it to the brief's items).

4. **What the orchestrator owns stays where lanes can't reach it.** (L)
   - Proposed: the orchestrator's scripts (the checks script, the gate script) are kept in its scratchpad, and every
     brief says "write only in your worktree and your own temporary space". The orchestrator stages by path and never
     commits with `-a` or `.`.
   - Why: RX3 overwrote the orchestrator's `checks.sh`, and FX24 wrote scratch files into the orchestrator's scratchpad.
     The orchestrator's own `commit -a` swept a temporary `.claude/launch.json` entry into `3e5345f46`, reverted in
     `08e4f4e27`.
   - Where: L (Lanes, the brief's "Where"; the merge step's commit).

5. **The engine, not the game: an application's work is labelled as one.** (P; D)
   - Proposed: a design or brief that composes views, navigation or optional analyses (the ring and its square, a
     showcase, a host's menus) says in its Intent that it is application-level work on the public surface. It names the
     library's parts it uses, and anything the library had to gain for it (RX: `section.center`, `PairsReader`,
     `settingOf`). The library never ships that composition as its only way.
   - Why: the author, 2026-11-23: "the base library is like a "game engine", not the "game itself". We should be
     explicit about this". Already said in the design README and ring-extension's Intent, and leading the storybook's
     Docs (DX).
   - Where: P (the instance's adaptations), D (the design template's Intent).

6. **A novice finding is confirmed before it is placed.** (L)
   - Proposed: the orchestrator reproduces every bug a novice pass reports in a visible pane before it goes to
     `fixes.md`, and splits the pass by story groups from the start, so one budget never covers everything.
   - Why: of the two passes' 5 bugs, 3 were not bugs: a browser tool's click, a wheel that pans by design, and a
     square checked on the SVG copy, which has none. The first pass ran out of budget before a third of its script, and
     a second pass was needed.
   - Where: L (Gates, the novice pass).

7. **The heavy lane's timing bounds hold only in isolation.** (P: amends B5; C)
   - Proposed: a heavy test with a wall-clock bound under 50 ms runs isolated, alone in its worker (`vitest.config.ts`'s
     `Isolated`, or `itWithin`'s own quiet guard). Its bound is never raised to pass under contention.
   - Why: at load 1.2, the lookup over 1,000,000 ranges (0.71 against 0.62 ms) and the pairs pyramid's slices (22–51
     against 20 ms) failed in the full heavy lane, and passed alone (M5). The fourth slice's B5 made the machine quiet,
     but the lane's own workers are still contention.
   - Where: P (B5's rule), C (the heavy lane's config).

8. **A plan's size line is a budget, not "not grown".** (L)
   - Proposed: a plan states each bundle's expected growth per stage, from its designs, as a number. "Not grown" is
     written only where nothing is added.
   - Why: the plan said "the minimal viewer not grown", while every stage added to it (+150.7 kB over the slice, each
     raise with its reason, per the soft-size rule). The line could not be met, and it went unchecked until M5.
   - Where: L (the plan template's measurement line).

9. **A test that awaits a lazy part says so.** (S: the storybook tests' rule)
   - Proposed: a test over a view with lazy chunks awaits them (`headers.loaded()`, the card's chunk) before the gesture
     it measures. A test under fake timers that still flakes has a late chunk or a real clock as its first suspect.
   - Why: the hover card story's flake (stage 4's gate) was a lazy header chunk landing after the card and relaying the
     headers (FX24).
   - Where: S (the storybook package's test notes, or the viewer `_spec`'s lazy-chunk section).

10. **The picture is agreed before it is built: what the author imagines and what the orchestrator imagines are made
    to meet early.** (P; D; L)
    - Proposed, for a new visible part or a new gesture (not a fix), the author's refinement (2026-11-23: "how about
      instead we mocked the interactions too during the design stage?"):
      1. **An interactive mock at the design stage.** Beside the design draft, a scoped lane (Sonnet) builds a mock the
         author can use, one story per gesture, in `Seq* (next)/Design/Mocks/`:
         - it is crude on purpose: fake data, plain SVG or DOM, no performance work, no tests beyond a smoke render;
         - it works where the design speaks: what is seen, what a press, drag, hover, wheel or key does, what a menu
           offers;
         - where the library's public surface already does the job (a host composition, an extension's layer), the mock
           uses it.

         The author plays with it at the design gate, and the design is agreed with the mock, not only the text. It is
         cheap at that point because nothing is built on it yet.
      2. **Name a precedent.** The design names the tool or view whose behaviour it borrows ("a window like IGV's
         locus box") and where it departs from it, and the mock follows that precedent.
      3. **The mock is the look's reference.** The lane's brief links the mock. The orchestrator's look and the novice
         pass compare the built view against it, gesture by gesture, not against their own picture.
      4. **The mock is thrown away, never grown into the code.** It stays as a picture of the agreement until the built
         story replaces it in the same merge, then it is deleted. Its code is never imported by a package; it is
         the storybook's alone, so a crude mock never becomes the implementation.
      5. **A miss corrects the mock.** When the author's use still contradicts something the mock agreed, the
         orchestrator corrects the design and notes it on the mock's story (as long as the mock exists), not only in
         `fixes.md`. A second miss of the same kind on a later part is a class (B2).
    - Why: the slice's largest findings were differences of imagination, not bugs. Each design text was met; the author's
      picture was another:
      - the square: cells clicked twice in the design; a draggable window in the author's;
      - letters centred on the visible cell rather than the whole one;
      - "+n" as a count rather than a button;
      - the ring's top labels;
      - what a link between chromosomes means.

      Reviews, tests and the orchestrator's looks all checked against the orchestrator's picture, so none could find
      them; only the author's use did, at the end.
    - Where: P (the instance's adaptations: a design for a visible part comes with a mock and a precedent), D (the
      design template: a "Mock" line linking its stories), L (the design stage's mock lane; the look and the novice pass
      against the mock; the mock deleted at the merge).

## C. The earlier fold-backs' rules, checked against this slice

From the fourth slice (B1–B11, "All 11"):

1. **A look is a gesture** — *held, and still not enough.* Every look listed its gestures. The author's ring findings
   (an undraggable box, hover on thin arcs) came after looks that had dragged and hovered. A lane's look gestures in the
   way its own story expects; the author used the views as an unfamiliar user would. No change; B6 and the use gate
   carry it.
2. **The orchestrator looks itself, every merge** — *held* (every merge that draws, in a visible pane).
3. **The lane limit counts every agent** — *held.*
4. **A pane's screenshot can be stale** — *held* (re-shot before every claim).
5. **Timed checks need a quiet machine** — *held for the machine, amended* (B7: the lane's own contention).
6. **Defaults checked** (`defaults-designed`) — *held*: green at every gate.
7. **`features.md` checked** (`features-planned`) — *held for names; amended* (B3: a brief's items, not only the
   plan's).
8. **Edge cases named** — *held* (`ec-ids` green; the close's worklist in the reconcile).
9. **Reviews in every lane** — *held*: 427 bugs in 53 of 56 lanes.
10. **The novice pass after the last drawing lane** — *held; amended* (B6: confirm and split).
11. **No unticked items at finalize** — *held*: the use gate carried by the author's word, M5 and the novice pass run
    before the close.

From the earlier slices:
- Sizes are bounds: *held, amended by B8*.
- §module-size: *held*. Four files are over 800 lines (record.ts, views.ts, navigator.ts at 1,290, viewer.ts at 860),
  all in the TODO, none over 1,500.
- A changed test expectation is explained: *held after review*. FX24's 17 card digests were decoded old against new
  in its fix round.
- The fold-back is a step of every close: *held*, this file.

## D. Not proposed

- The 14 design questions from the novice passes: they are the next plan's (the use gate carried with them), not rules.
- The Prism swap and the old Seq\*'s removal: deferred by the author, separate work.
- The sixth slice (scaffolds): drafted, starts after this close.

## Self-check (§self-check), on the claim "these proposals make the loop truer to what was built"

1. **Assume it holds.**
   - B3 closes the one gap that cost a lane's scope (the reported but unlanded items).
   - B4 and B7 stop the two kinds of interference this slice met.
   - B1, B2, B5 and B6 write down what the author already decided or the slice already did.
   - B8 and B9 make two lines true.
2. **What it rests on, negated.** It rests on the loop's misses being procedural: a missing check, a missing rule.
   Negated: the slice's most important findings (the ring's interaction) were not procedural. No rule here would have
   found them, because they came from the author using a new kind of view with intent.
3. **The negated world.** It explains why 427 review bugs and two novice passes still left the author's ring findings.
   It predicts that the next new view (the next plan's) will again be found wanting in use, whatever the loop's rules.
4. **Both held.**
   - *Integrates*: the procedural rules stand (they caught or would have caught their cases).
   - *Carried*: the use gate is carried to the next plan, as the author decided, with the 14 design questions. The
     next plan schedules the author's use of each new view as a stage, not as a gate at the end.

## Questions for the author left by this close

1. The 14 design questions from the novice passes ([novice-findings.md](novice-findings.md)): carry them all to the next
   plan with the use gate, or decide some now?
2. The minimal viewer's +150.7 kB (`measurements.md`): accept it as the slice's budget, or plan a trim in the next one?

## Decided by the author

Not yet.
