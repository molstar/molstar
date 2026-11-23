# seq: the principles here

`packages/seq` is a new implementation of Seq\*: SeqViewSpec, the language for sequence views, with its model and its
viewer, in the Mol\* v6 repository. `.seq/` is its instance of the principles package, and this README is the
manifest the `principles-router` skill reads. It governs `packages/seq/` and `.seq/` only; the rest of the repository
follows its root `AGENTS.md` and `.v6/`, which win where both speak about the same thing.

**`.seq/` is a development scaffold**, deleted once the package stands on its own. Nothing outside `.seq/` links into
it except the one pointer in `packages/seq/AGENTS.md`, and code never cites it. What must outlive it (the rules the
code follows, the decisions, what the format means) moves into the packages' `_spec/` folders before then.

## Provenance

- The package is `.seq/oz/`: the principles (`.principles/` in the upstream exploration repository), copied at commit
  `7d03f60` (2026-11-17). An update from it is a merge from that base. Edited once: its `workspace` and `layout`
  modules removed (not adopted here) and their citations reworded, at the review's recommendation.
- `queries/` and `agent-interface/` here: from the meso\* instance of the same package (`.meso/principles/` in the
  molstar-ludo repository), copied at commit `84c41cd84` (2026-11-17). Their links were repointed at `.seq/oz/`, and
  queries' rule 14 (§bounded-work) became a pointer, the package's representations module now holding it.

## Folder map

| The package says | Here |
|---|---|
| `.principles/` | `.seq/oz/` |
| `.design/`, `.sources/`, `.records/` | `.seq/design/`, `.seq/sources/`, `.seq/records/` |
| `.execution/`, `.tmp/` | `.seq/execution/` (committed and kept, never cited from outside), `.seq/tmp/` (gitignored) |
| `_spec/` | each package's own `_spec/`, as everywhere in Mol\* |
| skills linked from `.agents/skills/` | not linked: the skill files are read directly ([.seq/AGENTS.md](../AGENTS.md)), so removing `.seq/` leaves nothing behind |
| the checks: id, path, merge | `node .seq/oz/_checks/<check>.ts --package .seq/oz --instance .seq` ([.seq/AGENTS.md](../AGENTS.md), Checks) |

## Modules

The router picks from the "Read before" column. A module of this project's own named like a package module extends
it: read the package's first.

### From the package

| Module | Kind | Read before | Covers |
|---|---|---|---|
| [_core](../oz/_core/README.md) | rule set | any non-trivial work (the `principles` skill loads it) | how the principles work and how work is run |
| [architecture](../oz/architecture/README.md) | rule set | packages, ports, wiring, errors, commands, performance | ports, API surfaces (sync, async, streaming), composition root, error codes, commands, performance |
| [representations](../oz/representations/README.md) | rule set | the language, documents, ids, patches, derived forms | declarative state as recursive, streamable trees; committed vs transient state; addresses, hashing; derived forms; encodings |
| [ui](../oz/ui/README.md) | rule set | the viewer's models and views, Storybook | frontend state and logic, design system and Storybook, views as representations |
| [testing](../oz/testing/README.md) | rule set | tests, conformance suites | test tiers (here the repository's lanes), conformance library |

### The project's own

| Module | Kind | Read before | Covers |
|---|---|---|---|
| [queries](queries/README.md) | rule set | the query language, locators, coordinate spaces, mappings | query DSLs, the semantic layer (concepts, spaces, mappings), embedded languages (from meso\*) |
| [agent-interface](agent-interface/README.md) | rule set | exposing the language to agents | skills, the `<domain>_read`/`_write` tool surface taking expressions, outline-first reads, affordances in results (from meso\*) |

## Adaptations

- **workspace**, **layout**: not adopted, and removed from `.seq/oz/`; citations of their ids read "not adopted here".
  The repository's toolchain, tree and package rules (root `AGENTS.md`, `check:workspace`) win.
- **architecture**: its layout (`lib/`, `packages/<product>/*`, `apps/`) maps onto Mol\*'s. The libraries
  `packages/{core,io,model,graphics}` and `packages/ui` play `lib/`, and Seq\* is the product (`packages/seq/*`). The
  root package rules win where they differ. Kept for ports, commands, errors and performance.
- **testing**: its tiers are the repository's three lanes (quick, medium, heavy; root `AGENTS.md`). Allocation and
  heap guards only for the draw and hover hot loops, once they exist. Conformance suites (§conformance-suites) are
  kept.
- **ui**: the runtime and design system are `@molstar/core/reactive` and `@molstar/ui` (no rxjs; React only as an
  optional adapter); §frontend-state is taken as written.
- **representations**: one tree, encoded as JSON (the default) or MessagePack, with an S-expression printer (and
  reader) for inspection, never a file type of its own (the author, 2026-11-17). This is §one-recursive-form's "the spelling is the
  project's choice", made here.
- **_core, `.execution/`**: committed and kept after the close, as in meso\*; §execution-unreferenced stands.
- **_core, followups**: open items go to `packages/seq/TODO.md` (the root `AGENTS.md` sends them to a package's
  TODO); an entry says its context, why it was deferred and what done looks like.
- **Headers**: every source file carries the Mol\* header with the year and `@author` lines only (the author,
  2026-11-17); otherwise the root `AGENTS.md` header rule.
- **Bootstrap patterns**: none taken in yet. The reactive-models pattern is replaced by `@molstar/core/reactive`
  (`ReactiveModel`, `State`, `derive`), which the host already provides.
- **The drift check, edge-case ids**: `.seq/_checks/ec-ids.ts` (the instance's own check, beside the package's three;
  [.seq/AGENTS.md](../AGENTS.md), Checks) resolves every `EC-<n>` cited in `packages/seq` against the catalog,
  [sources/seqstar/edge-cases.md](../sources/seqstar/edge-cases.md) (§mechanical-drift-checks: a registry reference
  that must resolve). Why: the ids are numbers only, so code cites them without a `.seq/` path, and a wrong reading of
  where the catalog lives (the design only refers to it) once recorded an item as missing that was not.

## Folded back from the first slice (the author, 2026-11-18)

Rules for how work runs, learned in the first slice (eleven lanes, the
[record](../records/first-slice/2026-11-18-first-slice.md)) and approved as proposed. Each extends the package rule it
cites, and stands here until the package's home takes it ("For the package's home", below). The first slice's loop
files are frozen with its work folder; these rules govern the next loops, which write them in.

- **The look is part of the lane's review, before the merge** (extends §stories-are-tests, §used-before-done). The
  reviewer, or the orchestrator before the review, opens the lane's stories in light and dark and reports what it saw;
  the gate's screenshots stay as the record, and a loop file puts the look before the merge, not after it. Why: every
  lane passed its tests and an independent review, yet five bugs of what is drawn (values on a third of the width, the
  dark theme ignored, a base canvas not redrawn, nested domains drawn late, an SVG of height 0) were found only by
  opening a story after the merge, each with a one-line cause.
- **A fix's test is shown failing without the fix** (extends the mutation check of §whole-composition-tested). A bug
  fix lands with a test, and its report names the command that shows that test failing on the fix's parent; the
  reviewer reruns it. Why: a fix's first test passed with and without the fix (its settle step waited out the frames
  that hid the bug), and was caught only because the run was repeated on the parent.
- **An observation a fix rests on states its conditions** (extends §used-before-done). It names where it was made (a
  visible or hidden page, the browser, the device pixel ratio) and is reproduced headless on a fake clock before the
  fix starts; a hidden pane is inconclusive for anything that moves. Why: a hidden page runs no animation frames, and
  time went into chasing that artifact before a fake-clock test separated it from the real bug.
- **Defaults are confirmed from tables, the unsure ones first** (extends §gaps-queue). Each package's `_spec` keeps one
  Configuration table (item, default, why, alternatives); the gaps list points at the tables and lists only the
  defaults a lane was unsure of; a gate asks the author for those, not for all. Why: "configurable over fixed" grew
  about sixty documented, tested defaults in one slice, none confirmed, and listed in prose they could not be decided
  in one pass (the tension stays, below).
- **A merge names the review findings refuted or left unfixed, and why** (adopts the proposal of the package's
  standing tension "proportional review"; §merge-one-at-a-time). The loop file says so. Why: the first slice did not
  record them either, so it could not settle that tension; recorded, a later fold-back can.
- **A performance target names its path and its test when it is written, and is measured early** (extends
  §ui-performance, §performance-tested). The stage that builds a path runs a quick representative of its target (one
  scenario, one bound). Why: the targets were first measured four stages after the code they bound, and "a frame
  allocates nothing" was read two ways until the author decided ([design/performance.md](../design/performance.md)).
- **A state derived each frame publishes on change only** (extends §relevant-rerenders): it is set only when its value
  differs (structural equality for lists), and the package's `_spec` says which signals compare how. Why: a list set
  every frame with an equal value made a story re-export every frame, and the churn masked a bug.
- **A story that composes more than one consumer has a headless twin** (extends §stories-are-tests,
  §whole-composition-tested): a test next to the story composes it as the story does (a canvas and an export, two
  views), on a fake clock. Why: one such test found three bugs that the views' own tests could not see
  ([design/testing.md](../design/testing.md)).

## From the author during the second slice (2026-11-18)

- **The fold-back is a step of every close** (extends §reconcile-before-merge, §self-check, §permanent-record). The
  author, 2026-11-18: "we do the measure step and finalization (again we should include and reconcile the foldback.md
  as part of that step; make the fact that foldback is the step part of the principles themselves)". Finalizing a work
  folder writes `foldback.md` in it, always, with:
  - what the work taught, each proposal with its evidence (a commit, a record, a measurement) and where it goes
    (principles, a design, a `_spec`, a check);
  - every rule the earlier fold-backs added, checked against this work's evidence: held, amended or dropped;
  - the self-check on the proposals (the `principles-self-check` skill);
  - a work that taught nothing says so, with its self-check.
  Writing it is required; merging it is not. The author decides, and the decision is a line of its own:
  `Decided by the author, <date>: <merged | merged in part | not merged>: "<their words>"`. Work whose lessons should
  not spread (an exploration, a prototype, a spike) is `not merged`: its fold-back stays in its archive as evidence for
  a later one. Nothing in it is applied outside the folder before that line, and what was decided is applied before
  the record and the closing line of `state.md`. Checked: `_checks/closed-work.ts` (a closed folder has `foldback.md`
  with a decision line). The author, 2026-11-18: "all work should have the foldback file, not every foldback should be
  "merged" .. some exploration work for example should not merge its lessons". Why: the fold-back is how the
  principles learn from the work; as only one category of reconciliation it could be skipped or folded into the
  record, never decided; and lessons from work meant to be thrown away would bend the principles if merged by
  default.

## Folded back from the second slice (the author, 2026-11-18)

Rules learned in the second slice (sixteen lanes; the record in `records/second-slice/`), merged by the author ("Yes,
all of it"). Each extends the package rule it cites, and stands here until the package's home takes it.

- **The self-check runs on every design when it is agreed** (extends §design-first, §self-check). The claim is the
  design's own ("this design covers what a user does with X"); step 2 negates what it rests on; step 4 builds what
  integrates or plans it with a slice, and keeps the rest in the design as a standing tension (§surface-tensions). For
  an action, the negation is usually its inverse or its symmetric case (focus and back, leader and follower, one range
  and several). The author: "should we apply the recurrence to the designs too? the "4 phase cycle" is generally
  quite good at identifying residuals". Why: five of the six things the author's use found missing were the negation
  of a design's own premise, and none was in any design.
- **Interaction is a matrix** (extends §design-first): the UX design keeps a matrix of what a user acts on against the
  verbs, each cell a gesture and command or "no, because" ([design/ux-foundations.md](../design/ux-foundations.md)); a
  slice plan names the cells it fills, its gate checks them. Why: the gaps were empty cells nobody had drawn.
- **Precedents cover interactions, and a novice pass precedes the author's use** (extends §precedents,
  §used-before-done). What a user can do in the precedents (the current Seq\*, IGV, Jalview, the genome browsers) that
  the slice's views cannot is listed, each built, planned or refused; before the author's use gate a light agent plays
  a user in a visible browser with a short task script per view and reports friction. A slice that builds interaction
  has a short use gate after its first drawn stage, not only at its end; the finalize step counts what each gate found,
  by kind (bug, missing design, look). Why: the author's one use gate, at the end, found six missing designs at once and
  cost two extra lanes.
- **What is drawn is checked against the data** (extends §whole-composition-tested): an example's tests check, for each
  ranges layer at a few viewports (the fit view, a zoom across a tile edge), that every item of the data in view has
  its marks and every drawn mark is hit; not only one consumer against another. Why: a transcript's second exon was
  never drawn and merged runs had no hit box, while canvas, SVG and hits agreed with each other.
- **A restriction in a `_spec` states its reason** (extends the repository's "question the specs"): a rule that refuses
  or limits (an issue code, "only", "never") names why, and a review asks for it on the lines a lane adds. Why: "letters
  on their own system only" stood a slice with no reason and blocked a figure the author wanted.
- **Bounds are the best of three, and heavy runs never overlap** (extends §performance-tested): `itWithin`'s heavy
  test runs up to three times, setup and all (`BoundedTries`, `tries: 1` for a test about one run's time); the loop
  never starts a heavy run while another runs. Why: single-sample bounds failed four times under load and passed on
  rerun.
- **Gaps are asked at every gate; a default is confirmed by use** (amends "Defaults are confirmed from tables", above;
  settles the standing tension "configurable over fixed"). Each gate's message to the author ends with the gaps raised
  since the last gate, each a yes or no with the recommendation; finalize carries no unasked gap. A default is
  confirmed when the author uses it at a use gate, or stands as documented and tested; the gaps list carries only the
  defaults a lane was unsure of that change a picture or a gesture. Why: seven gaps went unasked from stage 1 to the
  close, and the author confirmed defaults only by using them ("the gestures are configurable anyway").
- **Sizes are measured at every gate** (amends "A performance target names its path and its test", above): the bundle
  check runs with each gate's checks. Why: the read bundle grew 26 kB at stage 2 and was found at stage 5.
- **A fix of an unreachable state says so** (amends "A fix's test is shown failing without the fix", above): a bound
  on a state no input can reach lands without a failing test, and its report and commit say why. Why: one such bound
  in lane T could have none.

## From the author at the third slice's Gate 0 (2026-11-18)

- **The Storybook's Design pages follow the designs** (extends §design-in-place): a change to a design in `design/`
  that a reader of the package would need (a new design, a rule, a decision, a table) is folded into the matching page
  under the Storybook's `Seq* (next)/Design` in the same change, or a new page is added; the pages never cite `.seq/`.
  The author: "we should be folding the things we are adding to the designs section in the storybook too". Why: the
  scaffold is deleted once the package stands on its own, and the pages are what readers see meanwhile.

## Folded back from the third slice (the author, 2026-11-20)

Rules learned in the third slice (over sixty lanes, figures; its fold-back in `execution/third-slice/foldback.md`),
merged by the author ("and then you can apply the foldback"). Each extends the package rule it cites, and stands here
until the package's home takes it. B10 (module and bundle sizes) went to the package as §module-size.

- **Defaults are designed with a picture** {#defaults-designed} (extends §design-first; amends "Defaults are confirmed
  from tables" and the self-check on every design, above): a design that adds a layer, a layout or a part has a
  Defaults section, what a reader sees with no setting (heights, labels, colors, legends; empty, crowded, at the fit
  view and the closest zoom, in dark), each with a sketch or a reference picture, and its gate reviews it. Why: eleven
  of the author's thirteen missing designs and looks in the third slice were defaults nobody had drawn ("feel empty",
  "really good defaults for everything"), costing an audit and four lanes; the self-check's negation found inverses,
  not defaults.
- **Configuration tables are checked** {#config-tables-checked} (extends "Defaults are confirmed from tables", above,
  and §mechanical-drift-checks): each package's test compares its default configuration's items with its `_spec`'s
  Configuration table, both ways (`seqx-core`'s `src/_test/config-table.test.ts`, `seqx-viewer`'s likewise); a design's
  Configuration row names an item the code has, or says "not built". Why: the close found seventeen rows missing and
  six items the designs named and the code never had.
- **A drift note is a gap** {#drift-notes-are-gaps} (amends "Gaps are asked at every gate", above; extends
  §gaps-queue): each gate's message to the author lists the open items of `drift/` beside the gaps, each a yes or no
  with its recommendation; an item leaves its note when asked or resolved, never parked for the close. Why: two lanes'
  design questions waited from their report to the close unasked.
- **An edge case is tested by name** {#ec-tested-by-name} (extends §whole-composition-tested): a lane's test of a
  catalog item names it (`(EC-43)` in its title or a comment of its file); `.seq/_checks/ec-ids.ts` fails a feature
  marked `done` whose items have neither a test naming them nor an outcome, and `--worklist` lists a built row's items
  no test names. Why: the close reconstructed 56 outcomes by reading tests for their behavior, and two test-only lanes
  then found an edge case not wired.
- **The Features page moves with the merges** {#features-move-with-merges} (amends "The Storybook's Design pages follow
  the designs", above): a lane that builds a feature moves its row in `design/features.md` and its line on the Design
  page "Features" to built, with its story's link, in its own change. Why: at the close the page still listed under
  Planned most of what the slice built.
- **A seeded generator is checked against its published algorithm** {#seeded-generators-checked} (extends
  §conformance-suites): reference outputs and an independent oracle (BigInt for a 64-bit generator) over many seeds,
  every carry reached; a generated dataset's tests state the statistics they rely on. Why: core's `PCG` was biased, and
  its fix changed Mol\*'s distinct colors, a stored snapshot and every seeded fixture at once.
- **A changed test expectation is explained** {#test-expectation-explained} (amends "A fix's test is shown failing
  without the fix", above): a lane that changes an existing test's expected value says why in its report, and the
  reviewer checks it against the design. Why: a test was changed to fit the code (tracks not wrapping on a cyclic
  axis) and only a review restored it.
- **A matrix cell is built when a test names its gesture** {#matrix-cell-by-test} (amends "Interaction is a matrix",
  above): a gate checks a cell by the test that exercises it, not by the design's word. Why: the cursor's rows were
  marked built while the cursor was column-only.
- **Use is continuous** {#use-continuous} (amends "Precedents cover interactions, and a novice pass precedes the
  author's use", above): the author's use as lanes merge is the use gate, its findings counted by kind as they come,
  with a look at what changed before the close; no fixed early gate. The novice pass stays, its task script fitted to
  the slice (a host's integration task where the slice builds for hosts). Why: the early gate was deferred twice and the
  author's use as lanes merged found what it would have; the novice pass found 21 before the author.
- **The heavy lane runs alone, at two workers** {#heavy-lane-alone} (amends "Bounds are the best of three, and heavy
  runs never overlap", above): no other test run beside it, `maxWorkers: 2` in heavy mode. Why: bounds failed in full
  runs under load and passed alone; two workers measured 1/3 to 5/6 passing in the same wall time.
- **A hot loop's quick representative is in every package that composes it** {#hot-loop-representatives} (amends "A
  performance target names its path and its test", above): the viewer's warm emphasis frame has its own case in the
  Storybook's perf tier at the stage that builds it, not only the core's draw loop. Why: stage 7 found the viewer's
  emphasis frame at 800 B while every core case passed.

## Folded back from the fifth slice (the author, 2026-11-23)

Rules learned in the fifth slice (56 lanes, the Genomes milestone; its fold-back in `execution/fifth-slice/foldback.md`),
merged by the author ("Apply all ten"). Each extends the rule it cites, and stands here until the package's home takes
it. The loop's own changes (B3, B4, B6, B8) are in the fifth slice's `loop.md`, "Folded back at the close", for the
method's loop.

- **A lane's model follows its scope** {#model-by-scope} (extends "Lanes at once" in the loop): a lane with a closed
  item list and named stories (a fix lane, a storybook-only lane, a pure split, a fast review) runs on the fast model;
  a lane that builds a design's core (the language, resolve, layouts, runtimes on untrusted input) and its independent
  review keep the strong one. The brief names the model and the findings row records it; the orchestrator's look and
  checks gate both alike. Why: 21 fast lanes this slice, their reviews finding bugs at the strong lanes' rate per item,
  none reverted.
- **A repeated finding is a class** {#findings-are-classes} (extends §whole-composition-tested): a finding that matches
  one already fixed is not logged again per story; its rule goes into the design or `_spec` it belongs to, a sweep test
  covers every story and layer under the rule, and the open instances fold into that one item, briefed to a strong
  lane. Why: the author, "we should not be fixing this one-by-one"; the hover contract (FX9) closed four instances, and
  the card over blocks (FX24) seventeen stories' cards at once.
- **The engine, not the game** {#engine-not-game} (extends §design-first): a design or brief that composes views,
  navigation or optional analyses says in its Intent that it is application-level work on the public surface, names
  the library parts it uses and what the library had to gain for it, and never makes that composition the library's
  only way. Why: the author, "the base library is like a "game engine", not the "game itself"" (the design README).
- **The heavy lane's tight bounds run isolated** {#heavy-bounds-isolated} (amends #heavy-lane-alone): a heavy test
  with a wall-clock bound under 50 ms runs alone in its worker, and its bound is never raised to pass under contention.
  Why: at load 1.2 two such bounds failed in the full heavy lane and passed alone (M5): a quiet machine is not enough
  when the lane's own workers contend.
- **The picture is agreed before it is built** {#mocks-at-design} (extends #defaults-designed; amends
  #use-continuous): a design for a new visible part or gesture comes with an interactive mock and a named precedent.
  - The mock is crude on purpose (fake data, plain SVG or DOM, no performance work), one story per gesture in
    `Seq* (next)/Design/Mocks/`, built by a fast lane beside the draft. It works where the design speaks: what is seen,
    what a press, drag, hover, wheel or key does, what a menu offers.
  - The author plays with it at the design gate; the design is agreed with the mock.
  - The lane's brief links it; the orchestrator's look and the novice pass compare the built view with it, gesture by
    gesture.
  - No package imports it, and the merge that brings the built story deletes it.
  - A use finding that contradicts an agreed mock corrects the design; a second of its kind is a class.

  Why: the slice's largest findings were imaginations that differed, not bugs: the square, letters centred on the
  whole cell, "+n", what a link means. Each design's text was met and only the author's use found them, at the end
  (the author: "how about instead we mocked the interactions too during the design stage?").

## For the package's home

Fold-back proposals for the principles package, from this instance (2026-11-18); until merged there they stand as the
adaptations above.

- §reconcile-before-merge: the fold-back is a step of every close, not a category of reconciliation: a `foldback.md`
  always, the earlier fold-backs' rules checked against the new evidence; the author decides merged, in part or not
  (an exploration's lessons stay in its archive), the decision line checked before the folder closes.
- §stories-are-tests, §used-before-done: the look is part of a lane's review, before the merge, in light and dark.
- §whole-composition-tested (the mutation check): a fix's report names the command that shows its test failing on
  the fix's parent, rerun by the reviewer. Optional at gates: a script that reverts the fix's non-test files on a
  scratch worktree and expects the new test to fail.
- §used-before-done: an observation a fix rests on states its conditions and is reproduced headless on a fake clock;
  a hidden page is inconclusive for anything that moves.
- §gaps-queue: one Configuration table per `_spec`; the gaps list names only the unsure defaults.
- Standing tension "proportional review": its proposal adopted (merges name refuted and unfixed findings); this
  instance's first slice adds no evidence, having not recorded them.
- §ui-performance, §performance-tested: a target names its path and its test, and the stage that builds the path runs
  a quick representative.
- §relevant-rerenders: a state derived each frame is set only when it differs.
- §stories-are-tests, §whole-composition-tested: a story composing more than one consumer has a headless twin on a fake
  clock.
- **Verdict, §orchestrator-owns-lockfile and the hot files of §commit-often** (no change): eleven lanes ran with no
  conflict on the lockfile, the package inventory or the root tsconfig; every hot-file change came as a request in a
  lane's report and was applied by the orchestrator after the merge.
- **Standing tension: configurable over fixed, against a config surface the author can confirm.** Making a choice a
  documented, tested, easy-to-change default instead of a fixed one keeps lanes moving; in one slice it grew about
  sixty defaults, and §gaps-queue cannot be decided at that size. Tables and asking only for the unsure defaults ease
  it but do not settle it: does a default need confirming at all once it is documented, tested and easy to change?
  Carried, not decided (from an adopting project, 2026-11-18). *Settled by the second slice's evidence (2026-11-18)*:
  the author confirmed defaults only by using them, never from a table; a default is confirmed by use at a use gate, or
  stands as documented and tested, and the gaps list carries only unsure defaults that change a picture or a gesture.
  Still open: defaults no story shows (store sizes, bounds) are never confirmed by use.
- **Standing tension: look before merge, against the cost of a look per lane.** Five visual bugs passed tests and
  review and were seen only in a story after the merge, all in the first stages of a new renderer and viewer; a look in
  every lane's review costs every lane. If later slices find no visual bug before the merge, the look was a first-slice
  cost. Carried, not decided (from an adopting project, 2026-11-18). *Released (2026-11-18)*: in the second slice the
  look before the merge found bugs in four of the five lanes that draw, and a renderer bug from the first slice; the
  look pays.
- §design-first, §self-check: the self-check runs on every design when it is agreed (recurrence.md "Checking the loop
  against itself" and the `principles-self-check` skill name the design gate); interaction is a matrix of objects and
  verbs; §precedents covers interactions, and a novice pass precedes the author's use, with a use gate after a slice's
  first drawn stage.
- §whole-composition-tested: what is drawn is checked against the data, not only one consumer against another.
- §performance-tested: a bound is the best of three runs of the whole test; heavy runs never overlap; sizes at every
  gate.
- §gaps-queue: gaps are asked at every gate; a default is confirmed by use.
- Standing tension "proportional review", later evidence (an adopting project, 2026-11-18): three reviews with one
  reviewer each (3, 5 and 8 findings), every finding citing a line and every fixed one confirmed by a test failing on
  the parent; none refuted. One reviewer whose findings each come with a failing test needed no skeptic.

- From the third slice (2026-11-20): §design-first, a Defaults section with a picture in every design that adds
  something drawn, reviewed at its gate; §gaps-queue, open drift notes are gaps asked at the next gate;
  §whole-composition-tested, a test names the catalog item it tests and a check fails a done feature with an item
  neither tested by name nor given an outcome; §whole-composition-tested, a changed test expectation is explained and
  checked against the design; §used-before-done, the author's use is continuous and counted by kind, the novice pass's
  task fitted to the slice.

## Open items

`packages/seq/TODO.md` once the package exists; until then [.seq/status.md](../status.md).
