# seq: the principles here

`packages/seq` is a new implementation of Seq\*: SeqViewSpec, the language for sequence views, with its model and its
viewer, in the Mol\* v6 repository. `.seq/` is its instance of the oz principles package, and this README is the
manifest the `principles-router` skill reads. It governs `packages/seq/` and `.seq/` only; the rest of the repository
follows its root `AGENTS.md` and `.v6/`, which win where both speak about the same thing.

**`.seq/` is a development scaffold**, deleted once the package stands on its own. Nothing outside `.seq/` links into
it except the one pointer in `packages/seq/AGENTS.md`, and code never cites it. What must outlive it (the rules the
code follows, the decisions, what the format means) moves into the packages' `_spec/` folders before then.

## Provenance

- The package is `.seq/oz/`: the oz principles (`.principles/` in the oz exploration repository), copied at commit
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

## For the package's home

Fold-back proposals for the oz package, from this instance (2026-11-18); until merged there they stand as the
adaptations above.

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
  Carried, not decided (from an adopting project, 2026-11-18).
- **Standing tension: look before merge, against the cost of a look per lane.** Five visual bugs passed tests and
  review and were seen only in a story after the merge, all in the first stages of a new renderer and viewer; a look in
  every lane's review costs every lane. If later slices find no visual bug before the merge, the look was a first-slice
  cost. Carried, not decided (from an adopting project, 2026-11-18).

## Open items

`packages/seq/TODO.md` once the package exists; until then [.seq/status.md](../status.md).
