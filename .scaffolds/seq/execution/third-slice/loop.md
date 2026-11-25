# Loop: the third slice

Governs the implementation of `plan.md` (§loop-file). The orchestrator is the session that runs this loop; lanes are
subagents it starts. Read this file and `state.md` before every step. The second slice's loop, with what it taught
folded in (the manifest, "Folded back from the second slice" and "From the author at the third slice's Gate 0").

## Each step

1. **Read** `state.md`, `gaps.md` and `drift/`. If `state.md` is missing, rebuild it from `plan.md` and the git log.
2. **Pick** the next unchecked items whose dependencies are merged (plan, Order and parallelism). Start at most two
   strong lanes (long work) and two fast ones (fast model, reviews, docs) at once (§parallelism).
3. **Brief**: before a lane starts, write its brief in `briefs/<lane>-<name>.md` from the template
   (`.scaffolds/principles/_core/patterns/workflow-files.md`, Lane brief), with the second slice's briefs as examples. Every brief
   has:
   - the goal and the designs to read: the agreed design the lane builds (plan, its item) and the sections of
     `ux-foundations.md` whose matrix cells it fills;
   - the model, the base commit;
   - the files it touches and what is out of scope;
   - the hot files it must not stage;
   - its checks, when to stop, what to report.

   Each also names:
   - **the agreed surface**, for a change in the base library (`@molstar/core`, `@molstar/io`): the additions its
     design lists (themes' color additions, the hydrophobicity table's move, the shared 2D reduction if it lands in
     core). Anything past them stops the lane (Stop and ask);
   - **the two consumers** its shared pieces are tested with, live in one test (canvas and SVG; the viewer and the
     headless export; a layout and hit-testing; a navigator and the view it moves);
   - **its performance representative**: the target the path it builds carries (performance.md and the design), as one
     quick case, and the heavy bound under `itWithin` if the target is a time;
   - **whose each choice is**: the design's Configuration table (spec param, `RenderConfig`, `SpecBounds`,
     `ViewerConfig`, port, call option, or fixed with its reason) copied into the package `_spec`'s table; nothing
     configurable without a second use; plain defaults, never inferred from what a section draws;
   - **minimal code**: reuse `@molstar/core` and the package's own helpers; no private copies; nothing without a use
     (the author: "try to minimize the amount of code"); the sizes measured in its report;
   - **drawing against the data**: the example views it touches checked (every item in view has its marks, every
     drawn mark is hit) at the fit view and across a tile edge.
4. **Start the lane** as a subagent in its own git worktree (isolation: worktree), on a branch `seq/<lane>` cut from
   `packages-seq`, with the brief as its whole instructions and the model the brief names. It commits on its branch only
   at each small green step (§commit-often). It never merges, and never edits the work folder or a hot file
   (§lanes-report-not-write, §orchestrator-owns-lockfile).
5. **On a lane's report**, never trust it as proof:
   - **check**: check out its branch and run its checks yourself (Checks, below). A red check goes back to the lane once
     with the failure; twice red after a fix stops the lane and comes to the author.
   - **look**, when something is drawn: point the `storybook-lane` launch configuration (port 6008, added to
     `.claude/launch.json` for the look and reverted after) at the lane's worktree, open its stories in a **visible**
     pane, in light and dark (`&globals=theme:dark`), and take screenshots. Before trusting a frame, confirm with
     `tabs_context` that the pane is shown. A hidden pane runs no animation frames and is never evidence. While it is
     hidden, ask the author once to show it, and carry on with what does not need the look.
   - **review**: a separate strong-model reviewer gets the diff, the designs and the screenshots, never the lane's
     conclusions. Beyond correctness, it checks:
     - the design followed: its rules, its Configuration table's "whose", its matrix cells; a departure is a finding;
     - **host neutrality**: no Seq\* package imports a host's model outside a composition root's adapter; legends and
       navigators are parts a host places, the animation loop a port;
     - minimal code: a private copy of something in `@molstar/core` or the package is a finding;
     - the two-consumer tests and drawing against the data are there;
     - each fix's test fails on the fix's parent: the report names the command, and the reviewer reruns it.

     The lane fixes the findings in one round, continued with SendMessage (§reuse-agents), or the orchestrator does. A
     finding that is a design question becomes a drift note.
   - **merge** the branch into `packages-seq` with `--no-ff`. The merge commit names the lane and the review's refuted
     or unfixed findings, with why. Then:
     - apply the hot-file changes the lane asked for (dependencies, inventory, tsconfig references, exports' `_test`
       nulls, Storybook registration, a `.v6/changelogs/` entry from its report, the `cli/` package for OP);
     - install, and run the checks again on the merged tree;
     - remove the lane's worktree.
   - **record**: copy the report's drift, gaps and followups into `drift/`, `gaps.md` and `followups.md`.
6. **Drift** (§drift-notes): where code and a design disagree, write `drift/<topic>.md` (what, where, the options, a
   recommendation). Continue with what it does not block, and bring it to the author at the next gate. A decision the
   author made is in the designs, and in its Storybook Design page, before the lane that builds it starts.
7. **What is drawn** lands with its story and the story's headless twin, composed as the story composes it, on a fake
   clock. It is looked at before the review (step 5); the gate's screenshots are the record.
8. **Designs and their pages move together**: a change to an agreed design (a drift resolved, a gap decided) updates
   its `.scaffolds/seq/design/` file and its Storybook Design page (`packages/seq/storybook/src/_stories/design/`) in one commit
   (principles, "From the author at the third slice's Gate 0").
9. **Update** `state.md` (done, in progress, next, decided, blockers, how to verify) and tick `plan.md` after every
   merged item, before starting a new lane, and before stopping. Commit the work folder (orchestrator only).

## Gate 0

Decided (plan, Gate 0), and the slice's designs agreed: the loop starts at stage 0. A lane waits only for its plan
dependencies.

## Running while the author is away

When the author says to continue without them, the loop runs through the gates without waiting:
- a decision a gate would ask for is made as recommended and listed in `gaps.md` for the author's review;
- the stop conditions (below) still stop it;
- a public API change beyond what a design agreed is never decided this way;
- the use gates are deferred, not passed (plan, "Use gates"; the author: "yes, update the plan to defer the use
  gates"): the loop runs on past them and leaves them unticked; the author's use comes when they are back, on the views
  as they then stand. Finalize waits for both (§used-before-done): nothing is called done before the author has used
  it. The novice pass runs in its place in stage 6 (a visible pane).

## Gates

At every `GATE` in `plan.md`:

- all of Checks on the merged tree, green by exit code with a non-zero test count (§green-is-exit-code);
- the quick lane's time for each touched folder (under 2 s) and in total (under 30 s, [testing](../../design/testing.md)).
  A breach is fixed (waste) or moved with a quick representative (nature) before the gate passes;
- the sizes (`pnpm -C packages/seq/core size`), against the targets as guides, not bounds (plan, Gate 0); a growth
  with no feature behind it is a finding;
- drawing checked against the data in every example view the stage touched;
- the interaction matrix's cells the stage fills, checked against `ux-foundations.md` (a cell built, or "no, because"
  stated);
- the reviews the gate names, findings resolved or named in the merge;
- the gaps raised since the last gate, each a yes or no question with a recommendation (only the unsure defaults);
  drift notes listed;
- a short report to the author: what was merged, the checks, screenshots from a visible pane where something is drawn,
  the performance representatives' numbers, the sizes, open drift and decisions needed. The loop continues to the next
  stage unless a decision is needed; it waits then.

At a **use gate** (deferred: plan) and the **novice pass**: each finding is fixed or placed (a lane, `TODO.md`, a drift note) and
counted by kind (bug, missing design, look) in `state.md`; the fold-back reads the count (the carried question: does an
early use gate pay for figures?). The novice pass is a fast agent playing a user in a visible browser with a task
script per view, before the author's use gate.

## Checks

Each its own command, from the repository root:

```sh
pnpm test:unit packages/seq/                 # and packages/core/src/<touched> when core changed
npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples   # and packages/core when it changed
pnpm -C packages/seq/core typecheck:test     # and the same for viewer and examples
pnpm lint
pnpm format:check
node .scaffolds/principles/_checks/principles-ids.ts --package .scaffolds/principles --instance .scaffolds/seq --allow '^(?!\.scaffolds/)'
node .scaffolds/principles/_checks/paths-check.ts --allow '^(?!\.scaffolds/|packages/seq/)'
node .scaffolds/principles/_checks/merge-gate.ts --package .scaffolds/principles --instance .scaffolds/seq
node .scaffolds/seq/_checks/ec-ids.ts
node .scaffolds/seq/_checks/closed-work.ts
node .scaffolds/seq/_checks/file-size.ts
```

The `.seq` checks read the git index: stage first. When each also runs:
- **at each gate:** the Storybook typecheck (`pnpm -C packages/seq/storybook typecheck`), `pnpm test:medium` for the
  touched packages, `pnpm test:heavy packages/seq/` (never two heavy runs at once), the size check;
- **at a gate after a base-library change (TH; RW if the reduction lands in core) and at finalize:** `pnpm build`,
  then `pnpm check:workspace`. Its known failure (`@molstar/ui`'s `glass-wordmark` export) is shown on master each
  time, not counted;
- **in stage 7:** the opt-in tiers (`pnpm -C packages/seq/core test:perf`, in the background) and the browser
  measurements in a visible pane at dpr 2.

## Stop and ask the author

- a permission refused (report the action and why; never work around it, §refusals-to-author);
- a check red twice after a fix;
- a design question a drift note cannot settle by the designs;
- a public API change beyond what a design agreed ([base-library](../../design/base-library.md));
- a test or performance budget that cannot be met without changing a design;
- a piece that can only work if Seq\* owns what a host owns (the undo stack, the keymap, the chrome, the animation
  loop).

## Finalize

When every item is ticked:

1. **Reconcile**:
   - designs and their Storybook Design pages edited where the code taught otherwise;
   - each drift note resolved, or turned into a followup in `packages/seq/TODO.md`, and `drift/` emptied
     (§mechanical-drift-checks);
   - the `_spec`s with what is built, their Configuration tables complete;
   - `status.md`, the features rows, the edge cases' outcomes for what was built;
   - the size targets reset to the measured sizes with room (plan, Gate 0).
2. **Fold back** (principles, "The fold-back is a step of every close"): write `foldback.md`:
   - what the slice taught, each proposal with its evidence and where it goes;
   - the earlier fold-backs' rules (principles, "Folded back from the second slice") each checked against this
     slice's evidence: held, amended or dropped;
   - the count of what each gate found, by kind;
   - the self-check.
   The author decides (merged, merged in part, or not merged); their words go on its `Decided by the author` line;
   then what was decided is applied (principles, designs, `_spec`s, checks) and committed, before the record.
3. Run the gate once more.
4. Write the record in `.scaffolds/seq/records/third-slice/`, end `state.md` with "Closed <date>: <record file>", and edit
   nothing in this folder after that.
5. Write `fourth-slice/plan.md` with what this slice taught.
6. **Clean up the worktrees** (the author, 2026-11-20: "part of the slice finalize should be cleanup of the
   worktrees"): every lane's worktree and branch whose work is merged is removed (`git worktree remove`, then the
   branch deleted when `git branch --merged` lists it); a parked lane's (CE) is kept and named in the record; nothing
   unmerged is deleted without the author's word.
