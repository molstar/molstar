# Loop: the second slice

Governs the implementation of `plan.md` (§loop-file). The orchestrator is the session that runs this loop; lanes are
subagents it starts. Read this file and `state.md` before every step. The first slice's loop, with what it taught
folded in (the manifest, "Folded back from the first slice").

## Each step

1. **Read** `state.md`, `gaps.md` and `drift/`. If `state.md` is missing, rebuild it from `plan.md` and the git log.
2. **Pick** the next unchecked items whose dependencies are merged and whose Gate 0 items are decided (Gate 0,
   below). Start at most two heavy lanes (strong model, long work) and two light ones (fast model, reviews, docs) at
   once (§parallelism).
3. **Brief**: before a lane starts, write its brief in `briefs/<lane>-<name>.md` from the template
   (`.seq/oz/_core/patterns/workflow-files.md`, Lane brief), and the first slice's briefs as examples. Every brief
   has:
   - the goal and the designs to read;
   - the model, the base commit;
   - the files it touches and what is out of scope;
   - the hot files it must not stage;
   - its checks, when to stop, what to report.

   Each also names:
   - **the agreed surface**, for a lane in the base library (`@molstar/core`, `@molstar/io`): the exports and
     signatures Gate 0 agreed. Anything past them stops the lane (Stop and ask);
   - **the two consumers** its shared pieces are tested with, live in one test (the plan's list: the viewer and an
     export; two views on one store; an agent and the user on one history; a Seq\* spec and a stand-in structure on one
     host stack);
   - **its performance representative**: the target the path it builds carries (performance.md), as one quick case,
     and the heavy bound if the target is a time;
   - **configurable over fixed**: choices that could go more than one way become config items with the recommendation
     as default, in the package `_spec`'s Configuration table.
4. **Start the lane** as a subagent in its own git worktree (isolation: worktree), on a branch `seq/<lane>` cut from
   `packages-seq`, with the brief as its whole instructions and the model the brief names. It commits on its branch only
   at each small green step (§commit-often). It never merges, and never edits the work folder or a hot file
   (§lanes-report-not-write, §orchestrator-owns-lockfile).
5. **On a lane's report**, never trust it as proof:
   - **check**: check out its branch and run its checks yourself (Checks, below). A red check goes back to the lane once
     with the failure; twice red after a fix stops the lane and comes to the author.
   - **look**, when something is drawn: open the lane's stories in a **visible** pane, in light and dark, and take
     screenshots. Before trusting a frame, confirm with `tabs_context` that the pane is shown. A hidden pane runs no
     animation frames and is never evidence. While it is hidden, ask the author once to show it, and carry on with
     what does not need the look.
   - **review**: a separate strong-model reviewer gets the diff, the designs and the screenshots, never the lane's
     conclusions. Beyond correctness, it checks:
     - **host neutrality**: no Seq\* package imports a host's model (`@molstar/ui/history`, a plugin's state) outside
       a composition root's adapter;
     - the two-consumer tests are there;
     - each fix's test fails on the fix's parent: the report names the command, and the reviewer reruns it.

     The lane fixes the findings in one round, continued with SendMessage (§reuse-agents), or the orchestrator does. A
     finding that is a design question becomes a drift note.
   - **merge** the branch into `packages-seq` with `--no-ff`. The merge commit names the lane and the review's refuted
     or unfixed findings, with why. Then:
     - apply the hot-file changes the lane asked for (dependencies, inventory, tsconfig references, exports' `_test`
       nulls, Storybook registration, a `.v6/changelogs/` entry for a base-library change);
     - install, and run the checks again on the merged tree;
     - remove the lane's worktree.
   - **record**: copy the report's drift, gaps and followups into `drift/`, `gaps.md` and `followups.md`.
6. **Drift** (§drift-notes): where code and a design disagree, write `drift/<topic>.md` (what, where, the options, a
   recommendation). Continue with what it does not block, and bring it to the author at the next gate. A decision the
   author already made (`state.md`, Decided) is applied to the design when its lane starts: changes.md's history line
   with lane N.
7. **What is drawn** lands with its story and the story's headless twin, composed as the story composes it, on a fake
   clock. It is looked at before the review (step 5); the gate's screenshots are the record.
8. **Update** `state.md` (done, in progress, next, decided, blockers, how to verify) and tick `plan.md` after every
   merged item, before starting a new lane, and before stopping. Commit the work folder (orchestrator only).

## Gate 0

The plan's Gate 0 items go to the author in one message, each with its recommendation and the evidence, as soon as the
first slice is closed. Stage 0 runs beside it. A lane waits only for the items it builds on:

| Lane | Waits for |
|---|---|
| K (patches) | `repr`'s `patch` and `diff` |
| L (MessagePack) | the codec in io |
| M (asset manager) | the asset manager change; without it, M is dropped and security.md keeps "what waits" |
| N (editing) | the history port: decided, 2026-11-18 |
| O, P, Q | the first slice's open decisions they touch (the selection tint for P) |

## Running while the author is away

When the author says to continue without them, the loop runs through the gates without waiting:
- a decision a gate would ask for is made as recommended and listed in gaps.md for the author's review;
- the stop conditions (below) still stop it;
- Gate 0's base-library items are never decided this way: they are public API changes;
- §used-before-done is never passed this way: stage 4's gate waits for the author's use.

## Gates

At every `GATE` in `plan.md`:

- all of Checks on the merged tree, green by exit code with a non-zero test count (§green-is-exit-code);
- the quick lane's time for each touched folder (under 2 s) and in total (under 30 s, [testing](../../design/testing.md)).
  A breach is fixed (waste) or moved with a quick representative (nature) before the gate passes;
- the reviews the gate names, findings resolved or named in the merge;
- `gaps.md` shown (the unsure defaults only), drift notes listed;
- a short report to the author: what was merged, the checks, screenshots from a visible pane where something is drawn,
  the performance representatives' numbers, open drift and decisions needed. The loop continues to the next stage
  unless a decision is needed; it waits then.

## Checks

Each its own command, from the repository root:

```sh
pnpm test:unit packages/core/src/repr packages/io/src packages/seq/   # the touched folders' quick lane
npx tsc -b packages/core packages/io packages/seq/core packages/seq/viewer packages/seq/examples
pnpm -C packages/seq/core typecheck:test     # and the same for viewer and examples
pnpm lint
pnpm format:check
node .seq/oz/_checks/principles-ids.ts --package .seq/oz --instance .seq --allow '^(?!\.seq/)'
node .seq/oz/_checks/paths-check.ts --allow '^(?!\.seq/|packages/seq/)'
node .seq/oz/_checks/merge-gate.ts --package .seq/oz --instance .seq
node .seq/_checks/ec-ids.ts
```

When each also runs:
- **at each gate:** the Storybook typecheck (`pnpm -C packages/seq/storybook typecheck`) and `pnpm test:medium` for the
  touched packages;
- **at stage 1's gate (core and io changed) and at finalize:** `pnpm build`, then `pnpm check:workspace`. Its known
  failure (`@molstar/ui`'s `glass-wordmark` export) is shown on master each time, not counted;
- **in stage 5:** the heavy lane (`pnpm test:heavy packages/seq/`), the opt-in tiers (`pnpm -C packages/seq/core
  test:perf`, in the background), the bundle sizes (`pnpm -C packages/seq/core size`) and the browser measurements in
  a visible pane.

## Stop and ask the author

- a permission refused (report the action and why; never work around it, §refusals-to-author);
- a check red twice after a fix;
- a design question a drift note cannot settle by the designs;
- a public API change beyond what Gate 0 agreed ([base-library](../../design/base-library.md));
- a test or performance budget that cannot be met without changing a design;
- a piece that can only work if Seq\* owns what a host owns (the undo stack, the keymap, the chrome).

## Finalize

When every item is ticked:

1. **Reconcile**:
   - designs edited where the code taught otherwise;
   - each drift note resolved, or turned into a followup in `packages/seq/TODO.md`, and `drift/` emptied
     (§mechanical-drift-checks);
   - the `_spec`s with what is built, their Configuration tables complete;
   - `status.md`, the features rows, the edge cases' outcomes for what was built.
2. **Fold back**: write `foldback.md` (what the slice taught, each proposal with its evidence and where it goes, and
   the self-check). It is applied and committed only on the author's approval.
3. Run the gate once more.
4. Write the record in `.seq/records/second-slice/`, end `state.md` with "Closed <date>: <record file>", and edit
   nothing in this folder after that.
5. Rewrite `third-slice/plan.md` from the draft with what this slice taught.
