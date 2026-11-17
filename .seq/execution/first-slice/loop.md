# Loop: the first slice

Governs the implementation of `plan.md` (§loop-file). The orchestrator is the session that runs this loop; lanes are
subagents it starts. Read this file and `state.md` before every step.

## Each step

1. **Read** `state.md`, `gaps.md` and `drift/`. If `state.md` is missing, rebuild it from `plan.md` and the git log.
2. **Pick** the next unchecked items whose dependencies are merged (`plan.md`, Order and parallelism). Start at most two
   heavy lanes (strong model, long work) and two light ones (fast model, reviews, docs) at once (§parallelism).
3. **Brief**: before a lane starts, write its brief in `briefs/<lane>-<name>.md` from the template
   (`.seq/oz/_core/patterns/workflow-files.md`, Lane brief): goal and the designs to read, model, base commit, the files
   it touches and what is out of scope, the hot files it must not stage, its checks, when to stop, what to report.
4. **Start the lane** as a subagent in its own git worktree (isolation: worktree), on a branch `seq/<lane>` cut from
   `packages-seq`, with the brief as its whole instructions and the model the brief names. It commits on its branch only
   at each small green step (§commit-often), never merges, never edits the work folder or a hot file
   (§lanes-report-not-write, §orchestrator-owns-lockfile).
5. **On a lane's report**, never trust it as proof:
   - check out its branch and run its checks yourself (Checks, below); a red check goes back to the lane once with the
     failure; twice red after a fix stops the lane and comes to the author;
   - **review** with a separate strong-model reviewer given the diff and the designs, not the lane's conclusions; the
     findings are fixed by the lane (one round) or by the orchestrator; a finding that is a design question is a drift
     note;
   - **merge** the branch into `packages-seq` (a merge commit naming the lane), then apply the hot-file changes the
     lane asked for (dependencies, inventory, tsconfig references, exports' `_test` nulls, Storybook registration),
     install, and run the checks again on the merged tree;
   - copy the report's drift, gaps and followups into `drift/`, `gaps.md` and `followups.md`.
6. **Drift** (§drift-notes): where code and a design disagree, write `drift/<topic>.md` (what, where, the options, a
   recommendation), continue with what it does not block, and bring it to the author at the next gate.
7. **What is drawn** lands with its story; at the gate the orchestrator opens it in Storybook (the browser pane), in
   light and dark, and takes a screenshot for the report.
8. **Update** `state.md` (done, in progress, next, decisions, blockers, how to verify) and tick `plan.md` after every
   merged item, before starting a new lane, and before stopping; commit the work folder (orchestrator only).

## Gates

At every `GATE` in `plan.md`:

- all of Checks on the merged tree, green by exit code with a non-zero test count (§green-is-exit-code);
- the quick lane's time for each touched folder (under 2 s) and in total (under 30 s, [testing](../../design/testing.md));
  a breach is fixed (waste) or moved with a quick representative (nature) before the gate passes;
- the reviews the gate names, findings resolved;
- `gaps.md` shown, drift notes listed;
- a short report to the author: what was merged, the checks, screenshots where something is drawn, open drift and
  decisions needed. The loop continues to the next stage unless a decision is needed; it waits then.

## Checks

Each its own command, from the repository root:

```sh
pnpm test:unit packages/core/src/repr packages/seq/   # the touched folders' quick lane
npx tsc -b packages/core packages/seq/core packages/seq/viewer   # and the examples package once it exists
pnpm lint
pnpm format:check
node .seq/oz/_checks/principles-ids.ts --package .seq/oz --instance .seq --allow '^(?!\.seq/)'
node .seq/oz/_checks/paths-check.ts --allow '^(?!\.seq/|packages/seq/)'
node .seq/oz/_checks/merge-gate.ts --package .seq/oz --instance .seq
```

At the integration gate (stage 5) also `pnpm test:medium` for the touched packages and the Storybook typecheck; in
stage 7 the heavy lane and the opt-in tiers. `pnpm check:workspace` after `pnpm build`, at the integration gate and at
finalize.

## Stop and ask the author

- a permission refused (report the action and why; never work around it, §refusals-to-author);
- a check red twice after a fix;
- a design question a drift note cannot settle by the designs;
- a public API change beyond what is agreed ([base-library](../../design/base-library.md): each change the author's
  call when its work starts);
- a test or performance budget that cannot be met without changing a design.

## Finalize

When every item is ticked: reconcile (designs edited where the code taught otherwise, drift notes resolved or turned
into followups in `packages/seq/TODO.md`, the `_spec`s with what is built, `status.md`, the features rows, the edge
cases' outcomes for what was built); run the gate once more; write the record in `.seq/records/first-slice/`; end
`state.md` with "Closed <date>: <record file>"; edit nothing in this folder after that.
