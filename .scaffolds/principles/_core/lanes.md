# Lanes {#lanes}

Running work in parallel: lanes in worktrees, an orchestrator that merges, briefs that carry everything a lane needs.
Read this only when fanning out; what every agent needs (small green commits, green as an exit code, reusing agents,
refusals) is in [workflow.md](workflow.md#agents-at-work). How each step runs is in the workflow-files pattern: the
[brief](patterns/workflow-files.md#lane-brief), the [merge step and findings row](patterns/workflow-files.md#loopmd).

1. **Parallelism**: {#parallelism} decide what runs at once, sized to the machine: by default 1-2 *heavy* lanes (builds,
   gates, long agent work) and 1-2 *light* ones (reviews, quick tests), so heavy work never starves quick feedback;
   remote offload is more lanes, not another model. The limit is the author's and counts every running agent (lanes,
   reviewers, fix rounds, other work's lanes); `state.md` lists what runs, with start times, and nothing starts past the
   limit. Each writing lane and each reviewer that runs checks has its own worktree, install, ports and scratch under
   `.tmp/`.
2. **Match the model to the task, by scope.** {#model-per-task} A lane with a closed item list (a fix lane, a story-only
   lane, a pure split, a fast review, scaffolding, mechanical edits) runs a fast model at low effort; a lane that builds
   a design's core (a language, resolution, layout, runtimes on untrusted input), design, reconciliation, anything
   ambiguous and the independent review of a core lane run the strongest at high effort. The plan and each brief name
   model and effort, the findings row records them, the orchestrator's look and checks gate both alike, and a task that
   turns out ambiguous is escalated, not guessed through.
3. **Name the root nouns before fanning out.** {#name-the-roots} The gate that agrees a design's map names its root
   concepts (the top of the state, the core kinds) before lanes start, and asks the author when a name is not obvious
   (§design-first). Why: a name changed after lanes spread it is a rename across every lane and merge.
4. **Parallel writers work in worktrees; the orchestrator writes shared files.** {#commit-often} Agents commit often
   (§small-green-commits). Parallel writers each work in their own worktree on a lane branch; a lane that only reads, or
   the only writer at the time, may work in the main checkout. Lanes commit only on their branch; the orchestrator does
   every merge, both ways, and removes finished worktrees. Before fanning out it names the files each lane touches:
   lanes sharing a **hot file** (a registry, a status document, a version line) are serialized, the file is split so
   each feature registers itself, or each lane writes an entry of its own (stating any version bump it needs) that the
   orchestrator folds in at the merge, in a commit of its own. A shared helper is named first, owned by the first to
   merge. A lane stops only processes it started, writes only in its worktree and its scratch, and shares no tool state
   (the stash); the orchestrator's scripts stay out of its reach. Why: shared files and state are where parallel work
   collides.
5. **Merge one lane at a time, and gate between.** {#merge-one-at-a-time} The orchestrator merges one finished lane and
   runs the full gate before the next, so a failure has one cause; it brings the work branch into a lane only at a stop
   the lane reported. A clean textual merge is not proof: when two lanes changed what a shared function means, the merge
   runs the other lane's checks too. The merge message names the review findings refuted or left unfixed, and why, so a
   later fold-back can weigh the review. The orchestrator keeps reports and exit codes in its context, not diffs.
6. **The orchestrator owns the lockfile and other generated shared files.** {#orchestrator-owns-lockfile} A lane
   declares the dependencies it needs in its package's manifest and may install locally, but never stages the lockfile
   or a file generated from many packages (codegen output, export indexes, aggregated registries); the orchestrator
   regenerates them once at each merge and commits that alone, and the gate installs without changing them. Why: every
   lane that adds a dependency rewrites the same file.
7. **A lane doesn't loosen its own checks.** {#lanes-dont-edit-checks} Tests, snapshots, baselines, guard thresholds,
   and lint, formatter and gate configuration are hot files for every lane: a lane that needs one changed stops and
   reports why, and the orchestrator decides with the author (§refusals-to-author). The report lists every assertion
   deleted, skipped or weakened. Why: an exit code is only as honest as the checks behind it.
8. **A lane reports; it doesn't write the work folder.** {#lanes-report-not-write} The work folder is the
   orchestrator's; a lane gets its brief in its instructions, and its drift, gaps and followups come back in its report,
   which the orchestrator files. Why: lanes writing `state.md` and the gap lists would make them hot files, and a lane's
   copy is stale from the moment its branch is cut.
9. **A lane's brief names every check its area requires.** {#briefs-name-checks} Its package's tests, the guards the
   area's rules demand (scaling for a bulk path, allocation for a hot one) and the repository's gate, one committed
   script the lane runs rather than a list it transcribes. The report says each check's exit code and what the lane left
   out, refused or fell back on; a report that omits one is not green. At the merge the orchestrator ticks each item of
   the brief against a test or story by name; an item without one is not landed and goes back to the lane or into the
   followups. Why: what the brief doesn't name the lane doesn't run, and what a report says done may not exist.
