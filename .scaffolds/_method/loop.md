# Loop: <work id>

<!-- The template for a work folder's loop file (§loop-file), written from the slices that built Seq*: each rule below
     carries a fold-back the author merged. Copy it to `.scaffolds/<name>/execution/<work-id>/loop.md`, replace every
     <placeholder>, and delete a section's comment once it is filled. Cut a section only where the work has no such
     step (a work that draws nothing has no look, no mock and no novice pass); say so in one line. -->

Governs the implementation of `plan.md` (§loop-file). The orchestrator is the session that runs this loop; lanes are
subagents it starts. Read this file and `state.md` before every step. <!-- Name the loop this one descends from and the
fold-back it applies ("the <previous> slice's loop with its fold-back, B1–Bn, applied"), so the rules have a source. -->

## Lanes at once

<!-- How many lanes run together, and who set it: the author changes this line, the loop never does. -->

**<n> lanes at once** (the author, <date>: "<their words>"). What runs together is limited by files, not by kind
(plan, Order and parallelism): never two lanes in one file, never two on <the parts one lane owns at a time>, never
two heavy test runs at once.

- **Every agent counts**: a lane's reviewer and fix round count as that lane, and so does every other agent the
  orchestrator runs (a review of designs, a measurement, a draft). Before starting one, the orchestrator lists what
  runs, with start times, in `state.md`'s In progress.
- **A lane's model follows its scope**: a closed item list and named outputs (a fix lane, a showcase,
  a pure split, a fast review) runs on the fast model; a lane that builds a design's core (<the language, the model,
  runtimes on untrusted input>) and its independent review keep the strong one. The brief names the model, the
  findings row records it; the orchestrator's look and checks gate both alike.

## Each step

<!-- The orchestrator's cycle for one lane, from pick to record. Every step is a rule a slice paid for; keep its
     order: the look and the review come before the merge, never after. -->

1. **Read** `state.md`, `gaps.md`, `drift/` and `fixes.md`. If `state.md` is missing, rebuild it from `plan.md` and the
   git log.
2. **Pick** the next unchecked items whose dependencies are merged (plan, Order and parallelism), within the limit
   above, counted from the list of what runs. At a fix point a <fix lane> takes the open items of `fixes.md`;
   **classify before logging**: a finding that matches one already fixed is the class's item
   (its rule in the design or `_spec` it belongs to, a sweep test over every case the rule covers, briefed to a strong
   lane), never a new line per instance.
3. **Brief**: every lane has its brief in `briefs/` ([brief.md](brief.md)), its base set to <the integration branch>'s
   head and the files it names checked to still exist.
4. **Start the lane** as a subagent in its own git worktree, on `<prefix>/<lane>` cut from <the integration branch>,
   the brief its whole instructions and the model it names. Its first step is `git reset --hard <the integration
   branch>` (a worktree can start on a stale base). It commits small green steps on its branch (§commit-often), never
   merges, never edits the work folder or a hot file (§lanes-report-not-write, §orchestrator-owns-lockfile), and
   writes only in its worktree and its own temporary space.
5. **On a lane's report**, never trust it as proof:
   - **check**: run its checks yourself on its branch (Checks, below), with its own quick-lane time for the folders
     it touched (<budget> a folder; a lane over it moves its slow tests, a quick representative kept, before its
     review). A red check goes back to the lane once; red twice after a fix stops the lane and comes to the author.
   - **items landed, not reported**: tick each numbered item of the brief against a test or an output by name in the
     findings row; an unticked item goes back to the lane or into `fixes.md`. A design's "Built" note lists ticked
     items only.
   - **look**, when something is drawn: the orchestrator looks itself, every time, in a **visible** pane, light and
     dark (a hidden pane is never evidence). For an interactive part it drives each gesture of the part's matrix row
     with real pointer events (click, a pointer path for hover, a press that could start a drag, keys) and says what
     each did, compared with the design's mock, gesture by gesture. A screenshot is evidence only
     with a second one after the pointer moved off, or with the DOM or the pixels agreeing.
   - **review**: a separate reviewer **in a worktree of its own**, given the diff, the designs and the screenshots,
     never the lane's conclusions; strong for a design's core, fast for a pure split (with the "no behaviour change"
     check). It checks:
     - the design followed: its rules, its Configuration table, its matrix cells;
     - host neutrality: nothing the host owns taken (<the undo stack, the keymap, the chrome, the animation loop>);
     - minimal code and module size (§module-size): reuse before new helpers; a file past 800 lines is a question in
       the report, past 1500 a failure;
     - a shared piece tested with its two consumers, live in one test;
     - what is drawn against the data; each edge case the brief names in a test title;
     - a test over lazy parts awaits them before the gesture it measures (a flake under fake timers suspects a late
       part or a real clock first);
     - a changed test expectation explained against the design;
     - each fix's test failing on the fix's parent (a fix of an unreachable state says so instead).
   - **fix round**: one round (SendMessage to the lane), or the orchestrator's; a finding that is a design question
     becomes a drift note.
   - **merge** with `--no-ff`, after the size check shows no failing row on the merged tree (read its rows, never only
     its exit code). The merge commit names the lane and the review's refuted or unfixed findings, and why. Then apply
     the hot-file changes the lane asked for; lower the file-size baseline (`file-size.ts --update`) when a listed file
     shrank; delete the design's mock when this merge lands the built view; install, run the checks on the merged
     tree, remove the lane's worktree and branch. The orchestrator stages by path, never `commit -a` or `.`, and keeps its own scripts in its scratchpad.
   - **record**: the lane's row in `state.md`'s findings table: the model, the items ticked, the review's and the
     look's findings by kind (bug, missing design, look), who looked and at what; its drift, gaps and followups into
     `drift/`, `gaps.md` and `followups.md`.
   - **features**: the lane's built rows in <the features table> moved to built, with their links, in the merge's
     commit or the next.
6. **Drift** (§drift-notes): where code and a design disagree, write `drift/<topic>.md` ([drift-note.md](drift-note.md)).
   A drift note is a gap: asked at the next gate, never left for the close.
7. **The author's findings**: a look or a bug that needs no design goes into `fixes.md` (what is seen, where, the kind,
   the author's words) for the next fix point; one that blocks their use jumps the queue; one that needs a design is a
   drift note asked at once.
8. **What is drawn** lands with its <example view> and a headless twin composed as the <example view> composes it, on a
   fake clock.
9. **Designs and their pages move together**: a change to an agreed design updates its file and <its published page>
   in one commit.
10. **Update** `state.md` and tick `plan.md` after every merged item, before starting a new lane, and before stopping.
    Commit the work folder (orchestrator only).

## Pausing

<!-- What happens when the author pauses or usage runs short, so a cold session can pick the work up. -->

No new lane starts; running lanes finish or are stopped with their branch committed; `resume.md` is written (where it
stands, what is merged, what runs, the next steps, the open questions), and `state.md` points at it.

## Running while the author is away

<!-- What the loop may decide alone. Keep the limits: an API change and a use gate are never decided this way. -->

When the author says to continue without them, the loop runs through the gates without waiting: a decision a gate
would ask for is made as recommended and listed in `gaps.md`; the stop conditions still stop it; a public API change
beyond what a design agreed is never decided this way; use gates are deferred, not passed (§used-before-done).

## Gates

<!-- What every GATE in the plan runs and reports. Timing is read on a quiet machine only: every timing miss the slices
     recorded was load. -->

At every `GATE` in `plan.md`:

- all of Checks on the merged tree, green by exit code with a non-zero test count (§green-is-exit-code);
- the quick lane's time per touched folder and in total, on a quiet machine (`uptime` under <n>) **with no lane
  running**; a test that fails once under load is run again alone and listed as a flake to watch, a second failure is
  a fix lane's;
- the heavy lane alone (no other run beside it); a bound under <n> ms runs isolated in its worker, and no bound is
  raised to pass under contention;
- the allocation tiers: a warm frame allocates nothing, a regression fixed before the gate passes;
- the sizes against the plan's budgets (a number a stage, from its designs); a growth past one is trimmed first, then lazy-loaded
  only where it is rare, then raised with its cause. Timing and allocation bounds are hard: never loosened to land a
  lane;
- drawing checked against the data in every view the stage touched;
- the interaction matrix's cells the stage fills: a cell built (a test names its gesture), or "no, because" stated;
- the gaps and drift notes raised since the last gate, each a yes or no with a recommendation;
- a short report to the author: what merged, the checks, screenshots from a visible pane, the representatives' numbers,
  the sizes, the decisions needed. The loop continues unless a decision is needed.

**The novice pass** runs after the last lane that draws: a light agent plays a user in a visible browser with a task
script per new view (a pointer path onto each menu and submenu), split by groups from the start. The orchestrator
compares what it finds with the design's mock, gesture by gesture, and reproduces every bug it reports in a
visible pane before placing it. **Use is continuous**: the author's use as lanes
merge is the use gate. At both, each finding is fixed or placed (`fixes.md`, a lane, the TODO, a drift note) and
counted by kind in `state.md`.

## Checks

<!-- The exact commands, each its own line, from the repository root: the repository's, the principles package's three,
     the method's (each with `--instance .scaffolds/<name>`) and the instance's own. Say which read the git index (stage
     first) and which run only at gates or finalize. -->

```sh
<the repository's tests for the touched folders>
<the type check, lint, format check>
node .scaffolds/principles/_checks/principles-ids.ts --package .scaffolds/principles --instance .scaffolds/<name> --allow '<regex>'
node .scaffolds/principles/_checks/paths-check.ts --allow '<regex>'
node .scaffolds/principles/_checks/merge-gate.ts --package .scaffolds/principles --instance .scaffolds/<name>
node .scaffolds/_method/_checks/<check>.ts --instance .scaffolds/<name>
node .scaffolds/<name>/_checks/<check>.ts
```

## Stop and ask the author

<!-- What stops the loop. The first four held in every slice; add the work's own (a dependency it may not change, a
     download, what the host owns). -->

- a permission refused (report the action and why; never work around it, §refusals-to-author);
- a check red twice after a fix;
- a design question a drift note cannot settle by the designs;
- a public API change beyond what a design agreed;
- the network (a fixture recorded from a service is asked before the first fetch; tests never reach it);
- a test, performance or size budget that cannot be met without changing a design;
- <what the host owns that the work may never take>.

## Finalize

<!-- The close, in order. The fold-back is a step of every close; the work folder stays. -->

The work folder stays in `.scaffolds/<name>/execution/<work-id>/`, closed and kept, never moved or deleted unless the
author says so; `records/<work-id>/` holds only the record file. It starts with the unticked plan items listed in
`state.md`, each run, cut by the author's word, or carried to the next plan. Then:

1. **Reconcile**: designs and their pages; `drift/` emptied (each note resolved and deleted, or turned into a
   followup); `followups.md` filed (each an issue or a line in <the TODO>) and cleared; the `_spec`s and
   their Configuration tables; `status.md`, the features rows, the edge cases' outcomes; `fixes.md`'s open items fixed
   or placed.
2. **Fold back**: `foldback.md` ([foldback.md](foldback.md)). The author decides on its own line; what was decided is
   applied and committed before the record.
3. Run the gate once more.
4. Write the record ([record.md](record.md)) in `.scaffolds/<name>/records/<work-id>/`, end `state.md` with
   "Closed <date>: <record file>", and edit nothing in this folder after that.
5. Write the next plan with what this work taught.
6. **Clean up**: every merged lane's worktree and branch removed; a parked lane's kept and named in the record;
   nothing unmerged deleted without the author's word.
