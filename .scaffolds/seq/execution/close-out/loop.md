# Loop: version 1's close-out

Governs the implementation of `plan.md` (§loop-file), item by item of `checklist.md` (§plan-and-checklist). The
orchestrator is the session that runs this loop; lanes are subagents it starts. Read this file and `state.md` before
every step.

This is the method's loop ([../../../_method/loop.md](../../../_method/loop.md)) filled in for Seq\*. It descends from
the fifth slice's loop with its fold-back applied (B1–B10, "Apply all ten") and the sixth slice's (B1–B5, "ok apply all
5"; [../sixth-slice/foldback.md](../sixth-slice/foldback.md)). The close-out draws: the mocks, the look, the use gate
and the novice pass all run.

**Stage 0 is the author's.** No lane starts before GATE 0 is ticked: version 1's rows decided, the 14 design questions
answered, the minimal viewer's budget, F3, F9, F11, F14, and the lanes written into `plan.md` and `checklist.md`. Before
it, the orchestrator only prepares: lists the rows, drafts the answers, and asks. The 14 questions have no mocks (the
re-plan of 2026-11-23): each answer is a sentence or a short list of gestures, agreed in text, and the built story is
looked at gesture by gesture before the merge. A new view or a larger interaction still gets its mock
(§mocks-at-design).

## Lanes at once

**Three lanes at once** (the author, 2026-11-22: "we can run up to 3 lanes in parallel now"). What runs together is
limited by files, not by kind (plan, Order and parallelism): never two lanes in one file, never two on one design's
view, on the language, or on the viewer's input layer, never two heavy test runs at once.

- **Every agent counts**: a lane's reviewer and fix round count as that lane, and so does every other agent the
  orchestrator runs (a review of designs, a mock, a measurement, a draft). Before starting one, the orchestrator lists
  what runs, with start times, in `state.md`'s In progress, and takes a finished one off it.
- **The orchestrator stays in its own directory**: every command it runs on a lane's tree names it, by absolute path or
  `git -C <worktree>`, `pnpm -C <worktree>`; never a bare `cd` into a lane's worktree (the sixth fold-back's B5).
- **A lane's model follows its scope** (§model-by-scope): a closed item list and named outputs (a fix lane, a mock, a
  showcase, a pure split, a fast review) runs on Sonnet; a lane that builds a design's core (the language, the model,
  the drag rule, a runtime on untrusted input: the parsers) and its independent review keep the strong one. The brief
  names the model, the findings row records it; the orchestrator's look and checks gate both alike.

## What a lane may touch

- Only the files its brief names, in its own worktree and its own temporary space.
- Never the work folder (`state.md`, `plan.md`, `checklist.md`, `findings.md`, `fixes.md`, `drift/`, this file).
- Never a hot file: the `package.json`s, `tsconfig*.json`, the lockfile, `biome.json`, `.v6/changelogs/`,
  `storybook/packages.ts`, `packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`. It names the change it needs in
  its report; the orchestrator applies it at the merge (§orchestrator-owns-lockfile).
- No network: `pnpm install --offline` only.
- `@molstar/core` and `@molstar/io`: only within the surface its design agreed; anything past it stops the lane.
- The language (SeqViewSpec) changes only in its own lane, with a version bump (now 0.10.0) and its changelog line named
  in the report.

## Each step

1. **Read** `state.md`, `checklist.md`, `gaps.md`, `drift/` and `fixes.md`. If `state.md` is missing, rebuild it from
   `plan.md`, `checklist.md` and the git log.
2. **Pick** the next unticked items of `checklist.md` whose dependencies are merged (plan, Order and parallelism),
   within the limit above, counted from the list of what runs. At a fix point a fix lane (FX25 onward) takes the open
   items of `fixes.md`; **classify before logging**: a finding that matches one already fixed is the class's item (its
   rule in the design or `_spec` it belongs to, a sweep test over every case the rule covers, briefed to a strong lane),
   never a new line per instance.
3. **Brief**: every lane has its brief in `briefs/` ([brief.md](../../../_method/brief.md)), its base set to
   `packages-seq`'s head and the files it names checked to still exist. Its numbered items each name their proof (a test
   title, an output); `features-planned.ts --plan plan.md` passes over them.
4. **Start the lane** as a subagent in its own git worktree (isolation: worktree), on `seq/<lane>` cut from
   `packages-seq`, the brief its whole instructions and the model it names. Its first step is
   `git reset --hard packages-seq` (a worktree can start on a stale base). It commits small green steps on its branch
   (§commit-often), never merges (§lanes-report-not-write).
5. **On a lane's report**, never trust it as proof:
   - **check**: run Checks (below) on its branch yourself, with its own quick-lane time for the folders it touched
     (under 2 s a folder; a lane over it moves its slow tests, a quick representative kept, before its review). A red
     check goes back to the lane once; red twice after a fix stops the lane and comes to the author.
   - **items landed, not reported**: tick each numbered item of the brief against a test or an output by name in the
     findings row; an unticked item goes back to the lane or into `fixes.md`. A design's "Built" note lists ticked items
     only.
   - **look**, when something is drawn: the orchestrator looks itself, every time, in a **visible** pane, light and dark
     (a hidden pane is never evidence). The lane's Storybook runs from its worktree on a port of its own (a temporary
     `.claude/launch.json` entry, `<worktree>/storybook`, `./node_modules/.bin/storybook dev -p <port> --ci --no-open`,
     reverted after, never committed). For an interactive part it drives each gesture of the part's matrix row
     ([design/ux-foundations.md](../../design/ux-foundations.md)) with real pointer events (click, a pointer path for
     hover, a press that could start a drag, keys) and says what each did, compared with the design's mock under
     `Seq* (next)/Design/Mocks/`, or with its agreed gestures where it has none, gesture by gesture. A screenshot is
     evidence only with a second one after the pointer moved off, or with the DOM or the pixels agreeing.
   - **review**: a separate reviewer **in a worktree of its own**, given the diff, the designs, the mock and the
     screenshots, never the lane's conclusions; strong for a design's core, Sonnet for a pure split (with the "no
     behaviour change" check). It checks:
     - the design followed: its rules, its Configuration table, its matrix cells; the `_spec` updated with it;
     - host neutrality: nothing the host owns taken (the undo stack, the keymap, the chrome, the animation loop, a
       structure viewer's state; "engine, not the game", §application-vs-base);
     - minimal code and module size (§module-size): reuse from `@molstar/core` before new helpers; a file past 800 lines
       is a question in the report, past 1500 a failure;
     - a shared piece tested with its two consumers, live in one test;
     - what is drawn against the data; each edge case the brief names in a test title (`EC-<n>`, §ec-tested-by-name);
     - a test over lazy parts awaits them before the gesture it measures (a flake under fake timers suspects a late part
       or a real clock first);
     - a changed test expectation explained against the design;
     - each fix's test failing on the fix's parent (a fix of an unreachable state says so instead);
     - a lane that rewrites text across many files (citations, ids, links): the identifiers it touched (`EC-<n>`, rule
       ids, spec sections), before and after, compared (`ec-ids.ts --compare packages-seq`); one gone with no reason is
       a finding.
   - **fix round**: one round (SendMessage to the lane), or the orchestrator's; a finding that is a design question
     becomes a drift note.
   - **merge** into `packages-seq` with `--no-ff`, after `pnpm -s -C packages/seq/core size` shows no failing row on the
     merged tree (read its rows, never only its exit code). The merge commit names the lane and the review's refuted or
     unfixed findings, and why. Then apply the hot-file changes the lane asked for; lower the file-size baseline
     (`file-size.ts --instance .scaffolds/seq --update`) when a listed file shrank; delete the design's mock when this
     merge lands the built view; `pnpm install --offline`, run the checks on the merged tree, remove the lane's worktree
     and branch. The orchestrator stages by path, never `commit -a` or `.`, and keeps its own scripts in its scratchpad.
   - **record**: the lane's row appended to `findings.md` ([findings.md](../../../_method/findings.md)): the model, the
     items ticked, the review's and the look's findings by kind (bug, missing design, look), who looked and at what; its
     drift, gaps and followups into `drift/`, `gaps.md` and `followups.md`; its line of `checklist.md` ticked with the
     merge commit.
   - **features**: the lane's built rows moved to built, with their links, in
     [design/features.md](../../design/features.md) and the Storybook's `design/features.mdx`, in the merge's commit or
     the next (§features-move-with-merges).
6. **Drift** (§drift-notes): where code and a design disagree, write `drift/<topic>.md`
   ([drift-note.md](../../../_method/drift-note.md)). A drift note is a gap: asked at the next gate, never left for the
   close.
7. **The author's findings**: a look or a bug that needs no design goes into `fixes.md` (what is seen, where, the kind,
   the author's words) for the next fix point; one that blocks their use jumps the queue; one that needs a design is a
   drift note asked at once.
8. **What is drawn** lands with its Storybook story and a headless twin composed as the story composes it, on a fake
   clock.
9. **Designs and their pages move together**: a change to an agreed design updates its file, its page under
   `Seq* (next)/Design` (`packages/seq/storybook/src/_stories/design/`) and the package's `_spec` in one commit.
10. **Tick and rewrite**: tick `checklist.md` with the commit of every item that landed, and rewrite `state.md` in place
    to where the work now stands (never appended), after every merged item, before starting a new lane, and before
    stopping. `plan.md` is not touched; a re-plan is a dated line at its end, with the author's word. Commit the work
    folder (orchestrator only).

## Pausing

No new lane starts; running lanes finish or are stopped with their branch committed; `state.md` is rewritten for a cold
pickup (where it stands, what runs, the next steps, the open questions) and committed.

## Running while the author is away

When the author says to continue without them, the loop runs through the gates without waiting: a decision a gate would
ask for is made as recommended and listed in `gaps.md`; the stop conditions still stop it. Never decided this way: stage
0 (it is the author's), a public API change beyond what a design agreed, a language change, and the use gate, which is
deferred, not passed (§used-before-done).

## Gates

At every `GATE` in `checklist.md`:

- all of Checks on the merged tree, green by exit code with a non-zero test count (§green-is-exit-code);
- the quick lane's time per touched folder (under 2 s) and in total (under 30 s), on a quiet machine (`uptime` under 3)
  **with no lane running**; a test that fails once under load is run again alone and listed as a flake to watch, a
  second failure is a fix lane's;
- the heavy lane alone (§heavy-lane-alone): no other run beside it, and a bound under 50 ms runs isolated in its worker,
  never raised to pass under contention;
- the allocation tiers: a warm frame allocates nothing, a regression fixed before the gate passes;
- the sizes against the plan's budgets (`pnpm -s -C packages/seq/core size`; default, minimal and ring bundles); a
  growth past one is trimmed first, then lazy-loaded only where it is rare, then raised with its cause. Timing and
  allocation bounds are hard: never loosened to land a lane;
- drawing checked against the data in every view the stage touched;
- the interaction matrix's cells the stage fills: a cell built (a test names its gesture), or "no, because" stated;
- the gaps and drift notes raised since the last gate, each a yes or no with a recommendation;
- a short report to the author: what merged, the checks, screenshots from a visible pane, the representatives' numbers,
  the sizes, the decisions needed. The gate's row goes into `findings.md` and its line of `checklist.md` is ticked. The
  loop continues unless a decision is needed.

**The use gate** is the author's: they try the genome views and the ring, each story against its agreed picture (the
fifth fold-back's B10). The orchestrator prepares the stories and a list of what to try; it never ticks the item itself.

**The novice pass** runs after the last lane that draws (standing here: the manifest's Adaptations): a light agent plays
a user in a visible browser with a task script per new view (a pointer path onto each menu and submenu), split by groups
from the start. The orchestrator compares what it finds with the design's mock, gesture by gesture, and reproduces every
bug it reports in a visible pane before placing it. **Use is continuous**: the author's use as lanes merge counts too.
At both, each finding is fixed or placed (`fixes.md`, a lane, `packages/seq/TODO.md`, a drift note) and counted by kind
in a row of `findings.md`.

## Checks

From the repository root. The scaffold checks read the git index: stage first (`git add .scaffolds packages/seq`). The
medium, heavy, perf and build rows run at gates and on a lane that touches what they cover.

```sh
pnpm test:unit packages/seq/ cli/seqx-render/ packages/core/src/repr packages/io/src/common
pnpm test:medium packages/seq/
pnpm test:heavy packages/seq/
pnpm -s -C packages/seq/core test:perf
pnpm -s -C packages/seq/storybook test:perf
pnpm -s -C packages/seq/core size
npx tsc -b
pnpm lint
pnpm format:check
pnpm build
pnpm check:workspace
node .scaffolds/principles/_checks/principles-ids.ts --package .scaffolds/principles --instance .scaffolds/seq --allow '^(?!\.scaffolds/)'
node .scaffolds/principles/_checks/paths-check.ts --allow '^(?!\.scaffolds/|packages/seq/)'
node .scaffolds/principles/_checks/merge-gate.ts --package .scaffolds/principles --instance .scaffolds/seq
node .scaffolds/_method/_checks/closed-work.ts --instance .scaffolds/seq
node .scaffolds/_method/_checks/defaults-designed.ts --instance .scaffolds/seq
node .scaffolds/_method/_checks/designs-in-specs.ts --instance .scaffolds/seq
node .scaffolds/_method/_checks/no-design-citations.ts --instance .scaffolds/seq
node .scaffolds/_method/_checks/file-size.ts --instance .scaffolds/seq
node .scaffolds/_method/_checks/features-planned.ts --instance .scaffolds/seq --plan .scaffolds/seq/execution/close-out/plan.md
node .scaffolds/_method/_checks/links.ts
node .scaffolds/seq/_checks/ec-ids.ts
```

`pnpm lint` has five warnings from before; a sixth is a finding.

## Stop and ask the author

- a permission refused (report the action and why; never work around it, never hand it to another lane,
  §refusals-to-author);
- a check red twice after a fix;
- a design question a drift note cannot settle by the designs;
- a public API change beyond what a design agreed, or a change to `@molstar/core` or `@molstar/io` past the agreed
  surface;
- a language change outside its own lane, or without its version bump;
- the network (a fixture recorded from a service is asked before the first fetch; tests never reach it);
- a test, performance or size budget that cannot be met without changing a design;
- anything the host owns taken: the undo stack, the keymap, the chrome, the animation loop, a structure viewer's state.

## Finalize

The work folder stays in `.scaffolds/seq/execution/close-out/`, closed and kept; `records/close-out/` holds only the
record file. It starts from `checklist.md`'s unticked items, each run, cut by the author's word or carried to the next
plan, and ticked so ("cut (the author, <date>)", "carried"). Then:

1. **Reconcile**: designs and their Storybook pages; the packages' `_spec`s and their Configuration tables; `drift/`
   emptied (each note resolved and deleted, or turned into a followup); `followups.md` filed (each a line in
   `packages/seq/TODO.md`) and cleared; `status.md`, the features rows, the edge cases' outcomes (`ec-ids.ts --worklist`
   empty or each item given an outcome); `fixes.md`'s open items fixed or placed; version 1's rows each built, cut or
   deferred.
2. **Fold back**: `foldback.md` ([foldback.md](../../../_method/foldback.md)). The author decides on its own line; what
   was decided is applied and committed before the record.
3. Run the gate once more.
4. Write the record ([record.md](../../../_method/record.md)) in `.scaffolds/seq/records/close-out/`, tick FINALIZE, end
   `state.md` with "Closed <date>: <record file>" in the same commit, and edit nothing in this folder after that.
5. Write the next plan with what this work taught.
6. **Clean up**: every merged lane's worktree and branch removed; a parked lane's kept and named in the record; nothing
   unmerged deleted without the author's word. Nothing is pushed.
