# Loop: the fifth slice

Governs the implementation of `plan.md` (§loop-file). The orchestrator is the session that runs this loop; lanes are
subagents it starts. Read this file and `state.md` before every step. This is the fourth slice's loop with its fold-back
applied (the fourth slice's `foldback.md`, B1–B11, merged by the author, 2026-11-21: "All 11") and the author's decisions
since.

## Lanes at once

**Three lanes at once** (the author, 2026-11-22: "we can run up to 3 lanes in parallel now"; one from 2026-11-21
while usage was tight, two before it). Each lane strong or fast by its work. Strong (the strongest model, high
effort) for design, a model or anything ambiguous; fast (the fast model) for specified work: fixes, stories, tests,
read-only. The author changes this line; the loop never does. What runs together is limited by files, not by kind
(plan, Order and parallelism): never two lanes on the language or on core's store and resolve, never two in one file,
never two heavy test runs at once. A lane's reviewer and fix round count as that lane, and so does every other agent
the orchestrator runs (a review of designs, a measurement, another slice's draft): before starting one, the
orchestrator lists what runs, with start times, in `state.md`'s In progress (B3: the fourth slice ran three against a
limit of two without anyone deciding it).

## Each step

1. **Read** `state.md`, `gaps.md`, `drift/` and `fixes.md`. If `state.md` is missing, rebuild it from `plan.md` and
   the git log.
2. **Pick** the next unchecked items whose dependencies are merged (plan, Order and parallelism), within the limit
   above, counted from the list of what runs (B3). At a fix point (plan, Fix points) an FX lane takes the open items of `fixes.md`.
3. **Brief**: every lane has its brief in `briefs/` (named in the plan). Before it starts, set the brief's base to
   `packages-seq`'s head, and check the files it names still exist (the third slice's tail and the module split move
   them). Every brief has:
   - the goal and the designs to read, the sections of `ux-foundations.md` whose matrix cells it fills;
   - the model, the base commit;
   - the files it touches, what is out of scope, the hot files it must not stage (plan, Hot files);
   - **the agreed surface** for a change in the base library (`@molstar/core`, `@molstar/io`): anything past it stops
     the lane (Stop and ask);
   - **the two consumers** its shared pieces are tested with, live in one test (for structure syncing: Seq\* and the
     fake viewer);
   - **its performance representative**: one quick case, a heavy bound under `itWithin` for a time, and a case in the
     allocation tier (core's or Storybook's) when it adds to a frame's path;
   - **whose each choice is**: the design's Configuration table copied into the package `_spec`'s table, every key
     the defaults hold present there (the config-keys check);
   - **minimal code** and **module size**: reuse `@molstar/core` and the package's helpers; a file it grows past 800
     lines is a question in its report, past 1500 a failure (§module-size);
   - **drawing against the data**: the example views it touches checked at the fit view and across a tile edge;
   - **the edge cases it covers, named by its tests** (`EC-<n>` in the test's name), so the catalog's coverage is read
     from the tests, not claimed.
4. **Start the lane** as a subagent in its own git worktree (isolation: worktree), on a branch `seq/<lane>` cut from
   `packages-seq`, with the brief as its whole instructions and the model the brief names. A worktree can start on a
   stale base: the lane's first step is `git reset --hard packages-seq` before it branches. It commits on its branch
   only at each small green step (§commit-often), never merges, and never edits the work folder or a hot file
   (§lanes-report-not-write, §orchestrator-owns-lockfile).
5. **On a lane's report**, never trust it as proof:
   - **check**: check out its branch and run its checks yourself (Checks, below; a lane that touches a frame's path
     also runs the allocation tiers, B5), with **its own quick-lane time** for
     the folders it touched (under 2 s a folder; a lane over it moves its own slow tests, a quick representative
     kept, before its review). A red check goes back to the lane once; twice red after a fix stops the lane and comes
     to the author.
   - **look**, when something is drawn: point the `storybook-lane` launch configuration (port 6008, set for the look,
     reverted after, never committed) at the lane's worktree, open its stories in a **visible** pane in light and dark
     (`&globals=theme:dark`), and take screenshots. Confirm with `tabs_context` that the pane is shown; a hidden pane
     is never evidence. **The orchestrator looks itself, every time** (B2), never on the lane's own browser check: for
     a part that is interactive, it drives each gesture of the part's matrix row with real pointer events (click,
     hover with a pointer path, a press that could start a drag, keys) and says what each did (B1). A screenshot is
     evidence only with a second one after the pointer moved off, or with the DOM or the canvas's pixels agreeing
     (B4: the pane's screenshots went stale). A part that mounts DOM over the canvas has a jsdom test that the view
     receives none of its pointer events (`isViewEvent`; B1).
   - **review**: a separate strong-model reviewer, **in a worktree of its own** (never the lane's), gets the diff, the
     designs and the screenshots, never the lane's conclusions. It checks the design followed (its rules, its
     Configuration table, its matrix cells), host neutrality, minimal code and module size, the two-consumer tests,
     drawing against the data, each `EC-<n>` the brief names in a test title (a missing one is a finding; B8), and
     that each fix's test fails on the fix's parent. A pure split (no behaviour change) gets the fast model with the
     "no behaviour change" check (B9).
   - **fix round**: the lane fixes the findings in one round (SendMessage, §reuse-agents), or the orchestrator does. A
     finding that is a design question becomes a drift note.
   - **merge** with `--no-ff`, only after `pnpm -s -C packages/seq/core size` shows no `FAIL` row on the merged tree
     (the script exits 0 on a failing row: read its rows, never its exit code; LA merged at 930.5 of 930 kB so, 2026-11-22).
     **Size bounds are soft** (the author, 2026-11-22: "don't ask me to raise the size limits constantly. try some basic
     stuff to reduce size, lazy load less frequent features, and otherwise keep it growing. the important parts are the
     performance and low allocations"): a row past its bound sends the lane to trim and lazy-load first; what remains, the
     orchestrator raises in `sizes.mjs` and `performance.md` with its cause, without asking. Timing and allocation
     bounds are the hard gates: never loosened to land a lane;
     the merge commit names the lane and the review's refuted or unfixed findings. Then:
     apply the hot-file changes the lane asked for (dependencies, inventory, tsconfig references, exports' `_test`
     nulls, Storybook registration, its `.v6/changelogs/` entry); lower the file-size baseline when a listed file
     shrank (`node .seq/_checks/file-size.ts --update`); install; run the checks on the merged tree; remove the lane's
     worktree and delete its branch.
   - **record**: the lane's row in `state.md`'s findings table (review findings and look findings, by kind: bug,
     missing design, look; who looked, with the stories and themes seen, B2), kept for every lane, not only at gates; its drift, gaps and followups into `drift/`,
     `gaps.md` and `followups.md`.
   - **features**: the lane's built rows in `.seq/design/features.md` and its entries on the Storybook Features page
     (`design/features.mdx`), with the story links, in the merge's commit or the next; a check fails a row still
     "planned" for a merged lane (B7).
6. **Drift** (§drift-notes): where code and a design disagree, write `drift/<topic>.md` (what, where, the options, a
   recommendation). A drift note is a gap: it is asked at the next gate, never left for the close.
7. **The author's findings** (use, a look at a story): a look or a bug that needs no design goes into `fixes.md` (what
   is seen, the story, the kind, the author's words) for the next fix point; one that blocks their use jumps the
   queue; one that needs a design becomes a drift note asked at once. The third slice's findings from the author's
   check of it go into `fixes.md` as they come (the author, 2026-11-20).
8. **What is drawn** lands with its story and the story's headless twin, composed as the story composes it, on a fake
   clock.
9. **Designs and their pages move together**: a change to an agreed design updates its `.seq/design/` file and its
   Storybook Design page in one commit.
10. **Update** `state.md` (done, in progress, next, decided, blockers, how to verify) and tick `plan.md` after every
    merged item, before starting a new lane, and before stopping. Commit the work folder (orchestrator only).

## Pausing

When the author pauses the loop or usage runs short: no new lane starts; running lanes finish or are stopped with their
branch committed; `resume.md` is written in the third slice's shape (where it stands, what is merged, what runs, the
next steps, the open questions), and `state.md` points at it.

## Running while the author is away

When the author says to continue without them, the loop runs through the gates without waiting: a decision a gate
would ask for is made as recommended and listed in `gaps.md`; the stop conditions still stop it; a public API change
beyond what a design agreed is never decided this way; the use gates are deferred, not passed (§used-before-done).

## Gates

At every `GATE` in `plan.md`:

- all of Checks on the merged tree, green by exit code with a non-zero test count (§green-is-exit-code);
- the quick lane's time for each touched folder (under 2 s) and in total (under 30 s), measured on a quiet machine
  (`uptime` under 3) **with no lane running** (B5: every timing miss of the fourth slice was load); a test that fails
  once under load is run again alone and listed as a flake to watch in the TODO, a second failure is a fix lane's;
- the allocation tiers (`pnpm -C packages/seq/core test:perf`, `pnpm -C packages/seq/storybook test:perf`): a warm
  frame allocates nothing; a regression is fixed before the gate passes (the third slice's emphasis frame grew to
  800 B across lanes between gates);
- the sizes (`pnpm -C packages/seq/core size`) against the targets as **bounds**, reset after the third slice's RR (the
  author, 2026-11-20: "reduce render first, then move the size targets"); a growth past one is fixed or brought to the
  author with the feature behind it;
- drawing checked against the data in every example view the stage touched;
- the interaction matrix's cells the stage fills (a cell built, or "no, because" stated);
- the gaps and drift notes raised since the last gate, each a yes or no question with a recommendation;
- a short report to the author: what was merged, the checks, screenshots from a visible pane, the representatives'
  numbers, the sizes, decisions needed. The loop continues unless a decision is needed.

The **novice pass** runs after the last lane that draws, its script with a pointer path onto each menu and submenu
(B10). At a **use gate** and the novice pass: each finding fixed or placed (`fixes.md`, a lane, `TODO.md`, a drift note)
and counted by kind in `state.md`.

## Checks

Each its own command, from the repository root:

```sh
pnpm test:unit packages/seq/                 # and packages/core/src/<touched> when core changed
npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples   # and packages/core when it changed
pnpm -C packages/seq/core typecheck:test     # and the same for viewer and examples
pnpm lint
pnpm format:check
node .seq/oz/_checks/principles-ids.ts --package .seq/oz --instance .seq --allow '^(?!\.seq/)'
node .seq/oz/_checks/paths-check.ts --allow '^(?!\.seq/|packages/seq/)'
node .seq/oz/_checks/merge-gate.ts --package .seq/oz --instance .seq
node .seq/_checks/ec-ids.ts
node .seq/_checks/closed-work.ts
node .seq/_checks/file-size.ts
node .seq/_checks/defaults-designed.ts
node .seq/_checks/features-planned.ts --plan .seq/execution/fifth-slice/plan.md
```

The third slice's fold-back added two of these: the config keys against the `_spec` tables are quick tests
(`packages/seq/{core,viewer}/src/_test/config-table.test.ts`), and `ec-ids.ts` fails a `done` feature whose edge cases
no test names (`--worklist` prints the close's list). The `.seq` checks read the git index: stage first. When each also runs:

- **at each gate:** the Storybook typecheck (`pnpm -C packages/seq/storybook typecheck`), `pnpm test:medium` for the
  touched packages, `pnpm test:heavy packages/seq/` (alone), the allocation tiers, the size check;
- **at a gate after a base-library change (BW's range reads, the standard components) and at finalize:** `pnpm build`,
  then `pnpm check:workspace`;
- **at M4:** the browser measurements in a visible pane at dpr 2 (the plan, stage 5), recorded in `measurements.md`.

Two checks the fourth slice's fold-back added are built by the fifth slice's first fix lane, then run with the rest:
a design whose Consequences name a layer, a layout or a part has a `Defaults` heading (B6, §defaults-designed; for the
designs promoted from the fifth slice on), and a `features.md` row still "planned" for a lane `plan.md` marks merged
fails (B7).

## Stop and ask the author

- a permission refused (report the action and why; never work around it, §refusals-to-author);
- a check red twice after a fix;
- a design question a drift note cannot settle by the designs;
- a public API change beyond what a design agreed ([base-library](../../design/base-library.md)): BW's range reads
  in the asset manager and the standard components stop at their surface;
- the network: recording a service's fixtures (UniProt, PDBe SIFTS, AlphaFold DB) fetches from it once; asked before
  the first fetch, tests never reach it;
- a test, performance or size budget that cannot be met without changing a design;
- a piece that can only work if Seq\* owns what a host owns (the undo stack, the keymap, the chrome, the animation loop,
  a structure viewer's state: structure syncing is a library, never a port).

## Finalize

**The work folder stays** in `.seq/execution/fifth-slice/`, closed and kept: never moved or deleted unless the author
says so (the author, 2026-11-21: "the execution folder stays unless i say otherwise"); `.seq/records/fifth-slice/`
holds only the record file.

Starts with the unticked plan items listed in `state.md`, each run, cut by the author's word, or carried to the next
plan (B11: the fourth slice reached it with two unticked). Then, when every item is ticked or placed:

1. **Reconcile**: designs and their Design pages; `drift/` emptied (each note resolved or a line in
   `packages/seq/TODO.md`); the `_spec`s and their Configuration tables; `status.md`, the features rows, the edge
   cases' outcomes; `fixes.md`'s open items fixed or placed.
2. **Fold back**: write `foldback.md` (what the slice taught with evidence and where it goes; the earlier fold-backs'
   rules held, amended or dropped; the count of what each gate, lane and use found, by kind; the self-check). The
   author decides on its own line; what was decided is applied and committed before the record.
3. Run the gate once more.
4. Write the record in `.seq/records/fifth-slice/` (the record file only), end `state.md` with "Closed <date>: <record
   file>", and edit
   nothing in this folder after that.
5. Write the next plan (version 1's close-out) with what this slice taught.
6. **Clean up the worktrees**: every merged lane's worktree and branch removed (`git worktree remove`, then the branch
   deleted when `git branch --merged` lists it); a parked lane's kept and named in the record; nothing unmerged deleted
   without the author's word.
