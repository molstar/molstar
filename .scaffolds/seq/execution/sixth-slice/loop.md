# Loop: the sixth slice (scaffolds)

Governs the implementation of `plan.md` (§loop-file). The orchestrator is the session that runs this loop; lanes are
subagents it starts. Read this file and `state.md` before every step.

This is the fifth slice's loop with its fold-back applied (B1–B10; the fifth slice's `loop.md`, "Folded back at the
close", and `principles/README.md`, "Folded back from the fifth slice"; the author, 2026-11-23: "Apply all ten"). It is
fitted to a slice whose subject is the scaffold itself: it moves `.scaffolds/seq/` into `.scaffolds/`, renames the principles
package, reconciles the designs into the packages' `_spec`s, and writes the method's templates. There is no new view,
so no mock (B10), no novice pass and no use gate. Its stage 0 is decided (`plan.md`, "Decided at stage 0").

## The slice runs alone

**No other lane, of any slice, runs while this one runs** (`plan.md`: it moves every path the loop, the briefs and the
checks use). Before stage 2 the orchestrator confirms with `git worktree list` that `packages-seq` is the only worktree.
The fifth slice is closed; the close-out plan waits for this slice's record.

## Lanes at once

At most **two lanes at once** in stage 1 (one per package `_spec`; never two on one file), and **one at a time** from
stage 2 on. Stages 2, 3 and 5 are the orchestrator's own work (renames and moves across the whole tree must be one
hand). Each brief names its model by scope (B1, #model-by-scope):
- the fast model (Sonnet) for a spec pass over one package against named designs, a link check, a template copied
  from a named file;
- the strong model for a spec pass that must decide what is a rule and what is history, and for the method's
  templates (stage 4).

The orchestrator lists what runs, with start times, in `state.md`'s In progress, every agent counted (reviewers and
measurements too).

## What a lane may touch here

In this slice, lanes **do** edit scaffold files (this slice's subject), but only those their brief names. They never
edit `state.md`, `plan.md` or this file, and never a hot file (`package.json`s, `tsconfig*.json`, the lockfile,
`biome.json`, `.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`).
**No code change**: a lane that finds code disagreeing with a design writes a drift note, never a fix (the fix is the
close-out's). A lane writes only in its worktree and its own temporary space (B4); the orchestrator keeps its scripts in
its scratchpad and stages by path, never `commit -a` or `.` (B4).

## Each step

1. **Read** `state.md`, `drift/` and the plan's stage. If `state.md` is missing, rebuild it from `plan.md` and the git
   log.
2. **Pick** the next unchecked stage item whose dependencies are done, within the limit above.
3. **Brief**: every lane has its brief in `briefs/`, with:
   - the goal, and the designs and the `_spec` it reads;
   - the model and why (B1), the base commit (`packages-seq`'s head);
   - the files it may edit, and what is out of scope;
   - **numbered items**, each with what proves it landed (a section of the spec, a check's line, a link) (B3);
   - for a spec pass: the package's `_spec` stays a rule and its reason, terse, no history, no `.scaffolds/seq/` path and no
     principle `§` id (the `.seq` rule; the root AGENTS.md, "Keep the docs true"); a decision that is only history stays
     in the design, marked design-only with its reason.
4. **Start the lane** as a subagent in its own git worktree (isolation: worktree), on `seq/<lane>`, first `git reset
   --hard packages-seq`. It commits small green steps on its branch, never merges.
5. **On a lane's report**, never trust it as proof:
   - **check**: run Checks (below) on its branch yourself.
   - **items landed** (B3): tick each numbered item against its proof by name in the findings row; an unticked item goes
     back to the lane.
   - **review**: a separate reviewer in a worktree of its own, given the diff, the designs and the spec, never the
     lane's conclusions. Strong for a spec pass that decided rule against history; fast for a copied template or a
     link fix. It checks:
     - that the spec still says what the code does: a spec claim not in the code is a finding;
     - that nothing in the spec cites `.scaffolds/seq/` or a principle id;
     - that every design decision the brief named is in the spec or marked design-only with its reason.
   - **fix round**: one round by SendMessage; a design question becomes a drift note.
   - **merge** with `--no-ff`, then the record, then remove the lane's worktree and branch.
   - **record**: the lane's row in `state.md`'s findings table (by kind: bug, missing design, look, drift; the model;
     the items ticked).
6. **Classify before logging** (B2): a finding that matches one already found (the same kind of drift in another
   design) becomes the class's item, its rule written once.
7. **Update** `state.md` and tick `plan.md` after every merged item and every stage the orchestrator does itself.
   Commit the work folder (orchestrator only).

## The orchestrator's stages (2, 3, 5)

- **Stage 2, the rename**: `oz` → the principles package, in folder names, text, skill descriptions, the manifest's
  folder map and the checks' `--package` argument. One provenance line keeps the upstream repository and commit, in
  words (§machine-invariant). Done when `git grep -n -i '\boz\b'` over the scaffold finds only that line.
- **Stage 3, the move**: `git mv` (history kept) into the layout of `plan.md`, one commit for the moves and one for the
  path rewrites, so `git log --follow` reads through them. Every path that names `.scaffolds/seq/` is rewritten: the instance's
  links, the checks and their arguments, the skill, `packages/seq/AGENTS.md`'s one pointer, `biome.json`,
  `.prettierignore`, `.gitignore`, the root `AGENTS.md`'s table (a `.scaffolds/` row). Then confirm that nothing in code,
  comments, tests or specs cites the scaffold (the rule; `cli/seqx-render` and `cli/structure-info` match `.seq` as a
  file extension: check each hit).
- **Stage 5, verify**: every check run from its new place on the Seq\* instance; a link check over `.scaffolds/` (no
  dangling relative link; a small script, kept in the method's checks); `.scaffolds/seq/` gone; the closed fifth slice's files
  unchanged but for their moved paths.

Each of them gets a fast reviewer of its own on the diff (moved paths, rewritten links, nothing dropped), before the
next stage starts.

## Gates

After stage 1, after stage 3, and at stage 5:
- all of Checks, green by exit code with a non-zero test count (§green-is-exit-code);
- the scaffold's own checks from their current place;
- the link check (from stage 3 on);
- drift notes and questions, each a yes or no with a recommendation;
- a short report to the author: what moved or changed, the checks, the decisions needed.

The loop continues unless a decision is needed. There is no timing gate: nothing here touches a frame's path, so the
quiet-machine rule does not apply (a lane that changes code was out of scope).

## Checks

From the repository root. Before stage 3 the scaffold's checks run from `.scaffolds/seq/`, after it from `.scaffolds/`
(their new arguments are part of stage 3's commit):

```sh
pnpm test:unit packages/seq/ cli/seqx-render/ packages/core/src/repr   # the specs' config-table tests and repr
pnpm lint
pnpm format:check
node .scaffolds/principles/_checks/principles-ids.ts --package .scaffolds/principles --instance .scaffolds/seq --allow '^(?!\.scaffolds/)'
node .scaffolds/principles/_checks/paths-check.ts --allow '^(?!\.scaffolds/|packages/seq/)'
node .scaffolds/principles/_checks/merge-gate.ts --package .scaffolds/principles --instance .scaffolds/seq
node .scaffolds/_method/_checks/closed-work.ts --instance .scaffolds/seq
node .scaffolds/_method/_checks/defaults-designed.ts --instance .scaffolds/seq
node .scaffolds/_method/_checks/designs-in-specs.ts --instance .scaffolds/seq
node .scaffolds/_method/_checks/no-design-citations.ts --instance .scaffolds/seq
node .scaffolds/_method/_checks/file-size.ts --instance .scaffolds/seq
node .scaffolds/seq/_checks/ec-ids.ts
```

From stage 4 on the generic checks are the method's (`.scaffolds/_method/_checks/`, each given the instance), their
Seq\* parts in `.scaffolds/seq/_checks/*.json`. The scaffold's checks read the git index: stage first. At stage 5 also `pnpm build` and `pnpm check:workspace` (the
moved `.prettierignore`, `biome.json` and `.gitignore` entries) and the Storybook build (the Design pages link to
nothing in the scaffold, checked).

## Stop and ask the author

- a permission refused (report the action and why; never work around it, §refusals-to-author);
- a check red twice after a fix;
- a design decision with no clear home: rule or history;
- a link, a check or a page outside the scaffold that would have to cite `.scaffolds/` to keep working (the rule says
  none may);
- anything in code: this slice changes none;
- `.v6/` (it moves later: the author, 2026-11-23, "Moves later"), the meso\* modules (kept in the instance: "keep
  inside for now").

## Finalize

1. **Reconcile**: `status.md` points at the new places; `drift/` emptied (each note resolved or a line in the
   close-out plan).
2. **Fold back**: `foldback.md` (what moving the scaffold taught, the fifth slice's rules checked, the self-check), the
   author's decision on its own line, applied before the record.
3. Run the gate once more.
4. The record in `.scaffolds/seq/records/sixth-slice/` (the record file only); `state.md` ends with "Closed <date>:
   <record file>".
5. The close-out plan (`close-out/plan.md`, drafted) moves to its place in the new layout with its paths fixed; it
   becomes the next plan.
6. **Clean up**: every merged lane's worktree and branch removed; nothing unmerged deleted without the author's word.
