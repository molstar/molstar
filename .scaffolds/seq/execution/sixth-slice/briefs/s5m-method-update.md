# Lane S5M: the method follows upstream's workflow files, and the sixth slice's fold-back applied

**Model**: strong. **Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch
`seq/s5m`.

**Why**:
- The author, 2026-11-23: "do we need to update .scaffolds/_method too with the checklist.md, updated plan, etc
  templates; since this changed in the upstream version?". The re-synced package (`c1ccd14`) changed its workflow
  files: §plan-and-checklist and §state-file in `.scaffolds/principles/_core/workflow.md`, the templates in
  `_core/patterns/workflow-files.md` (plan, checklist, state, findings, an optional status, the loop reading the
  checklist, the finalize script).
- The author, on the sixth slice's fold-back (`../foldback.md`, B1–B5): "apply all 5, and with questions go as
  recommended".

**Read first**: the package's `_core/workflow.md` (§plan-and-checklist, §state-file), `_core/patterns/workflow-files.md`
(whole), `.scaffolds/AGENTS.md`, `.scaffolds/_method/` (all), `../foldback.md`.

**Task**, a commit at each step, each numbered item with its proof in your report (B3 of the fifth slice):

1. **The workflow files**, in `.scaffolds/_method/`, following the package's templates and Seq\*'s tightenings
   (the instance manifest's Adaptations):
   - `plan.md`: what and why, written once; a re-plan a dated line; no status, no ticks;
   - `checklist.md` (new): one line per item, ticked with its commit; the gates, the fix lanes, the novice pass, the
     measurements and FINALIZE are items;
   - `state.md`: the current state only, rewritten in place, about a screen;
   - `findings.md` (new, append-only): the lane rows and gate findings, by kind (what `state.md`'s findings table was);
   - `status.md` (optional, for an instance): built and open now, edited in place;
   - `loop.md`: reads `checklist.md` and ticks it with the commit; rewrites `state.md` (never appends); findings to
     `findings.md`; Finalize starts from the checklist's unticked items;
   - `README.md` and `.scaffolds/AGENTS.md` say so.
2. **The finalize check** (`_method/_checks/closed-work.ts`): a closed folder with a `checklist.md` fails on an unticked
   item; a note (never a failure) when an open folder's `state.md` is past 600 words, as the package's does. Folders
   closed before checklists existed (no `checklist.md`) pass as before. A test or a `--self-check` case for each.
3. **The fold-back** (`../foldback.md`):
   - B1: `designs-in-specs` fails a pointer whose text starts more than one heading or bold lead in its spec, naming
     the candidates. Fix the pointers it then finds in `.scaffolds/seq/design/` (the closed slices' files aside);
   - B2: `brief.md` and the loop's review list ask a lane that rewrites text across many files for its identifiers
     before and after; `ec-ids` gets `--compare <ref>`, failing an `EC-<n>` cited at `<ref>` and cited nowhere now;
   - B3: `_method/README.md`: re-syncing an instance, each dropped rule names the package id and quotes the clause that
     covers it; a weaker clause leaves a one-line tightening;
   - B4: `_method/README.md` and `.scaffolds/AGENTS.md`: an instance keeps only what is its own and how it applies a
     package rule, within the rule budget; updates re-sync from upstream, never fork;
   - B5: `loop.md`: the orchestrator runs commands with absolute paths or `-C`, never a `cd` into a lane's worktree;
   - Q1: `.scaffolds/AGENTS.md`: rule ids are checked per scaffold; an instance's rule moves upstream when a second
     scaffold needs it.
4. **The close-out plan** (`.scaffolds/seq/execution/close-out/plan.md`, a draft, not started): split into `plan.md`
   (what and why) and `checklist.md` (its items), the next slice's first use of the new shape.
5. **A proposal for upstream**: `.scaffolds/seq/execution/sixth-slice/drafts/upstream-proposal.md`, B3 and B4 written
   as the package would state them (a rule or two in `adopting.md`, within its rule budget), for the author to take
   to the upstream repository. Never edit that repository.

**Where**: a worktree on `seq/s5m`. Edit only `.scaffolds/_method/`, `.scaffolds/AGENTS.md`,
`.scaffolds/seq/_checks/ec-ids.ts`, the designs' pointers that B1 finds, the close-out folder, and the draft in step 5.
Never the package (`.scaffolds/principles/`), no code, no hot file, never this slice's `state.md`, `plan.md` or
`loop.md`. Write only in your worktree and your own temporary space. No network.

**Checks**: every command in `.scaffolds/seq/AGENTS.md`'s Checks (after `git add .scaffolds`), each new check's
self-check, `pnpm format:check`, `pnpm lint`.

**Report** (under 300 words): commits; each numbered item with its proof; the pointers B1 found; anything you could not
tell.
