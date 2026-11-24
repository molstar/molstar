# Lane S4B: Seq\*'s instance re-adapted to the re-synced package

**Model**: strong (what the package now says, what stays Seq\*'s, what is history). **Base**: `packages-seq` at the
orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/s4b`.

**Why** (`../plan.md`, stage 4b; the author, 2026-11-23, on compressing the principles: "the base repo ... did this and
it saved quite a few tokens. let's do something similar here too? or maybe re-adapt for the principles we use in this
repo?"): the package is re-synced at upstream `c1ccd14`, which already has a rule budget (`_core/system.md`,
§rule-budget: a numbered rule states what and why in at most about 180 words, its evidence in a record), a lighter
`_core`, and Seq\*'s slices 3 to 5 folded back. The instance still carries its own copy of those lessons. An agent's
reading path (the guide's steps 0 and 1 and one module) is about 34k tokens; the manifest
(`.scaffolds/seq/principles/README.md`, 32k characters) and `status.md` (15k) are a third of it.

**Read first**: `.scaffolds/principles/_core/` (all of it: §rule-budget, §context-budget, the state and plan rules in
`workflow.md` and `patterns/workflow-files.md`, which upstream's last commit made "a short plan, a checklist beside it,
and the current state only"), `.scaffolds/principles/README.md`, then the instance: `.scaffolds/seq/principles/README.md`,
`.scaffolds/seq/status.md`, `.scaffolds/seq/AGENTS.md`, `.scaffolds/seq/_skills/seq/SKILL.md`, and the records in
`.scaffolds/seq/records/`.

**Task**, a commit at each step:

1. **The manifest** (`.scaffolds/seq/principles/README.md`): keep Provenance, the folder map and the modules table.
   Its rule sections (Adaptations, "Folded back from the first slice" … "from the fifth slice", "From the author …",
   "For the package's home", Open items) become one `## Adaptations` list of what is Seq\*'s and in force now:
   - a rule the package now states (search the package for it: e.g. `defaults-designed`, `module-size`, the design
     self-check, model by scope, classes, mocks, the heavy lane) is dropped from the instance, or reduced to one line
     saying how Seq\* applies it (a number, a path), citing the package's id;
   - a rule that is Seq\*'s and in force keeps its id and is written within the rule budget: what, and why in one
     sentence; the evidence (which slice, which lane, the author's words) goes out, it is in the records already
     (check; where a record lacks it, keep one short "(the author, <date>)" with the quote, nothing more);
   - a rule superseded by a later one is dropped;
   - "For the package's home" and Open items: whatever upstream took is dropped; the rest stays as a short list.
   Every id the instance defined and that anything cites (`principles-ids` tells you) still resolves: in the package,
   or kept here.
2. **`status.md`**: the current state only. What is built: one short paragraph per package (`packages/seq/*`, the
   shared ones) of what it holds now, not per slice; each slice one line with its record's link. Open, Next, Which
   design to read and Language version kept current and short. Nothing that only a record needs.
3. **`AGENTS.md` and the skill**: read them again against the re-synced package; fix what is stale; keep them short.
4. **Measure**: the reading path before and after (characters and about chars/4 tokens): the instance's `AGENTS.md`,
   `status.md`, the seq skill, the package's `principles` and `principles-router` skills, `_core/*.md`, the manifest,
   and `ui/README.md` plus `ui/rules.md`. State the numbers.

**Where**: a worktree on `seq/s4b`. Edit only `.scaffolds/seq/principles/README.md`, `.scaffolds/seq/status.md`,
`.scaffolds/seq/AGENTS.md`, `.scaffolds/seq/_skills/seq/SKILL.md`, and, where an id moved, the files citing it inside
`.scaffolds/seq/` (not the closed slices' folders or records). Never the package (`.scaffolds/principles/`: it is
upstream's), no code, no hot file, never this slice's `state.md`, `plan.md` or `loop.md`. Write only in your worktree and
your own temporary space. No network.

**Checks**: every command in `.scaffolds/seq/AGENTS.md`'s Checks (after `git add .scaffolds`), `pnpm format:check`.

**Report** (under 300 words): commits; the manifest's rules before and after (counts: dropped as the package's,
kept, superseded); status.md before and after; the measured reading path; any rule you could not place.
