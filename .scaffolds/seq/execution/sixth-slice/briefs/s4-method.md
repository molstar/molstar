# Lane S4: the method's templates and `.scaffolds/AGENTS.md`

**Model**: strong (what is generic and what is Seq\*'s is a judgment). **Base**: `packages-seq` at the orchestrator's
head (`git reset --hard packages-seq` first). Branch `seq/s4`.

**Why** (`../plan.md`, stage 4; the author, 2026-11-23: "put it in some common place, where next similar "scaffolds"
will live"): a next scaffold starts as `.scaffolds/<name>/` from `.scaffolds/_method/`, so the method must be written
down once, without Seq\* in it.

**Read first**: `.scaffolds/seq/AGENTS.md`, `.scaffolds/principles/README.md` and its `_core/` (the package's own rules
on loops, briefs, records: `_core/patterns/workflow-files.md`), `.scaffolds/seq/principles/README.md` (the instance's
adaptations, including "Folded back from the fifth slice"), and the fifth slice's files:
`.scaffolds/seq/execution/fifth-slice/{loop.md,plan.md,state.md,fixes.md,foldback.md}`, a brief or two from its
`briefs/`, and its record (`.scaffolds/seq/records/fifth-slice/`).

**Task**, a commit at each step:

1. **`.scaffolds/_method/`**, templates with Seq\* taken out, each short, with `<placeholders>` and one line a section
   saying what goes there and why:
   - `loop.md`: Lanes at once, Each step (brief, start, check, items landed, look, review, fix round, merge, record),
     Pausing, Running while the author is away, Gates, Checks (as placeholders), Stop and ask, Finalize. It carries the
     fifth slice's fold-back (B1–B10) and the earlier ones as rules, not history.
   - `brief.md`, `plan.md` (stages, gates, measurements as budgets: B8), `state.md`, `fixes.md`, `foldback.md`,
     `record.md` (the record's sections), `drift-note.md`, `design.md` (Intent, Precedents, Decision, Defaults with a
     Mock and a precedent: B10; Configuration; Consequences; Tensions; `## In the specs`).
   - `README.md`: how a scaffold is started from these, run, closed and kept.
2. **`.scaffolds/_method/_checks/`**: the generic checks, each taking the instance folder as an argument
   (`--instance .scaffolds/<name>`) instead of a hard-coded path, and its spec map (designs-in-specs) or baseline as the
   instance's file. Generic candidates: `closed-work`, `defaults-designed`, `designs-in-specs`, `no-design-citations`,
   `file-size`, `features-planned`; Seq\*'s own (`ec-ids`, its baselines and spec map) stay in the instance. Decide each
   and say why in the check's header. The instance then runs the method's checks with its arguments: update
   `.scaffolds/seq/AGENTS.md`'s Checks, the sixth slice's `loop.md` Checks, and remove the instance's copies that the
   method now holds (`git mv` where a file moves, so history follows). Every check passes on the Seq\* instance.
3. **`.scaffolds/AGENTS.md`**: what a scaffold is (a development scaffold's agent context: principles, sources,
   designs, work, records), the layout (`principles/`, `_method/`, one folder per scaffold), how one is started from
   `_method/`, run, closed and kept; the rules: kept after its close and editable, never cited by code, comments, tests
   or specs (the author, 2026-11-23); a module of an instance (such as Seq\*'s `queries` and `agent-interface`, from
   meso\*) moves into the principles package only when a second scaffold needs it (the author, 2026-11-23: "keep inside
   for now, once we have more scaffolds we generalize"); `.v6/` is not a scaffold yet ("Moves later").
4. **Link and command check**: every relative link and every `node .scaffolds/...` command in `.scaffolds/` (outside
   closed slices' folders and records) resolves; add the link check you used as `_method/_checks/links.ts` (closed
   history reported, not failed).

**Where**: a worktree on `seq/s4`. Edit only `.scaffolds/` (no code, no package, no hot file; never this slice's
`state.md`, `plan.md` or `loop.md` except the Checks section of `loop.md` in step 2). Write only in your worktree and your
own temporary space. No network. A rule you find that contradicts another, or a Seq\*-specific rule you cannot tell
from a generic one, goes in your report as a question.

**Checks**: every check from its new place on the Seq\* instance; `pnpm format:check`; `pnpm lint`; the link check.

**Report** (under 350 words): commits; the templates (names, a line each); each check's decision (method or instance)
and why; what you could not tell generic from Seq\*'s.
