# Workflow files {#workflow-files}

Rules §design-plan-implement, §drift-notes: designs are durable and global (`.design/`), execution state is working memory (`.execution/`): committed with the work, never cited, reconciled into durable places (design, spec, code, issues) before the close and then kept as the work's archive. Each execution session leaves a permanent record of how it went (`.records/`, §permanent-record), the part that is cited.

## Layout

Shown without an instance; in a project with one, `.design/`, `.records/` and `.execution/` are
`.<instance>/design/`, `.<instance>/records/` and `.<instance>/execution/` (adopting, step 0).

```
.design/                       # tracked, global only; a folder per area, one file per topic in it
  repr/
    stream-parsing.md
  sync/
    sync-port.md
  _prototypes/                 # designs of prototypes: what each explored; the product will differ (§iteration)
    sync-sketch/
      sync-sketch.md
_spec/                         # tracked, mirrors code (see spec-mirror.md)
.records/                      # tracked, a record per session, written at the end of each sitting, compaction or handoff (§permanent-record)
  streaming-parser/
    2026-03-02-parser.md
.tmp/                          # GITIGNORED local scratch: databases, run outputs, eval results (§machine-invariant)
.execution/                    # tracked, per piece of work, never referenced from outside; kept after the close as its archive
  streaming-parser/
    plan.md
    loop.md
    state.md
    research.md                # precedents found while designing; kept ones are cited in the design (§precedents)
    gaps.md                    # the three gap lists shown at every gate: built-not-agreed, agreed-not-built, built-not-adopted (§gaps-queue)
    briefs/                    # one brief per lane (see Lane brief)
    followups.md               # draft; each filed as an issue before merge
    drafts/                    # design drafts and the alternatives they weighed, until agreed into .design/
    drift/                     # open drift notes (§drift-triage); each deleted once resolved or turned into a followup
    notes/                     # working notes; generated output (screens, run logs) goes to .tmp/
```

```gitignore
# .gitignore: only the scratch; the work folder is committed
.tmp/
```

A prototype's own `_spec/` holds `residuals.md`: what the next iteration should pick up (open drift, followups, known
bugs, limits), one line each, grouped and cited (§residuals). It replaces filing the prototype's followups as issues.

## Design file

```md
# Streaming parser for representations

## Intent
Consume large state before it is complete.

## Precedents
- Incremental JSON parsers (SAX-style events): taken, top-level forms as the unit; not taken, event callbacks.
- This codebase's whole-document parser (`lib/repr/core`): the printed form and its round trip.

## Tensions
- Streaming vs one canonical parse: settled, the stream yields the same forms the
  whole-document parse would.

## Decision
Incremental top-level-form parser fed with chunks; patches for structural updates.

## Rejected alternatives
- Whole-document parse: blocks on large state.

## Consequences / invariants
parse(print(x)) == x for any chunking.
```

No status, no front matter: the file being in `.design/` means it is agreed. Edit it in place when it changes (git has the history), reflect current code behavior in `_spec/`, and never link to the plan. A draft lives in `.execution/<work-id>/` until agreed.

## plan.md

```md
# Plan: streaming parser
Design: .design/repr/stream-parsing.md

## Steps / checklist
- [ ] parseStream in repr-core  (heavy lane)
- [ ] property tests: any chunking == whole parse (quick)
- [ ] encodingSuite case (quick)
- [ ] spec entry + design updated
- [ ] GATE: unfiltered gate tests green
- [ ] GATE: review
- [ ] FINALIZE: reconcile, write the session record in .records/streaming-parser/, then close .execution/streaming-parser (state.md's last line)

## Testing strategy
Targeted: lib/repr/core quick. Widen to gate after each step group.
Guard: add scaling test (ratio) for parse on 10k vs 40k nodes.
```

## loop.md

```md
# Loop: streaming parser
1. Read state.md. If missing, create it from plan.md.
2. Pick the next unchecked item; prefer parallel lanes for independent items (heavy ≤2, light ≤2).
3. Implement. Run targeted tests first, then the tiers the plan names.
4. If code and design disagree: stop, write a drift note in drift/ (§drift-notes), continue with unaffected items.
5. Update state.md (and tick plan.md) before moving on, and always before stopping or running low on budget.
6. A lane is done when its gate is green and, if it changes what is drawn (a look, a layout, a geometry), it lands
   with a story or example showing it, and its report names the visual pass (light and dark); the reviewer opens it.
   Headless tests check the state a story makes, not what it draws.
7. When all items done: run FINALIZE (below).

FINALIZE: run the entry point's `finalize` task for streaming-parser (`mise run finalize -- streaming-parser` in the
examples); fix everything it reports; write the record; end state.md with "Closed <date>: <record file>"; edit
nothing in the folder after that.
```

## Lane brief

What the orchestrator hands each lane (§briefs-name-checks, §model-per-task, §commit-often). Short; with the rule
sections it names, it is the lane's whole context: a lane reads those, not all of `_core`.

```md
# Lane C: streaming parser, chunk boundaries
Goal: parse(print(x)) == x for any chunking (design: .design/repr/stream-parsing.md).
Read: the representation rules on encodings and streaming (by id), and the pattern README for the parser.
Model: fast model, low effort (well specified)   # or: strongest, high effort (design, review, ambiguous)
Base: work/parser at 3f9c2e1. Worktree: ../repo-wt-c on lane/c; commit on your branch only; never merge.
Environment: install in the worktree; port 5181; scratch in .tmp/lane-c/ (§parallelism).
Touches: lib/repr/core/src/stream.ts, its tests. Not in scope: the printer, the codec, other packages.
Hot files, never staged: the lockfile and generated shared files (declare new dependencies in the package manifest
  and say so in the report, §orchestrator-owns-lockfile); tests outside your files, snapshots, baselines, guard
  thresholds, lint and gate config (§lanes-dont-edit-checks).
Checks, each its own command, green only on exit 0 with a non-zero test count (§green-is-exit-code):
  - the package's tests at the tiers the plan names
  - the scaling guard for the bulk path (chunk count up to 10^5)
  - the gate, one committed script (`mise run gate`; run it, don't transcribe it): dependency boundaries,
    exports, principles ids
Already red at the base: name it in the report; don't fix it unless it's in scope.
Stop when: the checks are green, or the same check has failed twice after a fix, or a step is refused by the
  permission system (report the action and why; don't work around it, §refusals-to-author).
Report: commits; each check with its exit code and test count; assertions deleted, skipped or weakened; what you
  left out, refused or fell back on, each path with its code; drift, gaps and followups (you don't write the work
  folder, §lanes-report-not-write); anything left open.
```

## state.md (concise, written for a cold pickup)

```md
# State: streaming parser   (updated 2026-11-13 15:40, by agent A)
Branch: feat/streaming-parser   Worktree: ../wt-streaming

## Done
parseStream implemented + 3 property tests passing.

## In progress
encodingSuite case for chunked input (src/conformance.ts, half written).

## Next
spec entry; scaling guard test; design updated.

## Decisions taken
Emit forms eagerly at depth 0 (not buffered until complete). Rationale to go in the design on finalize.

## Blockers / open questions
drift/parse-sync.md: the design says parse() is sync only; streaming needs more. Proposed: edit the design (outcome 3); awaiting the author.

## Verify
the entry point's `test` task: `mise run test -- quick lib/repr/core` in the examples  (green at 15:38)
```

Rule of thumb: if a fresh agent cannot continue from this file alone in a few minutes, it is too thin; if it exceeds a page, it is too long.

## followups.md (draft, filed as issues before merge)

```md
- Binary encoding of the tree. Deferred: no consumer yet. Done = encodingSuite passes for it. -> issue
- Keyed diff for reordered lists. Deferred: needs design (identity semantics). -> issue (with the design notes in it)
```

## Session record (`.records/<work-id>/<date>-<name>.md`)

Written at the end of each sitting, compaction or handoff, and before `.execution/<work-id>/` is closed; never edited after (§permanent-record). It never links into the work folder: the folder shares its work id.

```markdown
# <work id>: <short name> (<start date> to <end date>)

## Asked
The owner's requests in order, one line each, key ones quoted.

## Decided
One line each: the decision, "don't reopen" where the owner said so, and the design section it now lives in.

## Delivered
| Lane | What | Commits | Verified by |
|---|---|---|---|

## Learned
What the work taught about the process or the principles, one line each, with what changed because of it (a rule,
a gate check, a pattern) or "proposed" with where the proposal is.

## Left open
Deferred items (who deferred them), drift still unresolved, and where the next session picks up.

## Self-check (when the session changed the principles)
The claim checked, what integrated (merged), and what is carried as a standing tension: the run's delta (§self-check).
```

## Finalize check

The project's own (adopting, step 4); a sketch. It takes the work folder's path, so it runs the same with or without
an instance.

```ts
// finalize <work folder>: lists what is not yet reconciled; exit 0 means "ready: write the record, close the folder"
import { existsSync, readFileSync, readdirSync } from "node:fs";

const dir = process.argv[2];                                              // e.g. .<instance>/execution/<work-id>
if (!dir) { console.error("usage: finalize <work folder>"); process.exit(2); }
if (!existsSync(dir)) { console.error(`${dir} does not exist`); process.exit(2); }
const read = (f: string) => (existsSync(`${dir}/${f}`) ? readFileSync(`${dir}/${f}`, "utf8") : "");
const problems: string[] = [];

if (/^\s*- \[ \](?! FINALIZE)/m.test(read("plan.md"))) problems.push("plan.md has unchecked items: finish them, or move them to followups.md");
if (read("followups.md").trim()) problems.push("followups.md not empty: file each as an issue, then clear it");
const drift = `${dir}/drift`;
for (const f of existsSync(drift) ? readdirSync(drift).filter((f) => f.endsWith(".md")) : [])
  problems.push(`drift/${f} is open: apply its outcome (fix code, update spec, edit the design) or file a followup issue, then delete it`);

if (problems.length) { console.error(problems.join("\n")); process.exit(1); }
console.log(`reconciled: decisions are in the designs, code facts in _spec/; write the record, then close ${dir}`);
```

The merge check is the package's own, [`_checks/merge-gate.ts`](../../_checks/merge-gate.ts): it rejects any tracked
file outside the work folder that names a path inside it or links into it (§execution-unreferenced); the path check
keeps local paths out of the folder as out of everything else (§machine-invariant).

## Why the work folder is kept but never cited

A tracked plans folder makes plans look like durable knowledge and invites references to them; a deleted one loses the
evidence a later fold-back needs (the alternatives a draft weighed, each lane's findings, what a plan expected against
what was built), which the short record cannot hold. So the folder is committed and kept, and the line is drawn by
citation instead of by deletion: designs, `_spec/` and code say what is true, the record says how we got there and is
cited, and the work folder is the archive behind the record, read as evidence and never as a source
(§execution-unreferenced, checked by the merge gate). `.records/` is not a plans folder either: a record is written
after each sitting, as history, and never edited afterwards (asked, decided, delivered, learned, left open).
