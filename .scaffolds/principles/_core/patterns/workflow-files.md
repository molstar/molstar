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
    plan.md                    # what and why, written once (§plan-and-checklist)
    checklist.md               # one line per item, ticked with its commit
    loop.md
    state.md                   # the current state only, rewritten in place (§state-file)
    findings.md                # append-only: a row per lane and gate, by kind (when work runs in lanes)
    research.md                # precedents found while designing; kept ones are cited in the design (§precedents)
    gaps.md                    # the three gap lists shown at every gate: built-not-agreed, agreed-not-built, built-not-adopted (§gaps-queue)
    briefs/                    # one brief per lane (see Lane brief)
    followups.md               # draft; each filed as an issue before merge
    foldback.md                # when the work taught something: proposals, earlier rules re-checked, the self-check, the author's decision
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

## Self-check
Claim: any chunking of a printed tree parses to the same forms. Premise negated: chunks arrive in order. Added: out
of order (a retried request) is refused with a coded error. Carried: a stream that never ends holds memory.

## Spelled while building (lane C, 2026-03-01; provisional, for the author)
- A chunk boundary inside a string escape is buffered, not split.
```

No status, no front matter: the file being in `.design/` means it is agreed. Edit it in place when it changes (git has the history), reflect current code behavior in `_spec/`, and never link to the plan. A draft lives in `.execution/<work-id>/` until agreed.

## plan.md

Written once at agreement; a re-plan adds a dated line, it doesn't rewrite. No status here.

```md
# Plan: streaming parser
Design: .design/repr/stream-parsing.md

## Stages
1. parseStream in repr-core (heavy lane; strongest model, high effort).
2. Property tests and the encodingSuite case (fast model, low effort), beside 1 once its API holds.
3. Spec entry and design update; review; finalize.

## Testing strategy
Targeted: lib/repr/core quick. Widen to gate after each stage.
Guard: add scaling test (ratio) for parse on 10k vs 40k nodes.

## Risks / out of scope
A stream that never ends (carried in the design). Binary encoding: out of scope.

Re-planned 2026-03-01: stage 2 split in two after the review (the author: "split it").
```

## checklist.md

One line per item, ticked with its commit; the loop and the finalize check read it.

```md
- [x] parseStream in repr-core (3f9c2e1)
- [x] property tests: any chunking == whole parse (5a1d0b4)
- [ ] encodingSuite case
- [ ] spec entry + design updated
- [ ] GATE: unfiltered gate tests green
- [ ] GATE: review
- [ ] FINALIZE: reconcile, the fold-back question, the session record, then close the folder (state.md's last line)
```

## loop.md

```md
# Loop: streaming parser
1. Read state.md. If missing, create it from plan.md and checklist.md.
2. Pick the next unticked item of checklist.md; prefer parallel lanes for independent items (heavy ≤2, light ≤2).
3. Implement. Run targeted tests first, then the tiers the plan names.
4. If code and design disagree: stop, write a drift note in drift/ (§drift-notes), continue with unaffected items.
5. Tick checklist.md with the commit, and rewrite state.md to where the work now stands (not a log of what happened),
   before moving on and always before stopping or running low on budget.
6. A lane is done when its gate is green and, if it changes what is drawn (a look, a layout, a geometry), it lands
   with a story or example showing it, and its report names the visual pass (light and dark); the reviewer opens it
   before the merge, not after. At the merge the orchestrator looks itself, every time: for an interactive part it
   drives each gesture of the part's matrix row with real pointer events and records what each did (a screenshot shows
   a result, never stands for a gesture; after an action, re-shoot or read the DOM, since a pane can paint late), and
   compares with the design's mock where there is one. It ticks the brief's items against tests by name (a design's
   "built" note lists only ticked items), and the lane's row in findings.md says who looked, what they saw, the model the
   lane ran on, and the findings by kind. A reviewer's brief lists the review's finding classes as its checklist.
   Use: the author uses each part as it merges (see "Use and the novice pass", below).
   Headless tests check the state a story makes, not what it draws.
7. When all items done: run FINALIZE (below).

FINALIZE: run the entry point's `finalize` task for streaming-parser (`mise run finalize -- streaming-parser` in the
examples), starting from checklist.md's unticked items (each run, cut by the author's word, or carried to the next plan,
in state.md); fix everything it reports; ask the fold-back question (anything learned about how work runs or code is
shaped? if so, write foldback.md and bring it to the author); write the record, its Learned section ending with
"Fold-back: none" or "Fold-back: proposed (foldback.md)"; end state.md with "Closed <date>: <record file>"; edit
nothing in the folder after that.
```

## Use and the novice pass

How §used-before-done runs. The author uses each part as it merges, as its people would, on real data, with the agent
standing by; before the close, a look at what changed since the last use. In a prototype run by one session, the
author's use is the last gate before the close.

- **A novice pass** (optional): a light agent plays a user in a visible browser with a short task script fitted to the
  work (a reader's task for a view, a host's for an integration), and compares what it sees with the design's mock
  where there is one. It runs after the last lane that draws, split into groups from the start, so one budget does
  not end a third of the way.
- **Its bugs are reproduced first**: the orchestrator reproduces each one in a visible pane before it is placed; a
  tool's click is not a person's.
- **Observations state their conditions** (a visible or hidden page, the browser, the pixel ratio) and are reproduced
  headless, on a fake clock where time matters, before a fix starts. A hidden page runs no animation frames, so it is
  inconclusive for anything that moves.
- **Findings are counted by kind** (a bug, a missing design, a look) in the lane's row in findings.md as they come, and at
  the close per gate, so a fold-back can tell which gate catches what. Why: a use gate only at the end finds the
  missing designs all at once, when each costs a new lane.

## Lane brief

What the orchestrator hands each lane (§briefs-name-checks, §model-per-task, §commit-often). Short; with the rule
sections it names, it is the lane's whole context: a lane reads those, not all of `_core`.

```md
# Lane C: streaming parser, chunk boundaries
Goal: parse(print(x)) == x for any chunking (design: .design/repr/stream-parsing.md).
Read: the representation rules on encodings and streaming (by id), and the pattern README for the parser.
Model: fast model, low effort (well specified)   # or: strongest, high effort (design, review, ambiguous)
Base: work/parser at 3f9c2e1. Worktree: ../repo-wt-c on lane/c; commit on your branch only; never merge.
Environment: install in the worktree; port 5181; scratch in .tmp/lane-c/ (§parallelism). Write only in your worktree
  and your own scratch.
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
Report: commits; each check with its exit code and test count; a version bump the change needs (as an entry, never
  the version line itself); assertions deleted, skipped or weakened; what you
  left out, refused or fell back on, each path with its code; drift, gaps and followups (you don't write the work
  folder, §lanes-report-not-write); anything left open.
```

## state.md (the current state only, about a screen)

Rewritten in place, never appended: what happened is in git, the checklist's ticks, `findings.md` and the record;
decisions are in their design. The fold-back check notes a `state.md` past 600 words.

```md
# State: streaming parser   (updated 2026-11-13 15:40, by agent A)
Branch: feat/streaming-parser   Worktree: ../wt-streaming

## Where it stands
parseStream and its property tests merged; the encodingSuite case half written (src/conformance.ts).

## Next
spec entry; scaling guard test; design updated.

## Blockers / open questions
drift/parse-sync.md: the design says parse() is sync only; streaming needs more. Proposed: edit the design (outcome 3); awaiting the author.

## Verify
the entry point's `test` task: `mise run test -- quick lib/repr/core` in the examples  (green at 15:38)
```

## findings.md (when work runs in lanes)

Append-only, one row per merged lane and per gate, so a fold-back can count what each catches (§used-before-done).

```md
| Lane or gate | Merged | Bug | Missing design | Look | Who looked | Model |
|---|---|---|---|---|---|---|
| C | 2026-03-01 | 2 (review: a split escape; an unbounded buffer) | 0 | 0 | orchestrator, light and dark | strongest |
```

## status.md (optional, at the instance's root)

Some projects keep one page that routes a new session: what is built and open **now**, and which design or spec to
read for each. It is edited in place at each merge (as an entry folded in by the orchestrator, §commit-often), never
appended: what was built when is in the records.


Rule of thumb: if a fresh agent cannot continue from this file alone in a few minutes, it is too thin; if it exceeds a page, it is too long.

## followups.md (draft, filed as issues before merge)

```md
- Binary encoding of the tree. Deferred: no consumer yet. Done = encodingSuite passes for it. -> issue
- Keyed diff for reordered lists. Deferred: needs design (identity semantics). -> issue (with the design notes in it)
```

## foldback.md (when the work taught something)

Written at the close when the fold-back question's answer is yes (§reconcile-before-merge): each proposal with its
evidence (a commit, a record, a measurement) and where it goes (the principles, a design, a `_spec`, a check). Writing
it is the work's; the author's decision is its last line, `Decided by the author, <date>: <merged | merged in part | not merged>: "<their words>"`, and nothing in
it is applied outside the folder before that line. An instance where learning is the point may require the file at
every close (its manifest's adaptations). A later session on closed work opens a new work id, or reopens the folder by
a new section in `state.md`, never by rewriting the old ones. Why the decision line: lessons from work meant to be
thrown away would bend the principles if merged by default. Why the earlier rules are re-checked: a rule re-checked at
each fold-back finds out when it was wrong.

```md
# Fold-back: streaming parser

## Proposals
1. **A parser's chunk tests cover a split inside every token kind.** (P: testing; C)
   - What happened: a split inside a string escape passed every test and failed in use (`3f9c2e1`).
   - Proposal: the chunking property test enumerates token kinds; where it goes: the conformance suite.

## Earlier rules, checked against this work
- The visual pass before the merge (from the last fold-back): held; no drawn change here, nothing to look at.

## Self-check
Claim, what it rests on and its negation, the negated world, what integrates and what is carried (§self-check).

Decided by the author, 2026-03-02: merged: "yes, take it".
```

## Session record (`.records/<work-id>/<date>-<name>.md`)

Written at the end of each sitting, compaction or handoff, and before `.execution/<work-id>/` is closed; never edited after (§permanent-record). It never links into the work folder: the folder shares its work id.

```markdown
# <work id>: <short name> (<start date> to <end date>)

## Asked
The owner's requests in order, one line each, key ones quoted.

## Decided
One line each: the decision, "don't reopen" where the owner said so, and the design section it now lives in; one
decided on the recommendation keeps the choice, the alternatives and the reason (§authors-words-decide).

## Delivered
| Lane | What | Commits | Verified by (tests, checks, review findings) |
|---|---|---|---|

## Learned
What the work taught about the process or the principles, one line each, with what changed because of it (a rule,
a gate check, a pattern) or "proposed" with where the proposal is. The last line answers the fold-back question
(§reconcile-before-merge): `Fold-back: none` or `Fold-back: proposed (foldback.md)`.

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

if (/^\s*- \[ \](?! FINALIZE)/m.test(read("checklist.md"))) problems.push("checklist.md has unticked items: run them, cut them by the author's word, or move them to followups.md");
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
