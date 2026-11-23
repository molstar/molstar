# Drift triage {#drift-triage}

Rules §recursive-loop: a disagreement between layers is a signal; triage it into one of four outcomes, never silently bend either side; agents fix code and stale specs; design and architecture changes are the call of the author the agent works for (§upward-authority).

## Drift note (a working note, deleted once resolved)

A drift note lives in `.execution/<work-id>/drift/<slug>.md` ([§workflow-files](workflow-files.md)). A note that exists is open; resolving it means applying the outcome and deleting the note. There is no status field.

```md
---
level: code-vs-design     # code-vs-design | design-vs-architecture | architecture-vs-principles | spec-vs-code
outcome: design-wrong     # proposed: code-wrong | spec-stale | design-wrong | undecided
refs: [.design/storage/document-store.md, lib/storage/src/memory.ts]
---
## What disagrees
The design says `watch()` emits only changes; every caller needs the current value first.

## Evidence
- 4 call sites call `get()` and then `watch()` (files...); the gap between the two caused issue #212.

## Proposal
`watch()` delivers the current value first. Update the 3 adapters and documentStoreSuite.

## Rejected alternatives
A parallel `watchWithCurrent()`: two ways to do one thing.
```

## Outcomes

| Outcome | When | Who may apply it |
|---|---|---|
| code-wrong | the design states intent we still hold; the code does not follow it | agent: fix the code |
| spec-stale | the code is right; `_spec/` lagged (level spec-vs-code) | agent: update the entry in the same change |
| design-wrong | following the design produces bad code; bending signals repeat (§bending-signals) | the author's decision, then edit `.design/` in place, then the code (§upward-authority) |
| undecided | the evidence does not settle it | nobody: stop building on the point; a followup issue if it outlives the work |

Above `code-vs-design`, every outcome edits an agreed document, so it is the author's call (§upward-authority). A drift is resolved at the lowest level that can absorb it (§recursive-scope).

## Finalize check

```ts
// finalize (excerpt; the full sketch is the finalize check in §workflow-files): any remaining note is unresolved
const drift = `${dir}/drift`;                                         // dir = the work folder it was given
for (const f of existsSync(drift) ? readdirSync(drift).filter((f) => f.endsWith(".md")) : [])
  problems.push(`drift/${f} is open: apply its outcome (fix code, update spec, edit the design) or file a followup issue, then delete it`);
```

Notes are committed with the work folder while they are open (§handoff), and the finalize check (§workflow-files) lists each one still there: at the close none is, so the archive holds no open note and git holds every resolved one.

## Design is wrong: edit the design in place

Edit the design directly once the human agrees. The evidence and migration go in the commit message; keep a line in the design only if future readers need the reason (§design-in-place):

```md
<!-- .design/storage/document-store.md (edited in place) -->
## Decision
`watch()` delivers the current value first, then changes.

## Rejected alternatives
- A parallel `watchWithCurrent()`: two ways to do one thing; every caller wanted the current value.
```

The spec entry then describes the new behavior; git history holds the rest. There is no revision list or separate log.
