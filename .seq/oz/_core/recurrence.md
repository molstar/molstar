# How work accumulates {#recurrence}

The frame the rest of `_core` sits in. Every piece of work here is one step of the same loop: **a session reads a
slice of `main`, does its work through that reading, and a merge folds the result back into `main`**. The next session
reads the `main` that merge left. The principles, the designs, the spec and the code are not described from outside
the loop; they are `main`, and they change only through it.

| In the loop | In a workspace built on these principles |
|---|---|
| `main` | the default branch: the code and its agent context (`.principles/`, `.design/`, `_spec/`, `.sources/`, `AGENTS.md`) |
| its record | git history and `.records/`: what each step asked, decided and delivered (§permanent-record) |
| a reading | what a session loads and works from: its skills, the designs and rules for its area, its plan and state in `.execution/`, centered on the author it works for (§shared-ownership) |
| a merge | a lane into its work branch, the work branch into `main`, each through the gate |
| a step | one merge; time here counts merges, not days |

## What follows

- **A merge can't be undone, only followed.** A later step can correct, summarize or supersede an earlier one, but
  not make it not have happened: records are never edited (§permanent-record), a removed rule leaves a tombstone
  ([system.md](system.md#ids)), and a design is edited in place with git keeping what it was (§design-in-place).
- **Every reading is one step behind.** A session works from the last merged state, so the spec trails the code,
  designs trail what prototypes learn, and agreement trails building. That lag is structural, not a defect: keep it
  visible rather than pretend it away (§spec-folders, §gaps-queue), and say which commit or record a reading is from.
- **Neither side moves alone.** The next reading depends on `main`, and the next `main` on the readings. A design
  fixed up front and structure left to emerge both fail; design first, prototype, and fold back what holds
  (§design-first, §iteration, §recursive-loop).
- **Every merge creates new work.** New structure makes new readings possible, so a healthy loop keeps producing
  gaps, open items and tensions it hasn't absorbed yet: residuals, the "Open" lists in patterns, "To decide" lists
  (§residuals, §surface-tensions). Their absence is a warning, not a goal.
- **Steps nest.** A lane's commits sit inside a session, a session's merge inside a pull request, a pull request
  inside the history of `main`. A prototype, a reference mock or a dry run is the same loop run on part of the
  context, outside `main`, whose result folds back only through a merge (§prototype-outside).
- **There is no view from outside.** No document stands beside the code to judge it from nowhere: the principles,
  the designs and the spec are part of `main`, read and rewritten by the same steps (§recursive-loop). A
  disagreement between them is drift to triage, not a verdict (§four-outcomes).
- **The principles are run on themselves.** `.principles/` is one more part of `main`, changed by the same steps,
  so every rule that governs a codebase governs the package too: it has a design, records, checks and drift like any
  other area. Applied to itself, the package should come out the same: that is its fixed point, and it is approached,
  never reached. Running the rules on the package is itself a step, and once its findings merge, `main` has moved by
  exactly what the run added: new rules, open items, a record of the run. That delta is not a defect to drive to zero;
  it is how the package keeps learning, and its size is what the check below measures. What comes closest to fixed is
  the checked structure, so a rule a tool can check becomes a check, here as anywhere (§self-check).
- **What doesn't fit stays visible.** Drift notes, the gaps queue and the open lists hold what a step could not
  absorb. The outside world enters only through readings: sources (§sources), signals from users and production,
  research for a design (§precedents).

## Checking the loop against itself {#self-check}

The check for the package (and for any claim the package makes about itself) is four steps, run on the evidence in
`main`, not described:

1. **Assume it holds.** Take the claim at face value and follow its consequences through the actual files: every rule
   it implies, applied to the package. Don't argue against it yet.
2. **Find what it rests on, and negate that.** Name the conditions the claim needs in order to hold, and turn the
   load-bearing one around. A real negation, not "on the other hand".
3. **Build the negated world on its own terms.** Take that negation as true and look at the same evidence through it:
   what does it explain that the claim doesn't?
4. **Hold both, and stop only here.** Name at least one thing that integrates (a fix, a new check, a rule restated)
   and at least one that doesn't (a tension that survives both readings). Merge the first; carry the second as a
   standing tension (§surface-tensions), visible until a later run dissolves it or shows it still holds. Don't smooth
   it away, and don't stop at step 3 because the two pictures look balanced.
   An integration can be a verdict: a claim that held against its negation, pinned by a check or by evidence a later
   run can repeat (the probe, the commands, what they showed), or a carried tension released with the reason it no
   longer holds. A bare "nothing changed" is
   not one: recording it is still a step, but it leaves nothing a later reading can rely on.

The run's delta (what it changed, what it carries) goes in the session record (§permanent-record). A run that
changes nothing and carries nothing is suspect: it is the warning that "every merge creates new work" names above. Run
it whenever the package changes, and on any tension that clarification alone doesn't dissolve.

