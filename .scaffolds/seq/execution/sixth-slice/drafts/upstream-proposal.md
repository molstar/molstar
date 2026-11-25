# Proposal for the principles package: an instance and its package (draft, for the author)

From Seq\*'s sixth slice, its fold-back B3 and B4 (`../foldback.md`). The author, 2026-11-23: they go into the method
here and are written up for the upstream package, not applied there. This file is the write-up; the author takes it to
the upstream repository. Nothing here edits `.scaffolds/principles/`.

Base: the package at `c1ccd14`. Both changes go into `_core/patterns/adopting.md`, a pattern (no rule id needed, so no
`§rule-budget` count); each is two bullets, within the budget a rule would have. The ids they cite all exist at the base.

## 1. Step 1, "Write the adaptations": what an instance keeps (B4)

Replace the bullet "Write the adaptations, one line each with why: …" with:

```md
- **Write the adaptations**, one line each with why: modules dropped; rules where the host wins (its layout, its test
  tiers, its package rules); where open items go if not to issues (§followups). The manifest keeps only what is the
  project's own and how it applies a package rule (a number, a path, a tool, a tightening), each within §rule-budget;
  the evidence behind a line goes into the session record, not the manifest. Why: a manifest that restates the
  package and carries its evidence grows past what a session reads (one instance's fell from 32k to 9k characters,
  its reading path from 34.3k to 26.5k tokens, once it kept only its own).
```

## 2. A new step 8, "Update from the package": re-sync, never fork (B3, B4)

After "## 7. Fold back":

```md
## 8. Update from the package

- **Re-sync, never fork.** An update replaces `.principles/` with the package at a newer commit, a merge from the base
  the manifest's provenance names, which it then records. A change the project needs in the package is proposed to
  the package's home (step 7), never made in the copy. Why: a forked copy stops receiving the fold-backs of other
  projects, and its next update is a conflict instead of a merge.
- **A dropped adaptation names what covers it.** Re-adapting the instance after an update, each adaptation dropped
  because the package now states it names the package's id and quotes the clause that covers it, in the update's
  record. Where the package states it more weakly ("may", "optional", none by default), the instance keeps one line
  that tightens it and names the id. Why: an instance dropped three rules as "the package's" where the package said
  less (a fold-back at every close, a standing novice pass, a mock's deletion), and the review, not a check, found it.
```

## Evidence

- B3: the sixth slice's lane S4B (re-adapting the instance after the re-sync to `c1ccd14`) dropped three rules as
  covered; its review found the package's wording weaker each time (§reconcile-before-merge allows "Fold-back: none";
  §used-before-done's novice pass "may"; §interactive-mocks names no deletion). The instance restored each as one
  tightening line.
- B4: after the same re-adaptation the instance's reading path went from about 34.3k to 26.5k tokens (its manifest 32k
  to 9k characters, its `status.md` 15k to 7k); the re-sync alone added about 1k.

## Seen while writing this, not proposed

- Step 7's last bullet, "Reconcile and remove the work folder before merge", contradicts §permanent-record and
  §reconcile-before-merge at the same base (the folder is committed and kept, and closes by `state.md`'s last line).
  Worth a one-line fix upstream: "Reconcile, then close the work folder (§reconcile-before-merge)".
- Step 1's `adopt-principles` folder names a plan, a checklist and a `state.md`, while "Done when" asks for "a plan, a
  state file and a gaps file": the checklist is missing from "Done when".
