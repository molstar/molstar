# Fold-back: what the sixth slice taught (draft, for the author)

The sixth slice moved Seq\*'s scaffold into a common place for scaffolds:
- stage 1: the designs' decisions reconciled into the packages' `_spec`s (S1A–S1D, the drift notes by S1F);
- stage 1b: the citations of designs repointed to spec sections (S1E);
- stages 2 and 3: `oz` renamed to the principles package, and `.seq/` moved to `.scaffolds/`;
- stage 4: the method written down (S4);
- stage 4b: the package re-synced from upstream, and the instance re-adapted under its rule budget (S4B);
- stage 5: verified.

It ran alone, with no code change but comments and one comment's rule id.

**Nothing here is applied before the author's line at the end.** Each proposal is one line: "Proposed: … Why: … Where: …".

Where to:

| Code | Where |
|---|---|
| **P** | the principles (the instance's adaptations, or upstream) |
| **M** | the method (`.scaffolds/_method/`) |
| **C** | a check |

## A. What each lane and review found

| Lane | Reviewer's findings, fixed in one round | Drift |
|---|---|---|
| S1A (core, repr) | drift notes against existing TODO entries, ambiguous pointers, a genomics example in core's spec, a wrong path | 5 |
| S1B (viewer) | a wrong rule written into the spec (the clipboard), a mispointed row, missing rows | 2 |
| S1C (data, drawing) | a missed rule (the hover outline's 3:1), five pointers landing on the wrong section, missing reasons | 5 |
| S1D (ring, catalogs) | a row citing a section that says nothing of it, draw order missing, a lazy chunk said absent, unmarked drift | 2 |
| S1E (citations) | an EC id dropped, the check's exemption too wide | — |
| S1F (drift resolved) | stale text left beside each amendment | — |
| S4 (method) | a template contradicting the package's `_core`, fold-back rules missing, a claim no check enforces | — |
| S4B (re-adaptation) | three rules weakened by being dropped as "the package's" when the package says less | — |
| The orchestrator's stages 2, 3, 5 | six stale lines, one link the move broke | — |

Read: **every lane's review found real errors**, eight of eight, though every lane's own checks were green. The checks
proved the shape (a row exists, a section exists, a link resolves); only the reviews read whether the content was
true.

## B. What this slice taught

1. **A pointer by text must be unambiguous.** (C)
   - Proposed: `designs-in-specs` fails a pointer whose text starts more than one heading or bold lead in its spec, and
     names the candidates.
   - Why: the check took the first match, and reviewers found about ten pointers landing on the wrong section ("Links",
     "Bounds", "Configuration", "The ranges layer") while the check passed.
   - Where: C (`_method/_checks/designs-in-specs.ts`).

2. **A mechanical rewrite keeps its identifiers.** (M; C)
   - Proposed: a lane that rewrites text across many files reports the set of identifiers it touched (`EC-<n>`, rule
     ids, spec sections) before and after, and the review compares the two sets.
   - Why: S1E dropped EC-169 from one comment while rewriting 240 citations. `ec-ids` passed, because the id was still
     cited elsewhere.
   - Where: M (`brief.md`, the review list), C (a `--compare <ref>` mode in `ec-ids`).

3. **Dropping a rule "because the package states it" needs the package's wording checked.** (M)
   - Proposed: re-adapting an instance to its package, each dropped rule names the package id and quotes the clause
     that covers it; a rule the package states more weakly stays as one line that tightens it.
   - Why: three of S4B's drops were weaker upstream ("may", "optional", "none"): the fold-back at every close, the
     standing novice pass, the mock's deletion.
   - Where: M (`README.md`, how an instance re-syncs), P (upstream's `adopting.md`, if upstream wants it).

4. **The instance keeps how it applies a rule, not the rule.** (P)
   - Proposed: an instance's manifest holds only what is its own or how it applies a package rule (a number, a path),
     each within the rule budget; the evidence lives in the records; an update re-syncs from upstream instead of
     forking.
   - Why: after the re-sync, the instance's re-adaptation cut the reading path from about 34.3k to 26.5k tokens (the
     manifest 32k to 9k characters, `status.md` 15k to 7k). The re-sync alone added about 1k tokens.
   - Where: P (the instance's Adaptations; upstream's `adopting.md`).

5. **The orchestrator never changes directory into a worktree.** (M)
   - Proposed: the orchestrator runs every command with absolute paths, or `git -C`, `pnpm -C`; never a bare `cd` into
     a lane's worktree.
   - Why: twice this slice a `cd` into a lane's worktree moved the session's working directory there, and the next
     commands ran against the wrong tree until it was noticed.
   - Where: M (`loop.md`, Each step).

## C. The earlier fold-backs' rules, checked against this slice

- **Reviews in every lane** (fifth B1–B2, earlier): *held, and decisive* (A above).
- **Items landed, not reported** (fifth B3): *held*: each lane's designs were ticked against the check.
- **What the orchestrator owns** (fifth B4): *held*; no lane wrote outside its worktree. One lane needed
  `pnpm install --offline` in its worktree: allowed, no network.
- **Classify before logging** (fifth B2): *held*: the PENDING conflicts and the stale-text findings were handled as one
  rule each.
- **Mocks, the novice pass, sizes, timing**: *no evidence*: nothing was drawn or timed.

## D. Not proposed

- How rule ids scope once a second scaffold exists: a question, below.
- The 102 dangling links in closed slices' drafts: history, reported by `links.ts`, not fixed.

## Self-check, on the claim "these proposals make the method truer"

1. **Assume it holds.** B1 and B2 turn two misses that a green check hid into checks. B3 and B4 keep an instance from
   drifting from its package. B5 removes a hazard of the orchestrator's own.
2. **What it rests on, negated.** It rests on the misses being mechanical. Negated: the reviews found mostly content
   errors (a rule written wrong, a drop that weakened a rule), which no check reads.
3. **Both held.** The checks catch the mechanical kind cheaply. The reviews stay mandatory for the content: this
   slice's evidence says every lane needs one, even a lane that only moves text.

## Questions for the author

1. **Rule ids across scaffolds** (S4's question): Seq\*'s own rule ids are defined in its manifest. When a second
   scaffold exists, should `principles-ids` run per scaffold (each instance's ids its own), or should rules that prove
   generic move upstream? Recommended: per scaffold, and an instance's rule moves upstream when a second scaffold
   needs it (the same rule as the meso\* modules').
2. **Upstream**: send B3 and B4 to the upstream package's `adopting.md`, or keep them in the method?

## Decided by the author

Decided by the author, 2026-11-23: merged: "ok apply all 5, and with questions go as recommended". Rule ids are checked
per scaffold, and an instance's rule moves upstream when a second scaffold needs it. B3 and B4 go into the method here,
and are written up as a proposal for the upstream package, not applied there (another repository; the author's).
Applied by lane S5M with the upstream workflow files (`c1ccd14`, §plan-and-checklist, §state-file): the author, the
same day, "do we need to update .scaffolds/_method too with the checklist.md, updated plan, etc templates; since this
changed in the upstream version?".
