# sixth-slice: the scaffolds (2026-11-23)

## Asked

- "i think we should reconcile it instead and put it in some common place, where next similar "scaffolds" will
  live"; "let's do the scaffolds as 6th slice"; "oz was a temporary working name".
- "let's actually remove the meso\* mentions from packages/core ... this is something that exists at a different
  "level" and should not leak into the common libraries. and afterwards write the loop file"; then "Yes, start it".
- On compressing the principles: "the base repo ... did this and it saved quite a few tokens. let's do something
  similar here too? or maybe re-adapt for the principles we use in this repo?"; then "i think i will do one more change
  to the upstream oz first and then get back here with it", and "ok, <the upstream repository> is updated" (the upstream at
  `c1ccd14`).
- "do we need to update .scaffolds/_method too with the checklist.md, updated plan, etc templates; since this changed
  in the upstream version?".

## Decided

- **Stage 0** (the author): the name `.scaffolds/` and the execution folders kept in the tree, both as recommended.
  meso\*'s modules stay in Seq\*'s instance until a second scaffold needs them ("keep inside for now, once we have more
  scaffolds we generalize"). `.v6/` moves later ("Moves later"). A closed scaffold is kept, editable and never cited by
  code ("i think we can still edit if needed").
- **Stage 1's fifteen drift notes**: "All as recommended". The designs were amended to the code's way, and six code
  items went to the close-out.
- **The fold-back**: "ok apply all 5, and with questions go as recommended". Rule ids are checked per scaffold, and an
  instance's rule moves upstream when a second scaffold needs it. B3 and B4 went into the method, with a proposal
  drafted for upstream.
- **Check scripts**: "yes, let it edit the check scripts", after the classifier refused a lane's edit, which was
  surfaced to the author and not worked around.

## Delivered

| Lane | What | Commits | Verified by |
|---|---|---|---|
| (orchestrator) | meso\* and `.seq/` links out of `packages/core` | `351d9e9fc` | repr tests, `tsc` |
| S1A, S1B, S1C, S1D | 46 designs' decisions in the packages' specs, a spec for the examples package, 15 drift notes | `936242251`, `81a82d453`, `3001f67d2`, `e79e53f32` | `designs-in-specs` (no design pending); an independent review per lane, each finding real errors |
| S1F | the drift notes resolved as decided | `eecdabccd` | review; the checks |
| S1E | about 240 comments cite spec sections, not designs; `no-design-citations` | `5b09abdaa` | the check and its self-check; review (an EC id restored) |
| (orchestrator) | `oz` renamed, `.seq/` moved to `.scaffolds/` | `0be663880`, `88f55a78c`, `e0f554e49` | 471 pure renames; every check from its new place; 612 links in the changed files resolve |
| S4 | the method: templates, six generic checks, a link check, `.scaffolds/AGENTS.md` | `35ded94db` | the moved checks identical on Seq\*'s instance; review |
| (orchestrator), S4B | the package re-synced at upstream `c1ccd14`; the instance re-adapted | `af15684a8`, `26dcae39e` | the reading path about 34.3k to 26.5k tokens; review (three weakened rules restored) |
| S5M | the method on upstream's plan, checklist, state and findings files; the fold-back's five | `f9ad57937` | `ec-ids --compare` failing on EC-169's real loss; every self-check |

**The gate at the close**:
- every scaffold check and self-check from its new place;
- the quick lane, 4,949 tests;
- lint, format, `pnpm build` and `check:workspace`.

## Learned

The fold-back's five (B1–B5, "apply all 5"), now in the method and its checks:
1. A pointer by text must be unambiguous: 126 were found and repointed.
2. A mechanical rewrite keeps its identifiers, checked per file.
3. A rule is dropped as "the package's" only with the package's clause quoted.
4. An instance keeps how it applies a rule, and re-syncs rather than forks.
5. The orchestrator never `cd`s into a worktree.

Every lane's review found real errors (nine of nine) while its own checks were green: checks prove the shape, reviews
the content.

## Left open

- The close-out plan (`execution/close-out/`, a plan and its checklist), the next work: the use gate, the novice
  passes' 14 design questions, the minimal viewer's growth, the drag rule, the drift notes' code items.
- The upstream proposal (B3, B4; the sixth slice's `drafts/upstream-proposal.md`), for the author to take to the
  upstream repository.
- The viewer spec's two `**Keys**` leads under Input, which no pointer can tell apart (a spec fix).
- 102 dangling links in closed slices' drafts, reported by `links.ts` as history.
- `.v6/` into `.scaffolds/`, later; the meso\* modules, when a second scaffold needs them.
