# Brief: drafting a third-slice design

The third slice (`plan.md` beside this folder) builds figures for Seq\* (next). Before any lane runs, each design it
needs is drafted, reviewed independently, and agreed by the author; then the slice runs autonomously on them. You draft
one design.

## Read first

- `.seq/AGENTS.md`, `.seq/principles/README.md` (all of it: the rules the first two slices folded back, above all
  "The self-check runs on every design when it is agreed", "Interaction is a matrix", "Precedents cover interactions",
  "What is drawn is checked against the data", "A restriction in a `_spec` states its reason", "Gaps are asked at every
  gate; a default is confirmed by use").
- `.seq/oz/_core/workflow.md` §design-first, §precedents, §surface-tensions; `.seq/oz/_core/recurrence.md` "Checking
  the loop against itself" (the four steps).
- `.seq/execution/third-slice/plan.md` (what the slice must deliver, its lanes, its Gate 0 decisions).
- `.seq/design/README.md` and the designs your topic touches (named in your task); `.seq/design/stories.md` is a good
  model of a recent design's shape.
- The code as built: `packages/seq/core/src`, `packages/seq/viewer/src`, and their `_spec/README.md`.
- Precedents: the current Seq\* implementation in `packages/sequences/` (read it freely; never name it in anything under
  `packages/seq`), and the outside tools your task names. `.seq/sources/` has material on them.

## What to write

One file, `.seq/execution/third-slice/drafts/<topic>.md`, in the shape of the existing designs (`stories.md`):

1. **Intent**: what the design is for, with the author's words where the plan or the record quotes them.
2. **Precedents**: how this was solved before (the current Seq\*, the outside tools, this codebase), one line each on
   what is taken and what is not. **Interactions too**: what a user can do with this in the precedents, each listed as
   built here, planned (which lane), or refused (why).
3. **Decision**: the rules, each with its reason. Name what is configurable (config item, default, why, alternatives),
   following "configurable over fixed": a choice that could go more than one way is a config item, the recommendation
   its default.
4. **The interaction matrix's cells** this design fills (rows: a section, a follower, an element, a range, the
   selection, the view, a track or layer, a bookmark or link, a navigator; columns: select, extend, add/remove, focus,
   back, zoom, pan, hide, edit), each the gesture or key and the command, or "no, because".
5. **Configuration**: the table.
6. **Tensions**: what pulls apart, with options and a recommendation; open ones stay visible.
7. **Self-check** (the four steps, on the claim "this design covers what a user does with <topic>"): assume it holds;
   name what it rests on and negate it; build the negated world (what does the negation surface: an inverse, a
   symmetric case, a plural, a missing view?); hold both: what integrates goes into the Decision above (or is planned
   with a named lane), what does not stays as a standing tension in the design.
8. **Consequences**: per package (`seqx-core`, `seqx-viewer`, `seqx-examples`, the storybook): what is built, the
   tests (the headless twin, drawing against the data, property tests, a quick representative of each performance
   target), and the `_spec` sections that change.
9. **Gaps for the author**: the decisions you could not settle from the designs, the plan and the precedents, each a
   yes/no question with your recommendation. Keep it short: only what changes a picture, a gesture or a public API.

Rules: terse, a rule and its reason, no history; cite edge cases as `EC-<n>` from
`.seq/sources/seqstar/edge-cases.md`; designs may link each other relatively; never link into `.seq/execution/` from a
design (your draft lives there until agreed, then moves to `design/`). Do not edit any other file. Do not write code.

## Report

Under 300 words: the file, the decisions you made that a reviewer should check hardest, the self-check's delta (what
it added to the decision, what it left as a tension), and the gaps for the author.
