# Lane CL13: the sixth slice's carried code items

**Model**: Sonnet (a closed list, each decided). **Base**: `packages-seq` at the orchestrator's head (`git reset --hard
packages-seq` first). Branch `seq/cl13`. Beside you: FX25 (`core/src/render/` ranges and genes, headers' canvas draw,
`region-hover.ts`, `ext-rondo/`, the region reads' hatch) and CL9 (the viewer's selection, input, navigator cursor and
keys, `core/src/lang/` validation, the Antibody example). Never edit their files; if an item needs one, stop and report.

**Why**: the sixth slice's drift notes, decided by the author as recommended (2026-11-23), their code items carried to
this plan (`plan.md`, "The sixth slice's drift").

**Read first**: `.scaffolds/seq/execution/close-out/plan.md` (that section's six items; the version-aware validate is
CL9's); `.scaffolds/seq/execution/sixth-slice/briefs/s1f-drift.md` and `../sixth-slice/state.md` for each item's note;
`packages/seq/TODO.md` (`bandTop`, Navigators' `rowsValue`); the specs you touch.

**Items**, a commit at each:

1. **The theme picker's BLOSUM62 text** (`core/src/theme/provenance.ts`): says what the scheme does: it compares each
   letter with the reference row. Proof: a test on the text.
2. **A folder-graph test in `packages/seq/core`**: each folder an entry, no import cycles between folders. Proof: the
   test, and it failing on a planted cycle.
3. **A lanes band's history step keyed by item, not `bandTop`**: back and forward restore the band's scroll by its top
   lane's item key, so a re-laid band returns to the same item. Proof: a test where lanes change between steps.
4. **`rowsValue` drawn by a navigator over a lazy alignment**: the rows' one-value-a-row fallback where a source serves
   no row summary. Proof: a test on a lazy alignment's navigator.
5. **The info issue `seq.render.pairs.filter`**: said when a data filter is not applied to bins. Proof: a test.
6. **The TODO's lines** for the above named in your report (you don't edit the TODO).

**How**: host neutral; reuse before new helpers; module size; nothing per frame (allocation tiers green); stories with
twins where something is drawn.

**Where**: a worktree on `seq/cl13`. Out of scope: FX25's and CL9's files. Headers per the root `AGENTS.md`. Hot files
never staged (the list in `fx25-fixes.md`'s Where). Never a bare `git stash`. No network.

**Checks**: as `fx25-fixes.md`'s.

**Stop when** green, or a check red twice after a fix, or a refusal (report it; don't work around it).

**Report** (under 350 words): commits; each item with its proof; the checks with exit codes and counts; the size rows;
what to look at; drift, gaps and followups; the changelog text.
