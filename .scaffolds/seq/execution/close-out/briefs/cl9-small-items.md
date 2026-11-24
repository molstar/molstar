# Lane CL9: small items version 1 shows

**Model**: Sonnet (a closed item list, each in a design already). **Base**: `packages-seq` at the orchestrator's head
(`git reset --hard packages-seq` first). Branch `seq/cl9`. Beside you: FX25 (`core/src/render/` ranges and genes,
headers' canvas draw, `region-hover.ts`, `ext-rondo/`, the region reads' hatch) and CL13 (`core/src/theme/provenance.ts`,
a folder-graph test, the lanes band's history, `rowsValue` on navigators, `seq.render.pairs.filter`). Never edit their
files; if an item needs one, stop and report it.

**Why**: version 1's rows decided by the author, 2026-11-23 (`drafts/version-1-rows.md`, "Small items", "take the rest
as recommended"); the version-aware validate carried from the sixth slice's drift (decided as recommended).

**Read first**: `.scaffolds/seq/execution/close-out/drafts/version-1-rows.md` (section 4); `.scaffolds/seq/design/`
`ux-foundations.md` (selection, the keyboard cursor), `format.md` (validation, versions), `extensibility.md`
(fallbacks), `examples.md` (the Antibody view); `.scaffolds/seq/sources/seqstar/edge-cases.md` for EC-149 and EC-154;
the core and viewer specs' sections you touch; `packages/seq/TODO.md`'s `repr` line (the version-aware validate).

**Items**, a commit at each:

1. **Selection snapped to codons (EC-149)**: a drag over nucleotides through a codon map snaps to whole codons where
   the design says. Proof: a test named for EC-149.
2. **The cursor across block gaps (EC-154)**: the keyboard cursor crosses a gap between blocks to the next block's
   first position (rows and lazy rows too). Proof: tests named for EC-154.
3. **The Antibody view's invalid labels shown**: the view's invalid numbering labels drawn and said as the validation
   design says. Proof: a test on the story.
4. **Built-in types' declared fallbacks, and the version-aware validate**: each built-in layer type declares its
   fallback; a spec of a newer minor version with an unknown kind gives a warning, not an error (an older reader
   draws what it knows, EC-235's way). Proof: tests per built-in type; a validate test both ways (newer minor: a
   warning; newer major: an error).
5. **The views' fixtures through MessagePack**: every example view's spec encodes to MessagePack and decodes equal.
   Proof: a sweep test (medium lane if slow, a quick representative kept).
6. **The features rows**: 178, 179/142, 37, 143, 39 moved in `features.md` and `features.mdx`. Proof: the rows.

**How**: tests name their EC ids; host neutral; reuse before new helpers; module size; nothing per frame (allocation
tiers green); stories with headless twins where something is drawn.

**Where**: a worktree on `seq/cl9`. Out of scope: FX25's and CL13's files. Headers per the root `AGENTS.md`. Hot files
never staged (the list in `fx25-fixes.md`'s Where). The language's grammar does not change (the validate's severity is
not a grammar change; if item 4 needs one, stop). Never a bare `git stash`. No network.

**Checks**: as `fx25-fixes.md`'s (unit, medium, perf tier, size rows, tsc, lint, format, the scaffold checks).

**Stop when** green, or a check red twice after a fix, or a refusal (report it; don't work around it).

**Report** (under 400 words): commits; each item with its proof; the checks with exit codes and counts; the size rows;
what to look at; drift, gaps and followups; the changelog text.
