# Lane TV: the text layout's viewer, and the rest of TX

Lane TX built the layout interface (merged), the text layout, its time targets, the p53 text figure and its stories,
and the wheel over a text section, on branch `seq/tx-text` (tip `5dd5548c2`). Its report left open what this lane
builds; its brief `tx-text.md` (beside this one) is still the goal, items 2 to 4, and [layouts](../../../design/layouts.md)
binds.

**Goal**:
1. Rebase `seq/tx-text` onto `packages-seq` at the base below (the interface is merged there, with TX's section-cache fix
   `dcefe7d76` applied as `9cfc6801f`; drop TX's interface commits as already merged) on your own branch
   `seq/tv-text-viewer`.
2. **The language version is 0.4.0 for the whole slice** (decided: one bump, unreleased): fold TX's 0.5.0 additions
   (`section.layout`, `group`, `numbers`, `group_labels`, `align`, `max_width`, `height`) into 0.4.0 everywhere (code,
   specs, examples, records, tests).
3. **The viewer for text** (TX's brief item 3): `ViewModel` through each section's `SectionNav` (TX's `TextNav`
   arithmetic is tested in Node; wire it: zoom, steps, reveal, clamp); `input` reading the placement (`at`, `offsetAt`,
   the section under a point); ↑↓ stepping lines; a drag selecting in reading order; a drag past a scrolling box's edge
   scrolling it, clock-timed; the hit target `window` with its `part` (labels hit as `item`, gutters as `row`); the
   `window` behavior provider (drags, slider semantics: `aria-valuetext` on settle, Esc cancels, arrows pan,
   Shift+arrows resize); the `reveal` group (`reveal.cursor` moved from `cursor.reveal`, `reveal.linked`);
   `selection.copyLetters` copying what is shown (reading order, the complement as drawn, shown positions only, rows as
   FASTA lines; an option for stored letters) with its key in `view.keys`; history steps restored by position across a
   layout change; `FrameInput.viewport` removed once the viewer goes through `SectionNav` (TX's drift).
4. **The rest of text** in TX's TODO that the design asks of this slice: the frames layer, `pieceBreak: 'line'`, numbers
   by a layer's numbering, the story `view` channel's fields for a text anchor; chips wait for RB's blocks (leave a
   seam, say so).
5. TX's drift to settle in your report with a recommendation: the interface additions past the design's signatures
   (`chrome(…, band)`, `fillBands`, `heights?`, `hitChrome?`, `followWindow?`, `Placement.clip?/window?/warpOf?`,
   `SectionNav` taking a `NavContext`): list them for the design.

**Model**: strong (opus). **Base**: `packages-seq` at `9cfc6801f`: `git checkout -B seq/tv-text-viewer 9cfc6801f`, then
bring TX's text commits over (`git cherry-pick 6fb54c720^..5dd5548c2`, skipping what is already merged), and check
`git merge-base --is-ancestor 9cfc6801f HEAD`.

**Coordination with RB** (blocks and rings, beside you, on the same base): keep your edits to text and the shared viewer
pieces this brief names; RB keeps to rings and blocks. Where the viewer's nav and input are shared, build the general
mechanism (`SectionNav` per layout, input through the placement) so RB's layouts plug into it; report what RB must
adopt. The orchestrator merges you first.

Everything else as in `tx-text.md`: Touches, Two consumers, Tests (the matrix cells as binding tests; back across a
layout patch; the nav in Node; twins on a fake clock), the performance representative (hover under 2 ms; a drag-select
step and a key step under 4 ms, `itWithin`), whose each choice is, minimal code, the Checks (with `pgrep -fl
"vitest.*heavy"` before a heavy run), Headers, fixes failing on their parent, Hot files, Stop when, Report.

The allocation tier: "p53 on its gene, whole protein" is over its 64 KB bound at the base (70.5 KB, a regression from
stage 1 the orchestrator is tracing); do not move that bound; report your numbers against the base's.
