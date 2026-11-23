# The novice pass: findings (stage 5, 2026-11-23)

The brief is [briefs/novice-pass.md](briefs/novice-pass.md). It ran on Sonnet in a visible Browser pane against
`packages-seq` at `678ce94fa`, Storybook on :6006.

**Coverage was partial.** The run ran out of budget before these: Hovering a stem, A pair, the scattered pick, the
made-up chr9–chr22 fusion, GenBank and BED, Dot plot and Liftover, a Hi-C cell, and the loading state. A second short
pass over them is a question for the author.

**Counted by kind**:

| Kind | Count |
|---|---|
| Bugs | 4 (3 confirmed by the orchestrator) |
| Missing design | 9 |
| Looks | 7 |

What the novice reported as working: finding a gene by name, the legend's "+n", dragging the square's window and
resizing it from an edge or a corner, and the translocation card's kind, strands and genes.

## Bugs → fixes.md (FX24)

1. **The ring's hover card shows a wrong position** (confirmed, screenshot). Hovering the t(1;17) ribbon gives
   "chr1 58,720,257", then "q32.2 chr1 207,100,001–211,300,000" twice. The band does not hold the position, and the band
   row is repeated. On a contact stem on chr14 the card read "chr14 58,720,257" with the same chr1 band. The position and
   band rows belong to the pointer's place on the ring and should name the block under it. They don't appear to.
   (Kinds of links; Hovering a stem.)
2. **The square's moves are not in the view's history** (confirmed). On The square beside, a click on chr4's cell moves
   the window, but Back (⌥←) stays `aria-disabled`. The novice saw the same after drags and resizes. A move made by the
   square should be one history step, like a click on a block's context (`navigator.ts`, "a jump by the block").
3. **A submenu doesn't open on click** (not confirmed). In the stem menu, a click on "Show with a partner" did nothing,
   while the Right arrow opened it. Check this by hand first: the browser tool's clicks may be the cause.
4. **The wheel over the ring is captured** (not confirmed). The page never scrolls while the pointer is over the ring.
   If the wheel is meant to zoom or turn the ring, this is a look, not a bug. Check what the wheel does there.

## Missing design → questions for the author

5. **A legend entry that hides its layer**. The kinds of links are hidden only from the track's Layers menu (the
   story's own text). The novice could not find that, and the Contacts entry has no swatch. Proposed: a click on a
   layer's legend entry toggles `layer.hidden`, as value categories already do. (Kinds of links)
6. **The square says nothing about its gestures and has no "whole genome" control**. Nothing on the page says the
   window drags or resizes. A double-click on an empty spot moved the window to chr12–chr16 instead of fitting the whole
   genome, and 0 (fit) did nothing. Proposed: a one-line hint under the square, and a double-click on the margin, or a
   button, that goes back to the whole genome.
7. **Synteny has no legend**. A bridge's colour and width are never stated, and inversion and translocation read alike
   on chr3. (Views › Synteny)
8. **The stem menu's partners**. "chr9 · 375 links", "chr9 · 65 links" and "chr9 · 24 links, reversed" are listed
   with no word on what separates them (bundles by strand, or kinds). chr22 is offered as its own partner. Ring-odd items
   are offered too: Move up, Reverse, Hide track.
9. **No key for the next feature**. In GFF3 › Find by name, Right moves the cursor a cell, and the keys panel lists no
   "next gene".
10. **No canonical transcript**. 40+ transcripts are labelled TP53-2xx and nothing marks which one is canonical.
11. **The VCF's 40 sample rows are unnamed**. Only hovering gives the name, as "State S11 / Value 0/1". A rows rail of
    names is the obvious answer (the alignments have one).
12. **Zoom in to read**. The hatched band says "zoom in to read", but the Links track is empty with no hint. The novice
    found no zoom gesture (+ and = did nothing, possibly the tool's keys). A zoom control, or the hint naming the
    gesture, is the question.
13. **No way to tell which layer is hidden**. This is a consequence of 5: the hidden-layer stories exist, but the page
    gives no route to them.

## Looks → packages/seq/TODO.md

14. Kinds of links in the dark: the faint contacts nearly vanish, so the legend's "faint in the chromosome hues" is not
    met (the dark layer's opacity).
15. The ring's labels: "chr1 … chr20" then bare "19", "21", "22". The labels are shortened where space is tight, and
    reading the two forms side by side is inconsistent. Either shorten all of them or none.
16. The expanded chromosome legend lists chrY while the ring draws 23 chromosomes with no chrY (the hues' table, not
    the blocks').
17. The card shows raw field names (`gene_a`, `q`, `Strands`). These are the file's columns and are honest, but a
    label map per field is a host's choice the ring could offer.
18. The keys panel scrolls past its last entry into blank grey, and its bottom rows are cut by the story frame.
19. The square dragged: at load the ring is clipped by the story's bottom, and a lone stub label "chr22 · 403" sits
    inside the ring unexplained.
20. Cards appear about 1.5 s after the pointer rests. That is the rest delay plus the tool's nudges, so check it by
    hand before acting.

## Pointer paths (B10)

| What | How the novice reached it |
|---|---|
| Keys panel | the keyboard button above a view (by sight); Escape closes it |
| Stem menu | pre-opened by the story. The right-click route was not tried, and neither was Shift+F10. "Show with a partner" opened with Right only |
| Hover cards | rest on a ribbon or stem |
| Square window | drag inside the window to move it, drag an edge or corner to resize it (found by guessing) |
