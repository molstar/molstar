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

## Second short pass (2026-11-23, on FX24's Storybook)

The author asked for it ("Ok run the other short novice pass"). It covered every item the first pass missed, at least
briefly. Not reached: the chr9–chr22 fusion's card, a double-click on a ribbon, a folded stub (none seen), the matrix
story's Hi-C cell, BED's cytobands, and dark on the stem stories.

**Counted by kind**:

| Kind | Count |
|---|---|
| Bugs | 1 (to check) |
| Missing design | 5 |
| Looks | 8 |
| Unverified | 2 |

What the novice reported as working:
- the square's Back from a scattered pick, in one step;
- GenBank's card;
- BED12's whole-gene highlight;
- a pair opened from the stem menu's submenu by a click.

### Bugs → fixes.md

27. **Not a bug** (the orchestrator, on the Canvas copy: the square reads "chr9:1 – chr22:50,818,468" with the pair's cells outlined). **After a stem menu opens a pair, the square isn't checked.** The novice saw the whole genome on the square after
    the ring opened chr22–chr9. Its screenshot shows the SVG copy, which has no square of its own, so check this on the
    Canvas copy. (The stem menu)

### Missing design → questions for the author

21. **The stem card's "Links" row reads against the rows under it**: "Links 24 links to 1 chromosome" is the hovered
    stem's own links (to one partner), while "Contacts 638" and "Translocations 24" are the chromosome's totals per
    layer (confirmed). The orchestrator's proposal is to name the first row by what it is ("This stem 24 links to
    chr9"), with the totals under "chr22, all links".
25. **A pair's ribbon card names one end only** (chr4 and its band); chr7's position and band are missing. The bins row
    "zoom in to see them" is the only hint.
26. **The partner submenu repeats chromosomes** (chr9 three times, chr22 twice) without saying which stem piece each
    row is (by kind, by strand). This is first-pass 8, seen again.
28. **A right-click on a stem opens nothing**. The menu exists only pre-opened in its story, so there is no route a
    user finds.
33. **A Hi-C card gives bin indices** ("chr1:177 × chr1:173, x 177, y 173"), with no kb or Mb on the card or the ruler.
35. **The links-pending state says nothing**: the ring's hole is empty, the square is an empty grid, and there is no
    hatch and no "loading" text.

### Looks → packages/seq/TODO.md

22. The translocation card: raw keys (`gene_a`, `event`, `q`), "Strands – –" unclear, and the chr12 and band rows
    repeated (check after FX24).
24. The Contacts legend entry has no swatch, and its label doubles ("Contacts: by chromosome").
30. The scattered pick in dark: bright blue grid lines over the cells, and a stray "chr9 · 470" stub label inside the
    ring.
31. The dot plot in dark: empty cells glaring white, and the axes named "chr1" against "1" (assembly B's names), each
    ruler labelled "ruler".
32. Liftover: overlapping gene labels, jargon in the description (A′, `seq`), no legend for bridges, and no card on a
    lifted gene.
34. The methylation matrix in dark: a missing cell is as dark as the ramp's low end; the samples are unnamed and called
    "State"; the chromosome is called "chr".
36. BED12: no names on the tracks, and the strand "–" is easy to misread.

### Unverified

23. The fusion's card was not found on the Canvas copy, and the SVG copy gives no hover card (it is a picture, by
    design).
29. Shift+F10 after a click opens nothing. The novice's click may not have been on a stem.

**The totals over both passes**:

| Kind | Count |
|---|---|
| Bugs | 5 (2 fixed by FX24, 3 not bugs) |
| Missing design | 14 |
| Looks | 15 |
