# The use gate: what to try (prepared by the orchestrator, 2026-11-24, at `6a533d536`)

The author's gate (plan, Carried; the fifth fold-back's B10): each view's story against its agreed picture, the
design's **Defaults, with a picture** and its **interaction matrix** (Storybook `Seq* (next)/Design/<page>`). Ticked by
the author only. Storybook from the main checkout (`storybook`, port 6006). A finding goes to `fixes.md` (a look or a
bug), `TODO.md` or a drift note (a design question); one that blocks use jumps the queue.

## 1. A genome in a row (Design › Genome Coordinates)

Stories: Features › Genome Coordinates › One Assembly In A Row, B In As Order, A Chromosome At A Time, Three Thousand
Contigs, A Bed Over The Assembly, Cytobands Over The Genome; Examples › Genomes › Views › Transcripts, Signal, Samples,
Methylation, Manhattan.

Try: click a chromosome's label (the whole block), Shift+click another (every block between), double-click (fit to it,
again the whole); 0 (the whole genome); a drag across a block's edge (a range a block a part); Shift+arrows across an
edge (announced); Alt+PageDown / Alt+PageUp (the next block); the navigator: click a chromosome (a jump); Alt+Left
(back). Ctrl+wheel and drag in Transcripts and Signal; the find box in Transcripts (a gene by name).

## 2. The ring (Design › Genome Ring)

Stories: Features › Genome Ring › The Genome At The Fit, One Chromosome, A Pair, Two Assemblies On One Ring, Crowded:
3,000 Contigs, Kinds Of Links (and its dark), Links Pending, Links Failed, No Links, Without The Extension; Features ›
Genome Coordinates › A Ring That Switches.

Try: Cmd/Ctrl+wheel or a pinch over a chromosome (the cut opens; one undo step; at the fit, option (c)); click a ribbon
(both ends lit), a stem (its range), a stub; double-click a ribbon (focus); right-click a stem (the menu: "Select the
stem and its partners", "Show with a partner ▸"); Shift+F10 at the cursor's stem; `.` and `,` (the next partner,
announced); a label's menu "Show only"; Alt+Left and undo (the blocks shown come back); the legend's eyes (contacts,
translocations hidden).

## 3. The square (Design › Genome Ring, The square navigator)

Stories: Features › Genome Ring › The Square Beside, The Square Dragged, The Square Resized, The Square Clicked, The
Square Scattered Dark; Features › Genome Coordinates › The Square.

Try: click a cell (that pair shown at once), its title or Home (all); drag the window's edge or corner over the ring (a
run of blocks grows or shrinks); undo.

## 4. Across sections (Design › Across Sections)

Stories: Features › Across sections › Two Genomes At The Fit, Zoomed Out Bins, One Side Zoomed In (and the other way), A
Stub, Pending A Side At A Time, Three Genomes Chained, Two Genomes As A Dot Plot, Liftover, Contacts Between Two
Chromosomes; Examples › Genomes › Views › Synteny, Dot Plot, Liftover.

Try: drag and wheel over the bridge's band (both sides), over one side (alone; with `lock: px` both); Shift+A at a side
(align); click a ribbon (both ends), double-click (each side to its hull); ↓ / ↑ across a bridge (the image of the
position, announced); the menu's "Select across"; canvas against SVG after each move (FX28).

## 5. Region reads and formats (Design › Region Reads, Genome Formats)

Stories: Features › Region reads › Zoom in, Pairs, Hi-C (each); Features › Genome formats › BED, GFF3, GTF, VCF,
GenBank (Tp53, Find By Name, Misread File); Examples › Genomes › Views › Contact Map.

Try: zoom in to read (the hatch, "zoom in to read", then the items); hover a bin (its count) and an item (its card, the
variant's title); a misread file's message; find a name.

## 6. What landed in the close-out (each with its look done; worth a try in use)

Search (CL5: Cmd/Ctrl+F, next and previous), track groups (CL6: Features › Track groups, fold and open), reference
lines (CL8: Manhattan), the legend row (FX29: Features › Structure link › Focus And Back), `height: 'content'` (FX30:
Features › Tooltips › By Data).
