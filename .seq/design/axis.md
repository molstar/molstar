# The axis: direction, letters through a path, collapsed introns

## Intent

A section's axis is one system laid out across the page. Three things are missing, all found by the author using the p53
views (2026-11-18):

- **Direction.** "when scrolling the protein, the gene slides the other way, is that expected?": a minus-strand gene
  following its protein is laid out in genomic order, so it moves against the protein ([maps.md](maps.md): the hull in
  genomic order).
- **One axis.** "having the both the protein and the gene on a same axis under each other": the protein's letters on the
  gene's axis, a residue a codon wide; refused today (`seq.render.layer.letters`, a restriction that had no reason).
- **Collapsed introns.** A gene following its protein drops its letters while an intron is in the window: the window
  spans the intron at full width (20 residues: 60 coding nucleotides at about 10.7 px, 140 with the intron at about
  4.6 px, below the letters' density step).

This design finishes the axis: a section's direction, another system's letters on it through a path, pieces that hide
introns, and what reads the axis (hit-testing, selection, follow both ways, focus, the view history, the export).

## Precedents

What is taken and what is not:

- **The current Seq\***: a block shows a piece of its system in a `direction` (1 or −1), a layout map, never new data
  (EC-31, EC-32); complementing is a separate `letter_map`, and reversing never complements by itself (EC-233); exons
  only as a position filter over an exon annotation (`exonsOnly`), the intron a break and a codon across it one residue;
  the protein's residues drawn over their codons, a split one in both parts, three-letter codes only with room in the
  narrowest part (EC-117). Taken: direction as a layout map, pieces over one system, the split codon, and (the fifth
  slice, [genome-coordinates](genome-coordinates.md)) a per-block direction. Not taken: a direction per piece inside a
  block (Tensions), complement off by default (the author's call).
- **IGV**: the sequence track's strand arrow shows the complement without reversing the axis; three-frame translation.
  Taken: complement as its own setting. Not taken: translation frames (a layer type, not the axis).
- **UCSC**: the "reverse" button lays the window out from the highest coordinate down, the ruler counting down;
  "complement bases" apart; multi-region mode shows exons only, with padding around each. Taken: all three (padding is
  `pad`).
- **Ensembl's transcript views**: exons with the intron's length between them, flanking intron bases on request. Taken:
  the hidden length on a break, `pad`.
- **SnapGene**: the translation under the DNA, one- or three-letter, a codon split by an intron drawn in each part; a
  reverse-strand feature's translation reads right to left; the whole map flipped to the reverse complement; a click on
  a residue selects its codon. Taken: the translation, the split codon, the flip, the click. Not taken: a
  reverse-strand translation read right to left on a forward axis; here a writer sets the axis `along` the protein.
- **JBrowse 2**: a linear view flipped horizontally; a sequence track with both strands and six frames. Taken: the flip.
  Not taken: six frames (a layer type, later).
- **This codebase**: `Follow.point` and `followRange` ([viewer-models.md](viewer-models.md)), followers that lead, offsets
  in display order with positions saved (the viewport rule), `projectBoundary` (EC-28, EC-29), closure (EC-30, EC-144).

**Interactions in the precedents**, each built here, planned or refused:

| What a user does there | Here |
|---|---|
| flip a view to the other strand (UCSC reverse, JBrowse flip, SnapGene) | built: `section.reverse` |
| show the complement without flipping (IGV, UCSC) | built: the `complement` param, by patch; no gesture of its own (no, because flipping covers reading the other strand) |
| show exons only, with padding (UCSC multi-region) | built: `pieces`, `pad`; `section.collapse` toggles |
| read an intron's length (Ensembl) | built: a break's hover carries the hidden range |
| expand one intron of a collapsed view | refused for now: the pieces are the spec's, a per-break override would make the export differ from the screen (Tensions) |
| click a residue to select its codon (SnapGene) | built: the hit's granularity is the residue, projected onto its codon |
| translation under the DNA, split codons in each part (SnapGene, the current Seq\*) | built: letters through a codon path |
| three or six translation frames (IGV, JBrowse) | refused here: a translation layer type, later ([features.md](features.md)) |
| the codon's letters in the hover label, split by `\|` (the current Seq\*) | planned: with the hover card (viewer-models, later behaviors) |
| the protein following the gene, waiting at the nearest coding residue in an intron (the current Seq\*) | built (as today) |

## Decision

### The author's decisions

The author, 2026-11-18, on the direction: "i think you are specializing on implementation details here instead of
letting the application drive the compositions. / rest as recommended"; and in general: "a lot of this is "configurable
behavior" and where appropriate, the code should allow the configuration via spec ports etc.; we already have a lot of
patterns for this, so use your own judgement what's appropriate to keep just as a single feature and what to make
configurable, with port, etc."; on sizes: "don't overdo it on the limits (they were sort of self-imposed and estimate),
but obviously try to minimize the amount of code". On the review of the decisions: "ok" to its recommendations.

1. The direction is the writer's: `genomic`, `reversed`, `follow`, `{ along }`, defaults `follow` for a follower and
   `genomic` otherwise, never inferred.
2. `pieces` defaults to `whole` for every section, followers included; the p53 views spell `coding`.
   `section.collapse` toggles to the spec param `collapsed` (default `coding`), with no memory of a form written before.
3. A collapsed intron is zero width on every layout, rings included, a hairline marking it.
4. A drag across a collapsed intron selects only the positions shown (`selection.acrossBreaks: shown`; `span` the
   setting).
5. Expanding a single intron is refused for now; Shift+C toggles the whole section.
6. A history step over a changed axis is restored by position, clamped to what is drawn, never silently dropped.

### Whose each choice is

Each choice below says whose it is (the same rule in every design of the axis, themes, layouts, grids and rows,
navigators):

- **a spec param** (the writer's: what the figure is), with a plain default the library never infers from data, from
  other nodes, from a column's presence or from an example;
- **a `RenderConfig` preset** (spacing, sizes, marks: how any figure is drawn);
- **a `SpecBounds` bound** (what a spec may size);
- **a `ViewerConfig` item** (the host's or the reader's: how this viewer behaves; its keys in `view.keys`);
- **a port**, **a registry entry**, or **an option of a call** (`renderSvg`'s, a builder's);
- **fixed**, with its reason.

A config item needs a second real use; else it is fixed. No new port: the axis is pure core over the resolved view; a
host drives it through the spec (patches, the two commands) and the config. Why: nothing here reaches outside the
library (no clock, no data, no measurement).

So, by the author's word on sizes, a size estimate cuts no feature here, and the code stays minimal:
the axis extends what exists (`AxisTransform`, `projectPoint` exported and the viewer's `pathPoint` removed, the
codon map's runs, `projectBoundary`, closure, the hit's `extent` and `hull` as ranges use them), with no private copy and
nothing without a use; the fixed behaviors below are fixed partly so.

### A section's direction

- **`direction`** on a section (spec): `genomic` (positions count up left to right), `reversed` (they count down),
  `follow` (the leader's direction times the follow path's: reads along its leader), or `{ along: <system or map> }`
  (reads along that system: the direction of the path from it to the section's system, so a minus-strand gene `along`
  its protein reads N to C left to right). Default: `follow` for a follower, `genomic` otherwise. The library resolves
  what the writer asked; it never infers a direction from the section's layers. Why: the composition is the
  application's (the author); a default read off the layers would flip the section when a layer is added.
- **A direction that cannot be resolved** (`follow` with no leader placed, `along` a system no path reaches) is
  `seq.view.direction.unresolved` (warning) and the section is `genomic`.
- **A path with runs in both directions** (an inversion in a `segments` map) resolves `follow` and `along` to the
  direction of the runs covering most of the path's length, with `seq.view.direction.mixed` (info). Fixed. Why: a
  direction that changed while panning would turn the axis over under the hand; no figure has asked for another rule
  (today's `pathDir`, the first run's, is replaced).
- **One direction per axis.** `Axis` gains `dir`. Pieces stay `[start, end, …]`, start below end, in display order (on a
  reversed axis the last piece first); a reversed piece shows position `end − 1 − k` at its offset `k` (EC-32). A
  direction per block (`block.direction`, an inversion drawn as its own block) is
  [genome-coordinates](genome-coordinates.md)'s; a direction per piece inside a block stays a tension, not built.
- **A reversed axis says so**: the ruler counts down (the ruler layer draws it from the axis, [themes.md](themes.md)), and the section's first label carries a
  direction mark (an arrow and "reversed"; with `along`, "along p53") drawn in the `labels` pass, so the export says it
  too. `axis.directionMark` (`RenderConfig`: `label`, or `none` for a host whose own chrome says it, as a Mol\* plugin's
  panel would).
- **A ruler's insertion labels** (decided by the author, 2026-11-20, after `112.1` never showed at the closest zoom): an
  insertion is labelled beside its residue's number wherever it clears what is drawn, a label wider than its position
  too: the round values first, then the insertions nearest a drawn major, the shorter, by position. One rule: every
  label is placed if it clears what is drawn, the order decides who wins; decided at layout, never on a hover, so
  nothing jitters, and the hit card names the rest. `ruler.labelGap` is the clearance.
- **A ruler's units** (decided, the author, 2026-11-19: "Yes, `units: none | bases`"; found by CS at 83 Mb): the ruler
  layer's `units` param, `none` (the number as its numbering labels it, the default) or `bases` (bp, kb, Mb, chosen by
  the tick step: 7,668,402 at a base tick, 7.67 Mb at a megabase one), formatted by `@molstar/core/util/number`'s
  genomic formatting ([base-library.md](base-library.md)). Why: a writer of a chromosome view says so; never inferred
  from the alphabet. A label in another numbering than ordinal ignores it (its labels are the source's).
- **Where a numbering does not reach** (decided, the author, 2026-11-19, the Antibody's light chain: Kabat and IMGT
  number only its variable domain, the rulers blank over the constant one with no word why: "customizable … 'not
  numbered' as default"): the ruler layer's `unnumbered` param. `mark` (the default): past the numbering's coverage the
  ruler's line is muted and dashed, one muted note at the start of each unnumbered run ("not numbered", with the
  numbering's name where the track header does not already say it), no ticks; the hover card and the live region there
  say the sequence position and "not numbered by Kabat". `fallback` with `fallbackNumbering` (a numbering's name;
  ordinal by default): the run is labelled by that numbering in the muted style, so the reader tells the two apart.
  Why: a blank is read as a bug; a fallback numbering is the writer's choice, never a silent default (its numbers are
  not the scheme's).
- **Complement** (spec, the `letters` layer's `complement`): `auto` (the default, decided) draws a nucleotide
  as its complement when it is read against its strand on screen: the axis's direction times its path run's direction is
  −1. So a reversed DNA axis shows the other strand read 5′ to 3′, and an inverted run of a one-to-one map shows the
  strand it lands on. `always` (IGV's and UCSC's complement without flipping) and `never` (EC-233's rule, recorded as
  changed) are the writer's. A protein is never complemented (the alphabet decides). Letters shown never change data:
  model, status and search read the stored letters; colors by letter read the shown one.
- **Navigation is in display order**, fixed: a drag, the wheel, the pan keys and the cursor's arrows move what is on
  screen the way the hand or key goes, whatever the direction. Why: a key that moves against the screen is never wanted.

### Pieces: collapsed introns

- **`pieces`** on a section (spec): `whole` (the system, or its `range`), `coding` (the coding runs of a codon map: the
  follow path's, else the one named), a ref to a `ranges` node or a codon `map`, or a list of them (their union: two
  isoforms' exons), intersected with `range`. Default: `whole` for every section, followers included (decided). Why:
  the application drives the composition; a default read off the follow path would change the picture when a follow is
  added, and introns are data a genome reader wants. The p53 views spell `coding`. `collapsed` (spec, default `coding`)
  is the form `section.collapse` switches to; the toggle keeps no memory of an earlier form. `pad` (spec, 0) keeps n
  positions each side of every piece (splice sites, UCSC's padding; overlapping pads merge).
- **The stop codon** is in the coding pieces when the map's coding runs hold it (a coding length counting the stop,
  EC-17); it draws the gene's letters and no residue letter.
- **A ref that reaches no piece** (`coding` with no codon map, a ranges node on another system) is
  `seq.view.section.pieces` and the section shows `whole`.
- **A break** is where two pieces leave positions hidden between them, drawn as a hairline in the break ink
  (`seq-break`, [themes.md](themes.md)'s token) across the section's band and a notch on the axis edge, zero units wide on
  every layout, rings included (no ring gap: decided). Fixed. Why: the axis stays one display unit a position, the cells
  either side touch and a codon follow on a coding axis is exact; a visible gap has no use yet. Two pieces contiguous
  across a cyclic origin have no break.
- **Hitting a break**: within 3 CSS px of its notch on the axis edge, and on its hidden-item ticks; on data the cells and
  items beside it win (they touch; the hairline has no width of its own). Fixed, as other hit boxes are. A break's hit is
  `{ target: 'break', position, system, extent }`: the section's position at its edge, the hidden range as `extent`,
  whose hover names its length and the hidden items' count. A double-click on a break does nothing; its hint names
  `section.collapse`. Why: the pointer aimed at the break, so a fit would surprise; showing introns is a spec change.
- **What a break hides is not dropped**, fixed: a `ranges` lane with items wholly in a hidden stretch draws one tick on
  the break, whatever their number, hit as the break (its extent holds them, its hover counts them); a Location there (a
  projected hover, a selection, a highlight) lights the break. Why: drawing is checked against the data, and an item in
  view with no mark is a lost item; one tick a lane keeps the marks bounded by the screen.
- **Marks cut by a break are square there** (the frame rule): a domain across an intron, a projected selection band.
- **An item split by a path is joined by a connector** (the author, 2026-11-21: "if there's no connector, add one"): an
  item projected through a path into runs with shown positions between them that the path does not reach (a domain
  through the codon map over a gene whose intron is shown) draws a hairline at its mid-height, in its colour, from one
  run's end to the next's start; the runs' ends there stay square. It is part of the item: hovering it is hovering the
  item, and its emphasis follows the item's. Across a break (zero units wide) nothing is drawn. Why: a gene model's
  convention, so the pieces read as one item, not as several. Configurable (Configuration: `layer.connector`, `none`
  turns it off; its width and ink a `RenderConfig` preset). With `pad`, a gap that holds shown positions on both sides of a collapsed break is one
  gap: the connector runs through it, the break's hairline drawn over it. Only a hover that lifts a residue through a
  codon map lights its own section (its codon on the section's letters tracks); an item through a path lights its
  outline and its connector, nothing else (ux-foundations: "nothing else changes").

### Letters through a path

- A `letters` layer whose data is on another system is drawn when its path from the data to the section **pairs
  positions one to one** (identity, offset, segments, pairwise, an alignment row: a letter its image's cell) **or a
  residue to a codon** (a codon step: a residue over its three nucleotides). Any other path keeps
  `seq.render.layer.letters`, its reason now in the message: "a nucleotide on a protein's axis is a third of a residue
  and an intron or UTR nucleotide has no place there; a many-to-one path gives a cell several letters".
- **A residue a codon wide**: its fine places `3k`, `3k+1`, `3k+2` carried along the path as runs; each run on the axis
  is a part. A codon split by an intron is two parts, the letter drawn in each part with room (fixed: each part names its
  residue, EC-117's room rule decides). `code` (spec, the layer's): `auto` (the default: three-letter codes where the
  narrowest part in view has room, EC-117), `one`, `three` (a figure's choice).
- **The layer's density follows its own cell**: a residue is three nucleotides wide, so the protein's letters show at
  zooms where the gene's do not.
- **One to one**: a section position with no letter mapped onto it is an empty cell (the gap mark); a data letter with no
  position (an insertion) is a caret at its boundary (`projectBoundary`, EC-29).
- **Status**, fixed: a residue its codon does not translate to is drawn in the conflict ink once the section's letters
  are held (`codonStatus`, EC-16); nothing is said before they are. Why: a conflict is data, shown as everywhere else.
- **The protein and its gene on one axis** is a gene section with the protein's `letters` (and its tracks) through the
  codon map; its writer spells `direction: { along: p53 }` and `pieces: coding`. The two-section view (the gene
  following the protein) stays.

### What reads the axis

- **Offsets and positions.** A viewport is offsets in display order, as today; `positionAt` and `offsetOf` go through the
  pieces and the direction (a binary search over the pieces, no allocation). A fractional offset `u` in a piece `[s, e)`
  at offset `o` is the continuous position `s + (u − o)` forward and `e − (u − o)` reversed, so a window's fractional end
  covers the part of its element nearest the start of the window. An offset on a break (a piece boundary) is the next
  piece's first edge as a window's start and the previous piece's last edge as its end, so a window ending on a break
  shows nothing hidden. A saved session, a bookmark and a link store positions.
- **Hit-testing** keeps the section's system and position in every hit (what the cursor, a drag's `positionAt`, the
  selection's `offsetOf` and the focus fallback read). Letters through a path add, as a ranges item does, the data's
  element on `extent` (the residue, a Location on p53) and its parts on the section's system as `hull` (its codon);
  target `cell`. Granularity comes from the hit (ux-foundations section 8).
- **Selection**: a click on a residue selects its `extent` (the residue), projected onto its codon (closure, EC-30); the
  anchor is the residue on p53, so Shift+click, Shift+arrows and Shift+drag extend by residues on p53, in display order.
  Plain arrows move the cursor on the section's system by the row's granularity: a codon on a residue row (to the next
  residue's first nucleotide in display order), a position elsewhere, across breaks; a step that crosses a break is
  announced with the hidden stretch it crossed ("crossed intron 3, 1,240 positions hidden"), as the cursor is
  announced. A drag or an extend across a break
  selects the positions shown between its ends (`selection.acrossBreaks`, `ViewerConfig`: `shown`, decided; `span` adds
  the hidden stretch, for a reader who selects a gene's whole span to export).
- **Follow** (both ways): the follower's window is the image of the leader's window edges in display order.
  **Exactly three units a residue** when the follower's pieces are `coding` with `pad` 0, it has no `range` cutting the
  coding sequence, and the window is within the protein's residues: the follower's offsets are then three times the
  leader's, so it moves smoothly, keeps its scale and never crosses a hidden intron. Otherwise (`whole`, `pad`, a ranges
  ref with UTRs, a `range`) `Follow.point` follows fine places: residue `x` is fine place `3x`, each nucleotide of a split
  codon its own edge (today a split codon runs straight from its first nucleotide into the intron). A follower leads
  through the inverse of the same map (`follow.followerMoves`, as built), so with `follow` a pan of the gene moves the
  protein the same way on screen; with an explicit `genomic` on a minus-strand follower the other way (as today).
  `projectPoint` is exported from `resolve`; the viewer's copy, `pathPoint` (`viewer/src/models/view.ts`), goes.
- **Scroll groups** need equal axes: system, length, pieces and direction.
- **Focus and reveal**: the hull of the address on the section's axis in display offsets; hidden stretches take no
  width. An address wholly in a hidden stretch moves to its break (centered, the zoom kept) and the outcome says it is
  hidden by collapsed introns and names `section.collapse`.
- **The view history** (one rule for a step over a changed axis or layout, in every design that touches history): a
  step keeps each section's offsets with its axis signature (direction, pieces, layout) and its centre position and px a
  position. Steps compare (`snapshotsEqual`, the restored step) by offsets only on equal signatures; onto a changed axis
  or layout a step is restored by position (the centre and the zoom, as an update does, EC-111), clamped to what is
  drawn (a centre now hidden goes to its break); never silently dropped. Changing the direction or the pieces is a
  patch: undo, not back, returns it.
- **Updates** (EC-111): reversing or collapsing keeps the center position and px a position; a center in a stretch now
  hidden goes to its break.
- **Commands**: `section.reverse` (flips what the section shows: `genomic` and `reversed` swap; `follow` or `along`
  becomes the explicit opposite of what it resolved to) and `section.collapse` (`pieces` between `whole` and the
  section's `collapsed`), each one patch and one undo step. Their keys are in `view.keys` (Shift+R, Shift+C), acting on
  the section of the keyboard cursor, else the first leading section; a host rebinds them there.
- **Export**: the draw lists carry everything above, so the SVG is the canvas's picture; the section's group has
  `data-direction` and `data-pieces`; `renderSvg`'s `viewport` option may be given in positions per section, so a
  caller (the render CLI) needs no knowledge of the axis's direction or breaks.
- **Stories** ([stories.md](stories.md)): a change of direction or pieces between steps is a `cut`; the `view` channel
  interpolates positions.

## The interaction matrix's cells

Cells this design fills or changes ([ux-foundations.md](ux-foundations.md), section 8); the others stand as built.

| | select | extend | add / remove | focus | back | zoom | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|---|
| **a leading section** | — | — | — | — | — | about the pointer, in display order | in display order (drag, wheel, keys) | no: hide its tracks (as built) | Shift+R `section.reverse`, Shift+C `section.collapse` (also the section header's icons and the context menu): patches, one undo step each |
| **a follower** | as a leading one | as a leading one | as a leading one | moves its leader | as a leading one | leads, through the inverse of its axis | leads; with `follow` the leader moves the same way on screen | as a leading one | Shift+R writes an explicit direction (the opposite of what it shows); Shift+C |
| **an element: a residue through a path** | click: the residue, its codon projected | Shift+click, Shift+arrows, by residues on its system | Cmd/Ctrl+click | double-click: its codon's parts | as built | — | — | `track.hide` | no: data is not edited in the viewer |
| **an element: a break** | no: what it hides is not shown, so it is not selected | no, as select | no, as select | no: a double-click does nothing (its hint names Shift+C) | — | — | — | no: a break is not a track; Shift+C shows or hides every break of its section | no |
| **a range** | drag: shown positions only across a break; whole residues on a residue row | Shift+drag, the same | Cmd/Ctrl+drag, Alt+drag, the same | Enter: its hull in display offsets | as built | — | — | — | — |
| **the selection** | — | Shift+arrows in display order, across breaks | Cmd/Ctrl+Shift+arrows | Enter | as built | — | — | — | — |
| **the view** | — | — | — | fit shows every piece | back lands on the same positions whether the axis changed or not | — | — | — | undo returns a direction or pieces change |
| **a track or layer** | as planned ([grids-and-rows.md](grids-and-rows.md)) | — | — | as planned | — | — | — | as built | a letters layer's `complement`, `code` by patch |
| **a bookmark or link** | — | — | — | runs on positions; into a hidden stretch, its break and the reason | as built | — | — | — | — |
| **a navigator** ([navigators.md](navigators.md)) | — | — | — | — | — | — | its strip shows the section's direction and breaks (the ruler rail, [navigators.md](navigators.md)) | — | — |

## Configuration

| Item | Whose | Default | Why the default; the alternatives |
|---|---|---|---|
| `section.direction` | spec param | `follow` for a follower, else `genomic` | decided defaults, by the section's own `follow` param, never inferred from layers; `genomic`, `reversed`, `{ along }` |
| `layer.units` (ruler) | spec param | `none` | plain, never inferred from the alphabet; `bases` (bp, kb, Mb by the tick step) |
| `layer.unnumbered` (ruler) | spec param | `mark` | a blank past a numbering reads as a bug; `fallback` with `layer.fallbackNumbering` (ordinal by default), muted |
| `section.pieces` | spec param | `whole` (every section, followers included: decided) | the application drives the composition; `coding`, refs (their union) |
| `section.collapsed` | spec param | `coding` | the form `section.collapse` switches to; a ref or a list |
| `section.pad` | spec param | 0 | the pieces as they are; n positions each side (splice sites, UCSC) |
| `letters.complement` | spec param | `auto` | read along the strand shown (decided); `always` (IGV, UCSC), `never` (EC-233's rule) |
| `letters.code` | spec param | `auto` | three letters where the narrowest part has room (EC-117); `one`, `three` |
| `axis.directionMark` | `RenderConfig` preset | `label` | the picture and the export say "reversed"; `none` for a host's own chrome |
| `layer.connector` (ranges) | spec param | `line` | a gene model's convention, the pieces one item (the author, 2026-11-21); `none` |
| `ranges.connector` (width, ink, its alpha) | `RenderConfig` preset | 1 CSS px, the item's colour, 0.6 | a hairline under the runs' weight; a host's own (a dashed or a muted ink) |
| `selection.acrossBreaks` | `ViewerConfig` (the reader's) | `shown` (decided) | what is selected is what is seen; `span` |
| keys of `section.reverse`, `section.collapse` | `ViewerConfig` (`view.keys`) | Shift+R, Shift+C | a keyboard way for each (ux-foundations section 8); any key the host prefers |
| a viewport in positions | option of `renderSvg` (`viewport`) | offsets, as today | a caller with no knowledge of the axis |

Fixed, with the reason in the Decision: the majority direction of a mixed path, navigation in display order, the break
as a zero-width hairline on every layout and its 3 px hit, one tick a lane for hidden items, a split codon's letter in
each part, codon status, a double-click on a break doing nothing, a history step over a changed axis restored by
position.

## In the specs

| Decision | Where |
|---|---|
| The direction the writer's, defaults by `follow`, never inferred (Decided 1) | core › Its direction |
| `pieces` default `whole`; `section.collapse` toggles to `collapsed` (Decided 2) | core › Its pieces |
| A collapsed intron zero width on every layout (Decided 3) | core › Breaks |
| A drag across a break selects what is shown (`selection.acrossBreaks`) (Decided 4) | viewer › Configuration |
| Expanding one intron refused; Shift+C the whole section (Decided 5) | viewer › The axis's commands |
| A history step over a changed axis restored by position (Decided 6) | viewer › History |
| Whose each choice is; a config item needs a second use | design-only: the rule designs place their choices by, applied in each spec's Configuration |
| `direction` (`genomic`, `reversed`, `follow`, `{ along }`); unresolved and mixed paths | core › Its direction |
| One direction per axis; a reversed piece's positions | core › An axis |
| A direction per piece inside a block | design-only: not built (a tension) |
| A reversed axis says so: the ruler counts down, the direction mark (`axis.directionMark`) | core › A reversed axis says so |
| A ruler's insertion labels, `units`, `unnumbered` | core › Units |
| Complement (`auto`, `always`, `never`) | core › Complement |
| Navigation in display order | core › Navigation |
| `pieces`, `collapsed`, `pad`; the stop codon; a ref with no piece | core › Its pieces |
| A break: a zero-width hairline; its hit | core › Breaks |
| What a break hides is not dropped | core › What a break hides is not dropped |
| Marks cut by a break square there | core › Every emphasis is a squircle |
| An item split by a path joined by a connector | core › The connector |
| Letters through a path: one to one or a residue to a codon; `code`; density; status | core › Layer types |
| Offsets and positions through pieces and direction | core › Offsets and positions |
| Hit-testing keeps the section's system and position; letters through a path add `extent` and `hull` | core › Lanes and Hits |
| Selection by residues; arrows by the row's granularity, a crossed break announced | viewer › The cursor |
| Follow both ways; three units a residue on a coding follower; `projectPoint` exported | core › A follow |
| Scroll groups need equal axes | core › Scroll groups' key |
| Focus and reveal through hidden stretches | viewer › Actions |
| The view history over a changed axis | viewer › History |
| Updates keep the centre and the zoom | viewer › Updates |
| Commands `section.reverse` and `section.collapse` | viewer › The axis's commands |
| Export: `data-direction`, `data-pieces`; a viewport in positions | core › A viewport in positions |
| Stories: a change of direction or pieces a `cut` | core › Stories |
| The spec params (`direction`, `pieces`, `collapsed`, `pad`, the ruler's, the letters') | core › The axis's params |
| Precedents, the interaction cells | design-only: intent |

## Tensions

- **Complement by default against EC-233** ("reversing never complements by itself", the JBrowse reply's "one strand of
  letters"). Complement is on by default (decided); the rule is explicit and visible (`complement: auto`, read against the
  strand on screen), and the stored letter is named wherever the letter is. EC-233 is recorded as changed.
- **One direction per axis**, narrowed (the fifth slice, [genome-coordinates](genome-coordinates.md), Decided 1): a
  direction per block is built (`block.direction`: an inversion between genomes drawn as its own block, with a gap); a
  direction per piece inside a block stays open, the pieces' shape (display order, start below end) not ruling it out.
- **Hidden is not gone**: a collapsed intron hides a splice-site variant. Ticks on the break, `pad`, focus naming the
  reason and the toggle answer it; expanding one break is refused for now (decided): a per-break session override
  would make the export differ from the screen; a per-break patch is a `pieces` ref edit, possible by hand. Standing.
- **The inverse placement**: the gene's letters on the protein's axis (a nucleotide a third of a residue, introns with no
  place). Refused with its reason; the one-axis view on the gene covers the figure. Standing until a figure wants it.
- **A reader cannot reverse or collapse**: a reader of a published, read-only figure can recolor (the color-blind
  setting) and fold clades (session state), but not reverse a section or expand an intron: direction and pieces are
  spec params, changed by patches. Standing: they are what the figure is; a session override would make the export
  differ from the screen.
- **Defaults the writer must spell**: a gene following its protein shows its introns, and the one-axis view reads
  against its protein, unless the spec says `coding` and `along`. Accepted (the author, decided): the composition is the
  application's; the examples and the builder's docs show the spelling.
- **Two proteins on opposite strands on one gene section**: one direction; the writer picks which it reads `along`, and
  the other reads right to left. Accepted: one direction per axis.

## Self-check

The claim: **this design covers what a user does with a section's axis.**

1. **Assume it holds.** A follower reads along its leader, and keeps its letters where its writer collapses its introns;
   a gene section carries its protein a codon wide, in the direction its writer asks; hit-testing, selection, follow,
   focus, history and export read positions through pieces and the direction.
2. **What it rests on, negated.** It rests on: (a) what the axis hides is not wanted; (b) the user reads along one
   reference, the leader; (c) letters go from the coarser system onto the finer axis; (d) one direction holds a whole
   axis; (e) one codon map decides the pieces. Negated: (a) the user wants what is in the intron; (b) the symmetric case:
   the reference is not a leader but another system (the protein a gene section shows), or the user wants genomic
   coordinates ascending; and the plural: two references asking opposite directions; (c) the gene's letters on the
   protein's axis; (d) an axis that turns over mid-way; (e) two isoforms on one gene section.
3. **The negated world.** (a) explains the splice-site variant lost in a collapsed view and a focus on an intronic
   address that finds nothing; (b) explains the one-axis figure reading C to N with a leader-only rule, and two
   references with no rule between them; (c) the symmetric case of letters through a path; (d) inversions, the Genomes
   milestone; (e) one map's coding pieces cutting the other isoform's exons.
4. **Both held.** Integrated into the Decision: ticks and hits on breaks, a break's hover with the hidden length, `pad`,
   a focus into a hidden stretch moving to its break with the reason, `section.collapse` both ways with `collapsed` as
   its plain target (from (a));
   `{ along }` as the writer's way to name any reference, the plural settled by the writer naming one (no inference),
   explicit `genomic`, `section.reverse` both ways, undo returning either (from (b)); `seq.view.direction.mixed` (from
   (d)); `pieces` as a list of refs, their union (from (e)). Kept as standing tensions: expanding one break; the gene's
   letters on the protein's axis (from (c)); a direction per piece (from (d)); defaults the writer must spell, and the
   reference not named read against itself (from (b)).

## Consequences

- **`seqx-core`**:
  - language: `section.direction` (`genomic`, `reversed`, `follow`, `{ along }`), `section.pieces`,
    `section.collapsed`, `section.pad`; the
    `letters` params `complement`, `code`; checks (a `pieces` or `along` ref's kinds); the builder's docs show the
    one-axis spelling.
  - resolve: `projectPoint` exported; fine-place points (each nucleotide of a split codon an edge); the coding runs of a
    codon map on its nucleotide system.
  - view: `Axis.dir` resolved from the param; pieces from `pieces` and `pad`, its breaks with their hidden ranges;
    `seq.view.section.pieces`, `seq.view.direction.unresolved`, `seq.view.direction.mixed`.
  - render: `AxisTransform` with the axis's direction (offset and position both ways, fractional ends, spans flipped, cut
    flags at breaks); breaks and their ticks; `letters` through one-to-one and codon paths (parts, codes, carets, status,
    complement); ranges mirrored; the direction mark (`axis.directionMark`); residue hits (`extent`, `hull`) and break
    hits; `seq.render.layer.letters` narrowed, with its reason; the SVG's `data-direction`, `data-pieces`;
    `renderSvg`'s viewport in positions.
- **`seqx-viewer`**: `ViewModel` offsets and positions through the axis; follow and lead on the new axis; history steps
  with axis signatures and centre positions; navigation and the cursor in display order across breaks; selection across
  breaks (`selection.acrossBreaks`) and by residue; focus into a hidden stretch; the `break` hit target; the commands
  `section.reverse` and `section.collapse` (change kind `patch`) with their keys in `view.keys`; `pathPoint` removed.
- **`seqx-examples`**: the application drives the composition: each p53 view spells its direction and pieces in its
  builder. The gene following its protein (`direction: follow`, the follower's default written out for the reader;
  `pieces: coding`); the protein and its gene on one axis (`direction: { along: p53 }`, `pieces: coding`); the genomic
  follower for contrast (`direction: genomic`); the reversed gene with its complement. Each three ways (Storybook canvas
  and SVG, a quick-lane test drawing against the data, a portable spec), with headless twins.
- **Storybook**: the stories above, and a Design page "The axis" from this design.
- **Tests**:
  - property tests over random sections (direction, pieces, `pad`, a cyclic origin): every shown position placed once;
    hit-testing the inverse of placement; offset and position round trips, fractional ends included; selection and
    export agreeing; a follow on a coding axis exactly three units a residue under its conditions; lead then follow
    returning the window; every nucleotide of a split codon reachable as a follower's edge.
  - the direction param: each value resolved, `along` a system no path reaches and `follow` with no leader reported
    unresolved, a section's direction unchanged when a layer is added.
  - hit-testing on a reversed collapsed axis with residue rows, through the readers of a hit: the cursor, a drag's
    `positionAt`, the selection's `offsetOf`, the focus fallback; a break's slop against the cells beside it.
  - the view history: back and forward across a reverse patch and a collapse patch landing on the same positions; a pan
    whose window ends on a break showing nothing hidden.
  - drawing against the data at the fit view and across a tile edge, on the reversed and collapsed p53 views: every
    item in view has its marks (a hidden one its break's tick), every mark is hit.
  - headless twins of the one-axis and following views (canvas and export, two sections) on a fake clock; records
    (`recordDrawLists`) of both p53 views.
  - performance: a quick representative under `itWithin` (a pan frame over a reversed axis of 150,000 single-position
    pieces, EC-33, under the pan target); draw and hover warm at 0 bytes on a reversed collapsed axis (the opt-in tier).
- **`_spec` sections that change**: core Language (section, letters), Resolve (Paths: points), Resolved view (axis,
  follow, issues), Render (Configuration, Layer types `letters`, Tracks layout, Lanes and Hits, Shapes, Export, Issues);
  viewer ViewModel (viewport, follow, history), SelectionModel, Commands, Bindings (the `break` target, the two keys),
  Configuration.
- **Other designs**: [features.md](features.md)'s three rows (the direction row gains `{ along }`); [edge-cases.md](edge-cases.md) (EC-233 changed;
  EC-31, EC-32, EC-117, EC-144 kept); [viewer-models.md](viewer-models.md) (collapsed introns built); [ux-foundations.md](ux-foundations.md) section 8's matrix
  (the cells above); [maps.md](maps.md) ("reversed pieces follow later" closed).
