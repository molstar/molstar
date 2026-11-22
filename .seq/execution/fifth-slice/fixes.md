# Fixes: the fifth slice's small findings

Small findings that need no design collect here for the slice's fix points, as the fourth slice's did (the author,
2026-11-21: "lets have a 'fix phase' for slice 5 too where we accumulate these smaller things?"). An FX lane takes the
open items at each fix point (plan, Fix points); a finding that needs a design becomes a drift note instead.

## Open

- **One hover contract for every mark, and a sweep that holds it** (the author, 2026-11-21, Examples › Plasmid › Ring
  at scale: the strands' density bands on the ring and in the detail show no hover box or card; "we should not be
  fixing this one-by-one btw. since we already fixed a similar issue, we need a general rule for this and global
  fix."). The same miss has come back layer by layer: a ring's own hover (fixed alone in FX3), the Regions density row
  (below), a ligand's score box (below), now density bands on a ring and a strip. **The rule** (to write into
  ux-foundations.md and the core `_spec`'s Emphasis): every drawn mark of every layer type, in every layout (row,
  ring, blocks, rows, grid) and at every level of detail (an item, a lane, a density bin or band, a summary, a chip,
  a reduced tile), answers a hover with the same three things: (1) a hit naming what is under the pointer (the item,
  else the bin's items or stretch, else the position), (2) an emphasis mark drawn in the hovered view and in every
  view linked to it, reaching every layer the covered positions reach (arcs, columns, the ring's position box), and
  (3) a card with that hit's fields; a layer type that can't is refused at registration, not left silent. A drag
  suspends (2)'s card (the drag item below). **The sweep**: one test (medium lane, canvas and SVG) over every story's
  spec and every layer instance, probing points on its drawn marks (taken from the draw list, not guessed), asserting
  the hit, at least one emphasis primitive in the hovered view and each linked view, and a non-empty card; a new
  layer type or layout is covered by adding nothing. The per-story items below that it covers (the density row, the
  score box, these bands) are checked by it and closed with it. A strong-model lane (it is a rule across layer types,
  not a closed list), its brief naming the items it closes.

- **A legend's hover works one way only** (the author, 2026-11-22, Features › Genome formats › VCF › TP53 ClinVar:
  "there is an asymmetry for many legends ... where the highlight only works one way. add to fix round? or at least
  don't show 'hand' icon for legends that don't have any action associated with it"; kind: bug, a class): hovering a
  variant outlines its legend entry ("deletion"), but hovering the "deletion" entry lights no variant, though
  themes.md lists "hovers a legend entry to see its items" as built (a highlight channel over its items in view). The
  rule, with the hover contract above and swept by it: every legend entry that names items (a category, a ramp's
  stretch) lights its items in view when hovered, in every layer type and layout, canvas and SVG alike; an entry with
  no action (no select, no filter, no highlight: a swatch only) shows the default cursor, never the pointer. The sweep
  probes each story's legend entries: hover → at least one emphasis primitive over that entry's items in view, or no
  pointer cursor. Taken by the hover-contract lane (strong), after RG merges (RG changes the hover texts and card).

- **No track drawn as an empty well, and a sweep that holds it** (the author, 2026-11-21: "Yes, group them", after
  the Promoters track, the stacked ATP block's tracks and the plasmid's strip each showed a blank band). Three causes,
  one symptom: (a) the data has nothing for a declared track (Promoters: a fixture gap, the story's fault); (b) a
  layout repeats tracks for a block or section that has no items there (stacked blocks: the library's); (c) a
  navigator or strip declared with no content (the plasmid's strip: the library's default). **The rule** (into
  ux-foundations.md and composition.md, with navigators.md for strips): a track with no items in what it covers never
  draws a full-height blank well. By default it collapses to a one-line row naming it with a quiet "no items" (so a
  host's "none found" stays visible); a stacked block's track with nothing for that block is left out of that block;
  a strip with no `content` draws its system's ticks. Configurable as the host's (`track.empty: 'note' | 'collapse' |
  'keep'`, `blocks.emptyTracks`), the defaults above. A track still waiting for its source is pending, never empty.
  **The sweep**: one test (medium lane, canvas and SVG) over every story's settled spec: every track, block row and
  strip draws at least one mark besides its chrome, or is collapsed by the rule; a story that means an empty track
  says so (`expectEmpty`), so a fixture gap like Promoters fails the sweep. The three items below (the Promoters
  track, the stacked block, the plasmid's strip) are its instances, closed with it; their story fixes (re-record the
  plasmid's GFF3 with promoters or a Sites track, place the strip by the detail) go in the same lane. A strong-model
  lane, like the hover contract.
  Also an instance (the author, 2026-11-21, the genome ring with its contigs: "feels weird it's empty"): the ring's
  block "2,976 contigs" is labelled with a leader line but draws no arc (its share of the genome a sliver at that
  size), so a named block shows nothing; the rule's form for a block too small to draw: drawn at a least arc
  (`blocks.minPx`) with a break mark on each side, its items as a density if any, never a label pointing at nothing.

- **A drag's work is bounded by the screen, not by the selection** (the author, 2026-11-21, tRNA › Views › Dot plot
  scale: "selecting a rect on the dotplot here is really laggy. both as a rect on the plot, and a range on the nav.
  are we doing too much work here? maybe related to the other perf issue? and similar lag to selecting on the genome
  ring with the contigs"). One class with the PAE drag and the 3,000-contig drag-select below. What I read (not yet
  profiled): every pointer move during a drag-select runs a full hit test (`viewer.hit`, all layers) and writes a new
  selection value (`selection.rectDragTo`, `viewer/src/models/selection.ts:586`; `select.ts:316` for a range), and
  each new value re-runs everything derived from the selection: its emphasis in every view and linked section (the
  grid's rect and its mirrored bands on both rulers, the ring's arcs), the hover card's lookup, the SVG twin. The
  fix lane profiles first (a perf-tier scenario a drag kind: a grid rect on the 4,000 × 4,000 dot plot, a range on its
  ruler, a sweep across 12 chromosomes on the contigs ring, a box on the PAE grid), then makes the rule hold: a drag
  move does work in proportion to what is on screen: the selection written at most once a frame (coalesced to the
  frame, not the event), its emphasis derived for the visible part only and incrementally (a rect's growth, not the
  whole rect again), no hover lookup or card while dragging, the SVG twin and linked panels updated on the frame, the
  full selection resolved once on release. Bound each drag kind's frame to the budget at 0 bytes warm in the perf
  tier. A strong-model lane; it closes the PAE drag item and the 3,000-contig drag-select part below.

- **The novice decisions, to build** (stage 1's gate, the author, 2026-11-21: "Yes, all five"; novice-findings.md):
  F3 the fake's "Within 5 Å" on a selection that is already the chain selects the residues of the other chains within
  5 Å (never the selection again), its card "n residues within 5 Å of chain A"; F9 the track the structure is coloured
  by shows its icon pressed (one at a time, the icon fixed in place, its tooltip naming it); F11 an empty section says
  its issue in its own row (with the empty-well rule above); F14 the "n filtered" card names the filter ("1 filtered:
  Species = Human") and says "click to clear"; F25 a navigator window at the fit is drawn as one quiet frame, no grips
  and no drag target (with the strip item), the Coverage strip story opening zoomed. Tests one a decision.

- **The canvas and the SVG of one spec lay out differently** (the author, 2026-11-21, tRNA › Views › Seed and chain:
  "seems like a layout issue. the canvas is big, svg small"). Same spec, two sizes: on the canvas the seed's rows band
  grows to every row (the page scrolls, the 1EHZ chain track pushed out of sight below it); on the SVG the band stops
  at about twelve rows with its own scroll and 1EHZ shows under it. In both, the first row is half hidden under the
  Conservation track (its top cut). The two renderers are twins over one layout (core `_spec`, Render); one of them
  takes a different height: likely the story frame giving the canvas its content's height and the SVG a fixed one,
  or the rows band's cap read on one side only. A class, by the rule above (fold-back B2): **the twin sweep**, one
  test (medium lane) over every story's spec, laying it out for the canvas and for the SVG at the same frame size and
  asserting every section's, track's and band's box equal (and the scroll state), so a story frame or a renderer that
  sizes on its own fails at once. With it the rows band's height follows the track-sizing design (content height to a
  cap, then its own scroll; the item below), the same on both, and the first row starts below the track above it.

- **The hover card: placement, wrapping, and the density row's missing highlight** (the highlight part closed by the hover contract above; the author, 2026-11-21, Files ›
  UniProt entry: "hovering the bottom 'regions' … is missing a highlight", "the hover placement is quite
  inconsistent"):
  - the Regions track's bottom row (its density row of overlapping regions) shows a card but draws no hover mark: the
    density row lights the hovered stretch (or the items under it), as a lane's item does (the top row's
    "Oligomerization" gets its box and a short card: Track, Extent); the density row's card lists the items under the
    pointer in that same short form;
  - the card jumps between below and above the pointer and covers the rows being read: one placement rule (beside
    the pointer on the side with room, never over the hovered row, the same side while the pointer stays in a track);
  - the card's width changes from one hover to the next and wraps values ("P04637 100–" / "370"): a value never
    wraps inside a range; the card sizes to its longest row up to a cap, then truncates the name, not the value;
  - unnamed features are labelled by a bare row number ("152", "153"): fall back to the feature's type and range.
  Tests: the density row's hover draws a mark; the card's box never intersects the hovered row; a range value on
  one line.

- **The navigator's window and zoom strip look inconsistent and rough** (the author, 2026-11-21, a ring that switches:
  "i don't understand why chr6 gets a 'box' and the other 'sub-section'"; "the 'subsection visual' is really ugly,
  with different line thicknesses"). Why it differs: a window wider than `navigator.zoomStrip.minPx` (48 px) is drawn
  as the box with grips (chr6, about 85 px), a narrower one (chr13) gets the zoom strip under the navigator, with its
  funnel lines. Needs a short look design before a fix lane: one rule a reader can see (the strip when the window is
  too small to grab, or always for a section of one block), one stroke weight for the window, the funnel and the
  strip, the strip's fill and border from the same tokens as the window's (ux-foundations: the loudest mark at rest),
  the funnel thin and muted or dropped; its picture in navigators.md, then the fix.

- **The vertical overview cut at the right on a wide page** (FX5 2026-11-22: not reproduced at 960, 1600, 2200 px; the grids made shrinkable, `minmax(0, 1fr)`, as hardening; open until seen at a wide page) (the author, 2026-11-21, Globins rows family-rail-pinned):
  at full width the story's canvas runs past the page's right edge, so the vertical overview is clipped (a narrow
  page, where canvas and SVG stack differently, shows it whole). The story's layout gives the canvas its measured width
  at every page width, the overview inside it; a story-layout test at a wide and a narrow width that the overview's box
  lies inside the canvas and the canvas inside the page.

- **A follower's rail turns a window over the cyclic origin into its complement** (the author, 2026-11-21, Files ›
  Plasmid GFF3, dragging the detail's rail: the bar leaves the pointer, grows to most of the axis, then jumps to the
  start; "maybe a symptom of a bigger problem"). `core/src/render/rail.ts:186`: a follower's window is read back from
  its two end positions and ordered with `min`/`max`; when the window crosses pUC19c's origin (end position < start
  position), that gives the rest of the ring (e.g. 60–2600) instead of 2600–60. The non-follower branch just below
  already keeps `end > total` over a cyclic origin and the rail already draws such a window as two pieces
  (`rail.ts:34`, `window.ts:43`); the follower branch should do the same (the end past the total when it runs over the
  origin), and a drag across the end wrap the bar through its two pieces, the pointer's offset kept modulo the total.
  **The bigger problem**: ordering two ends with `min`/`max` is wrong on any cyclic axis; audit every site that does it
  on positions or units that can be cyclic: `navigator.ts:517`, `:1163`, `navigator-parts.ts:255`, `:349`,
  `export.ts:400–405`, `items.ts:711`, `viewer/src/input.ts:1079` (a selection or export across the origin likely
  inverts the same way), fixing each through one shared helper (a cyclic interval from its start and end, the core
  `util` rule) and a twin test per site with a window across pUC19c's origin. Also in the same screenshot: a straight
  line from the ring's centre through the follower's window wedge, and the wedge wider than the detail's span; check
  with the fix.

- **The hover card follows a drag on the grid and makes it laggy** (the drag part an instance of the drag rule above; the author, 2026-11-21, Examples › Complex › PAE
  and contacts › CDK2 monomer: "showing the hover when dragging on the PAE plot makes the view feel really laggy"):
  while a drag sweeps a box on the PAE grid, the hover card (x, y, value) and the crosshair are re-resolved and redrawn
  on every move, on top of the selection box and its mirror on the section above. During any drag (a grid's box, a
  track's range, a rail, a navigator's window) the hover card is hidden and its lookup skipped, shown again on release
  where the pointer rests (the drag's own readout is what it shows meanwhile); then profile the remaining drag frames on
  this story (298 × 298 grid) and bound a box drag to the frame budget, 0 bytes warm, in the perf tier. Same family as
  the readout covering the active control above.

- **A selection's boxes on the grid don't meet the cells' edges** (the author, 2026-11-21, Examples › Complex › PAE and
  contacts › Over the complex, a range selected on cyclin A2's pLDDT: "mis-alignments on the boxes on the PAE plot (see
  the gaps at the edges)"). From the screenshot: the row band on the cyclin row runs from about x 597 to 1131 while the
  CDK2 cell spans 600–1128, so it overhangs a few px into the blocks' gap on both sides; the column band's ends and its
  grip run past the grid's top and bottom; the column and row bands overlap on the cyclin × cyclin cell (double
  shading). Likely the mirrored box is placed in the grid's whole extent (or with its stroke centred on the edge)
  rather than clipped to each cell the blocks' gaps cut it into (the story's text promises "a dragged rectangle split
  at the boundaries"). Each piece is clipped to its cell, its stroke inside, the intersection of a row and a column band
  drawn once; a twin test (canvas and SVG) that every selection piece's box lies inside its cell's box, at both a block
  boundary and the grid's outer edge. Not measured in the browser yet: the fix lane measures first.

- **The plasmid's strip is an empty bar at the top** (an instance of the empty-well rule above; the author, 2026-11-21, Examples › Plasmid › Ring › Detail with
  a strip: "feels quite weird to have the empty strip on the top"). `plasmidStrip` (`examples/src/navigators.ts:218`)
  declares `stripNavigator('/views/main/detail', { system })` with no `content`, so the strip draws only its window
  (here nearly the whole construct) over a blank well, no ticks or features; and it sits at the frame's top, above the
  ring, far from the detail it navigates. `tracksStrip` (`content: 'tracks'`) shows what it should look like. Fix: the
  story's strip gets content (the detail's feature tracks small, or the ruler's ticks at least) and sits next to the
  detail (between the ring and the detail, or the place navigators.md gives a strip over a section); and a library
  question for the design: a strip with no `content` draws the system's ticks by default rather than an empty well. A
  story test that the strip has marks besides its window.

- **Zooming a lazy alignment out past a pixel a row leaves a blank plot** (the author, 2026-11-21, Examples › Globins ›
  Rows › Lazy, the rows rail's grip dragged to widen the window: "trying to resize the vertical bar breaks the lazy
  views. does our streaming protocol correctly ask for a window to fill?"). The window is asked correctly while rows
  are a pixel or more: `readLazy` (`core/src/render/rows-read.ts:196`) reads the tiles of rows `[first, last]` at level
  1 only, the rest pending. Below a device pixel a row it asks nothing: it reports `seq.render.rows.reduced` and draws
  the neutral fill (core `_spec` Rows, Data), whether or not the source has coarse levels; the code's own comment says
  "a source's coarse row levels, later", and the message "serves no row bins" is said without asking. So the rail can
  widen the window into a state nothing can draw (a grey plot, the labels gone, a bare card "rows 2,752,442 / index
  2752441"). Two parts: (1) **now**: while a lazy source serves no row levels, the rows' zoom-out stops at a device
  pixel a row (the rail's grip, `[`, the wheel), so the window can't reach the blank state; (2) **a design** (with
  region-reads.md and protocol.md): ask the source for the row level that fits (`readTiles` at level `2^k`, k from rows
  per device pixel) and draw its bins by class like the held pyramid; the generated source serves them (cheap: sampled
  from its profile); the neutral fill and the issue only when the source declares no row levels. Tests: the rail can't
  pass the floor on a source without levels; with levels, a window of the whole 10 M rows reads tiles of one level
  only, bounded by the screen.
  **The limits, and making them the host's** (the author, same day: "what are the limits of this? seems like this
  needs to be configurable (e.g. disable the resizing for certain views)"). Today the rows' zoom bounds are fixed
  constants in core, `MinRowScale` 1/64 px and `MaxRowScale` 40 px a row (`core/src/render/rail.ts:53`), the same for
  every view and source, read by the zoom keys (`viewer/src/behaviors/rows.ts:180`), the rail's grips and its idle
  rule; nothing turns the zoom off (the columns' axis has `viewport.maxPxPerPosition`, 32, configurable). Proposed, in
  the same design: `rows.zoom` in `ViewerConfig` (`minPx`, `maxPx`, `enabled`, defaults today's) and a band's own
  override in the spec (a fixed row height, or zoom off for that band); zoom off keeps the rail as a pan-only bar (no
  grips, the keys and wheel-with-modifier inert, scroll kept); the effective floor the larger of `minPx` and the
  source's (a device pixel a row for a lazy source without row levels, part 1 above). The library rule (memory: the
  host configures optional behaviour) and grids-and-rows.md's zoom rows updated with it.

- **A tree opened from its toggle has no way back on screen** (the author, 2026-11-21, Examples › Globins › Rows ›
  Scale: "when i expand the tree, there doesn't seem to be a way to collapse it back. another port for app level
  control?"). The narrow-frame rule (grids-and-rows.md, "A narrow frame gives way in one order", step 3; core `_spec`
  Rows, `RowsLay.toggleX`) collapses the tree to its glyph; a click or Shift+T (`rows.tree`) opens it whole "and folds
  it back", but once open the glyph is gone, so only the key folds it. Keep the toggle at the band's top while the
  tree is kept open (the same glyph, pressed, `tree-toggle` hit, its tooltip naming Shift+T), and give the command an
  app-owned control the same way as the sort's clear: an entry "Show tree" / "Fold tree" in the band's menu, and
  `rows.tree` in the viewer's commands for a host's button. A twin test: open from the glyph, fold from the glyph, the
  layout back to its narrow form.

- **Hovering a ligand's score box lights its residues but not its arcs** (closed by the hover contract above; the author, 2026-11-21, Examples › Complex ›
  Views › As a complex: "mouse on empty patch highlights arcs, but on box only the columns"). Over the score track's
  empty band at PYZ's position, the hover resolves to the position (card "PYZ 1"), its contact arcs light dark and the
  linked CDK2 columns light. Over PYZ's score box itself (a ranges row over the molecule), the box gets its outline and
  the linked columns light, but no arc lights and no card shows. A ranges item's emphasis should reach every layer the
  positions it covers reach (the pairs layer's arcs with an end in it), the same as hovering those positions, and its
  card show (the item's fields: entry, resolution). A test on this story: hovering the box and hovering the empty band
  at the same position light the same arcs, and the box shows a card.

- **The blocks' Contacts track leaves ATP out** (the author, 2026-11-21, Examples › Complex › PAE and contacts ›
  Contacts: "does ATP really have 0 contacts here?"). It has: the recording (`complex/derived-data.ts`, `ligands`) holds
  ATP A's 43 heavy-atom contacts with 15 CDK2 residues and ATP C's 48 with 14, and the units stories draw them; but
  `crystalData` (`complex/blocks.ts:46`) builds the one pairs node `a: cdk2, b: cyclin`, so the ATP block shows as a
  chip with no arc. Add ATP's contacts to the blocks' Contacts track (a second pairs layer, or rows of ATP to CDK2 by
  residue, copies merged like the rest: an arc an ATP–residue contact of either copy), the story's text naming both
  kinds, and a story test that the ATP block has arcs.

- **A stacked block repeats tracks that hold nothing for it** (an instance of the empty-well rule above; the author, 2026-11-21, Examples › Complex › Blocks ›
  Stacked: "another view with empty tracks"): with `arrange: 'stack'` each block is a row of every track, so the ATP
  block (`complex/blocks.ts:69`, a one-position block drawn as a chip) gets empty Features, Interface and Secondary
  bands, three tall blank wells above one chip; its row also has no heading where the polymers have "CDK2 (A) ×2".
  In a `row` arrange the tracks are shared by all blocks, so nothing is empty there. Needs a short design before a fix
  lane (composition.md's blocks, molecules.md's ligands, the library rule that such a choice is the host's): a stacked
  block's track with no items for that block collapses to nothing by default (the letters track kept), configurable
  (`blocks.emptyTracks: 'collapse' | 'keep'`), and the block gets its heading like the others ("ATP ×2"). Together
  with the Promoters item just below: the same family (a declared track drawn empty), here a library rule, there a
  fixture.

- **The Promoters track is empty** (an instance of the empty-well rule above; the author, 2026-11-21, Files › Plasmid GFF3: "feels weird the promoters track is
  empty"): `files/plasmid-gff.ts`'s `GffTracks` declares a `promoter` track (on the ring and the detail), but the
  fixture `plasmid/pUC19c-rotated-features.gff3` has no promoter row: it is generated (`record.ts`, `gff3Text`) from
  `plasmid/derive.ts`, whose features are all found in the sequence (examples.md, Plasmid: none taken from another
  record), and no promoter is derived. Also the file's 14 `restriction_enzyme_recognition_site` rows have no track.
  Either derive the two promoters (lac: its -35/-10 boxes upstream of lacZα; AmpR's upstream of bla, across the origin
  like bla) and keep the track, or drop the track and add a Sites track for the rows the file holds (the better
  showcase of a type with many short items); re-record the fixture, and a story test that every declared track of the
  story has at least one item.

- **The 3,000-contig story could be faster** (the drag part an instance of the drag rule above; the author, 2026-11-21, on packages-seq's build: "better perf ... although
  still could be faster"): my measure, 40 hover steps along the row each drawn within two frames (about 33 ms, no
  drops). Take it at M5: profile hover and pan on three-thousand-contigs (ring and row), find the per-frame cost, and
  set a bound in the perf tier. **A drag-select across many chromosomes is slow** (the author, 2026-11-21: chr11 to
  chr22 selected on the row, mirrored on the ring): profile the drag's frames (the selection's blocks and their
  emphasis on both views), bound a drag across 12 blocks to the frame budget at 0 bytes warm.

- **Sizing a track of many lanes: content height, a cap, then its own scroll** (the author, 2026-11-21, gff3--tp-53:
  "can we do a better job at the track sizing for the transcripts in general? what if that section needs a scroll?";
  earlier: blank lanes above and below, the first row unlabelled). Needs a short design before a fix lane (Ranges,
  ux-foundations' heights): a ranges track sizes to its lanes as they are known (growing as a source's rows arrive,
  never a blank band), up to a cap (`track.maxPx` or a share of the view); past the cap its lanes scroll inside the
  track (the wheel over it scrolls lanes, a slim vertical rail like the rows', the header and labels fixed); the
  first row labelled. The fixed height of source-backed tracks (TODO, "Track heights … for a source's items") is the
  rule this revisits.
- **A bare number as a block name** (`:name 01`): GC1 refuses it with "quote it" (`seq.query.form.param`), since core's
  address reader turns `01` into 1 before Seq* sees it. Keeping the text would need a raw-token surface in core repr
  (allowed: repr's public API is open to this work); take it only if a real name set needs it.

- **A VCF allele's hover card names the record, not the allele** (VC's look, 2026-11-22, kind: look; for the hover
  contract's lane): hovering a ClinVar variant shows ClinVar's id, Track and Extent; it should name the allele
  (`pos:ref>alt`), its category and its INFO fields, as the design's matrix row "a variant allele" implies.

- **A click on one piece of a feature across a cyclic origin selects only that piece** (GB's look, 2026-11-22,
  kind: bug; likely shared with the GFF3 plasmid, the same renderer): Features › Genome formats › GenBank › Plasmid,
  bla (`complement(join(2312..2686,1..486))`): the hover names the whole item and lights both pieces, but a click on
  the piece at 2,312–2,686 selects only that stretch. A click selects the item's extent across the origin (both
  pieces), as the design's matrix row "an element: a transcript" says. Test on both plasmid views.

- **The VCF ClinVar story test fails under load** (2026-11-22: once in GB's run, once in RG's full quick run, each
  passing alone; the second failure makes it a fix lane's item, `TODO.md`'s watch rule; kind: test): find what is
  timing-sensitive in `vcf-stories.test.ts` (a real-turn wait, a timeout under load, a 30-read region on a real
  clock) and make it deterministic on a fake clock; run it 20 times under load.

## Fixed

- **"Try pUC19.gb" looks like it does nothing** (fixed by GB, merged 2026-11-22: the button loads the GenBank file) (the author, 2026-11-21, Examples › Services › Your own file:
  "pUC19.gb doesn't load anything"). By design today: GenBank is out of the slice's scope (plan, Out of scope), and
  `dropped-file.ts:71` refuses `.gb`/`.gbk` with "GenBank files are not read yet"; the story's "Try pUC19.gb" button
  (`drop-root.ts:30`, its sample in fact `sources.json`) exists to show that refusal. But the refusal is a small grey
  line above the buttons and the box stays empty, so it reads as broken. Drop the "Try" button for a format that isn't
  read (a Try button promises a load), keep the refusal for a dropped file of a format not read, shown as a notice in
  the box (the file's name, why, the formats that are); a story test that every Try button loads a view. GenBank
  itself comes into the slice (the author, 2026-11-21: "Add a GenBank lane"; plan, GB): after GB, "Try pUC19.gb"
  loads the recorded GenBank file and `.gb` leaves the refused list; until then the button goes.
- **The SVG export blank on a spec with an error** (fixed by FX6, merged 2026-11-22) (drift svg-on-error, decided as recommended by the author,
  2026-11-22: "go as recommended with the drifts"; kind: bug): `exportSvg` returns a blank picture when `checkSeqSpec`
  finds an error (`packages/seq/core/src/render/export.ts`, about 435 and 441, `if (hasErrors(issues)) return
  blank(true)`), while the canvas draws what resolves. Export what resolves, as the canvas draws it, the issues in the
  picture's `<desc>`; viewer-composition.md says so. Tests: Features › Builder › Check's twin asserts the SVG equal to
  the canvas's draw list (letters and an empty hits track) and the issue in `<desc>`; a spec with no errors exports
  unchanged.

- **A merged run of mixed categories drawn neutral** (fixed by FX6, merged 2026-11-22) (drift vcf-mixed-runs, decided as recommended by the author,
  2026-11-22; kind: look): `render/ranges.ts` (about 508, VC's `-2 - item` encoding) and `ranges-draw.ts` (about 79)
  paint a run whose items differ in fill in the neutral colour; ClinVar at TP53's fit is mostly grey. A mixed run takes
  the majority fill among its items (ties by the legend's order); the hover over it lists the categories under the
  pointer. Tests against the drawn lists: TP53's ClinVar band draws no neutral rect at the fit; a two-item tie picks
  the legend's first; a single-category run unchanged.


- **The navigator's readout covers the control in use** (fixed by FX5, merged 2026-11-22) (the author, 2026-11-21, BED cytobands: "the tooltip appears
  on top of currently active control. quite disruptive"): while the window is dragged, its readout
  ("chr17:36,190,896–62,395,127 of chr17") is drawn over the window and its grip. The readout never covers the window,
  its grips or the pointer: placed beside the window on the side with room, else above or below the navigator's
  strip. Test: during a drag the readout's box does not intersect the window's or the pointer's.

- **The square navigator's sweep has no feedback while dragging** (fixed by FX5, merged 2026-11-22) (the author, 2026-11-21, the square over blocks ×
  blocks): a drag across cells zooms both axes to the swept rectangle on release, but nothing is drawn during the
  drag. Draw the sweep's rectangle live (as the linear navigator's sweep does), snapped to the cells it will fit, and
  clear it on release or Escape. Test: a pointer down and move draws a sweep mark in the navigator's draw list at the
  swept cells; release fits both axes to them; Escape cancels with no change.

- **The palette control** (fixed by FX5, merged 2026-11-22) (the author, 2026-11-21, Structure link › Hover both ways, the 2OCJ track):
  - picking a palette does not always change the track (the residue letters keep their colours): the palette entry
    shows only on a track whose layer reads a categorical palette (a ranges or category theme), never on one that
    ignores it (a letters track coloured by chemistry);
  - it shows twice, as the header's palette icon and as "Color ▸" in the "…" menu: one place, the menu (the icon goes,
    or stays only where a host asks for it as a header part);
  - a submenu's "▸" sits right after its label: right-align it at the menu's edge, as menus do.
  - **a palette choice can break a track's colouring for good** (the author, Themes › Themes, the Variants track):
    picking a palette turns every lollipop grey, and picking "Categories" again does not bring the colours back. A
    palette change is a patch: undo restores it, and choosing the original entry restores the original colours exactly
    (the theme's own domain and mapping, not a palette applied over a field it cannot read). Test: change and change
    back, the draw list equal to the first; undo equal too.
  Tests: a letters track's menu has no Color entry; the Variants round trip above; a ranges track's Color change redraws with the new palette; one
  palette control a track; the arrow at the row's right edge (a jsdom layout check of the entry's parts).

- **The tree's leftmost vertical line is thinner** (fixed by FX5, merged 2026-11-22) (the author, 2026-11-21, Globins rows family): the root's vertical
  stroke at the tree column's left edge looks about half the width of the others, likely half of it clipped at the
  column's edge (a stroke centred on x = 0). Inset the tree by half a stroke (or snap strokes inside the column); a twin
  test that every vertical's stroke lies inside the tree column's box.

- **Three element classes exported from stories files show as stories** (fixed by FX5, merged 2026-11-22) (found by BS's guard, 2026-11-22, kind: look):
  `LegendPlacementsRoot`, `ThemesEditingRoot`, `SchemesGalleryRoot` are `export class` in `*.stories.ts`, so Storybook
  lists each as a story; their tests import them. Move each to a module of its own beside its story and drop it from
  `story-exports.test.ts`'s `Known` list.




- **The rail's place as a list, and a pinned place**: FX4, `railPlace` one place or a list; `pinned` at the
  frame's foot (story FamilyRailPinned).

- **The rows layout's horizontal rail**: FX3 (`ded361691`), `navigator.positionsRail` (`always`) and `navigator.railPlace`.
- **A ring's own hover draws no mark**: FX3 (`ded361691`), the hovered view draws its own position box.
- **An assembly override naming no sequence**: FX3 (`ded361691`), `seq.resolve.assembly.override`.
- **Block rulers unnumbered at the fit**: not reproduced in a real browser (FX3: numbered from about 150 px a block;
  unnumbered only at about 65 px, by the label rule).

- **GC2's leftovers from my second look** (2026-11-21, after its fix round; fixed by GC3, `61e71891b`):
  - a double-click on the 25-block row zooms to bases (161,976,790) instead of fitting the block;
  - B in A's order: the reversed middle block of chr2 shows no label at all, though "chr2 ←" fits;
  - the hover card calls a chromosome a "Chain" (`Chain chr3`, `Chain chr2 inv`): the row's name for a block's system;
  - the ring's group labels ("peach", "grape") still straddle the group arc and sit among the block labels.

- **A double-click on a long sequence jumping to bases**: FX2 (`e678d2dbc`), one zoom step about the pointer.

- **The hover card far from the pointer in a tall track**: FX2 (`e678d2dbc`), anchored by the pointer's lane.

- **The SVG not repainted on a cleared selection**: FX2 (`e678d2dbc`), the SVG twin follows the selection.

- **Clearing a rows sort**: FX2 (`e678d2dbc`), a third click runs `rows.unsort`; `viewer.rowSort` and a host "Clear sort".

- **A hidden track unnamed in the menu**: FX2 (`e678d2dbc`), listed by its layer's name.

- **The "…" actions button with no actions**: not reproduced on packages-seq (FX2, the orchestrator's look: the Genes
  menu holds its built-ins); the author, 2026-11-21: "might have been on older version and it got fixed already".

- **The navigator's readout; letters larger than their labels; a legend hover sticking**: FX1 (`2afb32e91`): a
  padded one-line readout whose number never truncates; letters capped by measured cap heights (`letters.maxPx`,
  `capHeight`); the legend's hover on its container.
