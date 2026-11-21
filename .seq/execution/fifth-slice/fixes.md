# Fixes: the fifth slice's small findings

Small findings that need no design collect here for the slice's fix points, as the fourth slice's did (the author,
2026-11-21: "lets have a 'fix phase' for slice 5 too where we accumulate these smaller things?"). An FX lane takes the
open items at each fix point (plan, Fix points); a finding that needs a design becomes a drift note instead.

## Open

- **The 3,000-contig story could be faster** (the author, 2026-11-21, on packages-seq's build: "better perf ... although
  still could be faster"): my measure, 40 hover steps along the row each drawn within two frames (about 33 ms, no
  drops). Take it at M5: profile hover and pan on three-thousand-contigs (ring and row), find the per-frame cost, and
  set a bound in the perf tier. **A drag-select across many chromosomes is slow** (the author, 2026-11-21: chr11 to
  chr22 selected on the row, mirrored on the ring): profile the drag's frames (the selection's blocks and their
  emphasis on both views), bound a drag across 12 blocks to the frame budget at 0 bytes warm.

- **A ring's own hover draws no mark** (the author, 2026-11-21, three-thousand-contigs): hovering the ring shows the
  card but no position box on the ring, while hovering the row draws it there; sent to FX2 as its item 9; not started.

- **Block rulers unnumbered at the fit** (the author, 2026-11-21, the-hosts-complex): chain blocks about 180 px wide
  show tick marks and no numbers at the fit (numbered before GC2 and GC3); sent to FX2 as its item 8; not reproduced headless (it draws "100"/"200"), likely the
  browser's text metrics: needs a look in the browser.
- **The rows layout has no horizontal rail** (the author, 2026-11-21, Globins rows, family: "it's weird this view is
  lacking the 'horizontal rail'. probably should be at the 'navigator' or at the bottom. in general, configurable"):
  the tracks layout has its default rail under the ruler; the rows layout has the vertical overview and no horizontal
  rail. Needs a short design before a fix lane: the rail's place as a setting (`rail.place`: under the ruler, with the
  navigator, at the bottom, none), the same for every layout that has a horizontal axis, its default per layout.

- **The Transcripts track's spacing** (FX2 `e678d2dbc`: lanes start at a tall band's top; left open: the
  first row unlabelled and the blank below the lanes, a source-backed track's fixed height) (my look at GG's fix round; the author, 2026-11-21: "the spacing around the
  transcripts is weird"): about 100 px of blank lanes above the first transcript in gff3--tp-53 and as many below the last (the track
  taller than its lanes, the lanes centred in it), the first row
  unlabelled, and lanes packed unevenly around the labels; a tall track's lanes start at its top, evenly spaced.
- **GC2's leftovers from my second look** (2026-11-21, after its fix round; fixed by GC3, `61e71891b`):
  - a double-click on the 25-block row zooms to bases (161,976,790) instead of fitting the block;
  - B in A's order: the reversed middle block of chr2 shows no label at all, though "chr2 ←" fits;
  - the hover card calls a chromosome a "Chain" (`Chain chr3`, `Chain chr2 inv`): the row's name for a block's system;
  - the ring's group labels ("peach", "grape") still straddle the group arc and sit among the block labels.

- **An assembly `sequence` override naming no sequence is dropped without a word** (GC1, found in its fix round): it
  should say `seq.resolve.assembly.*` like an alias that names nothing does.
- **A bare number as a block name** (`:name 01`): GC1 refuses it with "quote it" (`seq.query.form.param`), since core's
  address reader turns `01` into 1 before Seq* sees it. Keeping the text would need a raw-token surface in core repr
  (allowed: repr's public API is open to this work); take it only if a real name set needs it.

## Fixed

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
