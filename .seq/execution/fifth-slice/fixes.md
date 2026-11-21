# Fixes: the fifth slice's small findings

Small findings that need no design collect here for the slice's fix points, as the fourth slice's did (the author,
2026-11-21: "lets have a 'fix phase' for slice 5 too where we accumulate these smaller things?"). An FX lane takes the
open items at each fix point (plan, Fix points); a finding that needs a design becomes a drift note instead.

## Open

- **The rail's place as a list, and a pinned place** (the author, 2026-11-21; designed in navigators.md, "A list, and
  a pinned place"): `railPlace` takes one place or a list; `pinned` keeps a rail at the bottom edge of the visible
  frame at any scroll. For FX4, on FX3's `railRows`.

- **The 3,000-contig story could be faster** (the author, 2026-11-21, on packages-seq's build: "better perf ... although
  still could be faster"): my measure, 40 hover steps along the row each drawn within two frames (about 33 ms, no
  drops). Take it at M5: profile hover and pan on three-thousand-contigs (ring and row), find the per-frame cost, and
  set a bound in the perf tier. **A drag-select across many chromosomes is slow** (the author, 2026-11-21: chr11 to
  chr22 selected on the row, mirrored on the ring): profile the drag's frames (the selection's blocks and their
  emphasis on both views), bound a drag across 12 blocks to the frame budget at 0 bytes warm.

- **The Transcripts track's spacing** (FX2 `e678d2dbc`: lanes start at a tall band's top; left open: the
  first row unlabelled and the blank below the lanes, a source-backed track's fixed height) (my look at GG's fix round; the author, 2026-11-21: "the spacing around the
  transcripts is weird"): about 100 px of blank lanes above the first transcript in gff3--tp-53 and as many below the last (the track
  taller than its lanes, the lanes centred in it), the first row
  unlabelled, and lanes packed unevenly around the labels; a tall track's lanes start at its top, evenly spaced.
- **A bare number as a block name** (`:name 01`): GC1 refuses it with "quote it" (`seq.query.form.param`), since core's
  address reader turns `01` into 1 before Seq* sees it. Keeping the text would need a raw-token surface in core repr
  (allowed: repr's public API is open to this work); take it only if a real name set needs it.

## Fixed

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
