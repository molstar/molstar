## Seq* (next): version 1's close-out

- Seq\* viewer: a legend's eye hides or shows its layer, one undo step; a hidden layer's row stays in place, dimmed,
  with the eye that shows it again, also for a layer hidden from the start and a bridge's layer (`mountLegend` option
  `layers: false` drops the eye). The track header's menu and the eye name a layer the same way.
- Seq\* core: a bridge's ribbons have an `orientation` field (same or inverted; a data column of that name wins),
  keyed in a legend with the width's meaning; a continuous legend over an open domain says "in view".
- Seq\* genome formats: a transcript tagged `Ensembl_canonical` or `MANE_Select` in GFF3 or GTF is drawn first in its
  gene's stack with a bold label, and its card says so; a matrix's samples (a VCF's genotypes, a methylation matrix)
  are named beside its rows; a genotype card gives the call as the file wrote it (`0|1`, `1/2`, `1`) and what it is; a
  missing value has a colour of its own.
- Seq\* genome ring: a right-click, Shift+F10 or the Menu key on a stem opens its menu; partners grouped by chromosome,
  each with its pieces named by band and strand, the rest under "n more", links within the stem's own chromosome
  last; the stem card says what it counts, a pair's card names both ends; the hole says "Loading links…" or the
  failure.
- Seq\* viewer: a host's header entry can be shown pressed (`pressed`, `aria-pressed`), as the structure link's
  "colored by this track"; the "n filtered" chip's card names its filters ("Species ≠ Human") and says "click to
  clear"; a section the view cannot place says its issue in one row.
- Seq\* viewer: the square under the ring says its gestures; a double-click or 0 shows the whole genome, Back the
  window before; a Hi-C cell's card reads in kb or Mb; `.` and `,` step to the next and previous feature on a track
  (the ring's partners on a ring); the zoom-in hatch names the gesture; the square shows "Loading…" or the failure.
- Seq\*: a lanes band's back and forward return to the same item after its lanes are laid again; the theme picker's
  BLOSUM62 text says what the scheme does; a pairs layer's data filter that a source's bins do not apply is said
  (`seq.render.pairs.filter`, info).
- Seq\*: a drag over a gene's nucleotides snaps to whole codons, a codon split by an intron to its pieces; the keyboard
  cursor crosses the gaps between blocks, along a section and down a grid's y; every built-in layer type declares its
  fallback; a spec of a newer minor version with an unknown kind gives a warning naming both versions, not an error.
- Seq\*: a bookmark on a streamed node (a transcript read by region) survives its region's eviction; jumping to it says
  "Waiting for T2", asks the source where the node is and runs when it arrives; a key the source lacks is said once.
- Seq\* viewer: a drag-select's work is bounded by the screen: the selection is written once a frame and resolved on
  release, nothing is hovered while it runs, and the projections are indexed; a sweep across 12 chromosomes of a
  3,000-contig genome drops from about 300 ms to about 1 ms a frame.
- Seq\*: splitting many copies of a chain allocates about a quarter less: a ranges layer's layout object is built in one
  piece, its scratch arrays made on first use.
