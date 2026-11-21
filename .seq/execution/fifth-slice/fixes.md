# Fixes: the fifth slice's small findings

Small findings that need no design collect here for the slice's fix points, as the fourth slice's did (the author,
2026-11-21: "lets have a 'fix phase' for slice 5 too where we accumulate these smaller things?"). An FX lane takes the
open items at each fix point (plan, Fix points); a finding that needs a design becomes a drift note instead.

## Open

- **GC2's leftovers from my second look** (2026-11-21, after its fix round):
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

- **The navigator's readout; letters larger than their labels; a legend hover sticking**: FX1 (`2afb32e91`): a
  padded one-line readout whose number never truncates; letters capped by measured cap heights (`letters.maxPx`,
  `capHeight`); the legend's hover on its container.
