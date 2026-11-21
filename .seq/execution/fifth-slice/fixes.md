# Fixes: the fifth slice's small findings

Small findings that need no design collect here for the slice's fix points, as the fourth slice's did (the author,
2026-11-21: "lets have a 'fix phase' for slice 5 too where we accumulate these smaller things?"). An FX lane takes the
open items at each fix point (plan, Fix points); a finding that needs a design becomes a drift note instead.

## Open

- **A hidden navigator track has no name in the hidden-tracks menu** (the author, 2026-11-21, Globins rows, family):
  the entry is an eye with an empty label; every hidden track is listed by a name (its label, its layer's or
  navigator's title, else its kind). Sent to FX2 as its item 6.

- **The rows layout has no horizontal rail** (the author, 2026-11-21, Globins rows, family: "it's weird this view is
  lacking the 'horizontal rail'. probably should be at the 'navigator' or at the bottom. in general, configurable"):
  the tracks layout has its default rail under the ruler; the rows layout has the vertical overview and no horizontal
  rail. Needs a short design before a fix lane: the rail's place as a setting (`rail.place`: under the ruler, with the
  navigator, at the bottom, none), the same for every layout that has a horizontal axis, its default per layout.

- **The SVG not repainted when the selection is cleared** (the author, 2026-11-21, gff3--tp-53): the SVG twin keeps
  the cleared selection's line drawn while the canvas drops it. Every story's SVG redraws on every change the canvas
  does (selection set, cleared, hover where shown); test it in a twin: clear, then the SVG holds no selection mark.

- **The "…" actions button with no actions** (the author, 2026-11-21, on GG's GFF3 story: "could we stop showing the
  ... actions button if there are no actions registered?"): a track header shows "…" (tooltip "Actions") when its menu
  would be empty. The button is drawn only when at least one entry applies to that track (built-ins that apply, the
  host's, an extension's), decided per track and per state; the fourth slice fixed one case of this (the author: "some
  views show the 3 dot menu icon, but no menu is there"), so find why it is back and test every header kind.
- **The Transcripts track's spacing** (my look at GG's fix round; the author, 2026-11-21: "the spacing around the
  transcripts is weird"): about 100 px of blank lanes above the first transcript in gff3--tp-53 and as many below the last (the track
  taller than its lanes, the lanes centred in it), the first row
  unlabelled, and lanes packed unevenly around the labels; a tall track's lanes start at its top, evenly spaced.
- **The hover card far from the pointer** (the author, 2026-11-21, gff3--tp-53: "the hover card appears at a weird
  place"): hovering TP53-227 deep in the tall Transcripts track puts the card at the track's top, over the Genes track.
  The card anchors by the hovered item's row (or the pointer), not the track's top, in every tall track.

- **A double-click on a long sequence jumps to bases** (my look at BE's cytobands, 2026-11-21): chr17 at the fit (83
  Mb) to about 20 bases in one double-click; GC2's row did the same. A double-click should zoom one step (or fit the
  item under the pointer), not to the closest zoom. GC3 takes it on blocks; the plain section's case is the fix lane's.

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

- **The navigator's readout; letters larger than their labels; a legend hover sticking**: FX1 (`2afb32e91`): a
  padded one-line readout whose number never truncates; letters capped by measured cap heights (`letters.maxPx`,
  `capHeight`); the legend's hover on its container.
