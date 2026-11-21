# Fixes: the fifth slice's small findings

Small findings that need no design collect here for the slice's fix points, as the fourth slice's did (the author,
2026-11-21: "lets have a 'fix phase' for slice 5 too where we accumulate these smaller things?"). An FX lane takes the
open items at each fix point (plan, Fix points); a finding that needs a design becomes a drift note instead.

## Open

- **The hover card: placement, wrapping, and the density row's missing highlight** (the author, 2026-11-21, Files ›
  UniProt entry: "hovering the bottom 'regions' … is missing a highlight", "the hover placement is quite
  inconsistent"):
  - the Regions track's bottom row (its density row of overlapping regions) shows a card but draws no hover mark: the
    density row lights the hovered stretch (or the items under it);
  - the card jumps between below and above the pointer and covers the rows being read: one placement rule (beside
    the pointer on the side with room, never over the hovered row, the same side while the pointer stays in a track);
  - the card's width changes from one hover to the next and wraps values ("P04637 100–" / "370"): a value never
    wraps inside a range; the card sizes to its longest row up to a cap, then truncates the name, not the value;
  - unnamed features are labelled by a bare row number ("152", "153"): fall back to the feature's type and range.
  Tests: the density row's hover draws a mark; the card's box never intersects the hovered row; a range value on
  one line.

- **The navigator's readout covers the control in use** (the author, 2026-11-21, BED cytobands: "the tooltip appears
  on top of currently active control. quite disruptive"): while the window is dragged, its readout
  ("chr17:36,190,896–62,395,127 of chr17") is drawn over the window and its grip. The readout never covers the window,
  its grips or the pointer: placed beside the window on the side with room, else above or below the navigator's
  strip. Test: during a drag the readout's box does not intersect the window's or the pointer's.

- **The navigator's window and zoom strip look inconsistent and rough** (the author, 2026-11-21, a ring that switches:
  "i don't understand why chr6 gets a 'box' and the other 'sub-section'"; "the 'subsection visual' is really ugly,
  with different line thicknesses"). Why it differs: a window wider than `navigator.zoomStrip.minPx` (48 px) is drawn
  as the box with grips (chr6, about 85 px), a narrower one (chr13) gets the zoom strip under the navigator, with its
  funnel lines. Needs a short look design before a fix lane: one rule a reader can see (the strip when the window is
  too small to grab, or always for a section of one block), one stroke weight for the window, the funnel and the
  strip, the strip's fill and border from the same tokens as the window's (ux-foundations: the loudest mark at rest),
  the funnel thin and muted or dropped; its picture in navigators.md, then the fix.

- **The square navigator's sweep has no feedback while dragging** (the author, 2026-11-21, the square over blocks ×
  blocks): a drag across cells zooms both axes to the swept rectangle on release, but nothing is drawn during the
  drag. Draw the sweep's rectangle live (as the linear navigator's sweep does), snapped to the cells it will fit, and
  clear it on release or Escape. Test: a pointer down and move draws a sweep mark in the navigator's draw list at the
  swept cells; release fits both axes to them; Escape cancels with no change.

- **The palette control** (the author, 2026-11-21, Structure link › Hover both ways, the 2OCJ track):
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

- **The vertical overview cut at the right on a wide page** (the author, 2026-11-21, Globins rows family-rail-pinned):
  at full width the story's canvas runs past the page's right edge, so the vertical overview is clipped (a narrow
  page, where canvas and SVG stack differently, shows it whole). The story's layout gives the canvas its measured width
  at every page width, the overview inside it; a story-layout test at a wide and a narrow width that the overview's box
  lies inside the canvas and the canvas inside the page.

- **The tree's leftmost vertical line is thinner** (the author, 2026-11-21, Globins rows family): the root's vertical
  stroke at the tree column's left edge looks about half the width of the others, likely half of it clipped at the
  column's edge (a stroke centred on x = 0). Inset the tree by half a stroke (or snap strokes inside the column); a twin
  test that every vertical's stroke lies inside the tree column's box.

- **The 3,000-contig story could be faster** (the author, 2026-11-21, on packages-seq's build: "better perf ... although
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

## Fixed

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
