# Lane DB: tracks fit their content; labels inside or beside (strong)

Takes over lane LL (`ll-labels-and-extent.md`, read it: its two parts are in here).

**Goal** (the author, 2026-11-19, the defaults audit `drafts/defaults-audit.md` findings 1–5, "Accept as drafted",
"Yes, per layer"; [layouts.md](../../../design/layouts.md): "A track's height by default fits its content, in steps"
with its numbers, the section's `heights` param and the host's `trackSizing` port, "Labels inside if they fit, else
beside", "Defaults decided from the audit"; [themes.md](../../../design/themes.md): "A layer with no color takes the
next categorical color"):
1. Content heights: lanes 14 px (12 mark, 2 gap), 4 px padding, one lane 22, cap 5 lanes (78), packed over the whole
   data (never the view), marks centred; tracks of one kind in a section share the tallest; `heights: content |
   uniform | fixed` (validated; `ext:` names through the registry); `track.height` overrides; the `trackSizing` port
   (a factory on the composition, given each section's tracks and content heights and the named answer).
2. Labels: inside if the whole text fits, else beside (ring: outside with a leader, as now; linear: a label lane above
   the track at the item's start, a short leader when moved); never cut mid-word; one drop priority and text style; the
   label lane counts in the content height.
3. The origin-crossing extent bug (LL part 2): bla (AmpR) on pUC19c says "Extent pUC19c 1–2,686".
4. The label column measured (80–160 px, ellipsis); `sectionGap` 20; text `lineGap` 10, `bandGap` 3.
5. A colourless layer the next categorical colour by layer order.
Records (`SEQ_RECORDS=update`) and stories change: list every changed record and say why.

**Tests** (fail on the parent): each rule above, canvas and SVG alike; a pan never changes a track's height; the port
overriding; the label lane's hits; the extent's words. Perf: 0 bytes a frame kept (`test:perf`); layout cost of the
packing measured on p53 tracks and the chromosome. Look: Antibody › Schemes, p53 › Tracks and Gene by region, Plasmid ›
Ring, Complex › Confidence, light and dark at dpr 2 (headless Chrome on a free port, stop it after).

**Not in scope**: the legend, the rail, ring clusters, contrast (lane DS).

**Model**: strong (opus). **Base**: the commit the orchestrator gives. **Checks**: as
[hc-hover-card.md](hc-hover-card.md)'s, with size. **Headers**, **hot files**, **stop when**: as every lane's
([qt-quick-lane.md](qt-quick-lane.md)).

**Report**: commits; checks; tests failing on the parent; changed records; what you saw; a changelog entry's text.
