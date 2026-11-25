# Lane PZ: a plot's navigation and its chain boundaries

**Goal** (the author, 2026-11-19, using Examples › Complex › PAE and contacts: "still very unintuitive to scroll the PAE
plot and we could end up in very weird views … double click zooms to a spot, but double click again shifts things …
wheel scroll ends up scrolling the page and resizing the view … scrolling back up is even worse"; and on Over The
Complex: "the border between the regions here is very weird"). Build what
[grids-and-rows.md](../../../design/grids-and-rows.md) now says (decided, read the grid section's first bullets and
the matrix row "an element: a cell"):

1. **A plot's box is fixed**: set from the width and `aspect` (capped as today) and never changed by a zoom or a pan;
   cells are rectangles when x and y show different spans. Test: zoom one axis (the PAE's y following cyclin A's
   section), the plot's height and every track's place unchanged.
2. **Double-click on a grid cell zooms a step** about it on both axes, keeping their ratio; again further; Alt+Left
   back, `0` fits. The tracks' and rows' double-click unchanged. Tests: two double-clicks zoom twice about their cells
   (the second never re-centres to a residue window), back undoes both.
3. **The wheel over a plot is the plot's**: plain wheel pans (dx x, dy y), Cmd/Ctrl+wheel and pinch zoom both axes
   about the pointer; at an axis's end it stops and the page does not scroll (`preventDefault` on every wheel over the
   plot's box, in the viewer's input; the page's wheel outside it untouched). Tests: a wheel at y's end over the plot
   is consumed (no page scroll), one outside the plot is not.
4. **Chain boundaries** on a grid over blocks: the blocks' gap with no line in it, cells clipped exactly at each
   block's edge (today an off-centre `seq-break` hairline runs across the gap, and a one-pixel sliver of colour shows at
   blocks' left edges); a hairline only where blocks touch with no gap. Test: the draw lists at the fit of Over The
   Complex: no break line in a gap, no cell rect past its block's edge; canvas and SVG alike.

**Drawing against the data** and **the look**: Complex › PAE and contacts (Pae, Over The Complex, Contacts) and tRNA's
dot plot, in light and dark, at the fit, zoomed by double-click twice, and wheeled to an axis's end.

**Not in scope**: anything else. Starts after CF and MV merge (it touches core's grid drawing and the viewer's input).

**Model**: strong (opus). **Base**: `packages-seq` at the commit the orchestrator gives. Storybook port 6023.

**Checks**: as [hc-hover-card.md](hc-hover-card.md)'s.

**Headers**, **hot files** and **stop when** (a public API change): as every lane's (see
[qt-quick-lane.md](qt-quick-lane.md)).

**Report**: commits; checks; each item's test failing on its parent; what you saw; a `.v6/changelogs/` entry's text.
