# Lane RX: the rows' look: a row column's whole value only under the pointer, the tree's lines and wedges

**Goal** (the author, 2026-11-19, on Globins › Editing › Family: "the full labels should only display when you
actually hover on the species column ... otherwise it obscures the data"; [grids-and-rows.md](../../../design/grids-and-rows.md),
Row columns, "A cut value shows whole only under the pointer"): a row column's cell cut short is drawn whole in the
`tip` pass only while the pointer is on that cell; a hover on the row's data or label tints the row and leaves its
columns cut; the cursor on the row's header draws them whole. The core `_spec`'s "Emphasis along y" already says so
(the orchestrator changed it); the code in `core/src/render/rows.ts` (the tip pass, its hover test) and
`rows-columns.ts` follows. Tests (fail on the parent): a hover on a data cell of a row with a cut species draws no tip;
a hover on the species cell draws it whole to the view's edge; the header cursor draws them; canvas and SVG alike.
Look: Globins › Editing › Family and Rows › Family, light and dark.

**And the tree's look** (the author, 2026-11-19: "the (collapsed) clade visuals are quite ugly. and the different line
thickness too"; "Light wedge"; grids-and-rows.md, "The tree's look"): `branchRect` (`core/src/render/rows.ts:2002`)
centres a 1 px rect on fractional coordinates, so branches smear over one or two device pixels: snap them to device
pixels, one width (the renderer's device pixel ratio; reuse a snapping helper if the render code has one). The wedge
(`drawWedges`, `rows.ts:2020`, today a solid triangle in the branch ink to the tree's full width): a light tint of the
branch ink, outlined at the branches' width, from the clade's node to its deepest leaf's x, its far side about half the
row high with softened corners (the shape rule, ux-foundations §5), its fill a step darker on hover; its hit (the
toggle's 6 px rule) unchanged. Tests: every branch rect's edges on device pixels at dpr 1 and 2; a wedge's far x the
deepest leaf's; the SVG alike. Look: Globins › Rows › Folded, light and dark, dpr 2.

**And the hidden columns' chip** (the author, 2026-11-19: "the 2 columns hidden is very misleading, because it doesn't
show what columns where are when you hover it"; grids-and-rows.md, Row columns, Narrow): the chip names the hidden row
columns ("Species, Identity hidden"; "Species +1 hidden" when short of room), sits in the title band on the side those
columns live (never over the alignment's data), and its hover card lists each hidden column and why ("no room at this
width"); announced the same. Its placement (the author: "also the positioning of the label is wrong": today it
hangs at the alignment's top-right corner, over the data and past the numbers' gutter, since the columns it stands for
have no room): inside the title band's own area on that side (the right columns' band end, aligned to the numbers'
gutter; the left's beside the labels), never overlapping the data box, the rows rail or another title; when there is no
room for it on that side, it moves into the label column's title band. Tests: at a width that hides both, the chip's text and its card's rows; the SVG's text
alike. Look: Globins › Rows › Family at a narrow width.

**And the rows rail's idle state** (the author, 2026-11-19, Globins › Rows › Scale: "resizing vertical rail all the way
makes it impossible to make it smaller again"): the rows rail goes idle (no grips, no target) whenever every row is
shown, against [navigators.md](../../../design/navigators.md)'s idle rule (idle only when the view shows the whole
axis *and* zoom-in is at its limit); the core `_spec`'s rows-rail line now says so. Make the code follow it, for the
rows rail and check the ruler rail too (whole at the fit view but able to zoom in: grips shown, a grip draggable
inward). Tests (fail on the parent): drag the rows rail's grip to the band's end, then drag it back: the window
shrinks; a rail whole at its zoom-in limit is idle. Look: Rows › Scale, light and dark.

**Not in scope**: anything else. CF works in core (store, ruler, templates, accessibility, rows-model's measures), MV
in core's renderer tables and the viewer's composition: keep to the rows layer's tip.

**Model**: strong (opus). **Base**: the commit the orchestrator gives. Storybook port 6025.

**Checks**: as [hc-hover-card.md](hc-hover-card.md)'s, without size.

**Headers**, **hot files** and **stop when**: as every lane's (see [qt-quick-lane.md](qt-quick-lane.md)).

**Report**: commits; checks; tests failing on the parent; what you saw; a `.v6/changelogs/` entry's text.
