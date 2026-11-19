# Lane RX: a row column's whole value only under the pointer

**Goal** (the author, 2026-11-19, on Globins › Editing › Family: "the full labels should only display when you
actually hover on the species column ... otherwise it obscures the data"; [grids-and-rows.md](../../../design/grids-and-rows.md),
Row columns, "A cut value shows whole only under the pointer"): a row column's cell cut short is drawn whole in the
`tip` pass only while the pointer is on that cell; a hover on the row's data or label tints the row and leaves its
columns cut; the cursor on the row's header draws them whole. The core `_spec`'s "Emphasis along y" already says so
(the orchestrator changed it); the code in `core/src/render/rows.ts` (the tip pass, its hover test) and
`rows-columns.ts` follows. Tests (fail on the parent): a hover on a data cell of a row with a cut species draws no tip;
a hover on the species cell draws it whole to the view's edge; the header cursor draws them; canvas and SVG alike.
Look: Globins › Editing › Family and Rows › Family, light and dark.

**Not in scope**: anything else. CF works in core (store, ruler, templates, accessibility, rows-model's measures), MV
in core's renderer tables and the viewer's composition: keep to the rows layer's tip.

**Model**: fast (sonnet). **Base**: the commit the orchestrator gives. Storybook port 6025.

**Checks**: as [hc-hover-card.md](hc-hover-card.md)'s, without size.

**Headers**, **hot files** and **stop when**: as every lane's (see [qt-quick-lane.md](qt-quick-lane.md)).

**Report**: commits; checks; tests failing on the parent; what you saw; a `.v6/changelogs/` entry's text.
