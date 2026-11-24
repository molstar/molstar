# Lane CE: a pair's end on split copies

**Goal** (the author, 2026-11-19, on Globins › Hemoglobin › Copies Split: "all the contacts go to the 1st heme copy …
we need ability to address the individual copies too?"; "Caret, not first copy"; [grids-and-rows.md](../../../design/grids-and-rows.md),
the pairs layer's scope paragraph): a pair whose end is on a system shown as several split copies (`copies: each`) and
whose data names no copy is not placed on any copy (today: on the first block of that system, `core/src/render/pairs.ts`
with `ctx.place`): the known end draws a caret pointing at the copies, as an end not shown does; one info issue
`seq.render.pairs.copy` a layout naming the system and the count; merged copies (×n) unchanged. Tests (fail on the
parent): the split hemoglobin: no arc ends on a heme block, a caret at each contact's chain end, the issue once; merged:
arcs as before; canvas and SVG alike; the hover of a caret names the contact. Look: Hemoglobin and Copies Split, light
and dark.

**Not in scope**: units and per-copy pairs (molecules.md, the fourth slice). CF and RX work in core beside you
(CF: store, ruler, templates, accessibility, rows-model; RX: rows.ts, rows-columns.ts): keep to `pairs.ts` and its tests.

**Model**: fast (sonnet). **Base**: the commit the orchestrator gives.

**Checks**: as [hc-hover-card.md](hc-hover-card.md)'s, without size.

**Headers**, **hot files** and **stop when**: as every lane's (see [qt-quick-lane.md](qt-quick-lane.md)).

**Report**: commits; checks; tests failing on the parent; what you saw; a `.v6/changelogs/` entry's text.
