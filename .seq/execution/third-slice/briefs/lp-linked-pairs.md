# Lane LP: the `linkedPairs` behavior

**Goal** (the author, 2026-11-19: "Behavior, on by default"; [viewer-models.md](../../../design/viewer-models.md),
Behaviors, `linkedPairs`): on a hovered or selected range element, the pairs with an end in its extent are lit through
the highlight channel on every pairs layer shown (arcs, chords, dots), cleared when the hover or selection goes. Built
on the public API only (the query language, the highlight channel, the hover and selection models): if a piece it
needs is not public, stop and report it (a finding for the extension path, not something to reach into). Listed in
`DefaultSeqViewerSpec()`; params `on` (`hover`, `selection`, both by default) and the pair layers it reads (all). The
pairs layer's own rule (a hovered position lights its partners) unchanged.

**Tests** (fail without it): hover "cyclin box" on Complex › Contacts lights exactly the contacts with an end in
37–134 of cyclin A (both copies merged); a selected range lights them in the selection's pass; a viewer composed
without the behavior lights none (the opt-out); the tRNA's arcs and its dot plot alike; hover allocates nothing a frame
past its first (the hover loop's rule).

**Two consumers**: the canvas and the SVG export (an export of a selection draws its lit pairs).

**Look**: Complex › Contacts and tRNA in light and dark.

**Not in scope**: anything else. Starts after MV merges (it edits the default spec MV restructures).

**Model**: strong (opus). **Base**: the commit the orchestrator gives. Storybook port 6024.

**Checks**: as [hc-hover-card.md](hc-hover-card.md)'s.

**Headers**, **hot files** and **stop when** (a public API change past the behavior and its params): as every lane's
(see [qt-quick-lane.md](qt-quick-lane.md)).

**Report**: commits; checks; tests failing without it; what you saw; what the public API lacked, if anything; a
`.v6/changelogs/` entry's text.

**Round 2** (the author, 2026-11-19, after LP stopped on the two gaps; viewer-models.md, "What it needed"; and
locations-and-queries.md, "Pairs answer `items`"): build both, then the behavior. (a) A pairs annotation answers the
query model (`StaticModel.annotation`, the overlay's query model), a pair one item with a row per end, `QueryItem.system`
optional per row. (b) `EmphasisHighlight.elements`; `_updateHighlights` and `_svgEmphasis` pass a channel's elements;
the pairs layer draws the highlight pass for listed keys in the channel's colour (arcs, chords, dots, SVG). (c) The
behavior: hover and selection write to its own channel (the highlight pass, never the selection's), cleared when they
go; the hover's query cached per hovered element key. Tests also for (a) and (b) on their own. These public API
additions are decided; anything past them still stops you.
