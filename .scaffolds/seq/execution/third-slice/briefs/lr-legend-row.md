# Lane LR: the legend's channels row keeps its place (fast)

**Goal** (the author, 2026-11-19, p53 › Protein on its gene: hovering a legend entry adds a "Highlight: legend Clear"
row under the legend and the page below moves; "configurable? reserve by default"; [themes.md](../../../design/themes.md),
"Highlight channels in the legend part", and its Configuration row; the Design page `themes.mdx`): the legend part's
`channels` param, `reserve` (default: the row's height held whenever the part can show channels, empty with none) or
`collapse` (today's behavior). Validate the param.

**Tests** (fail on the parent): the legend part's height equal with and without an active channel under `reserve`;
`collapse` grows as today; the row's content and clear unchanged; headless (the part's model) and the mounted DOM
(jsdom) alike. Look: p53 › Protein on its gene, hovering a legend entry: nothing below moves, light and dark.

**Not in scope**: anything else in the legend. Keep to the legend part and its stories.

**Model**: fast (sonnet). **Base**: the commit the orchestrator gives.

**Checks**: as [hc-hover-card.md](hc-hover-card.md)'s, without size. **Headers**, **hot files**, **stop when**: as
every lane's ([qt-quick-lane.md](qt-quick-lane.md)).

**Report**: commits; checks; tests failing on the parent; what you saw; a `.v6/changelogs/` entry's text.
