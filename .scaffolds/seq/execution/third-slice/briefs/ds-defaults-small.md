# Lane DS: the legend, the rail, ring clusters, contrast (fast)

Takes over lane LR (`lr-legend-row.md`, read it: the channels row `reserve | collapse` is in here).

**Goal** (the author, 2026-11-19, the defaults audit `drafts/defaults-audit.md` findings 6–11, all four smaller
proposals chosen; [themes.md](../../../design/themes.md): "The legend's groups one a line", "Muted text holds 4.5:1",
"Highlight channels in the legend part" with `channels`; [navigators.md](../../../design/navigators.md): the rail 1 px
and muted at rest; [layouts.md](../../../design/layouts.md), "Defaults decided from the audit": ring labels grouped):
1. The legend part: groups one a line (title muted, 600, 11 px; 12 px between); the channels row `channels: reserve`
   (default) or `collapse`.
2. The rail at rest 1 px and muted, the accent only while hovered or dragged (its grab zone unchanged).
3. Rings: more than six labels within 12° drawn as one grouped label ("13 sites"), opened on hover or zoom.
4. Contrast: measure ruler tick labels and notes against their well in both themes; below 4.5:1, the text colour at
   75%; a test that asserts the ratio from the theme tokens. The logo's letters only above about 9 px glyph height, a
   flat scheme-coloured block below.

**Tests** (fail on the parent): each; the legend part's height unchanged by a channel under `reserve`; headless model
and jsdom DOM. Perf: 0 bytes a frame. Look: Antibody › Schemes, tRNA › Arcs, Plasmid › Ring, Globins › Rows › Family,
light and dark at dpr 2 (headless Chrome on a free port, stop it after).

**Not in scope**: track heights, item labels, label column, colours of layers (lane DB). DB works beside you in
config.ts and tracks.ts: keep your config changes to your keys.

**Model**: fast (sonnet). **Base**: the commit the orchestrator gives. **Checks**: as
[hc-hover-card.md](hc-hover-card.md)'s, without size. **Headers**, **hot files**, **stop when**: as every lane's.

**Report**: commits; checks; tests failing on the parent; what you saw; a changelog entry's text.
