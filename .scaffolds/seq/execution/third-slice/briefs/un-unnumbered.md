# Lane UN: where a numbering does not reach (fast)

**Goal** (the author, 2026-11-19, Antibody › light chain: the Kabat and IMGT rulers blank over the constant domain;
"customizable … 'not numbered' as default"; [axis.md](../../../design/axis.md), "Where a numbering does not reach", and
its Configuration row; the Design page `axis.mdx`): the ruler layer's `unnumbered` param. `mark` (default): past the
numbering's coverage the ruler's line muted and dashed, one muted "not numbered" note at the start of each unnumbered
run, no ticks; the hover card and live region there say the sequence position and "not numbered by <numbering>".
`fallback` with `fallbackNumbering` (a numbering name, ordinal by default): the run labelled by it in the muted style.
Canvas and SVG alike. Validate the params (an unknown numbering name: one issue, `mark` used).

**Tests** (fail on the parent): the Antibody light chain's Kabat ruler: no tick past the variable domain, the note
drawn once, muted dashed line; `fallback`: ordinal labels past it, muted; the hover card past it; an unknown
`fallbackNumbering`; the SVG export alike. Look: Examples › Antibody in light and dark at dpr 2.

**Not in scope**: anything else in the ruler (CF's insertion labels and rail are merged before you start: rebase if
not). Keep to the ruler layer, its params, the describe words, the antibody stories.

**Model**: fast (sonnet). **Base**: the commit the orchestrator gives.

**Checks**: as [hc-hover-card.md](hc-hover-card.md)'s, without size. **Headers**, **hot files**, **stop when**: as
every lane's ([qt-quick-lane.md](qt-quick-lane.md)).

**Report**: commits; checks; tests failing on the parent; what you saw; a `.v6/changelogs/` entry's text.
