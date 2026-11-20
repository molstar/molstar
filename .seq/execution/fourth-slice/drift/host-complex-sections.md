# Drift: no section can show a source's complex

**What**: lane HD's Host data story serves 2OCJ through the fake's adapter, but its view cannot take `of` the source's
complex or a dynamic unit: resolved views are built over the static model alone (`packages/seq/core/_spec/README.md`,
"no dynamic node"; `seq.render.axis.pending`), and a static `values.target` cannot name a dynamic unit. So the spec
uses static sequences `A` and `D` that read the source's units. After the fake drops chain D, section D stays with an
empty track. seq-data-ports' matrix row "a section of a complex after an update: reflows" has nothing that exercises it,
and the story never shows the host's own sequence or complex nodes in the view.

**Options**:
1. Keep it: host data reaches a view only through static nodes the spec names; a dropped unit leaves an empty section
   and an issue. Amend seq-data-ports' matrix row to say so. Nothing to build.
2. Sections over a source's complex (`of` a dynamic complex, resolved when the source opens, reflowed on `changes`),
   as the matrix row assumes. Resolve and the view model change; the static-model rule in the core spec gets an
   exception for opened sources.
3. The source brings `values` nodes too (the lane's sketch): the adapter serves the values tracks as nodes beside the
   units, so a view reads them by ref; sections still static.

**Recommendation**: 2, in the next slice's view-model work: a host's data is its complex, and the matrix row was
agreed. Until then, 1's issue so a dropped unit is not silent (HD's fix round shows it).
