# Lane RR: pairs and grids by region (stage 0, design only)

**Task**: draft `drafts/region-reads.md`: the protocol's `pairs`, `matrix`, `grid` and `rows` answer shapes
(protocol.md 72); `RegionRequest.y`, 2D tiles, the store's grid summaries (grids-and-rows.md, Data, 525-535); pairs
read by region (EC-208, EC-223) with a **pairs level of bin-pair counts** and its pyramid (Gate 0, 4: bundling is the
core's pairs level, which the ring extension reads); dense matrices on a ring (`seq.render.ring.image`, layouts.md
476); an alignment wider than 16,384 columns (MAF, grids-and-rows.md 280) only if it falls out. What the store keys,
evicts and summarizes; how a level is chosen at a zoom; what a hover over a bin says.

**Reads**: protocol.md, grids-and-rows.md, files.md, sources.md, performance.md, seq-data-ports.md, the store and
`render/layer-data.ts` (`readTiles`) as built, `execution/third-slice/research-genome-ring.md`; precedents: HiGlass
(its tiles), Juicebox (.hic levels), cooler (its pyramid), Rondo.

**Shape**: [design-drafts.md](../../third-slice/briefs/design-drafts.md) (read first: the principles, the self-check, the interaction matrix, precedents, Defaults with a picture, the Configuration table with whose each choice is, edge cases by id, gaps as yes or no questions with a recommendation). Gate 0 is decided ([plan.md](../plan.md), Gate 0): build on its answers.

**Model**: strong (opus), design only. **Where**: the main checkout, read only except your one draft file; no code, no commits. Other draft lanes write beside you in `drafts/`: touch only your file.

**Report**: under 350 words: the decision in brief, the public surface and language it adds, what it amends in which agreed design, the gaps with recommendations, the self-check.
