# Drift: region reads (a matrix's coarse form, a data filter on bins)

## A matrix's coarse form is `summary`, not `bins`

**What**: region-reads.md (The answer shapes) names a matrix's coarser-level form `bins` (`MatrixBins`). The code's is
`summary` (`MatrixSummary`, `data/protocol.ts`), and the core spec says why: a table's shape is told by its columns
(`regionShapeOf`), and pairs already use `bins`, so a matrix's coarse table under the same name would be read as pairs.

**Where**: `packages/seq/core/src/data/protocol.ts`, `tables-region.ts`; core spec, Data, Tables ("The region shapes'
forms": "`summary`: its own name, not `bins`").

**Options**: (a) the design follows the code; (b) rename the form `bins` and tell the shapes apart some other way.

**Recommendation**: (a). The name carries the reason; the design's text stays as history, its row pointing here.

## A data filter on an item's field and the pairs level

**What**: region-reads.md (Summaries, "A data filter on an item's field") says held pairs and a one-level source
summarized whole rebuild their bin pyramid for the filter (its version in the key), and a source's own levels filter
items only, saying so once (`seq.render.pairs.filter`, info: "the filter applies when zoomed in"). The code reads the
filter pair by pair on held pairs and on items read by region; bins read by region carry no fields and the filter is
not applied to them, no pyramid is rebuilt for it, and no `seq.render.pairs.filter` issue exists.

**Where**: `packages/seq/core/src/render/pairs*.ts`; core spec, Render, The pairs layer ("Themes and filter": "bins read
by region carry no fields, so neither is read on them").

**Options**: (a) the design follows the code, adding only the info issue so a reader knows the bins are unfiltered;
(b) build the design: summaries rebuilt per filter for held and one-level sources, the issue for a source's own levels.

**Recommendation**: (a) now, with the issue as a small fix lane's item; (b) when a filter on a numeric field over bins
is asked for (EC-225's Q-value case), since a rebuild per filter costs a whole read each time the filter changes.
