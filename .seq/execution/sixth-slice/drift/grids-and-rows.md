# Drift: grids and rows (a navigator's rows, a lazy alignment's order)

## What a navigator reads of the rows

**What**: grids-and-rows.md (The rows layer, "What a navigator reads") gives a navigator one value a row as its
fallback (`rowsValue`: identity to the reference row, else the non-gap share), computed with the column measures. The
code computes `RowsModel.rowsValue` (`render/rows-model.ts`) and nothing reads it: a navigator's picture lays out the
section's first rows track itself, fitted to the band, through the tracks layout (core spec, Navigators, "The
picture"), and a lazy alignment's rows are left out of it, said once (`seq.render.navigator.rows`).

**Where**: `packages/seq/core/src/render/rows-model.ts`; core spec, Navigators ("The picture", "What no summary serves
is left out, said once").

**Options**: (a) the design follows the code, and `rowsValue` goes (an unread member of a public model); (b) a
navigator draws `rowsValue` when the rows picture is not served (a lazy alignment), instead of leaving the rows out.

**The plan says (b)**: `packages/seq/TODO.md` (Navigators, after the first build) lists "the rows' one-value-a-row
fallback (`rowsValue`) and then the rows' structure where a source serves no row summary" as next work. So the member
is not dead code but a built half of planned work, and removing it would undo that plan.

**Recommendation**: keep `rowsValue` and the design's text; the row stays "not built as designed" until the TODO's
item lands, when a navigator over a lazy alignment draws it instead of leaving the rows out (and the core spec's
Navigators, "What no summary serves is left out", changes with it). Until then the full-rows picture is the rule for
held alignments: a miniature drawn by the view's own layout cannot drift from it.

## A lazy alignment in tree order

**What**: grids-and-rows.md (The rows layer, Order) orders a lazy alignment `data`, or `tree` when its source serves
rows in leaf order, and summarizes a lazy clade from the served row bins (approximate, said once). The code orders a
lazy alignment by data only (any other order `seq.render.rows.order`, info), so it has no clades and no clade summary.

**Where**: `packages/seq/core/src/render/rows-read.ts`, `rows.ts`; core spec, Render, The rows layer ("Data").

**Options**: (a) the design follows the code until a source serves rows in leaf order; (b) build tree order for such a
source.

**Recommendation**: (a): no runtime serves rows in leaf order (the alignment readers hold their rows), so the case has
no data yet.
