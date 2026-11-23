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

**Recommendation**: (a), with `rowsValue`'s removal a fix lane's item: the navigator's picture is the view's own layout
by design (navigators.md: a miniature drawn by other code drifts), and a lazy alignment's fallback would need its rows
whole, which is what lazy rows avoid. Keep (b) in mind if a lazy alignment's overview is asked for.

## A lazy alignment in tree order

**What**: grids-and-rows.md (The rows layer, Order) orders a lazy alignment `data`, or `tree` when its source serves
rows in leaf order, and summarizes a lazy clade from the served row bins (approximate, said once). The code orders a
lazy alignment by data only (any other order `seq.render.rows.order`, info), so it has no clades and no clade summary.

**Where**: `packages/seq/core/src/render/rows-read.ts`, `rows.ts`; core spec, Render, The rows layer ("Data").

**Options**: (a) the design follows the code until a source serves rows in leaf order; (b) build tree order for such a
source.

**Recommendation**: (a): no runtime serves rows in leaf order (the alignment readers hold their rows), so the case has
no data yet.
