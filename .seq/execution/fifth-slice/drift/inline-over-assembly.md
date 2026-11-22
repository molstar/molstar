# Drift: which column names a held pair's (or grid cell's) sequences over an assembly

**What**: genome-coordinates.md ("Annotations over an assembly") says `target` (ranges, values, pairs, grid) may name
an assembly, each row on the sequence its `target` field names. For ranges and values that is one column, `target`.
A pair has two ends and a grid cell two axes, and no design names their columns; maps have `from_name`/`to_name`,
BEDPE files `chrom1`/`chrom2`. GA2 reads sources over an assembly by sequence; inline pairs, grids and matrices over an
assembly read nothing, now said once (`seq.resolve.assembly.inline`, info).

**Options**: (a) `a_name` and `b_name` for pairs, `x_name` and `y_name` for grids, `name` for a matrix's items (parallel
to maps' `from_name`/`to_name`, the protocol's `a`/`b` and `x`/`y` sides); (b) reuse `target` as a list per row; (c)
inline data over an assembly stays refused: a host writes a source.

**Recommendation**: (a): the names follow the sides the protocol already uses, a row stays flat, and a builder can
type them. A language param (GC1's minor version holds the slice's params): a small bump or the next slice's.
