# Drift: the matrix layer's genotype codes borrow ramp tokens

**What**: RG's `matrix` layer colours codes 0, 1, 2, 3 and 255 (no copy, some copies, every copy, partly called,
uncalled) with fixed tokens borrowed from the value ramp (`seq-ramp-*`) and `seq-missing`, not a theme channel.

**Where it disagrees**: region-reads.md Defaults 9 ("tokens only") holds, but no design says which token a code takes,
and the ramp's steps are meant for values, so a theme that changes the ramp repaints genotypes by accident.

**Options**: (a) tokens of their own, `seq-genotype-0` to `seq-genotype-3` beside `seq-missing`, in both themes, the
design naming the mapping; (b) a categorical theme channel on the matrix layer (codes as categories), defaults the same
colours; (c) keep the ramp steps, the design saying so.

**Recommendation**: (a) now (a small token addition, the design's Defaults 6 naming them), (b) when a figure needs
another palette. Asked at stage 3's gate.

**Decided** (the author, 2026-11-23, at stage 3's gate): (a), tokens of their own. To a fix lane.
