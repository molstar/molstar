# Lane DP: the Plasmid dataset

Follow [datasets.md](datasets.md) for [examples](../../../design/examples.md) › Plasmid (pUC19c, INSDC L09137.2; the
features derived in code; restriction sites computed; the rotated and reverse-complement variants; the 200 kb scale
variant). Also, from the plan (the author, 2026-11-19: "should the plasmid rings have arrow on the ranges? usually saw
people draw it like that"): features drawn as arrows by strand on the ring and the linear detail (the builder sets
the `ranges` layer's `arrow` channel, EC-232), the arrowhead bent through the ring's polar warp, on the piece holding
the end it points to (a feature across the origin included), checked in canvas and SVG. The features may be served as
GFF3 (FS-A2's preset, the cyclic wrap) where the section's "derived in code" allows writing them so; say which.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/dp-plasmid
<that commit>`. **Storybook port**: 6012.
