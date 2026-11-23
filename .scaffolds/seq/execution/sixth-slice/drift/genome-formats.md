# Drift: a map's unequal rows, said by the runtime or by the layer

**What**: genome-formats.md (The `seq` runtime, "As a map") says a row of unequal ends in a pairs source read as a
map's runs is refused as a run by the runtime, counted, `seq.source.seq.run`. In the code the runtime serves every row;
the reader of runs (the liftover, `render/lift.ts`, `runsOf`) leaves out a row of unequal ends and says it by count,
`seq.render.lift.run` (info). No `seq.source.seq.run` code exists.

**Where**: `packages/seq/core/src/render/lift.ts`; core spec, Render, A liftover ("The runs by the near end").

**Options**: (a) the design follows the code: the source serves pairs as they are (the same rows feed a bridge's
ribbons, where an unequal block is a fine ribbon), and only a reader that needs runs refuses them; (b) the runtime
refuses them when the source is read as a map.

**Recommendation**: (a). One source serves both the ribbons and the liftover (the design's own reason for the pairs
shape), so the refusal belongs to the reader that needs equal lengths; the issue's code says where it was found.

**Resolved** (the author, 2026-11-23, as recommended): genome-formats.md (The `seq` runtime, "As a map") now says the runtime serves every row and the liftover's reader leaves out unequal ends, counted (`seq.render.lift.run`); its "In the specs" row points at core › A liftover.
