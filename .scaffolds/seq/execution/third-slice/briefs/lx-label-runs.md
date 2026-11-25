# Lane LX: leaders never cross; a run of same labels once (fast)

**Goal** (the author, 2026-11-20, Complex › Blocks › Row: the Interface track's labels "interface interface
interface", a leader drawn through a neighbouring label; "Label a run once"; [layouts.md](../../../design/layouts.md),
"Labels inside if they fit, else beside"): (1) a leader never crosses another label or another item's mark: a label
that cannot be placed without a crossing is dropped (or tries the next free place in its lane before dropping); (2)
neighbouring items of one layer with the same label, less than a label's width apart, get one label at the run's
start with no leaders to each; the hover card still names each item; the run's label is hit as its first item.

**Tests** (fail on the parent): the Row story's Interface track: one "interface" label per run, no leader crossing any
label's box or any mark (a geometric check over the draw lists), canvas and SVG alike; a dense mixed layer: no
crossing. Perf: 0 bytes a frame. Look: Complex › Blocks › Row and Antibody › Schemes, light and dark at dpr 2.

**Not in scope**: anything else. Keep to the label placement in `ranges.ts` (DB's) and its tests.

**Model**: fast (sonnet). **Base**: the commit the orchestrator gives. **Checks**: as
[hc-hover-card.md](hc-hover-card.md)'s, without size. **Headers**, **hot files**, **stop when**: as every lane's.

**Report**: commits; checks; tests failing on the parent; what you saw; a changelog line.
