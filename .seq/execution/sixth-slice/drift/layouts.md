# Drift: a lanes band's scroll in the history

**What**: layouts.md (Track sizing, 3, "History and sessions") keeps a lanes band's scroll as its top lane's first item
key and restores it by that key, clamped when the lanes change. The code keeps a band by its top row's key only for a
rows layer (`behaviors/rows.ts` gives `ViewModel.setBandRows`, which names rows); a lanes band has no namer, so a step
keeps it as its `bandTop` in lanes (viewer spec, ViewModel, "Bands": "with none, a band is kept as its `bandTop`").

**Where**: `packages/seq/viewer/src/behaviors/rows.ts`, `models/view.ts` (`KeptBand`); viewer spec, ViewModel, Bands.

**Options**: (a) the design follows the code: lanes are sticky (an item keeps its lane while in view) and grow only at
the bottom (EC-245), so a `bandTop` in lanes lands on the same items in most steps back; (b) build the namer for a lanes
band (the first item of the top lane in view as its key, found again by its lane).

**Recommendation**: (b) as a small fix lane's item when a filter or a new file is seen to move a band's step back; (a)
until then. The difference shows only when the lanes change between the step and the return (a filter, EC-246).
