# adopt-principles: the second review and the gap designs (2026-11-17)

## Asked

- Dynamic nodes in version 1 ("might be really hard to do later"); streaming in the features; "Streaming and patches"
  as one top-level area with what is undoable.
- A second independent review with Opus, looking for spec and design gaps; then "Yes, both in order": fix the
  contradictions, then write the missing designs in the review's order.

## Decided

- Dynamic nodes, the change model and the Streaming and patches section: [design/dynamic-nodes.md](../../design/dynamic-nodes.md),
  [design/changes.md](../../design/changes.md).
- From the review (recommended, provisional where marked): three tiers (static model and resolved view, overlay, draw
  lists); undo by `HistoryPort` checkpoints, one writer at a time; the asset policy wrapper; the designs
  [language-description](../../design/language-description.md), [locations-and-queries](../../design/locations-and-queries.md),
  [protocol](../../design/protocol.md), [layout-and-draw-lists](../../design/layout-and-draw-lists.md),
  [viewer-models](../../design/viewer-models.md), [examples-and-stories](../../design/examples-and-stories.md),
  [security](../../design/security.md); the slice's step 4 split in three; milestones after the slice.

## Delivered

Commits after `25f1edc9d`: dynamic nodes and streaming (`666e997a3`), the change model (`5f81a2512`), and this pass.
The scaffold checks and `pnpm format:check` pass.

## Learned

- The second review found what the first could not: gaps only visible once the contradictions were gone. Two reviews,
  one after the fixes of the other, were worth it.

## Left open

The provisional choices listed in `.seq/status.md`, item 5.
