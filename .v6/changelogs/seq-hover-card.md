## Seq* (next): the hover card and the keys panel

- The hover card: what the pointer rests on, in the screen reader's words (the position in every numbering, `–` where
  a numbering has no label, the letter, the item, a value in its unit, a pair's ends, a row, a grid cell's residues
  and value, a dragged rectangle's summary, a block's chain and copy), beside the pointer and never over what it
  describes; hidden on a press, Esc, a drag or a new layout; at the keyboard cursor when it moves. Drawn through the
  `hoverCard` behavior's `port` (the default spec's `hoverCardOverlay()`; a host's own adapter in its place,
  `disposeWithViewer: false` to share one); `hover.card` turns it off.
- The keys panel: `?` or the toolbar's new keyboard button ("Keys and gestures (?)", now in `DefaultToolbar`) lists the
  view's keys and gestures (`describeKeys`); non-modal, Esc closes it.
- A ruler in its sequence's own numbering is unlabelled in the examples; a field form declares a `unit`.
