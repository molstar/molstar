## Seq* (next): track header lines

- [Breaking] Track and section header lines in the label column: slots, state chips and a section's "n hidden" line,
  laid out in fixed columns, drawn, and hit as `header-slot`, `header-state`, `header-label` or `section-header`. The
  "n filtered" chip is the `seq:filtered` state; `TrackLayout.header` replaces `TrackLayout.filter`, and the `filter` hit
  target is gone. A click on the label column neither selects nor clears; a drag there pans.
