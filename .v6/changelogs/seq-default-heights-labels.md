## Seq* (next): tracks fit their content; labels inside or beside

- Tracks fit their content by default: lanes of 14 px, one lane 22 px up to five, tracks of one kind in a section
  equal, marks centred, heights decided by the whole data so a pan never resizes a track. A section's `heights` names
  the sizing (`content`, `uniform`, `fixed`, an extension's); a track's `height` overrides it; a host's `trackSizing`
  port has the last word.
- Item labels go inside their mark when they fit whole, else in one label lane above the track with a short leader;
  never cut mid-word.
- The label column is measured (80–160 px, an ellipsis past it); sections are 20 px apart; text lines 10 apart, bands
  3.
- A ranges layer with no colour takes the next categorical colour.
- An item across a cyclic origin is described by its two runs ("2,487–2,686, 1–175").
