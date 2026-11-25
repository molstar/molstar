## Seq* (next): declared numberings

- The `numbering` kind (a child of a `sequence`, `method: 'labels'`): a label a position, inline or from a records
  source's column, with an open lowercase `scheme` word (`kabat`, `imgt`, …). Queries (`(at :label 82B :numbering
  kabat)`, `(range … :numbering imgt)`), copied addresses, the ruler and hover read it, with insertion codes and IMGT's
  out-of-order labels; a numbering still loading resolves as waiting; a source's numbering is capped
  (`maxNumberingPositions`). Language 0.4.0 (unreleased, an addition).
