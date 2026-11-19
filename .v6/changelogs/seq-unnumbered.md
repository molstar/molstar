## Seq* (next): a ruler past its numbering's reach

- A ruler marks where its numbering does not reach (an antibody chain's constant domain under Kabat): the new ruler
  layer param `unnumbered` is `mark` by default, a muted dashed line and one "not numbered" note, no ticks, or
  `fallback`, which labels the run by `fallbackNumbering` (ordinal by default) in the muted style. The hover card and the
  live region say "not numbered by <numbering>" there. An unknown `fallbackNumbering` is
  `seq.render.numbering.fallback`.
