## Seq* (next): structure translation

- `@molstar/seqx-core/structure`: `structureFacts`, `toUnits`/`fromUnits` (the copy rule, unobserved positions cut, a
  size gate), `unitRows` and `colorRuns` (colors compiled as the track draws them) for a structure viewer next to a
  sequence view.
- [Breaking] `SeqPorts.structure` removed: syncing is a library, not a port.
