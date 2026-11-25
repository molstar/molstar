## Seq* (next): structure link

- `@molstar/seqx-viewer/link-structure` (`linkStructure`) wires a viewer to a structure viewer's side both ways (hover,
  selection, focus), without echoes or a port, `via` per binding; the viewer publishes `hover.setLocation`,
  `addresses.revealLocation`, `view.lastJump` (its kinds `focus`, `reveal`, `fit`, `zoom`, `move`) and `retiredSignal`.
- `toUnits`/`fromUnits` over an overlay holding no nodes read the model alone.
