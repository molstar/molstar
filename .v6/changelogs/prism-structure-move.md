## Prism: moving models by hand

- Prism: an opt-in Move section in the Selection tab (the `move` feature, `MoveFeature`; the Prism app includes it):
  the standard transform handles on a structure, a particle list, a shape or a volume (arrows move it along an axis,
  squares in a plane, rings turn it, the center moves it in the screen's plane), each drag one transform of its kind
  (`structure.transform`, `particles.transform`, `shape.transform`, `volume.transform`: one undo step), Reset back to
  where it was loaded.
- Prism: `shape.transform` and `shape.transform-reset`, `volume.transform` and `volume.transform-reset`, with the
  arguments of `structure.transform`.
- Prism: a move's release holds the canvas's frames instead of pausing it; only the model's own representations follow
  a drag; the Move section follows the models loaded and removed; a non-rigid placement keeps its scale or shear.
- `TransformShapeProvider` (plugin, `state/transforms/shape/transform`): a shape moved and turned as a whole, its
  instances' transforms premultiplied; a shape representation is made anew where its provider's transform changes,
  updated in place otherwise.
