## Volume Tools

Interactive volume tools for Mol*. Both are driven by polygons drawn over the viewport at a
given camera orientation, so they share the `ViewMask` type, its projection and the in-place
volume operations at the root of this extension:

- `types.ts` — `ViewMask`: a polygon plus the camera snapshot it was drawn in, which is
  everything needed to project any voxel into that view.
- `view-projection.ts` — `prepareMask` / `projectToNormInPlace`: per-view camera setup and the
  world-to-canvas projection used in the voxel loops.
- `volume-ops.ts` — `flipVolumeX`, `removeDust`: in-place edits of a volume's voxel data.

The UI for both tools lives in `src/examples/volume-tools/`, one page per tool.

```bash
npm run dev -- -e volume-tools       # watch + build
http-server -p 1338 -g               # serve
# open http://localhost:1338/build/examples/volume-tools/
```

### Mask Creator

Draw polygons from any camera angle to select regions of a volume, optionally combined with
atomic structure proximity, and export a single soft-edged mask as MRC/CCP4.

- `VolumeMaskBehavior` — add to the plugin spec.
- `MaskVolumeFromSource` — state transformer building the mask volume from the view polygons,
  structure proximity, or both.

```ts
import { VolumeMaskBehavior } from 'molstar/lib/extensions/volume-tools/mask';
```

### Segmentor

Interactive segmentation of a volume into several *bodies*. Each body is defined by polygons
drawn from any camera angle: voxels above the density threshold whose projection falls inside
all of a body's views belong to it (inverted views exclude). Bodies are resolved in list order,
an optional remainder body takes what is left, and dust can be removed from the source volume.
Views stay editable per body; labels are recomputed from the definitions after every change.
Export writes one soft-edged MRC mask per body (largest body first).

- `VolumeSegmentorBehavior` — add to the plugin spec. Registers the `body-label` color theme and
  creates a `VolumeSegmentorManager`.
- `VolumeSegmentorManager` (`VolumeSegmentorManager.get(plugin)`) — headless state and
  operations: target volume, threshold, body definitions (views, order, remainder), dust
  removal, handedness flip, undo, mask previews, export. State is observable via
  `manager.behaviors.state`.
- `BodyMaskFromLabels` — state transformer building the soft mask volume of one body (used for
  previews).
- `BodyLabels` — the per-voxel label store attached to the source `Volume` (`0` = unassigned,
  `1..255` = body id).
- `internal/` — pure compute: candidates, label operations, mask computation (extend + cosine
  soft edge from an exact Euclidean distance transform), export.

```ts
import { VolumeSegmentorBehavior, VolumeSegmentorManager } from 'molstar/lib/extensions/volume-tools/segmentor';

const plugin = await createPluginUI({ ..., spec: { ...spec, behaviors: [...spec.behaviors, PluginSpec.Behavior(VolumeSegmentorBehavior)] } });
const manager = VolumeSegmentorManager.get(plugin)!;

await manager.setTargetVolume(volumeRef);      // colors its isosurface by body label
const body = manager.addBody('Head')!;
await manager.addView(body.id, viewMask);      // polygon + camera snapshot, see `ViewMask`; labels recompute
await manager.addRemainderBody();
await manager.exportMasks({ bundleZip: true, compress: false });
```

Masks are written on the source grid with values in `[0, 1]`: `1` inside the body extended by
`extend` voxels, then a raised cosine falling to `0` over `softEdge + 1` voxels. Per-body
`extend` / `softEdge` override the defaults.
