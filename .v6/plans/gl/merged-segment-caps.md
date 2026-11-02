# Plan: Merged caps by segment

## Status

- **Done**: the solid-interior caps of a merged render object ([merged-render-objects.md](merged-render-objects.md))
  draw only the members the near plane or a clip object cuts, not every member where one is cut. From an issue report
  (WebGL2, GTX 980). No defines, no shader changes: the renderer lists the cut segments, the merged renderable draws
  them.

## Question

A merged surface cut by the near plane was slower than the same surface unmerged (the report: 4V6X Gaussian surface,
near at half the scene radius, 5.80 against 4.04 ms GPU time a frame). `cutsNearObject` (`renderer.ts`) tested the
segments' spheres only to decide whether the object has a near cap: then the mark and fill passes of the cap drew the
object's lists, every member the main pass's culling kept. Unmerged objects cap only those cut. The same held for the
caps of clip objects (`ClipInfo.capIndices` is tested against the object's sphere).

## Done

- **The renderer** (`getSolidInteriorCaps`): per cap the segments it cuts (`capSegments`), by the segments' spheres
  where the object draws several (`getSegmentSpheres`, `renderSegments`): the near plane as `cutsNear` per segment, a
  clip object by `Clip.canIntersectSphere` per segment; a cap cutting none is not drawn, one cutting all draws the
  object (`render`), else the mark and fill passes draw the cut segments (`drawObject` with segments). Recordings: the
  cap key adds each cap's segments (`walkCaps`), so a recording replays while the same segments are cut.
- **The merged renderable** (`renderSegments`): the listed segments, each whole at its levels (`fullSegment`), as lists
  of its own (the main lists stay), at the level view of the last lists; shifted to the identity region where culled
  on the GPU (`shiftLists`, as the CPU lists).
- **Measured** (`merging/4v6x-gaussian-merged-near` and `-unmerged-near`: the near plane at half the camera's radius,
  `cameraRadius`; RTX 4080 class, GPU ms a frame): WebGL2 merged 0.978 → 0.751 (draws 270 → 156, as unmerged), unmerged
  1.02-1.08; WebGPU merged 1.190 → 0.961, unmerged 0.92-0.95. Images: the renderer and merged visual cases
  (`renderer/**`, `geometry/merged/**`, the near cut of `geometry/merged/gaussian-surface-near-cut`) identical on
  WebGL2, WebGL1 and WebGPU, and on WebGPU with draw lists.

## Not done

- Segment ranges of the main lists (a cap drawing the main pass's culled entries of its segments): the cut segments are
  drawn whole at their levels instead, which costs nothing measurable for the few segments a plane passes through.
- A plane through thousands of segments (a slab clip across a mesoscale scene) lists them all per cap and frame;
  measured when such a scene is at hand (36ZA, 2888 segments, under a slab clip).
