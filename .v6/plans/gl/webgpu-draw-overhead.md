# Plan: CPU cost per draw call on WebGPU

## Status

- **Current**: done, steps 1–6. Continued by [webgpu-draw-submission.md](webgpu-draw-submission.md). Since: persistent
  slots are `UniformArena.createSnapshot` in the shared `uniformPages` (not a buffer per item), the writers' gap below
  closed (they skip unchanged blocks).
- Per draw of an unchanged object in direct mode (audited against `units-1000-spacefill-unmerged`, 2.31 ms CPU per
  still frame; draw lists replay without these), done: the global group (`global-bindings.ts`) the last per layout,
  compared by the uniform buffer, a generation of the textures list and the views, instead of its texture entries and
  key built every draw (as step 1 did for groups 1 and 2; 3.7% to 0.6%); group 2's storage buffers by cells resolved
  per program, one result object per cached group; `getDrawBlock` builds its layout and entries only for a new block;
  `getVertexBuffer` keeps its slot entry. 2.31 to 2.20 ms (`boxes-1000` 2.64 to 2.50). Then the culling of still
  frames: every gridded or leveled object culled again each frame (a frame token only); now skipped while the view and
  the culled cells' versions stay (`createCullView`, renderables and merged ones; not with an occlusion test, not in
  timing mode, invalidated by `uncull` and `cullSimple`; `cull-view` test). The view as a version (`LodView.version`):
  the renderer compares the camera plane, frustum and level scales once, at the first cull after an update, with those
  of the last cull (updates of other cameras between, as picking's, don't cull); compared per object first (28
  numbers, 7% of WebGL2's unmerged frame). 2.20 to 1.62 ms (rotating 2.43 to 2.21; merged 1000 units 0.31 to 0.20,
  100k CPU-culled instances about 1.0 to 0.12; WebGL2 the same for merged, its unmerged frame 5.8 to 5.5, bound by its
  draws). Not done: the `Uint32Array` offsets of `setBindGroup` (`PassState`
  copies the offsets to compare them, the native call is rare).
- Steps 1–4 done. `render/boxes-1000`: CPU time per frame 10.5 → 7.1 ms (frame 10.8 → 7.3 ms; WebGL2 4.3 and 5.4 ms);
  `render/boxes-100`: 0.68 → 0.40 ms (WebGL2 0.24 ms); the other render cases unchanged (GPU-bound or few draws). By
  step: 1–2 (caches, linking) 10.5 → 8.9 ms, 3 (persistent uniform slots) → 7.1 ms, 4 (redundant state) no change at
  1000 boxes (each has its own buffers and bind groups, only the pipeline and the global group repeat), 100 boxes 0.52 →
  0.40 ms.
- Steps 5–6 done:
  - **Shared material blocks**: render items of one material id share its block and bind group (as the WebGL2 material
    buffers; each item writes its values into the shared data, materials have no textures). The visuals of a
    representation share its material: `render/units-1000-spacefill-unmerged` (1000 units of one representation,
    `renderMerged` off) 7.5 → 5.3 ms CPU per frame, below WebGL2 (7.5 ms).
  - **Draw uniforms in the object's bind group**: `uDraw` is a binding of group 2 (`getWgslDrawBinding`, after the
    object textures) instead of group 3; the object block and the item's draw blocks (per global bindings and multi-draw
    level) are slots of one buffer per item (`UniformArena`, which also replaced the per-block buffers), so one bind
    group with two dynamic offsets serves all draws of an item. One `setBindGroup` less per draw (16% less time in
    them); `render/boxes-1000` unchanged within noise (7.1–7.4 ms), as each box is its own material.
  - Fixed on the way: after a change within a submission, the snapshot wrote to the ring on every read of that
    submission (now once per change).
- Remaining gap for separately drawn objects (`render/boxes-1000`, 7.2 ms against 4.4 ms on WebGL2): the uniform block
  writers (about 19%, the version checks of all members; shared with WebGL) and the native calls for the buffers and
  bind groups of each item. A cheaper unchanged check for the writers would help both backends.
- Uniform writers skip unchanged blocks: `ValueCell.observe` flags an observer when a cell changes (`update`, `set`),
  `createUniformBlockWriter` checks the versions of its cells only after one of them changed and stops observing on
  `dispose` (render items and global bindings of both backends, compute kernels). The writers no longer show in the
  profiles (were about 19% of `render/boxes-1000` on WebGPU, 17% on WebGL2). Alternating runs with and without (noisy
  machine, ±30% between rounds): WebGPU `render/boxes-1000` 6.0–6.4 against 6.7–8.9 ms CPU per frame,
  `render/units-1000-spacefill-unmerged` 4.6–4.9 against 5.5–7.8 ms; WebGL2 within noise (its unmerged units are
  dominated by the attribute offsets of the culled sphere levels, `changeOffset`).
- Sphere impostors on the GPU (`render/globule-300k-spacefill`, WebGPU 4.6 against 4.3 ms; `pnpm benchmark --gpu-timing`
  prints the GPU time per pass group): the whole gap is in the opaque pass (4.36 against 4.06 ms), postprocessing is the
  same. The frame is pixel-bound (0.5 ms at 320×180, 14.6 ms at 2560×1440) and mostly rasterization: the approximate
  spheres (no depth writes, so early depth tests, and a simpler shader) are only 26% faster. Not the cause: the depth
  format (`depth24plus-stencil8` instead of `depth32float-stencil8`, no change), robustness (`disable_robustness`, no
  change). Conservative depth was tried and dropped: the quad moved to the sphere's front with `depth_greater` (WebGL2,
  `EXT_conservative_depth`) measured no difference on NVIDIA/ANGLE D3D11; WGSL `@builtin(frag_depth, greater)` (language
  feature `fragment_depth`, Chrome 155+) fails to compile with `discard` on D3D12 in Chrome 156 (DXC E_FAIL; a fragment
  shader with only the depth output compiles), which the impostors need. The WebGL2 prefix had put
  `GL_EXT_conservative_depth` into the vertex shaders, now the fragment shaders. Generated shaders compared (Dawn's
  `dump_shaders` toggle for the WGSL, HLSL and DXIL; `WEBGL_debug_shaders` for ANGLE's HLSL): the WGSL sphere program
  mirrors the GLSL (same discards, lighting, fog, outputs; two color targets on both), Tint's HLSL is structurally alike
  (short-circuit operators as if/else chains, uniforms read from `uint4` constant buffers with `asfloat`, saturating
  float-to-int conversions, overflow-checked signed modulo); the DXIL is dominated by the inlined noise of bump mapping
  and animations (`sin`/`frac`, runtime-skipped branches), as in ANGLE's code. Same compiler and driver: WebGPU on
  Dawn's D3D11 backend (`--use-webgpu-adapter=d3d11`, Tint HLSL through FXC) takes 5.2 ms against 4.3 ms for WebGL
  (ANGLE HLSL through FXC), D3D12 (DXC) 4.6 ms: the gap is in the translation and its backend, not in the WGSL source.
  Unsigned index math in `readFromTexture` (no overflow checks) made no difference; no change of ours is indicated.
- Merged render items (`dSegmented`) are ported to WebGPU (they drew nothing before; representations of more than 100
  units merge them by default, `renderMerged: 'auto'`): the group and group-instance data offset per segment
  (`aSegment`), the invariant bounding sphere per segment (`invariantBoundingSphere(in)`, `segment.wgsl.ts`, for
  tumbling and instance clipping) and the per-segment lod scales of spheres (`aSegmentLod[uLodLevel]`). Visual cases
  `geometry/merged/*` (spacefill impostors and meshes, backbone cylinders, points, lines; colored by chain, and by
  element index for the group-instance data; chains of different lengths so that the group-instance offset is not zero);
  they fail with either offset removed from the WGSL. `render/units-1000-spacefill` on WebGPU: 0.29 ms CPU per frame for
  its 1000 segment draws (WebGL2 1.7 ms with one multi-draw). The cases also found merged sphere impostors blank on
  WebGL1: the merged values copied the members' emulated vertex ids (`aVertex`, never filled, as the members are not
  drawn), now left to the render item.
- Analysis (below). Found with the render benchmark ([backend-selection.md](backend-selection.md#default-backend)): at
  1000 separately drawn objects (`render/boxes-1000`) WebGPU takes about 10.5 ms of CPU time to issue a frame, WebGL2
  about 4.4 ms; GPU-bound frames are on par. Profiled with `pnpm benchmark --filter render/boxes-1000 --profile <dir>`
  (V8 CPU profiles per backend, open in the DevTools Performance panel).

## Profile

Samples under the frame loop of `render/boxes-1000` (1008 draws per frame), share of the CPU time and the time per frame
(about 11.4 ms on WebGPU, 5 ms on WebGL2):

| Cost                                                                                                      | WebGPU       | WebGL2         | Cause                                                                                                                                                                            |
| --------------------------------------------------------------------------------------------------------- | ------------ | -------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Program lookup per draw (`getWebGPUProgram` → `getLinearTextures`, `isWebGPUTexture`, `isLinearFiltered`) | 13% (1.5 ms) | –              | builds a key string from all textures of the item (a `WeakSet` check and the filter of each) on every draw, then a `Map` lookup by a template string                             |
| Bind groups 1 and 2 per draw (`getBindGroup` → `getTextureEntries`)                                       | 13% (1.5 ms) | –              | builds the texture entries and their key string on every draw, to compare with the cached group                                                                                  |
| Uniform snapshots (`getSlot` → ring `write`, `upload` in `submit`)                                        | 12% (1.4 ms) | –              | the ring is reset per submission: every block of every item (object, material, draw) is copied into it and uploaded again each frame, changed or not                             |
| `linkPrograms` in every commit                                                                            | 7% (0.8 ms)  | –              | calls `link` of all render items for each variant on each `commit` (WebGL only visits programs still pending)                                                                    |
| Vertex buffers per draw (`getVertexBuffer`)                                                               | 6% (0.7 ms)  | –              | maps the slot's names to cells (an array per call) and compares versions with `every`                                                                                            |
| WebGPU calls (self time of `render`)                                                                      | 23% (2.6 ms) | 11% (GL calls) | `setPipeline`, four `setBindGroup`, `setVertexBuffer` per slot, `setIndexBuffer`, `draw` per draw; repeated state is set again (one pipeline and one global group for all boxes) |
| Uniform block writers (version checks of all members)                                                     | 14%          | 17%            | shared, per block per draw                                                                                                                                                       |
| Renderer (`checkOpaque`, `isVisible`, `drawObject`)                                                       | 6%           | 9%             | shared                                                                                                                                                                           |

The rows above the shared ones are WebGPU's extra cost: about 5.9 ms per frame, the gap to WebGL2.

## Steps

Implemented as follows: the program per variant and the resolved textures of bind groups 1 and 2 are cached per render
item and invalidated (`invalidatePrograms`) with the defines, the transparency, replaced textures or a restored device,
bind groups compare their textures' views per draw; `WebGPULinkables` keeps the items pending per variant; render items
write their object, material and base draw blocks to buffers of their own (`createPersistentUniformSnapshot`, a ring
snapshot after a change within a submission that already read them), multi-draw levels stay on the ring; `PassState`
(`gl/webgpu/pass-state.ts`) skips repeated pipeline, bind group, buffer, viewport, scissor and stencil reference calls.

1. **Cache per render item, invalidate on change**:
   - linear textures and the program per variant: recompute only when the item's textures change (`update`, texture
     values or their filter), not per draw;
   - texture entries of bind groups 1 and 2: rebuilt with the program or the textures, the per-draw check is a buffer
     identity (the snapshot slot);
   - vertex buffer slots: the cells per slot resolved once, versions compared in a loop without allocations.
2. **Link only new items**: `linkPrograms` visits render items not yet linked for a variant (a pending set as on WebGL,
   emptied as items are linked), not all items per commit.
3. **Persistent uniform slots**: blocks keep their slot across submissions while unchanged (an arena of slots per block
   size, allocated per block and freed on destroy, written with `writeBuffer` when the data changed); a block changed
   after a draw of the current submission used it falls back to a ring snapshot for the rest of that submission
   (`writeBuffer` runs before the submission's commands). Uploads per frame shrink to the changed blocks; the global and
   draw blocks that change between passes stay on the ring.
4. **Redundant state**: the pass encoder wrapper skips `setPipeline`, `setBindGroup` (same group and offsets),
   `setVertexBuffer` and `setIndexBuffer` calls that repeat the current state of the pass.

Expected: steps 1–2 remove about 3.5 ms per frame of `render/boxes-1000`, step 3 about 1.4 ms, step 4 part of the 2.6 ms
of WebGPU calls; together close to WebGL2. Each step on its own commit, measured with
`pnpm benchmark --backends webgl2,webgpu --filter 'render/*'` and checked with the visual suite (`webgpu`,
`webgpu-core`) and the Viewer smoke test (device loss restores the caches: `restorables`).
