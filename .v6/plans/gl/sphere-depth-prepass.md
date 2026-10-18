# Plan: Depth pre-pass for sphere impostors (WebGPU)

## Status

- **Current**: dropped after step 1: the pre-pass alone costs about as much as the impostors without lighting, so it
  can't pay off (results below). Step 2's read-only depth passes are kept as a general facility.
- Step 1 done (WebGPU, opaque pass, 2560x1440): the spheres' `depth` variant in place of the opaque draw, the cost of
  a pre-pass, against impostors as today, without lighting (`ignoreLight`) and approximate spheres:

  | Scene                             | Impostors | Depth-only pre-pass | Without lighting | Approximate |
  | --------------------------------- | --------- | ------------------- | ---------------- | ----------- |
  | globule, 300k atoms               | 14.2 ms   | 10.4 ms             | 10.5 ms          | 11.0 ms     |
  | assembly, 60 copies of 7qpd       | 2.42 ms   | 1.70 ms             | 1.70 ms          | 1.00 ms     |
  | flat layer, 90k atoms, no overdraw | 0.18 ms   |                     | 0.12 ms          | 0.14 ms     |

  Most of an impostor's cost is the intersection and depth write without early depth testing, which the pre-pass
  repeats; the lighting it would save is the difference to `ignoreLight` (3.8 ms for the globule, 0.7 ms for the
  assembly), less than the shading pass's own rasterization (at least 4.3 ms for the globule's quads, about 1 ms for the
  assembly). What could help is early depth testing for the depth output itself (conservative depth, option 4 of
  sphere-impostors.md), blocked on WebGPU by Chrome's D3D12 compile; approximate spheres stay the faster option.
- Step 2 done: `RenderPassOptions.depthReadOnly` (and `describeTarget(target, { depthReadOnly })` for prepared
  pipelines): draws test against the depth but don't write it, it is not cleared (a clear throws). WebGPU begins the
  pass with `depthReadOnly` (and `stencilReadOnly` with stencil formats) and no load or store ops; the target gets its
  own key (`withDepthReadOnly`), whose pipelines mask depth and stencil writes (`getDepthStencilState(..., readOnly)`),
  so any render state draws in it; WebGL masks the writes in its render state applier (the attached texture must not
  be bound there). Test on Dawn (`webgpu-depth-read-only.test.ts`): the depth attachment copied into the color target
  in a read-only pass, the same pass with a writable attachment a validation error, a depth clear throwing; the
  masking in a unit test. Visual suite (all profiles) and smoke unchanged. Option 5 of [sphere-impostors.md](sphere-impostors.md), viable since option
  3 showed that WebGPU keeps early depth testing for fragment shaders without a depth output, also with `discard`.

## Question

Sphere impostors write their per-pixel depth (`frag_depth`), which turns early depth testing off: every fragment of
every sphere is shaded, including those hidden behind nearer spheres. In dense scenes most are hidden: drawing the 300k
atom globule a second time at equal depth costs 14.3 -> 28.5 ms on WebGPU (opaque pass, 2560x1440), while approximate
spheres, without a depth output, cost 11.0 -> 11.7 ms. Can a depth-only pass for the spheres, followed by a shading
pass without depth output, give impostors the cost of visible fragments only?

## Idea

1. **Depth pre-pass**: the spheres' `depth` variant (ray-sphere intersection, `discard`s, depth write, no shading)
   into the opaque depth attachment, no color writes.
2. **The other opaque objects** as today (with depth writes): what lies behind a sphere is rejected early, as the
   spheres' depth is already there.
3. **Shading pass for the spheres**, the depth attachment read-only (`depthReadOnly`), the `color` variant without a
   depth output:
   - the quad at the sphere's front-most depth instead of its center, so that it bounds the depths it covers: with
     `greater-equal`, early depth testing rejects spheres hidden as a whole (coarse culling);
   - in the shader, the exact depth (as in the pre-pass) compared with the depth attachment, bound as a texture as
     well (allowed while it is read-only): `discard` where a nearer sphere won (pixels of partially hidden spheres
     that pass the coarse test); the pixels that remain are those the pre-pass kept, shaded once.

Cost: the spheres' vertex and rasterization work twice (about 1.5 ms of the globule at 2560x1440) and their intersection
twice; saved: the shading of hidden fragments. It pays off with overdraw and costs a little without; measured per
scene (step 1), the default decided on the numbers (decision 2).

## Design

- **Renderer**: pass draws for the pre-pass (the `depth` variant of objects with depth-writing impostors, a render
  state without color writes) and the shading pass (their `color` variant with `dDepthPrepass`, `greater-equal`,
  depth read-only), next to the existing opaque draws, which leave those objects out while the pre-pass is on. As
  `PassDraws` entries, so that `renderer.collect` and `Canvas3D.prepare` cover their pipelines.
- **`DrawPass`** (opaque, not to the drawing buffer): `begin` the opaque target for the pre-pass, the other opaque
  objects, then a render pass with the depth read-only for the spheres' shading, as `PassBuilder` plans (prepare and
  draw share them). Transparent spheres, picking, marking, the transparent depth, illumination and multi-sample keep
  their paths (multi-sample repeats the draw pass, so it inherits the pre-pass).
- **Neutral API**: a read-only depth option for render passes (`RenderPassOptions`), on WebGPU `depthReadOnly` with
  the depth attachment also bindable as a texture in that pass; WebGL does not use the pre-pass (decision 3).
- **WGSL spheres**, `dDepthPrepass` in the `color` variant: no `frag_depth` member (its presence alone disables early
  depth testing); the quad at the front-most depth (the view-space center moved by the radius towards the camera,
  clamped to the near plane, as the near-plane handling does today); the depth comparison against the attachment
  with a small tolerance (the same intersection code as the pre-pass, so equal up to rounding). To check: the ray of
  `sphereImpostor` comes from the unprojected quad point (`vPoint`), which moves along the view ray when the quad
  depth changes (perspective) or keeps its direction (orthographic); the intersection must not depend on where on
  the ray the point lies.
- **Must match the pre-pass**: every `discard` of the `depth` variant (clipping, LOD fade, caps, interior) also in the
  shading pass's coverage, else pixels the pre-pass left empty would be shaded, or vice versa holes; both variants
  share the intersection and discard code, which the visual tests compare against the current image.

## Decisions

### 1. Scope

- **A. Spheres** first: the measured case, the largest impostor counts (spacefill).
- **B. All depth-writing impostors** (spheres, cylinders of ball-and-stick, text with depth): one mechanism, more
  shader work; cylinders follow A if it pays off.

Recommended: A, then B for cylinders with their own measurement.

### 2. When it is on

- **A. Always on WebGPU** if the cost without overdraw is negligible.
- **B. A renderer prop** (`depthPrepass: 'auto' | 'on' | 'off'`), `auto` by a heuristic (e.g. the number of
  impostor spheres, or their projected area against the viewport).
- **C. Off by default**, a prop for scenes that need it.

Recommended: decide after step 1's numbers; B if the no-overdraw cost is measurable.

### 3. WebGL

- **A. WebGPU only**: WebGL has no read-only depth attachments (the shading pass would read a copy of the depth), and
  ANGLE on D3D11 did not reject hidden fragments of the approximate spheres early (with `discard`), so a gain there is
  uncertain.
- **B. Both**, WebGL with a depth copy.

Recommended: A; B only after measuring WebGL with a minimal case.

### 4. How the shading pass reads the pre-pass depth

The shading pass uses the depth twice: the depth test (coarse culling of hidden spheres) and the shader (the exact
comparison of partially hidden ones). WebGPU forbids a texture as a writable attachment and a sampled texture in one
pass.

- **A. Read-only depth attachment** (`depthReadOnly`, `stencilReadOnly` with stencil formats): the attachment's usage
  becomes read-only, which allows binding the same texture for `textureLoad`; no load or store ops, pipelines without
  depth or stencil writes. Needs a read-only option in the neutral render pass API. WebGPU only: on WebGL sampling a
  texture attached to the bound framebuffer is a feedback loop, also without depth writes.
- **B. A copy of the depth** after the pre-pass, read by the shader, the attachment writable as today: a full-size copy
  per frame (about 14 MB at 2560x1440) and a second texture, no new API, and the way WebGL would go.

Chosen: A (no copy, no extra texture); B is the route for WebGL, should decision 3 change.

## Steps

1. **Measure the potential**: the opaque pass of `render/globule-300k-spacefill-large`, `render/7qpd-spacefill`,
   `render/assembly-60-spacefill-rotating` and a case without overdraw (a single layer of spheres facing the camera),
   each with impostors as today and with approximate spheres (a lower bound of the shading pass's cost), WebGPU.
2. **Read-only depth passes** (decision 4A) in the neutral API and the WebGPU backend (`depthReadOnly` and
   `stencilReadOnly`, the attachment bound as a texture), with a test on Dawn.
3. **The spheres' shading variant** (`dDepthPrepass`): quad at the front-most depth, no depth output, the comparison
   against the pre-pass depth; checked against the current image in the visual suite (spacefill cases, clipping, LOD,
   interior, near plane, orthographic and perspective).
4. **Renderer and `DrawPass`**: the pre-pass and shading draws, the opaque draws without those objects, the plans;
   `Canvas3D.prepare` covers the new pipelines; the visual suite unchanged within the WebGPU tolerances.
5. **Benchmark** the cases of step 1 with and without the pre-pass; decide decision 2; record the results.
6. **Cylinders** (decision 1B), if step 5 pays off.
