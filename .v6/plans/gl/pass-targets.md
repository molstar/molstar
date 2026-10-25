# Plan: Pass targets allocated on use

## Status

Done (steps 1, 2 and 4; step 3 measured, not needed for now). Pass targets per pixel of the drawing buffer, once the targets
the props and scene don't draw into are released (headless WebGPU; `webgpu-pass-target-bytes` checks them):

| Props, scene (defaults: WBOIT, FXAA, bloom emissive, temporal multi-sampling) | B/px | at 1872×1536 |
| ----------------------------------------------------------------------------- | ---- | ------------ |
| opaque                                                                        | 28   | 77 MiB       |
| transparent                                                                   | 48   | 132 MiB      |
| transparent, DPOIT                                                            | 76   | 208 MiB      |
| transparent, highlighted (marking stays once used)                            | 68   | 186 MiB      |
| illumination (marking kept)                                                   | 100  | 274 MiB      |
| illumination, multi-sample on                                                 | 120  | 329 MiB      |

Before: 932 MiB allocated up front by `Passes` at this size in every case; after steps 1–2, 110 MiB opaque, 208 MiB
transparent, 428 MiB with illumination. A screenshot's `ImagePass` holds its own set while rendering and a second
after. Follows the review item M5 ([../../webgpu/followup.md](../../webgpu/followup.md), Memory and performance 5).

## Measured

Headless WebGPU, drawing buffer 1872×1536 (2.9 Mpx), by allocation site (instrumented `createTexture`; "used" = a
render pass attachment, a bind group entry or a copy in three frames):

| Props, scene                                    | Allocated (`Passes`) | Used    |
| ----------------------------------------------- | -------------------- | ------- |
| defaults (WBOIT, FXAA, bloom, temporal MS), opaque | 932 MiB           | 181 MiB |
| defaults, transparent (alpha 0.5)               | 932 MiB              | 291 MiB |
| illumination on, transparent                    | 932 MiB              | 411 MiB |

- `Passes` allocates every target of every feature up front, about 330 B per pixel: both OIT modes (WBOIT 110 MiB,
  DPOIT 132 MiB), the tracing and illumination targets (290 MiB), bloom (95 MiB), marking (110 MiB), multi-sample,
  outline, shadow, DOF, SSAO, SMAA. Most are full-size rgba8 (4 B/px), rgba32f (16 B/px) or with a depth-stencil
  renderbuffer (WebGPU `depth32float-stencil8`, 8 B/px).
- An `ImagePass` (screenshots, thumbnails; the headless plugin keeps one) holds a second full set at its size: the
  headless plugin at this size holds 1.9 GiB of pass targets.
- Used but not needed: the shadow and outline targets are in the postprocessing compose bind group with both off (the
  shader reads them only under `dShadowEnable`/`dOutlineEnable`); render targets' depth renderbuffers are allocated
  where the pass attaches the shared depth texture instead (`colorTarget` with WBOIT).

- Cross-checked with the driver (Windows `\GPU Process Memory(pid_*)\Dedicated Usage` of the headless test process,
  Dawn on D3D12): its growth matches the bookkeeping (illumination +202 MiB against +197 MiB; a screenshot +247 against
  +242), and releases show there once Dawn processes its deferred deletions, which it does when work is submitted after
  the destroy: without frames (an idle canvas) memory stayed taken until the next submission. The WebGPU context
  therefore submits once (empty) after destroying resources outside a frame; then a release returns the memory without
  frames (481 → 305 MiB turning illumination off).

## Design

### Lazy storage (step 1)

- Textures and renderbuffers get their GPU storage when first written: a render pass attachment, a copy into them, an
  upload (`load`). `define(w, h)` only records the size; the object, its references in value cells and render
  targets stay as they are, so passes change little.
- Sampled before written: a placeholder (WebGPU: the 1×1 texture per sample type and dimension of loading images;
  WebGL: the texture's 1×1 image), not allocated on bind: bound-only textures take no memory, and the results are as
  before, since new storage reads zero on both backends (WebGPU zero-initializes, `texImage2D(null)` is zero) and a
  target once written keeps its storage. Reads to the CPU of an unwritten texture allocate it (zeros).
- Bindings only for what a program reads: a WGSL program's source may depend on the variant's defines
  (`WgslProgram.source` as a function of the define lookup), so a feature's helpers and the reads of its textures are
  left out of the source of variants that turn it off (the compose passes, the marking overlay, SSAO); defines inside
  the source are for small branches, as ported from GLSL. A headless check on Dawn (`webgpu-unused-bindings`) rebuilds
  every recorded pipeline without each of its bindings: one that still validates is not used by the shader. WebGL needs
  nothing: samplers of unused functions are inactive and get no texture unit.
- `getByteCount` counts storage only; `GraphicsStats` gains the allocated and reserved bytes.
- Why at the resource level, not per pass: the textures are referenced from many renderables' value cells at
  construction; creating them per pass on use means updating those cells and every pass's bookkeeping. Lazy storage
  keeps object identity, so the passes, the WebGPU bind group caches (keyed by texture version) and WebGL's attachments
  need no new references.

### Releasing (step 2)

- `PassTargets` (one per `Passes` and per `ImagePass`): the passes create their render targets, textures and
  renderbuffers through it (data textures, picking and Hi-Z are not), so it knows what it may release
  (`releaseStorage` on textures and renderbuffers: the storage goes, the size stays, the next write allocates again).
- What frames use is what their plan begins: a collector (`createTargetCollector`, a `PassBuilder` that only records
  the targets begun, by their attachments: `Framebuffer.getAttachments`) runs the plan of the current props and scene,
  hovering included (picking, marking), as `Canvas3D.prepare` plans it. Plans begin every target a frame draws into,
  also for a clear only (`PassBuilder`). The targets the pass set's render builders began since the change count too.
- When: a second after the last props change or scene warm-up (`scheduleRelease`), so that a slider or a sequence of
  changes releases once. Accumulating targets (multi-sample, tracing) of a mode turned off restart their accumulation
  when it is turned on again, as on any props change.
- No per-use stamps or time window: sampling through cached bind groups and redrawn WebGL items pass no hook, so a
  target only sampled (the image marking is redrawn over) would look unused; the plan names it.
- `ImagePass`: all its targets a second after its last render, so that kept image passes (the screenshot dialog's,
  the headless helper's) take no memory between renders and sequences of renders don't reallocate.

### Sharing intermediate targets (step 3: measured, not needed for now)

- Measured per frame (headless WebGPU, 1872×1536, written before read within the frame = shareable, the most of one
  descriptor alive at once = needed): defaults 252 MiB used, 22 MiB shareable (66 MiB if every pass that loads its
  attachment overwrote it instead of blending); with SSAO, outline, shadow and DOF 291 MiB, 22 (66); SMAA 252, 11 (55);
  illumination 384, 55 (99). That is a few full-size rgba8 intermediates (3–4 alive of 5–10), the draw pass's and
  WBOIT's depth-stencil buffers (never alive together) and, with illumination, rgba32f tracing targets.
- Not done: every pass would declare the last read of its intermediates (errors are images broken in subtle ways, worst
  for targets read in later frames: the image marking is redrawn over, accumulation, Hi-Z's depth), and WebGL can't
  alias storage between texture objects (its framebuffers attach texture objects), so sharing there means swapping the
  textures renderables hold per frame. The targeted items of step 4 save about as much on both backends. Revisit if
  memory at 4K (2.9× the measured) stays a problem.

### Fewer and smaller targets (step 4)

- One depth-stencil buffer, the transparent depth target's, for every pass that needs depth or stencil only while it
  draws: WBOIT and depth peeling (caps' stencil, impostors' depth), marking's depth and mask, illumination's thickness.
  Each clears depth when it begins and the caps clear stencil per cap; none is read after its pass, and the
  transparent depth pass that reads it within a frame clears it first.
- Accumulation in place: tracing blends each frame into its running average (color by the blend constant
  `1 / (iteration + 1)`, alpha replaced with that weight for the denoiser); the multi-sample compose targets
  (illumination's, temporal multi-sampling's) hold the output, `w * (s_1 + ... + s_k) + (1 - k * w) * s_0`, by adding
  `w * (s_k - s_0)` per sample: float targets, negative terms unclamped (an 8 bit compose target keeps the output
  target). On screen, the drawing buffer gets a copy of the compose target; an image pass blends marking in place, so
  it gets a copy of its own. No hold or output targets to recombine every frame.
- Passes only when they draw: bloom counts when something is emissive or in luminosity mode (`BloomPass.isActive`),
  the postprocessing pass when one of its effects draws (`PostprocessingPass.isEnabled(scene, props)`); bloom reads
  the transparent targets with blended transparency only (`dBlendTransparency`, as the compose).
- Only the channels and precision read: WBOIT's weights in r32f (read from `.r`), its accumulation in rgba16f as in
  the original method (half-float sums round faint layers: 400 overlapping spheres at alpha 0.1 differ by up to 6
  levels in 55 of 131k channels); WebGL1, without mixed attachment formats, keeps rgba32f. Depth peeling keeps one back
  color target for both ping-pong sides (the peel reads the previous depth and front color only, each back layer is
  blended back before the next peel), and geometry programs don't bind it. Marking's edges are computed in the overlay
  from the antialiased mask (the edge program read them at the same pixel) and its edges target is gone.
- Kept: tracing's average in rgba32f (fp16 loses precision over up to 2^16 iterations); the temporal hold target (the
  first sample, weighted until all are in); marking's layer in rgba16f.
- Dropped: depth without stencil where caps are absent (solid interiors are on by default for surfaces and impostors,
  caps also cover the near plane, impostors write depth, and switching formats by scene switches pipelines);
  `depth24plus-stencil8` (loses reversed-Z precision, Metal has no 24-bit depth).
- Illumination's normals octahedral-encoded in rgba8 (`packNormalToRGBA`: 2 × 12 bits over rgb, emissive in alpha;
  ≤ 0.06°) instead of rgba16f. Background pixels hold whatever the input pass cleared them to, so the trace gives
  background hits a zero normal (as a transparent background's clear): with an opaque background, the cleared color
  read as a normal lit bounces off the background (mean +0.1 levels once decoded differently). Against before: the
  noise of 8 iterations, converging (mean bias 0.02 levels at 128 iterations).

## Steps

1. Done: lazy storage in both backends; WebGPU programs bind only what their shaders read (`webgpu-lazy-storage`,
   `webgpu-unused-bindings`).
2. Done: `PassTargets`, the target collector, the release after changes (props, the canvas context's transparency) and
   after image renders (`webgpu-release-targets`, `webgpu-context-props`).
3. Measured, not needed for now (see step 3).
4. Done (see step 4), guarded by `webgpu-pass-target-bytes`.

## Decisions open

- None. Breaking changes are in [../breaking-v6-changes.md](../breaking-v6-changes.md) (Lazy texture storage, Pass
  targets).
