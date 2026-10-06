# WebGPU migration status

The native backend lives in `src/mol-gl/webgpu`. It uses WebGPU buffers,
textures, render pipelines and WGSL directly. It does not translate WebGL calls
or acquire a WebGL context. WebGPU is the default for the Viewer, plugin APIs,
examples and headless tools. The supported rendering paths and exports are
implemented and verified below. WebGL remains available only when explicitly
selected for legacy deployments.

## Running the current backend

```sh
bun install
bun run test:webgpu
bun run test:webgpu-viewer
bun run test:webgpu-headless
bun scripts/build.mjs -a viewer --prd
bun run serve
```

Open `http://localhost:1338/build/viewer/` in a WebGPU-capable
browser. The verification script uses an installed Google Chrome and runs actual
GPU rendering, picking and readback tests on localhost. It fails when WebGPU is
unavailable rather than silently skipping the GPU checks.

The headless check runs native Dawn/Metal rendering in Node and uses `ffprobe`
and `ffmpeg` to independently decode its MP4 output. Headless command-line tools
load the optional `webgpu`, `@napi-rs/canvas`, `pngjs` and `jpeg-js` modules;
install them with Bun when using the published library. `mvs-render` defaults to
WebGPU and accepts `--renderer webgl` for an explicit legacy override.

The Viewer and plugin APIs select WebGPU without configuration. To explicitly
select the legacy backend, use `?renderer=webgl`,
`Viewer.create(element, { renderingBackend: 'webgl' })`, or set
`PluginConfig.General.RenderingBackend` to `'webgl'` before initialization. `Canvas3DContext.fromCanvasWebGPU` is asynchronous;
`Canvas3D.create` accepts the initialized context. WebGPU Canvas3D instances expose
`webgpu` and have no `webgl` property. Molecular atom and bond representations
select sphere and cylinder data for the native backend. Gaussian density grids
and marching cubes use native WebGPU compute. Orbital and electron-density grids
also use native compute when the WebGPU context is supplied.

## Verified on this machine

- Theme dependency tracking skips per-vertex/per-instance expansion and GPU theme
  uploads for transform-only updates and unrelated rendering settings. Geometry,
  source themes, logical instance IDs, size, marker, material, alpha and overlay
  dependencies invalidate the cache independently. Only dependency and opacity
  metadata stays on the CPU, avoiding another retained copy of expanded theme
  arrays. GPU checks verify moving/restored geometry with unchanged theme build
  and upload counts; existing material, color, size, transparency, emissive,
  substance and image tests verify invalidation and opacity classification.
- Point/line exports retain native spatial colors, overpaint, transparency and
  instance transforms. GLB keeps compact point/line modes by default, with exactly
  the original positions and correct spatial strides. OBJ/USDZ/STL convert points
  to spheres and lines to cylinders instead of omitting them; GLB can request the
  same conversion through `setOptions`. Analytical gradient checks and independent
  archive/binary checks verify generated vertices, finite geometry, opacity and
  matching triangle counts in all four formats. Spatial lookup uses exported
  vertex indices independently of mappings used for source themes.
- Generated sphere/cylinder meshes export native spatial colors, overpaint and
  transparency in GLB, OBJ and USDZ; STL retains the same triangles. Interpolation
  uses generated vertices and each instance's transform, independently of the
  source-primitive mapping used for group themes. Analytical gradient references
  verify every generated GLB vertex and linear glTF color conversion. Independent
  archive/text/binary checks verify per-instance colors, opacity, finite geometry
  and matching triangle counts in all four formats. Instanced spheres now generate
  one copy of the source spheres per instance. OBJ/USDZ face transparency samples
  the indexed vertex, and color quantization counts faces rather than vertices.
- Mutable geometry, instance, theme, uniform, clipping and animation buffers use
  incremental queue writes when their existing capacity suffices. Validation
  completes before live buffers change, so rejected updates preserve their
  previous contents. Shrinking clears unused capacity, including packed byte
  channels whose shader bounds use the storage array length. GPU checks verify
  visible marker/material/transform changes without allocations, failed updates
  without committed writes, clipping-mask shrink/clear behavior, growth within
  capacity and replacement of only the buffer that outgrows its capacity.
- Geometry construction and vertex/index buffers survive marker, theme and
  material changes. Instance transform buffers update independently of geometry.
  Attribute/source revisions, cell identities and geometry-dependent flags
  invalidate cached geometry. GPU checks verify changed pixels, exact restoration,
  canonical picking and upload counts, including a replacement value box with
  an equal revision and recovery after a rejected spatial-theme update. Borrowed
  compute buffers retain their external ownership throughout rebuilds.
- Image/font, palette, group/value and spatial-theme textures update independently
  and survive unrelated geometry, transform, marker and material changes. Grid
  source revisions invalidate uploads even without a ValueCell update. GPU checks
  verify independent upload counts, green-to-blue grid changes, exact marker/frame
  restoration and canonical picking. A rejected rebuild after allocating a new
  palette leaves the previous textures usable; newly allocated resources are
  released, and borrowed compute textures remain externally owned.
- Density textures survive marker and material updates. Transfer and palette
  changes upload only their changed texture; density source revisions invalidate
  the grid even when the surrounding ValueCell is unchanged. GPU checks verify
  upload counts, changed selection and exact frame restoration after density
  and transfer reloads, and analytical palette colors. Value-box and source
  identities are tracked alongside revisions; replacing a transfer value box at
  an equal revision invalidates only that texture and restores the exact frame
  when reverted. Failed resource rebuilds
  release newly allocated textures while retaining the previous item's textures.
- Native stereo rendering retains asymmetric left/right camera projections and
  composes independent color and canonical picking outputs. Each eye owns its
  full or temporal multisample accumulation; GPU references cover both modes,
  odd/offset viewports, resizing, idle retention, focus changes and restoration
  to a single view. Picking snapshots the corresponding eye matrices before
  asynchronous readback and reconstructs molecular world positions correctly.
  The production settings panel verifies trusted stereo enable/disable clicks
  and two pickable protein views. Camera/clipping properties now report current
  camera state, and zero-radius clipping updates preserve the fitted view.
  Illumination and the image API retain their established single-camera behavior.
  Stereo currently recomputes AO per eye instead of sharing an eye's cached AO.
- Native camera orientation axes reuse the existing geometry, per-axis colors,
  labels and camera-axis loci. The overlay renders after molecular effects and
  before antialiasing, with its own depth and camera uniforms; it never enters
  illumination G-buffers. GPU checks cover all four corners, orientation, model
  scale independence, display pixel ratio changes, custom labels, highlighting,
  canonical picking and transparent screenshot exports. New native plugins show
  pickable axes even before a molecule loads. Canvas and image settings can
  disable axes independently.
  The production viewer also verifies a trusted X-axis click changes camera
  orientation through the existing plugin behavior.
- Illumination groundwork: native opaque G-buffers retain un-fogged shaded
  color, inward view-space normals, emission, raw albedo, density and nearest /
  farthest depth for thickness estimation. A WGSL compute port provides PCG
  hemisphere sampling, expanding ray steps, binary hit refinement, diffuse
  bounces, Russian roulette, directional soft shadows and float32 ping-pong
  accumulation. GPU tests verify isolated-surface lighting and emission,
  analytical accumulation weights, restart determinism, transparency exclusion,
  fog independence, resized / offset viewports, and green-to-red bounce-color
  transfer after a neighboring material update. Native normal-guided denoising
  matches independent Gaussian/color-weight references, including perpendicular
  normal boundaries and progressive threshold adjustment. Composition applies fog
  once, preserving transparent alpha and solid background colors. Canvas ticks
  now run progressive illumination, stop at the configured iteration count and
  restart after scene/settings updates. Screenshots inherit illumination settings
  and await every GPU iteration. Picking and the normal transparent/volume paths
  remain independent. Stable iterations reuse the native opaque and thickness
  inputs; fixed-thickness tracing skips the back-depth draw entirely. GPU pass
  counters verify one input refresh for four tracing dispatches and refreshes after
  material changes. The native frame-rate controller preserves the existing cheap
  first iteration, then adapts ray count, march steps and refinement within the
  configured bounds. Tests cover spare frame time, stalls, target changes, zero
  FPS targets and long pauses. Progressive supersampling covers all six jitter
  levels, canonical edge picking, camera offset restoration, cropped viewports
  and complete screenshot convergence against independently weighted GPU frames.
  Effect checks cover outlines and their suppression, bloom, DOF, FXAA/SMAA,
  sharpening and combined transparent screenshot exports. Opaque tracing skips
  duplicate SSAO and screen-space shadows while retaining transparent SSAO.
  The verified scenes cover the supported effect combinations; arbitrary custom
  combinations should still be checked when adding a new effect.
- Default Viewer and plugin initialization selects WebGPU without an override.
  The rebuilt production UI loads a protein and produces native GPU pixel
  readbacks without requesting any WebGL context. Explicit `?renderer=webgl`
  still initializes and loads the protein through the legacy renderer. Trusted
  UI clicks also verify the production illumination toggle, default 32-iteration
  convergence, supersampled convergence, molecular picking, visible protein pixels
  and zero GPU errors. The first camera fit waits for asynchronous molecular
  geometry; later temporarily empty geometry updates preserve the existing camera.
- Production trajectory controls verify next, previous and first model actions,
  forward/backward wraparound, visible geometry changes and picking into the
  displayed model. Returning to the first model restores the original pixels;
  stepping preserves the camera and state references. Trusted animation start
  and stop clicks cover every frame, and forced device loss during playback
  verifies recovery resumes the active trajectory. A stopped trajectory remains
  stationary through subsequent redraws. No WebGL context is requested.
- Production screenshot controls verify transparent PNG and opaque JPEG
  downloads at custom dimensions, including disabling auto-crop and changing
  the output format through the UI. Independently decoded images retain visible
  molecular pixels and PNG transparency. Reused image and preview passes refresh
  renderer settings and depth-peeling iterations; preview illumination also
  follows current settings. A second download verifies changed exposure without
  changing geometric coverage, with color-dependent antialiasing disabled for
  that comparison. Exports request no WebGL context and produce no GPU errors.
- Production CUBE loading and trusted density isovalue/visibility controls
  verify shrinking Gaussian surfaces, unchanged camera position/target, picking
  into the source volume, and exact restoration after returning to the original
  value or showing a hidden surface. Hidden surfaces retain no picking IDs.
  Numeric sliders now retain a pending typed value when Enter triggers a blur
  within the same React update batch; rapid isovalue entry no longer commits the
  previous value. The tests exercise that immediate entry path without delays.
- Trusted production controls switch the density from isosurface to direct volume,
  apply the representation update, set opacity to zero and restore it, and select
  byte, float and half-float density formats. Visible density cells remain pickable
  with a threshold appropriate for the faint default transfer ramp. Restoring
  opacity and the byte format restores the exact frame, and format changes preserve
  the camera. These checks request no WebGL context and produce no GPU errors.
- Production geometry-export controls save GLB, STL, OBJ and USDZ through
  trusted format-selection and Save clicks. Downloads are independently parsed:
  GLB headers and populated position/normal accessors, finite molecular vertices,
  packed STL records, OBJ triangles/materials and USDZ geometry/material scenes.
  ZIP archives are read from their central directories and decompressed with
  Node's zlib, independently of the exporter. GLB/STL/OBJ triangle counts agree.
  The molecular fixture includes surfaces and generated sphere/cylinder geometry;
  exports acquire no WebGL context and produce no GPU errors.
- Native asynchronous ray picking uses an independent narrow orthographic view
  and a small integer-ID/depth target. Concurrent queries retain separate camera
  snapshots and inspect padding in spiral order after a single readback. Geometry,
  volumes and world helpers use the native pipelines without replacing the
  displayed selection/color targets. Checks cover analytical straight/oblique
  triangle intersections, pickability, density-cell loci, camera scaling,
  forward molecular positions, misses, invalid zero-direction rays, concurrent
  requests, cached ray identify and unchanged main-camera/frame output. Pending
  or cached queries on a lost/disposed canvas cannot return stale hits.

- Native weighted mesh-color smoothing compute for Gaussian and molecular
  surfaces, at unit and whole-structure granularity. Ordered conservative bins
  preserve accumulation order; CPU-reference checks cover one-, three- and
  four-channel grids, logical instance reordering, transforms, sampling strides,
  empty cells and deterministic batches. Colors and synchronous overpaint,
  transparency, emissive and substance updates
  write native RGBA8 storage textures. The renderer and screenshot passes borrow
  those textures directly, preserving ownership across cache rebuilds and pass
  disposal, including allocation replacement without ValueCell updates; grid readback is needed only for verification. Rapid layer changes
  and clearing reproduce the original frame.
- Native spherical orbital and electron-density grids through angular momentum
  L=4, Gaussian/CCA/reversed-CCA coefficient ordering, contracted primitives,
  cutoff radii and fractional/zero occupancies. CPU comparisons cover isolated
  angular channels and their sums, single-primitive bases, batched readback,
  public grid tasks, signed isovalues and Bohr-to-Angstrom transforms. The orbital
  extension state transforms render both signed lobes and electron-density
  surfaces, update orbital indices, honor pickability, resolve volume loci and
  export screenshots without a WebGL context.
- Density vertex and vertex-instance color themes interpolate continuously in
  the position iterator's canonical grid order, independently of source tensor
  cell IDs. Vertex-instance overpaint uses the same trilinear sampling and honors
  logical instance IDs. Neighbors clamp at grid boundaries and exact voxel
  coordinates avoid division by zero. GPU checks compare a spatial color ramp
  against an analytical reference in blended, weighted and depth-peeling modes,
  retain density alpha, and verify transparent screenshot parity, independent
  instance colors, overpaint strides and picking after instance reordering.
- Transparent SSAO with separate depth normals and bilateral blur, alpha-weighted
  transparent occluders and scene transparency thresholds. Opaque and transparent
  color layers composite independently without changing alpha. Translucent
  protein checks cover disabling/enabling thresholds, multiple radii, blur,
  picking, screenshots and analytical colored-occlusion fog fading. Transparent
  color after SSAO is retained for subsequent outline and shadow composition.
- SSAO resolution scaling, including display pixel ratio, with reduced native
  occlusion targets and full-resolution alpha-preserving composition. Native
  compute builds opaque and transparent depth/alpha pyramids at full, half and
  quarter SSAO resolution; multi-scale sampling selects the appropriate depth
  level. Analytical GPU readbacks verify all mip values, odd dimensions,
  one-pixel dimensions and resize/reuse. Protein checks verify reduced-resolution
  shading and exact restoration when returning to full resolution.
- Postprocessing updates preserve camera snapshot fog and projection settings;
  camera properties update only when their corresponding canvas settings change.
- Native horizontal/radial gradients, image backgrounds and six-face cubemaps,
  including viewport/canvas coverage, aspect-cover cropping, ratio, opacity,
  saturation, lightness, rotation and GPU mip-level blur. URL and managed file
  assets load asynchronously and notify the canvas; screenshots await readiness.
  Native geometry and density fog fade into the environment while picking stays
  independent. Asset replacement cancels stale uploads and releases references.
  Analytical pixel checks cover gradients, image orientation/cropping, all cube
  directions, blur, opaque/transparent composition and source coverage. Shared
  asset tests cover redraw, screenshot disposal and clearing; invalid images and
  invalid cube dimensions fail without GPU validation errors. Colored SSAO over
  fogged environments is checked for both opaque and transparent proteins.
- Native marking outlines use separate unmarked-depth and marked-mask passes,
  with independent selection/highlight colors and strengths, hidden-edge opacity,
  inner contrast, display-scaled thickness and source-depth fog. Canvas and native
  screenshot exports share the pass, before depth of field and antialiasing.
  GPU checks cover analytical ghost/fog composition, marker clearing, highlight
  priority, offset viewports, screenshot resizing and unchanged picking.
- Separate opaque and nearest-transparent outline depth/opacity, with independent
  curvature checks, near-opaque suppression, layer dilation and fog. Transparent
  color is retained in a separate pass within the portable 32-byte MRT limit.
  Overlap checks verify opaque outlines behind transparent planes, nearest
  source alpha rather than accumulated alpha, and independence from selection
  thresholds, pickability and color-only flags. Ray-marched volume outlines and
  screenshot exports use the same separate depth pass.
- Native meshes, dynamic uniform colors, instanced transforms and logical instance IDs.
- Native canvas contexts use weighted blended transparency by default. Intersecting
  triangles and separate objects accumulate premultiplied color/emission and
  multiplicative revealage without depending on drawing order. Meshes, impostors,
  labels, images and density volumes share the portable 24-byte accumulation pass.
  Independent visible-color depth preserves volume bloom and postprocessing,
  while molecular selection retains its own opacity threshold. Checks verify
  analytic overlap color/alpha, opaque occlusion, nearest picking, live blended/
  weighted mode changes, supersampling, screenshots and device recovery.
- Native dual-depth peeling supports `dpoit` for geometry and ray-marched
  density volumes. Each iteration searches both nearest and farthest remaining
  depths with full-precision depth textures, accumulates front and back colors
  and emission, and resolves their coverage against opaque geometry. The existing
  1–10 iteration setting applies to canvas and image passes. Analytical checks
  cover six overlapping layers, reversed drawing order, coplanar ties, closely
  spaced surfaces, opaque occlusion, picking, supersampling and transparent
  exports. Density checks cover highlighting, emissive bloom and live mode
  changes. Headless unavailable-adapter recovery preserves the selected mode,
  iteration count, camera and byte-identical screenshot output.
- Physical GGX lighting with multiple colored directional lights, colored ambient
  illumination, exposure, metalness and roughness, including a numeric matte
  reflectance check. Cel shading and step counts, flat/flipped normals, front/back
  face culling, procedural bump frequency/amplitude and material bumpiness.
- Molecular surface material updates and native GGX screenshot exports; density
  gradients use the same physical material model while retaining their opacity.
- Normal-dependent x-ray and inverted x-ray opacity, edge falloff, interior
  colors and material strengths, and off/on/opaque transparent back-face modes.
  Real protein x-ray surfaces and custom-size transparent exports.
- Point sprites, spheres, thick lines, cylinders and SDF text labels.
- Molecular bond representations select native cylinder data on WebGPU. Dual
  endpoint themes support half-bond gradients, fixed dash blends and single-color
  modes; dynamic PDB representation updates and screenshot exports are verified.
- Molecular sphere visuals select native data for atoms, polymer backbone
  spheres and nucleotide atoms. Protein tests verify thickness-dependent alpha
  and exports, switching between explicit meshes and native spheres, and both
  per-unit and whole-structure atom loci.
- Solid sphere and cylinder interiors close camera near-plane cuts in perspective
  and orthographic views. Native patches preserve nonuniform atom scaling,
  interior colors/materials, transparency and independent picking. Tests verify
  near/far tracing depth, plane normals, screenshot exports and occlusion. Pixel
  clipping exposes the retained back intersection; uncapped bonds acquire solid
  interior ends without blending a second hidden layer.
- Plane, sphere, cube, cylinder and infinite-cone clip objects, inversion,
  transformed clipping coordinates, instance clipping and per-location masks.
- Native density volume clipping, including inverted and instance variants.
- Position/group wiggle, per-instance wiggle weights, instance tumble,
  animation disable, and continuous canvas drawing from its animation clock.
- Real protein cartoons, backbone, spacefill, Gaussian surfaces, molecular
  surfaces, lines and points, with molecular loci picking and visual inspection.
- Native Gaussian density and mesh extraction with analytical voxel/closest-atom validation,
  sparse subsets, radius offsets, smoothness, padded physical grids, empty
  selections and deterministic ties. Bounded output batches preserve all voxels
  and IDs. Protein surface generation dispatches native compute, honors the
  CPU/GPU setting and reproduces its frame after toggling back to WebGPU.
  Unit and whole-structure mesh/wireframe variants retain atom loci. Filtered
  density inputs exclude one-past-end indices, and zero-radius atoms contribute
  no density.
- Visibility, integer object/instance/group picking, and depth readback.
- Trusted Chrome mouse input at 1x and 2x display scaling: rotation, panning,
  wheel zoom, asynchronous atom hover, single-atom selection and toggling,
  deselection on empty backgrounds, visible marking and marker removal.
  Interaction screenshots are saved under `tmp/webgpu/interaction-1x.png` and
  `tmp/webgpu/interaction-2x.png` by `bun run test:webgpu`.
- Independent selection depth and opacity thresholds, including selection through
  faint, unpickable and color-only foreground objects; x-ray opacity and nearest
  transparent triangles within one object. Faint density volumes retain their
  colors while obeying the selection threshold.
- Representation opacity updates use the canonical alpha value and alpha factor.
  Per-fragment opaque/transparent surface classification prevents opaque regions
  in mixed-opacity objects from blending with hidden transparent fragments.
- Transparent spheres support radius-dependent alpha thickness, including dynamic
  updates, clamping and theme-alpha selection thresholds.
- Native substance overlays support instance, group-instance and vertex-instance
  material data, pre-mixing, interpolated strength, interior materials and GGX
  shading. Spatial substance grids use native RGBA textures and WGSL slice
  interpolation in transformed molecular coordinates. Mesh smoothing builds
  the grid on the CPU without WebGL, preserving reordered instance IDs;
  native unit/whole-structure surface updates, screenshots and layer clearing
  are verified. Group-instance grid accumulation uses the native WebGPU
  weighted compute path; the compatibility API may still read the resulting
  packed grid back when a CPU-owned texture is explicitly requested.
- Smoothed overpaint, transparency and emission use native spatial grids as
  well. Scalar grids are packed into RGBA alpha channels; overpaint retains
  pre-mixing and strength before fragment blending. Combined overlays match
  sampled vertex references at multiple camera scales, preserve selection at
  the configured opacity threshold and restore the frame when cleared.
  Spatial transparency participates in both opaque/transparent surface passes.
  Unit/whole-structure protein overlays and transparent exports are verified.
- Invariant and instance color grids use native WGSL interpolation, retaining
  molecular coordinates across camera scales. Group overpaint is pre-mixed
  after the color-grid lookup; spatial palettes preserve nearest/linear
  filtering in the fragment shader. Gaussian and molecular mesh surfaces build
  RGB grids without WebGL, render element colors and screenshots, and reproduce
  their frames after toggling smoothing. Grid accumulation uses the native
  weighted WebGPU compute path and preserves ordered sample bins and logical
  instance IDs; CPU accumulation remains only as the legacy WebGL fallback.
- Ordinary geometry retains encoded palette values through vertex interpolation
  and resolves nearest/linear palette colors in the fragment shader. Analytical
  barycentric pixel references verify the interpolation order. Group overpaint
  uses a separate raw overlay buffer, preserving pre-mixing before fragment
  blending, including partial alpha/strength and reordered instance IDs.
- Slice/image palette lookup honors nearest/linear filtering, with discrete
  opaque pixel checks, smoothly interpolated values, cell picking, exports and
  frame restoration. Direct volume palettes honor the same filtering contract;
  constant-density fixtures match analytical uniform-color references in both
  canvas rendering and screenshots, preserving volume-cell picking.
- Canvas resize, transparent backgrounds and RGBA image readback.
- Plugin initialization, PDB parsing and ball-and-stick representation generation.
- Camera-reset snapshot callbacks, current scene bounds, custom orientation,
  automatic fitting with partial snapshots, and configured transition defaults.
- Mapping picked pixels back to molecular loci.
- Custom-size transparent screenshot export, including native outlines and FXAA.
- Native FXAA edge smoothing, depth outlines, optional transparent outlining,
  disabling effects, and preservation of the original picking attachments.
  Transparent outlines use background-independent geometry coverage retained
  in the emissive attachment alpha, even when transparent bloom is excluded.
  Outline alpha doubles/clamps the source opacity, then composes over the
  configured background. Multiple opacity checks and molecular screenshot
  exports match premultiplied compositing; outline-only pixels remain unpickable.
  Outline edges retain source view distance and opaque/transparent classification.
  Fog fades their opacity using the camera smoothstep, with analytical pixel-alpha
  checks for both surface classes and exact restoration when fog is disabled.
- Native three-pass SMAA using the existing area/search lookup tables, color
  contrast thresholds, bounded searches, gamma-correct transparent coverage,
  offset viewports, molecular lines and custom-size screenshot exports.
- Native emissive attachments for meshes, images, labels and density volumes,
  with material lighting-independent emission and fog fading.
  Disabled geometry emission overlays ignore cached data while retaining uniform
  material emission. Enable/disable/re-enable checks preserve frames and picking
  without clearing the overlay texture.
- Five-level Gaussian bloom in luminosity and emissive modes, threshold,
  radius/strength controls, optional transparent emission, halo alpha coverage,
  picking preservation, and emissive volume screenshot export.
- Native contrast-adaptive sharpening (RCAS), including denoise controls.
- Native planar and spherical depth of field, camera-target and scene-center
  focus references, zero focus-range handling, and transparent image export.
- Opaque SSAO with the existing hemisphere sample generator, reconstructed view
  normals, depth-aware separable blur, multi-scale radii, and alpha preservation.
- Native directional screen-space shadows with multiple colored lights, ambient
  lighting, orthographic depth reconstruction, distance/tolerance/step controls,
  fog-aware composition, transparent backgrounds, picking and screenshot export.
- Native direct volume ray marching for byte, float and half-float density data.
- Volume-cell picking with original tensor indexing, highlighting, instancing,
  opaque geometry occlusion, and offset viewport ray reconstruction.
- Grid and arbitrary-plane slices, all four interpolation modes, palette decoding,
  per-cell picking and highlighting, grid-space trimming, and iso-value masking.
- Native texture-mesh positions/normals, packed byte/float group IDs, transformed
  position-theme locations, group colors, instance picking and geometry updates.
- Camera scale in perspective/orthographic geometry, text labels and density
  volumes; preserved picking IDs and molecular positions; scaled screenshots
  and head-rotated lighting exports.
- Native GPU marching cubes for float32 volume isosurfaces, including periodic wrapping, all 256
  cube cases at positive/negative iso-levels, interpolated CPU-reference normals,
  multi-cell topology, tensor axis orders, per-edge atom IDs against CPU conventions
  and exact batch equivalence. Missing/ignored IDs retain their sentinel behavior. GPU
  prefix scans compact vertices in order, including partial workgroups and
  empty cells; sparse-fixture readback shrinks by more than half. Volume
  surfaces retain loci picking, iso-value updates and CPU/GPU setting changes.
  Periodic boundaries preserve CPU positions/normals in all six axis orders and
  original-grid cell IDs. Floodfilled tensor views honor their custom sample
  getters; CPU fallback IDs use the original grid accessor. Empty surface updates
  preserve the camera clipping radius so returning geometry restores the view.
- Cropped marching-cubes regions preserve global positions, full-grid gradients
  and source-cell IDs. Empty regions dispatch no GPU work; invalid regions fail
  before access. Segmented volume surfaces use native extraction, with segment
  colors, original-grid voxel loci, selection controls and CPU/GPU switching.
  Segment masks zero-pad outside every source-grid face, avoiding row aliasing
  in cropped extraction. Boundary-touching segments close within half a voxel,
  retain valid source-cell IDs and match CPU fallback triangle counts.
- Native handle and pointer overlays reuse CPU geometry with native color/depth
  passes. Handles retain helper loci, transformed picking, highlighting and
  scene occlusion; translation survives rotation and display scaling. Pointer
  themes refresh on color changes, retain display sizing under model scaling,
  blend transparently and preserve scene picking. Plugin checks cover settings,
  global highlighting alongside camera axes, disabling and screenshot ownership.
  Helper overlays are excluded from illumination geometry.
- Native debug-helper registry and CPU scene ownership for all five extension
  helpers: scene/visible/object/instance bounding spheres, clipping shapes and
  indicators, mesh normals, image/trim edges and direct-volume boxes. Plugin
  GPU checks verify visible overlays, screenshot inclusion and borrowed scene
  ownership, clearing/rebuilding, disabling, unchanged picking and exclusion
  from illumination inputs. Normal geometry refreshes after transform changes
  and uses inverse-transpose directions; unbounded clipping helpers resize with
  scene bounds, and cached sphere geometry refreshes its bounds after updates.
- DOM-free native WebGPU context and asynchronous `HeadlessPluginContext.create`
  / `HeadlessScreenshotHelper.create` factories, defaulting to WebGPU with a
  caller-provided GPU implementation. Native Metal checks in Node cover PDB
  representations, raw RGBA/PNG/JPEG, custom and odd dimensions, transparent alpha,
  cropped images, full supersampling and SMAA/FXAA. SMAA uploads lossless CPU
  lookup bytes rather than using browser image decoding; a unit test compares
  every byte with the shared PNG tables. Headless canvas initialization resolves
  plugin readiness and exposes context ownership for cleanup.
- WebGPU defaults for the MVS-render CLI and image/GLB example tools. The MVS
  CLI is tested with local CIF input, native protein rendering, font-backed
  labels, PNG/JPEG output and Mol* state snapshots. Two native molecular
  snapshots export a seven-frame MP4; independent decoding confirms both colors
  and correct 160×128 cropping for a requested 161×129 frame. Headless animation
  checks use the registered MP4 transformer identifier and preserve the encoder's
  crop. CPU-owned spheres/cylinders export populated GLB meshes without WebGL;
  native texture meshes and spatial theme grids also export without WebGL,
  including readback of GPU-generated float geometry textures.
- Thirty-four geometry/volume/spatial-bin/segment/workload/stereo/debug/lookup unit tests, including reordered instances, sphere distance subsets, overpaint,
  transparency, packed point sizes, constant-density normalization and invalid
  input rejection.

The GPU suite requires zero uncaptured WebGPU errors. Screenshot artifacts are
written under `tmp/webgpu/`. Run `bun run jest --runInBand` for the existing library
tests; legacy headless GL tests are skipped if the optional `gl` module is absent.

Repository checks at this checkpoint: lint and TypeScript passed; the production
viewer and both library builds passed; 141 test suites / 1,524 tests passed, with
11 suites / 14 tests skipped because optional headless GL support is absent.

Native multi-sampling now supports `on` and progressive `temporal` modes at all
six sample levels (1–32 canonical jitter samples), premultiplied float16
accumulation, and full sampling for screenshot exports. Picking IDs and depth
remain from the unjittered frame, and camera view offsets are restored after each
batch. Temporal marker changes finish in one frame when `reduceFlicker` is enabled.
`reuseOcclusion` retains the baseline SSAO texture and shifts its sampling coordinates
for each jitter, avoiding repeated SSAO and depth-pyramid work. Browser checks
compare all sample levels with independent weighted GPU readbacks and verify
convergence, idle retention, marker flicker controls, cropped cameras, screenshot
alpha, and the actual number of SSAO passes with reuse enabled and disabled.

Native geometry exports now read CPU-owned texture meshes and GPU-owned RGBA8
spatial grids without a WebGL framebuffer. DOM-free tests parse GLB vertex/color
accessors and alpha materials, OBJ materials, USDZ scene text, and binary STL
triangle counts. They check packed byte and normalized float group IDs, padding,
per-instance spatial colors, overpaint and transparency. Grid copies discard
WebGPU row alignment padding. Transparency overlays enable GLB alpha blending;
STL records now count triangles and pack one 50-byte record per triangle.

Native Gaussian surfaces now produce RGBA32F position/normal textures, packed
RGBA8 groups and storage vertices in compute shaders. The renderer binds these
vertices directly, including after theme/overlay changes and in screenshot
renderers. Compatibility mirrors preserve existing synchronous theme and location
APIs. Browser checks verify empty and populated surfaces, affine positions,
inverse-transpose normals, texture/mirror equality, direct buffer identity,
CPU/GPU/parent setting transitions, molecular loci and all smoothed material
layers. Headless geometry exports independently read the GPU float textures.
Spatial export interpolation also preserves exact integer grid samples.

Headless background images and cubemaps use the configured native canvas module
to decode RGBA pixels and upload them directly to WebGPU. Browser assets retain
ImageBitmap uploads. DOM-free checks verify image opacity in the actual molecular
screenshot path, all six cubemap face directions, odd dimensions and mipmap blur.

Forced device-loss checks now verify automatic recovery in browser and headless
plugins. Recovery recreates native GPU representations and renderers, retains CPU
molecular data and state refs, restores selection/camera/settings/manual overlays,
and replaces viewport event subscriptions and debug helper registrations. Browser
checks verify the original frame and molecular loci; headless checks verify
byte-identical exported surfaces with spatial overpaint and transparency. Production
UI checks continue camera-axis, stereo and illumination interactions after recovery.
The animation clock resumes from its prior elapsed time. Screenshot calls await an
ongoing headless recovery. Initialization failures remain visible in the plugin log;
`recoverWebGPU()` allows an explicit retry. A controlled unavailable-adapter test
verifies that retry retains the animation clock, camera, molecular state and
spatial layers, and reproduces the pre-loss screenshot byte for byte.
Changing unrelated canvas settings preserves the trackball's adjusted distance
limits. Recovery restores the saved camera after restarting the controls, so
their initial distance enforcement cannot change the canonical recovered frame.

## Known boundaries

- WebGPU readback is asynchronous. `asyncIdentify` provides fresh hover/click
  results; synchronous `identify` deliberately returns only a matching completed
  result so it never blocks the browser event loop or returns a stale hit for a
  different target.
- WebXR is reported as unavailable by the native backend. It does not acquire a
  WebGL context as a fallback.
- The remaining limitations are performance or precision trade-offs: general occlusion culling and device-limit-aware sharing for very large
  assemblies still use their existing safe paths. Numeric scalar fields are accepted
  by native marching cubes through validated f32 conversion; values outside f32
  precision remain a documented precision trade-off.
- Unsupported theme granularities fail explicitly instead of silently falling back
  to WebGL. Directional lighting supports up to 1,024 lights per frame within the
  portable uniform-buffer budget.

## Current implementation limits

Representation distance bounds and overlap fading run in native geometry shaders.
Color and tracing passes use the molecular four-by-four ordered coverage mask;
picking and marking passes apply the configured opacity threshold. Camera scale,
offset viewports and signed sampling jitter are included. GPU checks compare
every covered triangle and sphere pixel with an independent mask; spheres use
their center distance rather than the tessellated surface distance. Checks test both overlap edges
and zero-width bounds, and verify selection thresholds and exact restoration
without rebuilding theme data.

Native spheres use twelve bounding/cap vertices and forty-two indices per sphere,
with WGSL ray intersections providing surface normals and exact fragment depth.
Color, picking, transparency, outlines, marking and illumination passes use that
depth. Spatial themes sample the intersection point. This replaces the previous
225-vertex / 1,158-index sphere meshes; geometry storage drops by approximately
95 percent. GPU checks unproject selected pixels onto the analytical sphere,
including nonuniform ellipsoids in perspective and orthographic views with
reordered instance colors, and retain coverage, distance fading, near-plane caps
and clipping checks.
Device-loss recovery also preserves canvas representation membership: state
representations previously removed from the canvas are rebuilt for future use
without becoming visible again, including after a failed-adapter retry.

Sphere distance levels now retain the prepared position/group prefixes and
radius scales. Native draw ranges cull out-of-range levels; vertex culling and
smooth radius overlap transitions run on the GPU. Browser checks compare near
and coarse subsets, preserve molecular picking, unproject scaled radii, verify
perspective/orthographic and runtime camera scaling, measure halved submitted
indices for a stride-two level, and check exact restoration after camera and
level-setting changes. Sphere instance-grid culling now coalesces physical draw
ranges while preserving logical molecular IDs, with safe per-instance fallback
for animated, sheared, stale or incomplete grid metadata. General occlusion culling remains performance work. Cylinders use native
triangle meshes with deterministic radius-aware tessellation.
It expands color/size/marker data per vertex and instance. Geometry construction,
vertex/index buffers and instance transform buffers are retained until their
inputs change; mutable geometry/theme/instance buffers update in place when
their capacity suffices. Theme data is expanded on the CPU when its dependencies
change, and retained GPU theme buffers skip expansion on unrelated updates. Direct volumes retain unchanged density,
transfer and palette textures, and other geometry retains unchanged image/font,
palette, group/value and spatial-theme textures. Large assemblies still benefit from buffer sharing, incremental updates and
device-limit-aware batching. Unsupported geometry and theme granularities fail
explicitly. Density uploads currently use RGBA32Float with manual trilinear
sampling; image overlay data is expanded per cell and instance.

Directional lighting currently supports up to 1,024 lights per frame and rejects
larger configurations explicitly. Light uniforms leave the minimum eight storage
buffer bindings available to native volume themes.

Gaussian density compute uses conservative spatial bins with atom lists in input
order, preserving density sums and closest-atom tie-breaking. Bins coarsen to fit
device storage limits. Browser checks compare every voxel and atom ID against a
full atom scan, including separated atoms, shifted coordinates and batch
boundaries; a dispersed unit fixture eliminates more than three quarters of the
atom evaluations. Gaussian meshes and wireframes use native GPU marching cubes and compaction,
including parent-inclusive surfaces. The packed Gaussian density buffer is consumed directly by
the marching-cubes shader when no flood-fill view is requested, avoiding a second scalar-field
upload; edge smoothing and line conversion keep their compatibility CPU steps. Shared mesh
topology across related representations remains a performance opportunity for larger assemblies. Its WGSL
exponential matches the existing GPU density kernel; the CPU path uses a fast
approximation. `WebGPUContext.stats.computeDispatches` counts submitted compute
dispatches.

Camera snapshot callbacks receive `Camera.SnapshotScene`, containing
`boundingSphere` and `boundingSphereVisible`, in both backends. Callbacks that
previously accessed WebGL renderables through the scene must use representation
data instead; camera fitting no longer depends on WebGL resources.

Native marching cubes uses a stable GPU prefix scan to compact triangle batches
before readback and mesh upload. Only a vertex count and populated vertex records
are downloaded. Volume isosurfaces use it for float32 grids, including wrapped and floodfilled views,
that fit device storage limits; other grids retain CPU extraction. Native Gaussian
texture meshes retain the compacted GPU records, convert them to storage vertices
and geometry textures, and bind those vertices directly for rendering. One final
CPU mirror supports synchronous themes and location APIs. Parent edge smoothing
and wireframe line conversion retain their compatibility CPU steps.
`WebGPUContext.stats.marchingCubesDispatches` distinguishes extraction and
compaction dispatches from Gaussian density computation.
