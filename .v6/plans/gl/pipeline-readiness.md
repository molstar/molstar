# Plan: Pipeline readiness without rendering

## Status

- **Current**: steps 1–5 done; 6 (snapshot replay with precompile lists) done (6a–6g), 6h (WebGL's first draws into
  several targets, an ANGLE compile once per browser) analysed and accepted. Step 7 (precompiling when a
  representation transform applies or updates): 7a measured (adding: 60–190 ms after the commit to hide, up to 1.1 s;
  updates: frames to hold rather than compiles to start earlier, done); 7b and 7c done: representations prepared with
  their data while they build; 7d (particle targets) done and measured: a mesoscale load compiles nothing after its
  commit; what it waits for there is the commit's own work.
- Step 6f done, measured in Chrome (RTX 4070 Ti, 1280x800, the Viewer, `kinase-story` and a fade story of three
  snapshots), loading the story and playing it right away:

  | backend, story          | precompile (background)                   | playing compiles            | playing without precompile   |
  | ----------------------- | ----------------------------------------- | --------------------------- | ---------------------------- |
  | WebGPU, fade            | 16 signatures, 65 programs, 1.6 s         | nothing (WBOIT and blended) | 8 programs, 10 pipelines     |
  | WebGPU, kinase          | 71 signatures, 211 programs, 7.2–8.0 s    | nothing (WBOIT and blended) | 13 programs, 37 pipelines    |
  | WebGL2, fade            | 16 signatures, 60 programs, 1.5 s         | 2 programs (`tracing`)      | 12 programs                  |
  | WebGL2, kinase          | 71 signatures, 131 programs, 3.5 s        | 4 programs (`tracing`, ...) | 13 programs                  |

  Loading takes longer (WebGPU kinase ready after 9.6 s instead of 2.9 s), the first snapshot shows meanwhile. Fixed on
  the way: precompile program items created their define values in the signature's order and the context defines
  last, which WebGL keys programs by (`resources.program`: hashes in order), so they were other programs than the
  objects' (now in schema order, unit test); the stand-in scenes take both sides of what the scene averages decide
  exclusively (opacity average 1 or below: postprocessing blends transparency; emissive average 0 or above: bloom in
  emissive mode composites the background, which selects the blended pass's state), before only the transparent and
  emissive side. **Not fixed**: WebGL with WBOIT still freezes the first time impostors draw (fade 694 ms, kinase #7
  689 ms; 657 and 660 ms without precompile); with blended transparency it doesn't (29 ms, also without precompile, 40
  ms). So linking ahead isn't what is missing: ANGLE (D3D11) builds a shader for the draw's render target layout when
  it first draws (WBOIT: two float targets), see 6h.
- Step 6e done: `precompileSnapshots(plugin, snapshots)` (`plugin/state/precompile.ts`) derives the snapshots' list
  and precompiles it for their props, each distinct props once (a story's snapshots usually share theirs), in the
  background; called when the snapshot manager loads a state (`setStateSnapshot`, which `.molj` and `.molx` files go
  through) before the current snapshot is applied, and when an MVS state or story is loaded (`loadMVS`, after the
  snapshots are converted, before the first is applied). A file of a single snapshot (no entries) is applied at once,
  which the warm-up after the commit covers. Test on Dawn (`webgpu-precompile-snapshots.test.ts`): a state of two
  snapshots (a cartoon; then surfaces, transparent ball-and-stick and overpainted spacefill), as read from a file,
  loaded: applying the second snapshot compiles no render pipeline; with the first applied without loading the state it
  does.
- Step 6d done: lists from params alone, generic rather than declared per representation type:
  - **Visuals** (`repr/signatures.ts`, `getGeometrySignatures`): the define values of an empty geometry of the visual's
    kind with its props (`createValuesSimple`), crossed with every value of what data decides: the color types the
    theme may give (its granularity, 'instance' also as 'group'), smoothing (a grid with a context, per-vertex colors
    without), palette, size types, marker types, the layers set on the representation (overpaint, transparency,
    emissive, substance, clipping, wiggle) and their types; linear textures by the schema (`getSchemaLinearTextures`,
    a texture spec's filter and `usedIf`); texture meshes in chunks where compute isosurfaces are supported, direct
    volume grids as 3d textures where supported. A visual adds what it decides (`signatureDefines`: surfaces' varying
    groups by flat interpolation, the Gaussian density's packed group, a volume's six axis orders).
  - **Representations** (`Representation.getSignatures(props, ctx)`): the visuals the props select, created without
    data (cheap: closures); factories that choose name the visuals they may create (`selectVisual`, `chooseVisual`:
    Gaussian surfaces and volume isosurfaces and segments, GPU or CPU by size, both when the props and context allow
    the GPU; impostors or meshes by the props and context, the one chosen); the quality presets as
    representations apply them (`getQualitySignatures`: 'auto' quality as medium and low, which differ in
    double-sidedness, i.e. pipelines only); units representations also merged where `renderMerged` is on, or
    `gpuCulling` (an object alone merged, `gpuCulled`: the culling's kernels prepared with it; `withMergedSignatures`,
    with the merged programs' `dSegmented`); 'auto' of either is taken as unmerged and culled on the CPU, by the
    objects' counts (large objects, which compile while they build).
  - **Themes**: created without data, as the theme legends are (`getColorThemeSignature`); themes whose data decides
    the granularity name what it may be (`ThemeProvider.getSignatureGranularities`: the MVS multilayer theme, its layer
    themes with any selection granularity).
  - **State trees** (`plugin/state/precompile.ts`): transformers register how their signatures follow from their
    params (`registerTransformSignatures`: structure, volume, particle and shape representations, the measurements,
    unit cells, bounding boxes, MVS primitives), transformers that set layers register them (`registerTransformLayer`)
    and shape providers their geometries (`registerShapeGeometry`: boxes, MVS primitives); `getTreePrecompileList`,
    `getSnapshotsPrecompileList` (each snapshot's tree, and its canvas props with the context's transparency, for
    `Canvas3D.precompile`). Transforms without registration (e.g. extensions' representations) compile when first
    drawn.
  - **Resolved, not over-approximated**: reflection in transforms as none (pipelines only). The instance granularity
    needs no resolution, a uniform since (fewer-defines.md, candidate 2).
  - Signatures leave out the context's defines (light count, color marker), which the preparing context sets, so that
    recorded and derived lists compare. `getRenderObjectSignature` signs render objects without a context.
  - **Tests**: on Dawn (`webgpu-precompile-derive.test.ts`), every structure representation type with its defaults
    and with props, themes and layers on a protein with a ligand and on DNA: the derived list covers the recorded one
    (17 signatures; 44 derived, 20 of them merged variants), and in a new context, precompiling the derived list first,
    building the scene and preparing its frames and first interaction compiles no render pipeline (the GPU surfaces'
    kernels are their own). Without a context (`mvs/runtime/_test/precompile.test.ts`), `kinase-story` with its
    structures read from a local file: the list derived from the 11 snapshots, before they are applied, covers what
    playing them draws (7 signatures; 30 derived). Found on the way: signature stand-ins set `uDoubleSided` for kinds
    without it, which the object state reads apart from false (points and images got other pipelines).
  - Measured on Dawn (Node, D3D12): precompiling the 44 derived signatures took 38–41 s for 501 pipelines, compiled
    one after the other in the background; Chrome compiles in parallel (26 pipelines in 440 ms when 1grm loads, step
    5's timing); to be measured in the browser (6f).
- Step 6c done: recording in `Canvas3D` (only renderables know their kind): the signatures of the scene's renderables
  when a commit finishes and after `update` (values may change defines: `update(repr)` records that representation's
  renderables again, `update()` the scene's), deduplicated (`PrecompileRecorder`),
  `getPrecompileList()` and `resetPrecompileList()`; `PrecompileList` (`objects`, not stored, so without a version).
  `Canvas3D.precompile(list, props[])`: per props (the context's transparency mode
  among them, the current by default) the mode as the device draws it (`DrawPass.resolveTransparency`: blended where
  WBOIT or DPOIT are not supported), the scene's globals as `setProps` sets them, program items of the signatures
  (kept until the next call, reused across props), then the passes as `prepare` with the interaction plans them, under
  `DrawPass.withTransparency` (plans for another mode than the current, restored after) and against a stand-in scene
  (groups with the counts of the list's primitives and volumes, averages taking every branch that draws objects:
  opaque and transparent, emissive); the prepare builder takes a `collect` that answers the scene's groups (also the
  real ones, which the pick helper plans with) with `renderer.collectSignatures`. Test on Dawn: a list recorded in one
  context, as JSON, precompiled in another for its mode (WBOIT) and blended; adding the same representations and
  preparing their frames and first interaction compiles nothing, without the list it does. Found on the way: adding
  objects to a second scene with other globals compiles their family's configurations for the new programs (family key
  reuse), which the first test's own scene did. Not covered by a test: switching to the other mode after precompiling
  for it (the headless context has no transparency setter; 6f plays it).
- Step 6b done: `ObjectSignature` (`gl/precompile.ts`: render object kind and whether merged, define values without
  the render variant, linearly filtered textures, double-sided, reflection), `getObjectSignature(renderable)`,
  `getSignatureKey`; each kind's shader code and composed schema registered (`getRenderableProgramSource`, the
  renderables share one schema object per kind, merged per type); `GraphicsContext.createProgramItem` (WebGL: programs
  as render items create them, released by `destroy`; WebGPU: a pipeline preparer with the render items' program
  creation, `createVariantProgram`) and `createSignatureProgramItem(graphics, signature, transparency, globals)` (the
  context's own defines, light count and color marker, from `globals`); `PipelineRequest.item` takes program items;
  `RenderItem.getLinearTextures`; `renderer.collectSignatures(pass, items, target)` through stand-ins of the objects
  (the values the pass draws and object states read). Test on Dawn (`webgpu-precompile.test.ts`): signatures of a scene
  (impostors, meshes, text; opaque, transparent, double-sided) prepare exactly the objects' programs and pipelines in
  the same context, and in a new context prepared from the signatures (as plain data) before the objects exist,
  preparing the objects compiles nothing more into the target, which without the signatures it does. WebGL's program
  items are not covered by a test yet (no headless WebGL in the native tests' runs here).
- Step 5: `ComputeContext.prepareKernels(requests)` (`KernelRequest`: code, schema, define values) compiles kernel
  pipelines in the background (`createComputePipelineAsync`, tracked as pending work, timing logged) into the cache
  `createKernel` uses; a kernel not prepared still creates its pipeline synchronously, which returns at once in Chrome but
  blocks the queue, and frames with it, until compiled. The modules name their kernels (`GaussianDensityKernels`,
  `PrefixSumKernels`, `getIsosurfaceKernels(input)`, from the same input defines the kernels are created with), and the
  asynchronous paths wait for them: `computeIsosurface` (Gaussian surface, volume isosurface, segments) and the GPU
  Gaussian surface before its density. Measured in the Viewer (WebGPU, 1grm, adding a GPU Gaussian surface): the
  longest frame gap while adding the first one from 90–150 ms to 34–38 ms (the kernels compiled 13–52 ms each, in
  parallel, while frames went on); the following ones unchanged (18–33 ms, run to run). Also prepared: the grid
  compute kernels (`createGrid3dComputeRenderable().run` waits for its kernel, `kernels` exposed), and the color
  smoothing kernels, which the visual's synchronous update dispatches: the GPU Gaussian surface prepares them with its
  own when color smoothing is not off (both color types, the type not known yet); turning smoothing on later creates
  them synchronously once. The Gaussian surface prepares all its kernels before computing the density: the density
  texture is shared by the surfaces, and an `await` between computing it and extracting the isosurface let another
  surface write it meanwhile (the WebGPU visual cases of the GPU Gaussian surface differed in 13,337 pixels);
  `computeIsosurface` awaits only kernels not ready (`ComputeContext.hasKernels`).
- Frames after a props change wait for what they draw with (`holdFramesUntilPrepared` in `Canvas3D.setProps`): the
  props are prepared, and when that started work (pipelines compiling, programs linking) frames are skipped until it is
  ready, the last frame stays shown. Before, switching to illumination for the first time showed frames of only the
  background while its pipelines compiled (the Viewer smoke test reads back every frame after the switch: 13 of 50 such
  frames without the hold, none with it).
- Timing (timing mode): WebGPU logs each batch of pipelines compiling in the background together (count, wall time from
  the first request to the last ready, the longest with its label; e.g. loading 1grm in the Viewer: 26 pipelines in
  440 ms, the longest 336 ms) and the synchronous creations (the fallback after a failed asynchronous compile, kernel
  pipelines), as WebGL logs finalizing programs.
- Step 4: `Canvas3D.prepare({ props?, interaction? })` makes what frames with the props draw ready without drawing: the
  variants the props require linked for every render item (`ShaderManager.getRequired`, as setting the props links them),
  then the frame's passes as `render` chooses them (illumination, else multi-sample or the draw pass, and Hi-Z) planned
  with the prepare builder; with `interaction` also the first hover: the pick and marking variants linked, the pick pass
  (`PickHelper.plan`) and marking as if objects were marked (`DrawPass.plan` with `asIfMarked`, the marking passes
  collected with `all`). The idle warm-up after scene changes is `prepare()`, then `prepare({ interaction: true })` once
  the frame's pipelines are ready: requested together, the pick and marking pipelines share the compiler's threads with
  the frame's, which the first frame waits for (WebGPU, Viewer loading 1cbs: first draw 260–280 ms after the state
  update instead of 350–390 ms, 20 pipelines compiling at once before; WebGL2 unchanged, its programs from the browser's
  cache). Precompiling a state update's tree (`getTreePrecompileList`) when the update starts was tried and dropped: the
  representations' transforms arrive in the update's last tree, about 30 ms before the commit for 1cbs, and the list's
  stand-in scenes compiled 50 pipelines instead of the 32 the load needs. Removed: collect mode
  (`resources.preparePrograms(run)`, WebGPU's preparing state), `DrawPass.prepareMarking`, `ShaderManager.warmUp`,
  `PickHelper.warmUp`. Pick, ray and Hi-Z passes are plans too (Hi-Z prepares its first level, from the scene depth, and
  one following level). Linking matters on WebGPU beyond the drawn objects: setting props that require a variant links
  it for every render item, including helpers that never draw with it, and family key reuse then compiles their
  pipelines; without linking ahead, the first hover compiled the camera helper's marking pipeline. Checks: the Viewer
  smoke test's first hover compiles no pipeline (as before), and new: props prepared with `prepare({ props })` (multi-
  sampling with a transparent background; illumination) compile no pipeline in the frames after setting them (with the
  props ignored by `prepare`: 1 and 12 pipelines).
- Step 3: a pass's control flow is written once, as a plan against a `PassBuilder` (`canvas3d/passes/pass-builder.ts`)
  that two builders interpret: the render builder draws, the prepare builder requests the pipelines of each draw
  (`renderer.collect` for renderer passes, the full-screen renderables with their state) for the target described at
  that point, without drawing. A plan's code that only matters for drawing (per-frame uniforms, viewports, clears, timers,
  flushes, jitter, sample and iteration counters, the remembered marking base) is behind `b.drawing`; code that selects
  programs or targets (props applied to renderables, texture bindings, attachment changes) runs under both. Requests
  are made as the plan runs, as a render item's program depends on its current defines and textures, which a plan may
  change between draws (the bloom blur's kernel size per mip level). Plans: the draw pass (blended, wboit, dpoit, bloom,
  helpers, output, copy), postprocessing (SSAO, outlines, shadows, background, compose), antialiasing (FXAA, SMAA,
  CAS), DOF, marking, multi-sample and illumination (tracing, compose, its multi-sample). Preparing covers what draws
  differ between repetitions: one dpoit layer and one multi-sample sample stand for all, temporal multi-sampling plans
  its first and a following frame, illumination its first and a following iteration; the tracing pass keeps the quality
  adapted to frame times (uniforms). `ImagePass.prepare` updates the background (loads), links the scene's variants
  and plans the render with the prepare builder; `render` prepares before each render. The re-render loop stays as a
  safety net and reports in debug mode what a render started that `prepare` missed (pending work labels name the
  pipeline: program, variant, state, target, draw mode). Over the WebGPU and WebGPU core visual profiles no render
  starts work (with the planning in `prepare` disabled: 57 of 142 cases report missed pipelines).
- Step 2: the renderer's passes are descriptions (`PassDraws` in `gl/renderer.ts`: variant, the objects a pass is for
  by geometry (`select`), the objects it draws in their current state (`check`), the draws per object with flag and pass
  state, solid interior caps inline or after all objects) used by both `renderPass` (every `render*` method: its uniform
  updates, then the pass) and `renderer.collect(pass, group, target, { all })`. The per-object state is one pure
  function (`getObjectState`: cull and front face from flip-sidedness, double-sidedness, reflection, geometry type and
  the pass flag), the solid interior mark and fill states are cached per pass state (were created per draw). `collect`
  ignores visibility and culling (out-of-view objects are collected); `all` also takes objects the state checks leave
  out (e.g. unmarked objects for the marking passes); caps only for targets with stencil (`TargetDescriptor.hasStencil`).
  Test (`plugin/headless/_test/webgpu-collect.test.ts`, on Dawn): a scene of cartoon, ball-and-stick, points, a
  double-sided transparent surface and an emissive transparent surface, every renderer pass with and without a
  transparent background: every draw the pass issues (render item, variant, render state) was collected, and after
  `preparePipelines` the draws start no work; collecting with flags ignored makes it fail. Solid interior caps (clip
  objects) are not in the test scene.
- Step 1: `TargetDescriptor` (attachment formats, compared by `key`) from `GraphicsContext.describeTarget(target)` as
  currently attached (WebGPU: the color and depth formats; the drawing buffer without taking the canvas texture; WebGL:
  one descriptor for all targets), and `resources.preparePipelines(requests)` with `PipelineRequest` (render item,
  program variant, render state, target). Requests take any render item: full-screen passes draw compute items with the
  `compute` variant, which on WebGPU are render pipelines too; renderables expose their item (`renderItem`). WebGPU
  render items register how they request a pipeline (the program for the variant, as drawing selects it), and drawing
  and preparing share one keyed path (`getPipelineFor`); WebGL links the variant's program. Test
  (`plugin/headless/_test/webgpu-prepare.test.ts`, on Dawn): a prepared request's draw starts no work, an unprepared
  render state's draw starts one compile. Found on the way: the headless WebGPU tests skipped silently since the
  optional `webgpu` package was no longer installed; with it, the existing test failed on two upstream changes
  (`DefaultPluginSpec` moved to `@molstar/plugin/default-spec`; the headless plugin has no `canvas3dContext`, the test
  waits for `Canvas3D.contextChanged` with `restore`), both fixed. `webgpu` is now an optional peer of
  `@molstar/plugin-headless`, installed with `pnpm native:install -- --webgpu`; tests without it show as skipped
  (`describe.skipIf`), as the WebGL context tests without `gl` now do. Clears on WebGPU draw with their own pipelines, created
  synchronously on first use per target and clear kind (`render-pass.ts`, logged in timing mode); not covered by
  requests, as measured not needed: only clears that are not load operations draw (a partial viewport's scissored
  clears, clears after drawing), in the Viewer seven pipelines over a session with a relative-frame viewport and
  illumination, each created in at most 0.2 ms with no frame gap above 17 ms (the shader is trivial).
- Compile time per pipeline (Chrome trace of Dawn's categories, `gpu.dawn`; Viewer loading 1cbs on WebGPU, NVIDIA,
  32 pipelines): DXC is about 90 % of it (Tint's HLSL about 9 %), and the vertex stage cost more than the fragment
  stage (75–157 against 16–93 ms, compiling in parallel). The wiggle and tumble noise (`fbm`, 9 calls inlined into every
  vertex stage, off by a uniform) was about half: alone on Dawn (D3D12, DXC) mesh 95 → 45 ms, spheres 133 → 73,
  cylinders 153 → 82 ms without it. Now compiled only where an object animates (`dWiggleAnimation`, `dTumbleAnimation`
  from the amplitudes; wiggle also with the wiggle layer, `dWiggle`): the load's pipelines 3.5 → 1.9 s of compile time,
  the longest 235 → 130 ms; first draw 170–185 ms after the state update (was 260–280 ms); WebGL2 with a cold cache
  1.2 → 1.0 s to the first draw. Turning an animation on compiles the object's pipelines once (its draws are skipped
  meanwhile, as for any define change). Images unchanged (new visual cases `renderer/animation/wiggle-*`). The bump
  noise of the fragment stage likewise (`dBump`, `BaseGeometry.hasBump`: frequency, amplitude and the bumpiness of the
  material or the interior's substance above 0; with a substance layer, `dSubstance`, too, its bumpiness per group): the
  representations' frequency and amplitude default above 0, their bumpiness to 0. Alone on Dawn, the color pipelines
  (before both, after both) mesh 104 → 37 ms, spheres 131 → 58, cylinders 146 → 69, pick pipelines mesh 69 → 23,
  spheres 97 → 46; the load's pipelines 1.76 s of compile time, the longest 113 ms; first draw 140–170 ms after the
  state update; WebGL2 with a cold cache 0.88 s to the first draw. Images unchanged (new visual cases
  `renderer/bump/*`, the material's and a substance layer's bumpiness). The other runtime-switched paths, each
  disabled alone (median of 11 interleaved compiles on Dawn): the lod fade of the fragment stage 4–6 ms (mesh, spheres,
  cylinders: 9–12 %), now behind `dLodFade` (the object's `lod` with a min or max distance; images unchanged, new visual
  cases `renderer/lod-fade/*`): mesh 36 → 33 ms, spheres 59 → 53, cylinders 70 → 64, the mesh pick pipeline 23 → 18; the
  load's pipelines 1.60 s, the longest 97 ms. Without it the mesh programs read no draw uniforms: programs bind
  `uDraw` only where they read it. The spheres' quadratic projection (the stereo and XR eye cameras' asymmetric
  projections, by a uniform) cost 7.6 ms: the bounding rectangle of `sphereProjection` now takes the off-axis terms
  (`uProjection[2].xy`, zero for symmetric projections) instead, so the quadratic projection and
  `uIsAsymmetricProjection` are removed: sphere color pipelines 57 → 49 ms, pick 42 → 33. Images unchanged; a new
  visual case `renderer/camera/stereo-spheres` (sphere impostors, a wide eye separation) differs from the quadratic
  projection in about 10 pixels by 1/255. Fog, the interior mix, the transparent branch, the vertex stage's lod
  levels and the color and size type switches cost 1–3 ms each, not worth a variant. Lighting alone is about 13 ms of
  the mesh's; the rest is a floor of Tint and DXC. Found on the way: tumble seeded its noise by the object's id
  (`uObjectId`), a creation counter that differs by backend, merging and each rebuild (WebGPU's images differed from
  WebGL2's, and a rebuilt representation jumped to another pose): now by the object's bounding sphere center placed by
  the instance's transform (the same merged, `aSegmentSphere`), visual cases `renderer/animation/tumble-unmerged` and
  `tumble-merged`; and the byte-count test on WebGPU (`webgpu-byte-counts.test.ts`) failed now and then: its peptide is
  drawn as ball-and-stick, whose spheres and cylinders hold their data in data arrays, so the attributes after the first
  frame were mostly the helpers', not all created yet (256 against 336 bytes settled); it now compares against counts
  settled by a first load and clear.
- Planned: replaces the collect mode (`resources.preparePrograms(run)`, `DrawPass.prepareMarking`) and the
  `ImagePass` re-render loop from [resource-readiness.md](resource-readiness.md). A dry run of the `ImagePass` passes in
  collect mode was tried: over the WebGPU visual profile renders that started work went from 63 of 142 cases to none
  (with the background loaded first); not kept, as it runs the draw path in a special mode instead of preparing from
  params.

## Goal

A `prepare` that, given a context and params, makes every pipeline, shader and resource a render with those params
needs ready, without drawing, recording commands or running pass code in a special mode. Later, readiness from params and
schemas alone (no data), so that replaying or animating state snapshots compiles nothing while it plays.

## Where things stand

What a draw needs, per backend:

- **WebGL**: a linked program per render item and variant (defines from the render object's values). Render state is
  applied per draw and costs no compilation. `linkPrograms(variants)` links them ahead of drawing, tracked as pending
  work (`KHR_parallel_shader_compile`).
- **WebGPU**: a render pipeline per program and _configuration_: render state (blend, depth, stencil, cull, front
  face), the pass target's formats and the draw mode (`getPipeline` in `gl/webgpu/context.ts`, key
  `${renderStateId}|${target.key}|${drawMode}`). The configuration is known only when a pass draws:
  - the pass state comes from the renderer method (`OpaqueState`, `BlendedTransparentState`, the solid interior stencil
    states, ...) and from props (`transparentBackground` picks another blended state);
  - the renderer adds cull and front face per object (`drawObject`: `dFlipSided`, `uDoubleSided`, `hasReflection`,
    `dGeometryType`, the pass flag);
  - the target formats come from the pass's render targets (capabilities and props, e.g. half float), which
    `RenderTarget` doesn't expose;
  - full-screen passes (postprocessing, SSAO, outlines, SMAA, bloom, DOF, compose, copy) set defines from props inside
    `render` (e.g. `ssao.update` paths), so their programs also change only when they draw.
- **Compute** (WebGPU): kernel pipelines are created synchronously on first dispatch (`createComputePipeline` in
  `gl/webgpu/compute.ts`), a hitch on first use of a GPU surface or other kernels.
- **Other resources**: pass render targets are allocated when a pass is first enabled or resized; lookup and background
  images load asynchronously (pending work, done).

What exists to work around this: family key reuse (configurations drawn by one program of a family are compiled for new
programs of that family; removed after 6g, see there), the idle warm-up (pick, marking through collect mode) and the `ImagePass` dry run (collect
mode around the pass sequence). Collect mode runs the real pass code with draws turned into pipeline requests; it works
but couples readiness to the draw path, needs every pass to tolerate a draw-free run (marking invalidation, multi-sample
indices, illumination iterations) and gives no way to prepare without a scene and a camera.

## Design

### 1. Pipeline requests (neutral)

```ts
interface PipelineRequest {
  item: GraphicsRenderItem; // or an object signature, see 4
  variant: GraphicsRenderVariant;
  state: RenderState; // including cull and front face
  target: TargetDescriptor; // color formats, depth/stencil format
  drawMode: DrawMode;
}
resources.preparePipelines(requests: Iterable<PipelineRequest>): Promise<void>;
resources.prepareKernels(requests: Iterable<KernelRequest>): Promise<void>;
```

- WebGPU: program (module) as `linkPrograms` creates it, then `createRenderPipelineAsync` for the configuration, keyed
  as `getPipeline` keys it, tracked as pending work; resolves when compiled. Kernels with
  `createComputePipelineAsync`, and dispatch uses the prepared pipeline.
- WebGL: links the item's program for the variant; state, target and draw mode are ignored.
- `TargetDescriptor` is what `getPassTarget` derives its key from today; render targets expose it (`RenderTarget.formats`),
  and passes can describe targets they haven't allocated yet.

### 2. One source for draw state (renderer)

The renderer's state selection becomes pure functions shared by drawing and preparing:

- `getObjectState(values, flag, passState)` → render state with cull and front face (today inline in `drawObject`);
- per renderer method a description of its draws: which variant, flag and pass state per object kind and props (opaque,
  blended back/front, solid interior mark/fill/cap, wboit, dpoit, emissive, depth, pick, marking).

`renderer.collect(scene, method, props)` yields the requests for every render object of the scene that the method would
draw (ignoring visibility and culling: a superset, so that moving the camera compiles nothing). The draw path calls the
same functions, so drawing and preparing can't diverge; a debug check (below) proves it.

### 3. Passes prepare from props

Each pass gets `prepare(ctx, props)` next to `render`:

- applies props to its renderables (defines) without drawing: the prop handling in `render` moves into an `update`
  step that both call;
- describes or allocates its targets for the props and the size;
- collects its requests: its own full-screen renderables with their targets, and `renderer.collect` for the scene parts
  it draws.

`DrawPass.prepare`, `PostprocessingPass.prepare` (and the passes it owns), `MarkingPass`, `IlluminationPass`
(`TracingPass`, all iteration-dependent draws), `MultiSamplePass`, `PickPass`, `HiZPass`, `BloomPass`. Composed:

- `ImagePass.prepare(props)` (what `render` will draw) — replaces the dry run;
- `Canvas3D.prepare(props?, { interaction })` — the frame's passes, with `interaction` also pick and marking (replaces
  `ShaderManager.warmUp` and `DrawPass.prepareMarking`); called after commits in idle time, and by applications.

Collect mode and `preparePrograms(run)` are removed when nothing uses them.

### 4. Readiness without data: precompile lists

A program is determined by a _signature_: the render object kind (its schema and program), its define values and the
render variant. A _precompile list_ holds object signatures (without variant), derived from a state's params when it is
loaded; the device derives the variants, render states and targets from each snapshot's canvas props (step 6). Defines
that depend on data (e.g. color and size granularity, varying groups) are over-approximated. A session also records
the signatures it drew, the reference the derivation is tested against.

## Steps

1. **Requests and targets**: `PipelineRequest`, `TargetDescriptor` (`RenderTarget.formats`), `preparePipelines` on both
   backends, WebGPU `getPipeline` keyed from the same descriptor. Unit tests: a prepared request makes the next draw
   with that configuration find its pipeline (no pending work).
2. **Renderer state functions**: `getObjectState` and the per-method draw descriptions, used by the draw path;
   `renderer.collect`. Check: WebGL byte-identical, WebGPU visual profile unchanged.
3. **Pass prepare**: the `update`/`render` split per pass, `prepare` per pass, `ImagePass.prepare` from props; the dry run
   removed. Check: over the WebGPU visual profile no render starts pending work (as with the dry run: 0 of 142), and a
   debug assertion that a draw never requests a pipeline that the last `prepare` didn't (reports the missing
   configuration).
4. **Canvas3D prepare**: `Canvas3D.prepare` with `interaction`, the idle warm-up through it; collect mode and
   `prepareMarking` removed. Check: the Viewer smoke test's first hover compiles nothing (as today), toggling marking,
   transparency modes and postprocessing options that were prepared compiles nothing.
5. **Kernels**: `prepareKernels`, async compute pipelines, the GPU surface visuals prepare their kernels when created.
6. **Snapshot replay**: precompile lists derived from a state's params when it is loaded, see Step 6 below.
7. **Representation transforms**: precompiling when a representation transform applies or updates, with its data, see
   Step 7 below.

## Step 6: snapshot replay

Replaying or animating snapshots (state files, MVS stories) should compile nothing while it plays. Steps 1–5 prepare
what the current scene and props draw (`Canvas3D.prepare`, frames held after props changes); what a later snapshot will
draw is not known until it is applied. Decided: a **precompile list** derived from the state's params when it is
loaded, compiled before playing.

### Measured (6a)

`kinase-story.mvsj` (11 snapshots) and a story of three snapshots (a cartoon fading to transparent, a transparent
ball-and-stick ligand appearing, the cartoon fading back to opaque), played in the Viewer (RTX 4070 Ti, 1280x800):

- **What compiles**: only when a snapshot brings an object kind or define combination not drawn before: 4 of 11
  kinase snapshots, 8–16 pipelines and up to 13 programs each (mesh and text variants, then spheres and cylinders of
  ball-and-stick, with their pick and marking pipelines from the warm-up after the commit). The fades compiled nothing:
  opaque and transparent variants were prepared when the story loaded.
- **WebGL2**: the first ball-and-stick of a sequence freezes the canvas for 740–860 ms (both stories); the main thread
  is idle meanwhile (longest busy stretch 30 ms, no WebGL call over 10 ms): ANGLE compiles the impostor shaders in the
  GPU process, as for illumination in [program-switch-frames.md](program-switch-frames.md).
- **WebGPU**: no frame over 55 ms; the new objects appear when their pipelines compiled in the background.
- **Not compilation**: frames of 80–110 ms when representations of a snapshot are built (CPU, both backends).

So the list matters for both backends: WebGL's freezes, WebGPU's late appearance.

### The precompile list

- **Entries are objects, not pipelines**: per render object kind (mesh, spheres, cylinders, points, lines, text, image,
  direct volume, texture mesh) its define values without the render variant, and the values the renderer's state
  selection reads (`getObjectState`: flip-sided, double-sided, x-ray, solid interior, ...). Backend- and device-neutral.
- **Pipelines are derived on the replaying device**: for each snapshot's canvas props, after the
  device's fallbacks are applied (e.g. no float targets: blended transparency instead of WBOIT or DPOIT), the
  renderer's pass descriptions (step 2) give the variants and render states per entry, the passes (step 3) the
  targets; WebGL creates the variants' programs, WebGPU their pipelines. A list prepares blended
  transparency where that is used, and the other way round, as do other device differences (formats, capabilities).
- **Derived from params**: the signatures a state's snapshots imply, from their representations' params, without data
  or building (6d).
- **Recording**: `Canvas3D` records the object signatures of the scene's renderables since a reset (when commits
  finish and after updates), deduplicated: what a session drew, the reference the derived lists are tested against.
- **Preparing**: when a state or story is loaded, `Canvas3D` prepares its list for every snapshot's props in the
  background (pending work), so that playing compiles nothing; entries of unknown kinds are skipped.
- **Not in the state file**: a stored list would add nothing a derived one lacks except exactness on data-dependent
  defines, and would need a recorded replay (or recording code when writing files), go stale with other Mol\* versions
  and add a versioned field to `.molj` and `.molx`. Representations that don't declare their signatures (e.g. from
  extensions) compile when first drawn, as before.

### Decisions

1. **Approach**: precompile lists (B), derived from the state's params when it is loaded, not saved in the state file
   (above). Decided.
2. **Name**: "precompile list" (`PrecompileList`, `Canvas3D.getPrecompileList()`, `Canvas3D.precompile(list,
   props[])`). Decided.
3. **Transparency and other fallbacks**: entries without render variant, pipelines derived per device and snapshot
   props (above), so fallbacks work in both directions. Decided.
4. **Lists from params alone**: without building representations (too expensive) or drawing: no replay is needed,
   and MVS stories need no extension (their list is derived when they are loaded). Decided.

### Steps

- **6a. Measure**: done (above).
- **6b. Object signatures**: the signature of a render object (kind, defines without variant, state inputs); programs
  created from a signature without a render item on both backends (`resources.createProgram`), render items getting
  theirs from the same caches; the renderer's pass descriptions taking signatures as well as render objects.
- **6c. Record and prepare**: done (above). The recording in `Canvas3D`, `Canvas3D.precompile(list, props[])`
  deriving the variants, states and targets per props (with the device's fallbacks), tracked as pending work.
- **6d. Lists from params alone**: done (above). The object signatures a state's params imply, without data or
  building: visuals and representations derive theirs generically from their geometry kinds, props, themes and layers,
  what data decides over-approximated; transformers register how their params give signatures; a test that the derived
  list covers the recorded one.
- **6e. Loading**: done (above). When a state or story is loaded, its list derived and `Canvas3D.precompile` called
  with every snapshot's canvas props, in the background.
- **6f. Check**: done (above). Load `kinase-story` and the fade story and play them right away (also with transparency forced to
  blended, and the other way round), and assert that playback compiles nothing (pending work, program and pipeline
  counts) on both backends; WebGL's 740–860 ms freeze gone. Measured in the browser (Chrome, both backends), not only
  on Dawn: the time to precompile the derived lists of the stories and the visual suite's scenes (programs, pipelines,
  wall time, how much compiles in parallel) and the memory they take, against what playing them compiles.
- **6g. Fewer variants**: done, see [6g](#6g-fewer-variants) below.

- **6h. WebGL first draws**: analysed, accepted. A Chrome trace (`gpu.angle`) of the fade story on WebGL with WBOIT: ANGLE
  (D3D11) compiles every pixel shader when the program links, on worker threads in parallel, for one render target
  only (its default output layout takes a shader's first output); a program drawing into several targets compiles
  again at its first such draw, `D3DCompile` on the GPU process' main thread, which stalls all of its work: the
  spheres' and cylinders' WBOIT color programs 447 and 277 ms (the freeze at #1), the cartoon's 185 ms while loading.
  It affects every draw into several targets: WBOIT, DPOIT and picking (four targets). Chrome caches the compiled
  shaders on disk: with a persistent profile only the first run froze (723 ms), the next ones not (27–32 ms), so it
  is a once per browser cost. Options: warm draws (each prepared program drawn once into its multi-target layout while
  loading: the stall moves to loading, where the D3D compile then also populates the cache; a draw that reaches
  ANGLE's state sync needs the program's uniform blocks bound and its samplers' texture types valid, i.e. close to a
  real draw), one-target WBOIT on ANGLE (accumulation and revealage in two passes, twice the transparent geometry
  work), or accepting it (first visit only). Decided: accepted, a once per browser cost; warm draws are the fix if
  first-visit stalls come to matter (e.g. embedded stories).

### 6g. Fewer variants

Measured in Chrome on `kinase-story` (the precompile against a run without it, which compiles only what loading and
playing draw):

- **Signatures**: 71 derived, 10 drawn. All drawn ones are derived; the rest are what the over-approximation adds:
  merged variants (half of them; nothing is merged in the story), both double-sidednesses ('auto' quality), the MVS
  multilayer theme's color types (group drawn; uniform, instance and groupInstance derived) and both solid interiors
  for spheres and cylinders (only `true` drawn; why `false` is derived is open, possibly how a representation's mesh
  and sphere params merge).
- **WebGPU pipelines**: 664 precompiled, 69 used. Per pass: emissive 174 (none used: bloom in emissive mode, nothing
  emissive in the story), pick 99 (21), marking 99 (20), color 95 (10), depth 95 (7), color into the WBOIT targets 83
  (7). Every signature gets every pass (`collectSignatures` draws with `all`), with the render states of every flag
  (opaque and transparent, front and back faces, solid interior caps).
- **WebGL programs**: 140 created while precompiling, 89 distinct; 47 of them a run without precompile creates too, 42
  never: mostly emissive variants, the over-approximated signatures' color variants, and postprocessing and bloom
  variants of the stand-in scenes' other branches.

Done:

1. **What a signature is drawn with**: `ObjectSignature.opaque`, `transparent` and `emissive` (how the object may draw;
   not part of the key, merged by OR where signatures of one key meet: `mergeSignatures`, the recorder, the lists). Derived
   from the props (alpha 1 or below, emission) and the layers (transparency, emissive); recorded from the values (the
   final alpha `uAlpha`, transparency and emissive averages). `renderer.collectSignatures` makes a stand-in per way the
   object may draw (alpha 1 or 0.5, emissive or not, with the values the passes check) and lets the passes' own checks
   select, `all` only where a plan asks for it (marking as if marked), as `collect` does for objects.
2. **What large data decides, resolved as small data does**: 'auto' merging as unmerged (`renderMerged: true` gives
   both, as groups of one stay unmerged; before, `true` dropped the unmerged ones), 'auto' quality as medium.
3. **Interaction later**: already so, `precompile` requests the frames' passes first, then marking and picking, and
   pipelines compile in request order.

Found on the way, then removed: WebGPU's family key reuse (from
[resource-readiness.md](resource-readiness.md) 4b, before preparing existed) compiled, for every new program, the
configurations its program family (name and variant) drew with. Frames link the variants the props require for every
object (`ShaderManager`), so opaque objects got `depth` programs and, by reuse, the transparent objects' depth
pipelines, which they never draw. Measured in the Viewer (WebGPU) with and without it: loading 1grm and the first hover
skip no draws either way (the warm-up and the interaction preparing cover them; 206 and 205 pipelines); the kinase
story with precompile compiled 110 pipelines (12%) more with it, none needed; without precompile, playing compiled 146
pipelines and skipped draws in 37 frames with it, 81 and 22 without (the extra compiles delay the needed ones). Removed:
`familyPipelines` and the compiles in `onProgramCreated` (pipeline configurations stay cached by key), and with it
`ProgramItem.link`, which precompile had used to mirror it. Checked: the headless precompile tests, the Viewer smoke
test (first hover, context loss restored) and the visual suites.

Measured in Chrome after 6g (the same runs as 6f):

| backend, story | signatures | precompile                                 | playing compiles             |
| -------------- | ---------- | ------------------------------------------ | ---------------------------- |
| WebGPU, kinase | 71 → 24    | 664 → about 150 pipelines, 7–8 → 1.6–1.7 s | nothing (WBOIT and blended)  |
| WebGPU, fade   | 16 → 8     | 1.6 → 0.7–0.8 s                            | nothing (WBOIT and blended)  |
| WebGL2, kinase | 71 → 24    | 131 → 78–86 programs, 3.5 → 1.9–2.2 s      | 3 programs (`tracing`, ...)  |
| WebGL2, fade   | 16 → 8     | 60 → 36 programs, 1.5 → 0.8–1.0 s          | 2 programs (`tracing`)       |

WebGL's WBOIT freeze on the first impostor draws is unchanged (6h).

Fewer defines (uniforms with uniform branches instead of defines, which shrinks the program space for drawing too) is
a separate plan: [fewer-defines.md](fewer-defines.md).

## Step 7: representation transforms

A representation's programs compile after its geometry is built: the state update commits, the first frame requests
the pipelines, and draws are skipped (WebGPU) or the frame waits (WebGL) until they are ready. A transform that
applies (`StructureRepresentation3D` and the other registered ones) or updates could start compiling first, so that
compiling overlaps building.

- **Why not the tree** (step 4, tried and dropped): precompiling a state update's tree when the update starts gained
  little, as the representations' transforms arrive in the update's last tree (about 30 ms before the commit for
  1cbs), and the params-only list compiled 50 pipelines where the load needed 32, competing for the compiler's threads
  with the ones needed.
- **What the transform adds**: the time its geometry takes to build, which the overlap can hide (nothing for 1cbs;
  seconds for large structures, surfaces and mesoscale scenes, where compiling took 440 ms in Chrome for 1grm), and its
  data, which resolves what the params-only derivation crosses: the color theme's granularity, 'auto' quality by size,
  the visuals' choices (GPU or CPU, impostors or meshes), merging; a narrower list than the tree's.
- **Updates**: a theme or props change of an existing representation compiles programs its new defines create;
  compiling them first keeps the old geometry drawn until they are ready, instead of skipping its draws meanwhile.

Design:

- **Additive precompile**: `Canvas3D.precompile` releases the previous call's items; transforms precompile alongside
  each other (a preset's representations, one update), so items are kept per owner (the transform) and released once
  its objects are committed and drawn (their own items then hold the programs).
- **Batched**: the requests of one update are merged before planning the passes against the stand-in scenes (CPU work
  per call), so that a preset of several representations plans once, not once each.
- **Not awaited**: the transform starts precompiling and builds right away; awaiting would serialize the two again.
- **With data**: `Representation.getSignatures(props, ctx, data)` beside the params-only derivation: the factories'
  choosers and the themes' granularity resolved from the data (themes created with it, cheap next to building).
- **Not where linking blocks**: WebGL without `KHR_parallel_shader_compile` links on the main thread, which would delay
  building; skipped there.

Steps:

- **7a Measure** (done): in Chrome (RTX 4070 Ti, 1280x800, the Viewer, a fresh browser per case), the time from
  adding a representation to its commit (its build), and from the commit to the first frame that draws it fully (WebGPU:
  skips no draw; WebGL: the first frame, which waits for linking); then the color theme set to `uniform` (another
  color type) and `alpha` 0.5 on it:

  | case (WebGPU / WebGL2)    | add: build      | add: commit to full | theme: commit to full      | alpha: to full |
  | ------------------------- | --------------- | ------------------- | -------------------------- | -------------- |
  | 1cbs cartoon              | 25 / 26 ms      | 57 / 67 ms          | 51 / 68 ms                 | 17 / 0 ms      |
  | 1cbs Gaussian surface     | 93 / 476 ms     | 82 / 80 ms          | 68 / 69 ms                 | 17 / 0 ms      |
  | 4v6x Gaussian surface     | 233 / 1183 ms   | 96 / 83 ms          | 69 / 86 ms                 | 18 / 0 ms      |
  | 3j3q cartoon (`medium`)   | 3209 / 4121 ms  | 1115 / not measured | 1063 / 61 ms               | 9 / 0 ms       |
  | 3j3q Gaussian surface     | 1428 / 1489 ms  | 72 / 81 ms          | 71 / 76 ms                 | 17 / 0 ms      |
  | PetWorld spacefill        | 1409 / 1497 ms  | 183 / 99 ms         | not measured               | not measured   |

  - **Adding**: what follows the commit (8 pipelines, 25 to 34 with a GPU surface's kernels; 2 programs on WebGL) is
    60–190 ms, 1.1 s for 3j3q's cartoon on WebGPU; the build is longer than that except for 1cbs's cartoon, so
    compiling during it would hide most of it: 60–190 ms per added representation, up to 1.1 s.
  - **Updating**: a theme change compiles 1–3 pipelines (WebGPU) or a program (WebGL), 50–90 ms (1.1 s for 3j3q's
    cartoon on WebGPU). Meanwhile WebGPU skips the representation's draws (4–5 frames without it), WebGL holds the frame.
    The update's own build is often shorter than the compile (20 ms against 51 ms for 1cbs), so starting the compile
    earlier would not keep it drawn: frames after an update need holding until its pipelines are ready, as after a
    props change (`holdFramesUntilPrepared`), more than precompiling.
  - **Updates held** (done): `Canvas3D` records the renderables a representation update changed in place (its
    `updated` re-adding it, or `update(repr)`) again, and where a signature changed (`PrecompileRecorder.recordAgain`),
    frames wait for what they need (`holdFramesUntilPrepared`, at most 1 s). Measured again (WebGPU): the theme change
    holds the last frame 36 ms (1cbs cartoon), 53 ms (1cbs Gaussian surface), 1.0 s (3j3q cartoon) and skips no draw,
    `alpha` none either. Test on Dawn (`webgpu-update-hold.test.ts`). Not held: updates that replace render objects (a
    rebuild), whose new objects compile after the commit as added ones do.
  - Not measured: WebGL2's 3j3q cartoon add (no frame caught after the commit) and PetWorld's updates (the theme
    change reached no canvas update; `alpha` rebuilt its geometry, 9 s). One run each.

- **7b With data** (done): `SignatureContext.data` (the representation's data, `repr/signatures.ts`): 'auto' quality
  resolved by it (`getQualitySignatures`), the visuals the data selects (`signatureVisuals` with data: per symmetry
  group for units representations, the structure for complex ones), merging and GPU culling by the objects' counts
  (`withMergedSignatures` with counts: 'auto' resolved, the draws per object by their levels of detail; where merging
  resolves on, unmerged as well, for values that don't merge); the themes are the
  transform's, created with the data (their granularity). Test on Dawn (`webgpu-prepare-representation.test.ts`): the
  lists with data cover what every representation type (and some props and themes) on 1cbs and DNA draws, and 32 chains
  merged by 'auto' merging, which the params-only list misses (it takes 'auto' as unmerged). Small structures: the
  lists equal the params-only ones (16 signatures each); the data decides for large ones.
- **7c Prepared while they build** (done): `Canvas3D.prepareObjects(owner, list)`: requests of one task planned
  together (`planSignatures`, shared with `precompile`) for the frame's passes with the current props, the stand-in
  scenes' averages from the objects and the scene (not every branch), without the interaction (the warm-up's, once
  drawn), with the culling's kernels of objects culled on the GPU; the program items kept per request until each owner
  (the representation) was added and committed and the work is ready, at most 60 s. `prepareRepresentation`
  (`plugin/state/precompile.ts`) in the apply and update of the structure, volume, particle and shape representation
  transforms (shapes colored and sized by shape group), before building, not awaited; not
  where compiling blocks the main thread (`GraphicsCapabilities.parallelCompile`: WebGL without
  `KHR_parallel_shader_compile`). An update's hold waits for work in flight too (`holdFramesUntilPrepared(true)`): its
  pipelines are requested while it builds. Test on Dawn: the first frame after adding a cartoon and ball-and-stick of
  1cbs compiles no pipeline, with `prepareObjects` stubbed out it does; a volume's isosurface (a GPU texture mesh) and
  particles' spacefill are covered by what their transforms prepare. Not prepared: the measurements, unit cells and
  bounding boxes (small, built at once), which compile when first drawn. Measured again in Chrome (as 7a): the commit to
  the full frame of an added representation goes to 0 on both backends where the build is longer than the compile
  (WebGPU: 3j3q cartoon 1.3–1.5 s to 0, its total 5.2–6.0 s to 4.8–5.5 s; 3j3q Gaussian surface 107–120 ms to 0, 1cbs
  Gaussian surface 82 ms to 0; WebGL2: 4v6x, 3j3q, PetWorld 80–100 ms to 0), unchanged for 1cbs's cartoon (built in
  25 ms); builds vary by up to 1.3 s between runs, three runs each with and without. PetWorld's spacefill (8 objects, a million
  instances) waited 230–250 ms after its commit for spheres derived as merged and culled on the GPU, drawn unmerged:
  the draws per object count their levels of detail, and the objects' values may not merge (`canMergeRenderObjects`,
  known once built), so where merging resolves on, both are prepared; now 0 (three runs).
- **7d Particle targets** (done): the particle target representation (`repr/particles/representation/target`:
  particles drawn as the structures, shapes or volumes their targets are, a visual per target, merged by type) has no
  `getSignatures`, so neither precompile lists (mesoscale states, stories) nor `prepareRepresentation` cover it: its
  objects compile when first drawn, after the build.
  - **With data** (the transform's apply and update, `SignatureContext.data` the particle list): per target of
    `targetMapping` with particles (as `createOrUpdate` walks `getParticleTargetGroups`), its props
    (`createParticleTargetPropsProvider`: the kind's group, quality resolved by the target's data) and the geometry
    `createTargetGeometry` would make, without making it: structures as sphere impostors or meshes
    (`structureTargetUseImpostor`) or a blob surface (mesh), volumes as an isosurface (mesh) or dots (spheres, the
    uniform dot size theme), shapes as their shape's geometry; its signatures from that geometry's utils with the
    props (`getGeometrySignatures`), colored by the theme's granularity (instances: a particle per instance, one
    group). Merging as `mergeRenderObjectsByType` decides, per geometry type: the targets of the type, their particles
    and levels of detail (`withMergedSignatures` with counts).
  - **Without data**: no signatures. The targets are not known, and taking every kind the params configure
    (structure impostors, meshes and blob surfaces, volume isosurfaces and dots, shape meshes, lines and spheres, each
    by color types and passes) would compile dozens of pipelines a particle list never draws. Not needed either: the
    representation prepares with its data when it is added and on every update (7c, updates held until ready), and a
    story's later snapshots usually keep the particle list its first build prepared.
  - **Done**: `getTargetGeometryChoice` (the geometry a target is drawn with, which `createTargetGeometry` follows) and
    `getTargetGeometryUtils`; the representation's `getSignatures` with data (sizes as the visual sets them: structure
    impostors per element, volume dots uniform; surfaces with varying groups where flat interpolation is supported, a
    shape's mesh as its geometry has them). Found on the way: `prepareObjects` planned on a timer, which a build that
    never yields to it (a small target) only reached after its commit and first frame; it plans at the caller's first
    await now (a microtask: less batching across representations, but never late). Test on Dawn
    (`webgpu-prepare-particle-targets.test.ts`): RELION particles with a structure (1crn), a volume (a DX grid) and a
    shape (an OBJ cube) target each: the lists with data cover what their first frame records, which compiles no
    pipeline (all three fail without the signatures).
  - **Measured** (Chrome, RTX 4070 Ti, 1280x800, the Viewer, a fresh browser per run, three runs): loading
    `mmcif_particles` files (its format's parse, then its visuals: the particle target representation), the time from
    the visuals' return (built and committed) to the first frame that skips no draw (WebGPU) or the first frame (WebGL,
    which waits for linking), prepared and with `prepareObjects` off:

    | file (parse, build)                     | WebGPU prepared | WebGPU not  | WebGL2 prepared | WebGL2 not |
    | --------------------------------------- | --------------- | ----------- | --------------- | ---------- |
    | PetWorld postsyn (1.4–1.8 s, 1.2–1.6 s) | 340–370 ms      | 650–730 ms  | 160–230 ms      | 290–310 ms |
    | MG_1184 cellpack (2.0–2.3 s, 2.1–2.4 s) | 980–1000 ms     | 1.7–1.8 s   | 360–370 ms      | 480–500 ms |

    Prepared, the 16 pipelines (5 sphere variants, postprocessing, SSAO, Hi-Z, the culling's and the raster's kernels)
    are requested as the build starts and none after the commit; they are known ready only when frames run (their
    promises settle once the main thread yields). What the wait is then: the commit's own work before the first frame
    (PetWorld 255 ms, MG_1184 880 ms, of it the scene commit about 480 ms with 137 ms of buffer writes), and a first
    frame of 70–90 ms. Not prepared, the pipelines are requested after that work and the wait grows by their compile
    (WebGPU 0.3–0.8 s). What is left belongs to the load's work, not to compiling.

## Open questions

- Whether `prepare` covers only the given props or also likely next ones (e.g. both transparency modes, illumination
  on), traded against compile time and memory for pipelines never used.
- Params-only derivation (6d): the over-approximation grows the lists about 2.5 times (merged variants the most);
  the compile time and memory for pipelines never used, measured in the browser on the stories (6f), may call for
  resolving 'auto' merging, and for fewer variants and defines (6g).
- Order and priority when preparing many pipelines (visible first, then interaction, then the precompile list).
