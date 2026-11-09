# Plan: Modular passes, effects and material layers

## Status

Proposed; analysed, decisions taken but one (shadow map application), nothing built. Goals: new scene features (transmission with IOR and thickness, shadow maps) and
screen effects (chromatic aberration, motion blur, screen-space reflections, custom ones) added as modules, each
optional to bundle and optional to compile and allocate; and the per-group and per-object material properties
generalized so that a new one (transmission) is not a 30–40 file change.

## Today

### Passes

- `DrawPass` owns and constructs every pass up front: WBOIT, DPOIT, marking, `PostprocessingPass` (SSAO, shadow,
  outline, bloom, background), `AntialiasingPass` (FXAA, SMAA, CAS), DOF, volume (`passes/draw.ts:47-160`). `Passes`
  adds pick, multi-sample and illumination; every `ImagePass` builds a second set (`passes/image.ts:105-115`).
- The order is hand-written: `planEye`, `_renderBlended`/`_renderWboit`/`_renderDpoit` and `_renderOutput`
  (`draw.ts:521-824`) choose targets by flags (`postprocessingEnabled`, `antialiasingEnabled`, `dofEnabled`, …), and
  `getColorTarget` repeats the chain's order. `IlluminationPass` has a second copy of the output chain (outline, SSAO,
  bloom, background, antialiasing, DOF, marking) that reaches into `drawPass.postprocessing.*`
  (`illumination.ts:185-290, 460-680`).
- One compose program (`postprocessing`) folds SSAO, screen-space shadows, outlines, bloom, fog and the background
  into the scene color, switched by defines (`postprocessing.ts:40-82`).
- Props are one static `PostprocessingParams` (`postprocessing.ts:142-215`); `limitToSupported` turns off the effects
  the backend lacks, effect by effect (`supported-props.ts`).
- Programs: GLSL by name in a module-level map filled at import of `gl/webgl/programs.ts` (`registerGlslProgram`
  exists); WGSL in a static object in `gl/webgpu/shader/programs.ts` (no registration). Every effect's sources ship
  with its backend; `check:backend-bundles` keeps the backends apart, nothing keeps effects apart.
- Compile: WebGL render items create their programs on construction (`webgl/render-item.ts:541-544`,
  `resources.program`), so building `DrawPass` starts compiling every effect's quad program, used or not. WebGPU
  creates pipelines on request (`pipeline-readiness.md`).
- Memory: pass targets already take storage only when written and are released when the plan stops using them
  ([pass-targets.md](pass-targets.md)). What is left is the JS objects, render items, uniform buffers and, on
  WebGL, the programs.
- The renderer's passes are a closed union (`Renderer.Pass`, `renderer.ts:576-593`) over a `PassDraws` table, and
  render variants a closed set (`color`, `pick`, `depth`, `marking`, `emissive`, `tracing`; `render-item.ts:23`): every
  geometry program is written for each variant on both backends.
- Lights are directions in view space (`renderer.ts:1171-1180`): they turn with the camera, so a shadow map would
  have to be drawn again on every rotation.
- The existing `shadow` effect is a screen-space ray march over the opaque depth ("Simplistic shadows").

### Material layers

Per-group layers (marker, overpaint, transparency, emissive, substance, clipping, wiggle) are each written out by
hand; nothing describes them as data. Adding one, with emissive as the template, touches:

- data: a `geo/geometry/<x>-data.ts`, a schema in `renderable/schema.ts` spread into `BaseSchema`, the WebGPU
  `DataLayers` list (`webgpu/data-layers.ts:18`, appended to keep the page bits), `createEmpty<X>` in the nine
  geometries' `createValues`, merged objects (`renderable/merged.ts` `BaseTexSpecs` and type checks,
  `repr/merged.ts` `MergeTypeKeys`);
- shaders: GLSL `color-vert/frag-params`, `assign-color-varying`, `assign-material-color`/`apply-light-color`, the
  image and direct-volume programs; WGSL `varyings.wgsl.ts` (an inter-stage location: 0–8, 14, 15 are taken of 16),
  `material.wgsl.ts`, `output.wgsl.ts`; `ignoreDefine` per variant (`shader-code.ts:66-85`) and the WGSL type list;
- repr: `repr/signatures.ts` layers, a near-copy `set<X>` in `repr/visual.ts` (setOverpaint … setWiggle are six
  copies), smoothing per layer, `theme/<x>.ts` (Loci layer: remap, merge, filter, script, bundle),
  `Representation.State`, and the setters and remaps of units, complex, volume, particles and shape representations;
- plugin: an effects transform with `registerTransformLayer`, a helper, the catalog, the default registry, the
  component manager, theme strength; Prism's style layers and the geo exporters.

Generic already: WebGL data textures from any `DataSpec` (generated `read_<name>` accessors), the WebGPU accessors
from the `DataLayers` list, and the plugin's state sync (`SyncStructureRepresentation3DState`, any field of the
state).

Per-object material: `Material { metalness, roughness, bumpiness }` (`core/util/material.ts`) as material uniforms;
the shading switches are per-geometry params and defines (`dXrayShaded`, `dCelShaded`, …). The substance layer's rgba8
is full (three channels and the weight). No transmission or IOR exists anywhere.

## Design

### 1. Screen effects as modules

An effect is a descriptor, not a field of `DrawPass`:

```ts
interface ScreenEffect<P extends PD.Params> {
  readonly id: string;                      // also its props key
  readonly params: P;
  readonly stage: 'compose' | 'image' | 'antialias' | 'final';
  readonly order: number;                   // within the stage
  readonly reads: readonly EffectInput[];   // 'depthOpaque' | 'depthTransparent' | 'normals' | 'previousView' | …
  readonly programs: readonly EffectProgram[]; // ShaderCode with the sources per backend (§2)
  isEnabled(props: PD.Values<P>, scene: Scene): boolean;
  create(ctx: EffectContext): EffectInstance<P>; // targets through `PassTargets`
}
interface EffectInstance<P> {
  plan(b: PassBuilder, io: EffectIO, props: PD.Values<P>): void; // reads io.input, writes io.output
  setSize(w: number, h: number): void;
  dispose(): void;
}
```

- **One chain for every output path.** The draw pass and illumination produce the scene color and depths; an
  `EffectChain` plans the enabled effects in stage order over them, ping-ponging between chain targets from
  `PassTargets` (the last writes the drawing buffer or the output target), then marking presents. `_renderOutput`,
  `getColorTarget` and illumination's copy of the chain become the chain's. Why: one place decides the targets, and a
  new effect works under illumination, multi-sampling and image passes without edits there.
- **Stages.** `compose`: producers whose result the compose program folds into the scene (SSAO, screen-space or
  mapped shadows, outlines, bloom); the compose program and its slots stay one program (they interact: occlusion
  before fog, outlines over transparency, bloom over the background) and stay in core: its slot code is small, the
  producers are the cost. `image`: full-screen passes on the composed image (SSR, motion blur, DOF, chromatic
  aberration, custom). `antialias`: FXAA or SMAA (one). `final`: CAS. Today's order is kept as the default orders.
- **One quad per effect**, no fusion of cheap effects into one program for now: a quad over an rgba8 target is a read
  and a write per pixel; fusing brings back define combinations. Revisit if measured.
- **Inputs on demand.** An effect declares what it reads; the chain makes it available: `normals` runs the tracing
  input pass (normals and emissive, already a G-buffer) unless illumination did; `previousView` keeps the last frame's
  view-projection (without jitter).
- **A registry.** `Canvas3D.effects` (`add`, `remove`, `has`, an event on change); its `ImagePass`es share it. What is
  registered is up to the app: the plugin registers its spec's effects and behaviors add or remove effects at runtime,
  as they do for other providers ([../../designs/plugin-composition.md](../../designs/plugin-composition.md));
  `DefaultEffects` (`@molstar/graphics/canvas3d/effects/default`) is today's set and the default specs list it. An
  effect not registered is neither bundled (unless imported elsewhere) nor run. Built-in and outside effects register
  the same way.
- **Props built from the registry.** `postprocessing` is a group generated from the registered effects
  (`getCanvas3DParams(registry)`), keyed by effect id, rebuilt on change; the registry's event tells UIs and the
  plugin's canvas params to rebuild. The built-in effects keep their ids and param shapes, so snapshots, MVS and the
  UIs keep their paths when the defaults are registered. Props of an effect not registered are kept as given and
  applied when it registers (snapshots load before behaviors run); `limitToSupported` turns off registered effects
  whose programs the backend lacks, from the descriptors.
- **Public API.** `ScreenEffect`, `EffectInstance`, `EffectChain`'s inputs (`EffectIO`, `EffectInput`), the registry,
  `PassBuilder`, `PassTargets` and `registerProgram` are public and stable: effects are written outside the workspace.
  They get docs (`apps/docs`, v6) and an example effect in `examples/`.

### 2. Optional to bundle

- Programs per context, not per module: `GraphicsContext.registerProgram(shaderCode, sources)`; the backend's
  static tables keep the geometry and core programs (copy, compose, OIT, marking, Hi-Z, compute), the effects'
  programs move to their modules. Registering an effect registers its programs with the context (and with a new
  context after a backend switch or a lost context). Unregistered names fail as today.
- Sources per backend, resolved at build time: an effect's GLSL in `<x>.glsl.ts`, WGSL in `<x>.wgsl.ts`, referenced
  from the descriptor (`programs: [{ code, glsl, wgsl }]`); the backend-set builds resolve the other backend's files
  to an empty stub (an esbuild resolve plugin next to `getBackendSetAlias`, by file suffix). A descriptor whose
  backend's sources are stubbed is unsupported there. Outside builds that don't use the plugin bundle both, as they
  bundle both backends today. Not `import()` on first enable: the app builds are single files, and enabling would
  turn async.
- `check:backend-bundles` gains cases: a Canvas3D with the default set minus SMAA contains no SMAA module (its
  inlined area and search textures are about 55 KB of `smaa.ts`); a set with only FXAA contains none of SSAO, bloom,
  DOF, outline, shadow, background, CAS.

### 3. Optional to compile and allocate

- Instances are created when the effect is first enabled (`create`), not with the draw pass: no render items, so on
  WebGL no programs compiled for effects never turned on. They are disposed when released, on the schedule
  `PassTargets` already releases on (a second after the last props change).
- `prepare` plans the chain like any pass (`PassBuilder`), creating instances for hypothetical props as needed and
  releasing them if the props are not applied; precompile lists cover effect programs as they cover the compose
  program now.
- The same holds for the non-effect passes that are optional: DPOIT and WBOIT, illumination and tracing,
  multi-sample, volume: created on first use.

### 4. Scene features: renderer passes and variants

Transmission and shadow maps change how geometry draws, so they are not screen effects. The render variants stay a
closed set (each variant is a program per geometry per backend); new renderer passes are added to the union.

**Shadow maps** (`shadow: 'screen-space' | 'map' | 'off'`; today's `'on'` read as `'screen-space'`).

- A depth map per shadowing light (the first light by default), an orthographic light camera fitted to the scene's
  bounding sphere (molecular scenes are bounded: no cascades), drawn with the `depth` variant (impostors ray-cast in
  it; orthographic through `uIsOrtho`) in a new renderer pass `shadowDepth`, frustum-culled on the CPU by the light
  camera (the GPU culling's Hi-Z is the main camera's).
- Applied deferred, in the compose program: the world position from depth into light space, PCF over the map,
  darkening as the screen-space shadows do. No geometry program changes; a slot of the compose stage. Forward
  application per light in the lit programs (correct per-light diffuse and specular, transparent surfaces shadowed)
  is a later option: a global texture with `usedIf` and a global define that doubles the lit programs.
- Cost: a depth pass of the scene per light whenever the scene or the light's world direction changes, which with
  view-space lights is every rotation; measured before ([sphere-depth-prepass.md](sphere-depth-prepass.md) found a
  depth-only pass about as costly as the impostors without lighting). A light option fixed to the world would make
  the map static while the camera moves.

**Transmission** (`transmission`, `ior`, `thickness`).

- Per object and per group: material uniforms `uTransmission`, `uIor`, `uThickness` (`Material` extended,
  `BaseGeometry` params), and a transmission layer kind (§5) mixing over them as substance does over the material.
- Drawn in a new renderer pass `transmission` after the opaque pass: a draw fact `Transmissive` keeps them out of the
  opaque pass; before it, the opaque color is copied into a backdrop target (with mips for rough transmission). The
  color variant under `dTransmission` refracts the view ray by the IOR at the surface normal, offsets the backdrop
  lookup by the refracted ray over the thickness, attenuates by the object color (Beer–Lambert) and mixes by Fresnel
  with the lit surface. Thickness: the uniform, or measured: analytic chords for sphere and cylinder impostors, the
  back-face depth for meshes (tracing's `depthOpaqueBack` pass already draws it).
- Transmissive surfaces write depth and count as opaque for the passes after (SSAO, outlines, transparency). What is
  behind them is opaque only: transparent objects behind glass are not seen through it (as in other real-time
  renderers). Illumination refracting in the trace is separate and later.
- Colored transmission in blended transparency without a backdrop is
  [dual-source-blending.md](dual-source-blending.md)'s candidate 1, a cheaper and lesser alternative.

**New screen effects** (modules of §1, in the `image` stage):

- Chromatic aberration: radial per-channel offsets, one quad; the first effect written outside the built-ins, the
  check that the interface suffices.
- Motion blur: camera motion only, reprojecting depth with the previous view-projection (no velocity target). Skipped
  while the camera is still, so accumulation (temporal multi-sampling, illumination) is unaffected. Object motion
  (animations, trajectories) needs a velocity output from the geometry programs: not planned.
- Screen-space reflections: reads `normals` and depth, marches in screen space, weighted by roughness and metalness.
  The tracing G-buffer has no roughness or metalness: it needs a channel there (both backends' tracing outputs). Off
  under illumination, which traces reflections already.

### 5. Material layers as kinds

- A `LayerKind` descriptor in graphics: name, data format (`r8`, `rgba8`), channels and their encoding, which
  layer types it supports (texture, vertex-instance, volume-instance), whether it is smoothed, which variants read
  it, its draw facts (as transparency and emissive select passes). Generated from the registered kinds: the
  schema cells (`t<X>`, `d<X>`, `d<X>Type`, the grid cells, `u<X>Strength`) spread into `BaseSchema`, `createEmpty`
  in every geometry, the WebGPU `DataLayers` (built-in kinds in today's order, then new ones: the page bits stay),
  merged specs and `MergeTypeKeys`, signature layers, and one `Visual.setLayer(kind, layer)` replacing the six
  `set<X>` copies.
- Kinds are a public registry like effects: the built-in kinds always registered (they are small and read by the
  core programs), new kinds registered before the representations that use them are created; a kind brings its
  shader code for the programs it supports.
- The per-kind API goes (breaking changes allowed): `theme/<x>.ts` become one `Layer<T>` with a value codec per kind;
  `Representation.State` holds `layers` keyed by kind instead of the named fields (`overpaint`, `transparency`, …,
  `themeStrength` by kind too); `Visual` has `setLayer` only; the plugin has one layer transform and one helper taking
  a kind, replacing the six effects transforms, their helpers and the component manager's per-kind code. Snapshots
  with the old transforms are mapped to the new one when loaded; MVS, Prism's style layers and the geo exporters move
  to the kinds.
- Shaders stay written by hand per kind (what a layer does to shading is not generic), reading through generated
  accessors. Inter-stage locations are the constraint on WebGPU (two left): a texture-type layer can be read in the
  fragment shader by the flat group index instead of a varying (measure:
  [vertex-inputs-and-varyings.md](vertex-inputs-and-varyings.md) packs varyings too).
- Transmission per group: a kind of its own (`rgba8`: transmission, thickness, IOR, weight), not a wider substance.
- Layer presence stays a define ([fewer-defines.md](fewer-defines.md) leaves it open): it keeps a kind's code and
  bindings out of programs that don't use it, which is what "optional to compile" means for layers.

## Steps

0. Measure: bytes per effect in the Viewer bundle (esbuild metafile), programs compiled at Canvas3D creation on WebGL
   with the default props, time to the first frame with and without the eager passes, a depth pass's cost per light
   on the benchmark scenes (`pnpm benchmark`).
1. Programs per context (`registerProgram` on both backends, the WGSL table registrable), the resolve plugin and the
   new `check:backend-bundles` cases.
2. `ScreenEffect`, `EffectChain`; CAS and FXAA moved first, then SMAA, DOF, bloom, and SSAO, outline and shadow as
   compose producers; illumination's output on the chain. Guarded by the visual suite's postprocessing subset and
   `webgpu-pass-target-bytes`.
3. Instances on first use and their release; the optional non-effect passes (§3); `prepare` and precompile over the
   chain.
4. The effect registry in `Canvas3D` (shared by its `ImagePass`es), props built from it, kept props of effects not
   registered; the plugin spec's effects and behaviors registering effects; `limitToSupported` from the descriptors;
   the docs and an example effect.
5. Chromatic aberration (as the example, outside `packages/graphics`), camera motion blur, SSR (with the G-buffer's
   roughness and metalness).
6. Shadow maps, deferred.
7. Transmission per object.
8. Layer kinds (§5) with their registry and the per-kind API replaced; then transmission per group.

Breaking changes (`DrawPass`'s public sub-passes, used by the headless plugin and tests; `PostprocessingParams`
built from the registry; the per-kind layer API) go into [../breaking-v6-changes.md](../breaking-v6-changes.md) as
they land.

## Decided

- Effect sources resolved per backend at build time, not imported on first enable.
- Effects in a registry, registered by the app (the plugin's spec and behaviors); props built from it.
- Transmission per object and per group.
- Layer kinds may change the public API as needed.
- Effects and layer kinds can be written outside the workspace: their interfaces are public, stable API.

## Decisions open

- Shadow maps applied deferred only, or forward per light as well; a light option fixed to the world.
