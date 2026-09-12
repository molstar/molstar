# Mol\* 6.0: rendering backends and a path to WebGPU

Design and blast-radius analysis against the current 5.x source tree (`molstar@5.11.0`, inspected September 2026). This complements the [v6 architecture](v6-architecture.md). Names below are proposed contracts, not implemented APIs.

## 1. Recommendation and scope

Introduce a **rendering-backend boundary around device resources, compiled scenes, render passes, and GPU operations**. Retain the existing WebGL implementation behind it. Keep molecular models, representation selection, CPU geometry, camera interaction, and plugin state above that boundary.

The useful starting point is the render object produced by a representation. A backend accepts its geometry/material data and owns turning it into GPU resources and drawing it. This boundary must also cover GPU-produced geometry and readback; wrapping only `Renderer.create()` would leave much of the application dependent on WebGL.

For 6.0, establish the boundary and preserve WebGL behavior. A production WebGPU renderer, WGSL shader migration, and feature parity are later work. A small second-backend experiment can validate the design without becoming a release feature.

Preserve a route to portable scene extraction for [offline renderers such as Blender](#71-offline-rendering-with-blender). They share geometry/material readback with interactive backends but use a separate render-job contract. A Blender extension and a complete scene-snapshot format are later work, not prerequisites for 6.0.

The interface can be small, but extraction is a substantial cross-cutting refactor. The current abstractions expose GPU resources in geometry and themes, and expose WebGL passes through Canvas3D. This work should have its own PR sequence alongside packaging, rather than being treated as a rename of `WebGLContext`.

The minimum useful 6.0 outcome is:

- Backend injection at canvas/context creation, with an asynchronous initialization path.
- A shared render-object contract with CPU data and opaque backend-owned resources.
- Backend-owned scenes, rendering passes, and asynchronous readback.
- Semantic capability checks and a narrow geometry-service boundary for existing GPU operations.
- Existing WebGL-specific integrations isolated behind explicit adapters where a portable implementation is deferred.
- A tested WebGL backend with unchanged default output and no unnecessary per-frame copying.

This makes another backend implementable without changing molecular providers, plugin state, or the general canvas controller again. It does not make every current accelerated feature portable automatically.

## 2. Current dependency surface

### 2.1 Inventory

A source scan counted files with direct relative imports into `mol-gl/webgl`, excluding `_spec/` and `src/tests`. It includes type-only imports because those also constrain published interfaces. These are dependency indicators, not counts of files that all require substantial edits.

| Area | Files directly importing `mol-gl/webgl` | Main reason |
| --- | ---: | --- |
| `mol-canvas3d` | 31 | Context, rendering passes, picking, helpers, XR |
| `mol-repr` | 30 | Context parameters, capability checks, GPU geometry |
| `mol-geo` | 9 | Texture-backed geometry and smoothing |
| `extensions` | 12 | Geometry export, orbitals, tunnels, debugging |
| `mol-math` | 2 | GPU Gaussian density |
| `mol-plugin` | 2 | Headless and viewport screenshots |
| `mol-theme` | 1 | GPU-backed color grids |
| `apps` | 1 | Direct WebGL use |
| **Outside `mol-gl`** | **88** | Direct edges only |

Additional consumers reach WebGL through properties such as `plugin.canvas3d.webgl`, or through `Scene`, render-object, and schema modules. They are not all represented in that count.

The implementation that ultimately needs an equivalent on a new backend is also sizable:

| Current folder | TypeScript files | Approximate source lines |
| --- | ---: | ---: |
| `mol-gl/webgl` | 16 | 5,300 |
| `mol-gl/renderable` | 11 | 1,100 |
| `mol-gl/shader` | 91 | 8,100 |
| `mol-gl/compute` | 8 | 1,100 |
| `mol-canvas3d/passes` | 21 | 7,500 |

Counts exclude tests and include comments. They omit other files such as `renderer.ts`, `scene.ts`, and GPU math. Much of this code can move behind the WebGL adapter unchanged in 6.0; the counts describe the eventual backend implementation surface, not a requirement to rewrite it now.

### 2.2 Why the existing Renderer interface is too low

[`Renderer`](../src/mol-gl/renderer.ts) already has an interface, but its methods accept concrete `Scene.Group` and WebGL `Texture` objects. Its surface includes individual transparency and depth passes such as `renderWboitTransparent` and `renderDpoitTransparent`.

Replacing that implementation would still leave:

- [`Scene.create`](../src/mol-gl/scene.ts) constructing WebGL renderables and sorting them by GL program identity.
- [`createRenderable`](../src/mol-gl/render-object.ts) dispatching directly to the WebGL implementation for each render-object kind.
- [`Renderable`](../src/mol-gl/renderable.ts) exposing `Program`, WebGL render variants, and draw machinery.
- [`DrawPass`](../src/mol-canvas3d/passes/draw.ts) and the other passes allocating framebuffers/textures and issuing backend-specific work.

These objects belong inside a backend. A future WebGPU implementation should choose its own pipelines, resource bindings, and pass organization.

### 2.3 Representation and geometry leaks

[`RepresentationContext`](../src/mol-repr/representation.ts) and [`VisualContext`](../src/mol-repr/visual.ts) expose `webgl?: WebGLContext`. That parameter also appears in visual constructors, `mustRecreate`, value processing, and smoothing operations.

Some uses only ask whether an algorithm is supported. For example, [sphere/cylinder impostor selection](../src/mol-repr/structure/visual/util/common.ts) checks WebGL extensions. Others create and own resources directly:

- [`TextureMesh`](../src/mol-geo/geometry/texture-mesh/texture-mesh.ts) contains vertex/group/normal WebGL textures, a double buffer, and context-restoration metadata.
- [`DirectVolume`](../src/mol-geo/geometry/direct-volume/direct-volume.ts) holds a WebGL grid texture. Its [representation](../src/mol-repr/volume/direct-volume.ts) chooses 2D atlas versus 3D texture storage and uploads it.
- [Color themes](../src/mol-theme/color.ts) can return GPU-backed grids; this is not just a representation issue.
- [Mesh smoothing](../src/mol-geo/geometry/mesh/color-smoothing.ts) and [texture-mesh smoothing](../src/mol-geo/geometry/texture-mesh/color-smoothing.ts) mix material updates with GPU work.

The existing `GraphicsRenderObject` shape is a useful starting point, but it is not yet backend-independent. Its value types derive from [schemas](../src/mol-gl/renderable/schema.ts) coupled to GL uniforms, attributes, and textures.

### 2.4 Canvas, interaction, and host integration

[`Canvas3DContext`](../src/mol-canvas3d/canvas3d.ts) owns WebGL and `Passes`; it explicitly supports multiple Canvas3D instances sharing a context. Canvas3D constructs a scene, renderer, picking/ray helpers, Hi-Z, shader manager, and XR manager, then coordinates concrete passes during drawing.

Its public surface also exposes `webgl`, `getImagePass()`, render objects, and WebGL-shaped statistics. These are migration points for downstream integrations, not implementation details hidden by renaming a class.

Picking has both synchronous and asynchronous paths already. [`PickHelper.asyncIdentify`](../src/mol-canvas3d/helper/pick-helper.ts) returns a polling ticket; [`interaction-events`](../src/mol-canvas3d/helper/interaction-events.ts) supports both modes. Build on that instead of forcing all picking into a synchronous contract.

[Plugin initialization](../src/mol-plugin/context.ts) has async-named entry points, but their internal viewer/context construction is currently synchronous. A backend factory requires that the creation path actually await initialization.

### 2.5 Readback, compute, and extensions

| Consumer | Existing coupling | Needed boundary |
| --- | --- | --- |
| [Headless screenshots](../src/mol-plugin/util/headless-screenshot.ts) | Creates headless GL and Passes; binds a color target and reads pixels | WebGL headless adapter plus a backend-neutral capture result |
| [Viewport screenshots](../src/mol-plugin/util/viewport-screenshot.ts) | Inspects GL limits/extensions and concrete image passes | Capture service, output limits, effective quality settings |
| [Geometry export](../src/extensions/geo-export/mesh-exporter.ts) | Reads texture-mesh vertices/normals/groups and sampled colors through framebuffers | Asynchronous geometry/material readback |
| [Gaussian surfaces](../src/mol-repr/structure/visual/gaussian-surface-mesh.ts), [volume isosurfaces](../src/mol-repr/volume/isosurface.ts) | Direct GPU density/marching-cubes calls and extension tests | Geometry operations with supported CPU fallbacks |
| [Alpha orbitals](../src/extensions/alpha-orbitals/gpu/compute.ts) | GLSL snippets passed into a generic grid-compute helper | An extension-owned backend implementation behind a typed operation |
| [Debug helpers](../src/extensions/debug-helpers) | Construct WebGL scenes/helpers | Explicit WebGL diagnostics adapter or portable render objects |
| [Tunnels](../src/extensions/sb-ncbr/tunnels/algorithm.ts) | Passes WebGL into geometry construction | Shared geometry services or CPU path |
| [XR](../src/mol-canvas3d/helper/xr-manager.ts) | WebGL session/layer handling | Optional backend presentation integration |

Fallback coverage is uneven. Isosurfaces and Gaussian surfaces have CPU paths. Alpha-orbital [single-orbital computation](../src/extensions/alpha-orbitals/orbitals.ts) has a CPU fallback, but [density computation](../src/extensions/alpha-orbitals/density.ts) currently throws without its GPU support. Do not claim that disabling acceleration makes every feature portable.

## 3. Alternatives and tradeoffs

| Approach | Initial change | Limitation | Decision |
| --- | --- | --- | --- |
| Rename/wrap `WebGLContext` | Small-looking, many call-site edits | Preserves GL state, extensions, textures, and synchronous assumptions in the contract | Reject |
| Replace only `Renderer` | Local interface work | Scenes, passes, compute, and readback still require WebGL | Insufficient |
| General buffer/texture/pipeline abstraction for every call | Large rewrite throughout GL, passes, and compute | Commits to a cross-API resource model before a second implementation validates it | Defer |
| Render-object/view boundary with geometry services | Small public surface; substantial extraction | Requires isolating texture-backed data and host integrations | Recommend |
| Replace Mol* rendering with another graphics engine | Broad implementation and behavior change | Custom impostors, picking, molecular materials, compute, and pass quality still need porting | Separate project |

WebGPU uses WGSL, explicit pipelines/command encoding, and different clip-depth and framebuffer conventions. A GL-shaped facade would transfer those mismatches to every caller. Keep those differences inside the backend. See [Chrome’s WebGL-to-WebGPU comparison](https://developer.chrome.com/docs/web-platform/webgpu/from-webgl-to-webgpu).

## 4. Proposed boundary

```text
Plugin state, representations, themes, CPU geometry
                     |
        render objects + geometry services
                     |
Canvas3D controller --+-- Rendering backend
(camera, input,          (compiled scene, passes,
 scheduling, loci)        resources, GPU operations,
                          picking/capture readback)
                              |
                       WebGL now / WebGPU later
```

Retain one `@molstar/graphics` package. Separate modules and enforce import direction; another npm package is not necessary for this boundary.

Suggested subpaths:

```text
@molstar/graphics/rendering/contracts
@molstar/graphics/rendering/render-object
@molstar/graphics/rendering/geometry-services
@molstar/graphics/backends/webgl
@molstar/graphics/backends/webgl/interop
@molstar/graphics/backends/webgpu           # later
```

Initially, the WebGL backend can delegate to existing `gl/` and pass modules. Move files only when their ownership is clear. Backend-neutral modules must not import WebGL implementation modules, including through type declarations or barrel exports. Keep full backend imports at the application/composition entry point.

### 4.1 Context ownership and views

Use a factory for asynchronous backend initialization and a view for each Canvas3D instance. The backend owns the shared device/context and surface resources; a view owns its compiled scene and view-specific rendering state. This preserves existing shared-context usage.

Sketch of the responsibilities, not a complete public API declaration:

```ts
interface RenderBackendFactory {
    readonly id: string;
    create(init: BackendInit): Promise<RenderBackend>;
}

interface RenderBackend {
    readonly capabilities: RenderCapabilities;
    readonly lifecycle: BackendLifecycle;
    readonly geometry: GeometryServices;
    createView(): RenderView;
    resize(size: DrawingSize): void;
    dispose(): void;
}

interface RenderView {
    sync(changes: RenderObjectChanges): void;
    commit(budgetMs?: number): CommitResult;
    readonly bounds: SceneBounds;
    render(frame: FrameRequest): FrameStatus;
    pick(request: PickRequest): Readback<PickResult | undefined>;
    capture(request: CaptureRequest): Promise<PixelData>;
    dispose(): void;
}
```

The named request/result types must contain neutral data, not aliases for current WebGL interfaces:

| Contract | Contents and semantics |
| --- | --- |
| `BackendInit` | Drawing surface or headless target, asset access, initial size, and backend-specific options owned by the factory |
| `DrawingSize` | Pixel dimensions and pixel ratio; no access to `gl.drawingBufferWidth` |
| `RenderObjectChanges` | Added/removed objects, updates, and visibility/marker changes referencing existing versioned data |
| `CommitResult` | CPU preparation progress and remaining work; no promise that the GPU has finished |
| `SceneBounds` | CPU-visible bounds and visibility information needed for camera fitting |
| `FrameRequest` | Camera matrices, viewport, time, rendering settings, and helper layers |
| `FrameStatus` | Whether work was submitted and whether progressive rendering needs another frame |
| `PickResult` | Object/instance/group IDs and position, using existing loci identity conventions |
| `CaptureRequest` / `PixelData` | Size, crop, quality, and normalized pixels; no framebuffer or render-target access |
| `BackendLifecycle` | Ready/lost/recovering/disposed state and invalidation notifications |

This sketch deliberately omits existing property setters and optional facilities. Extract them by responsibility, rather than mechanically copying every method of `Renderer`, `Scene`, or `WebGLContext` into a new interface.

`render()` submits work without waiting for GPU completion. Picking, capture, and geometry export carry their own completion semantics. Shared-context views must preserve viewport clearing, resource lifetime, and ordering; do not create a second device/context for each view as a shortcut.

### 4.2 Preserve render objects and ValueCell updates

Separate `GraphicsRenderObject` data definitions and `createRenderObject()` from the backend-specific `createRenderable()` factory. Retain IDs, geometry kinds, material IDs, visibility/picking state, bounds, and `ValueCell` versions.

For ordinary mesh/point/line/sphere/cylinder/text/image geometry, carry typed arrays and material data across the boundary. The backend decides how to bind/upload them and compile shader variants. Preserve existing field names initially if that avoids a broad producer rewrite; moving shared value definitions out of modules that also construct GL renderables is more important than renaming `aPosition` or `uAlpha`.

Keep shader programs, GL schema binding descriptors, framebuffer attachments, and draw-call objects private. A neutral data schema describes values; each backend maps those values to its binding layout. Existing WebGL GLSL and schema machinery remains usable inside its adapter.

Avoid copying geometry into a new generic scene representation every frame. Synchronize deltas by object identity and `ValueCell` versions. Preserve incremental commits and time budgets, including shader preparation, so large structures do not acquire new long blocking updates.

### 4.3 Texture-backed geometry and resources

A neutral resource handle needs ownership and lifetime, not GL binding methods. It may represent a volume field, a generated surface, or a sampled material grid. Give it a logical kind, backend-owner identity/generation, metadata needed by consumers, and explicit release/readback behavior. Its implementation may use textures on one backend and buffers on another.

For CPU-authored data, keep the CPU source and move upload/storage choices into the backend. For GPU-generated data, return an opaque result that the same backend can render without a GPU-to-CPU round trip.

Apply this to:

- `DirectVolume.gridTexture`: retain grid dimensions, transforms, statistics, and sampling intent; move atlas/3D-texture allocation behind the volume resource operation.
- `TextureMesh` vertex/group/normal textures: expose a generated-surface result. Keep texture layout and double buffering private to WebGL.
- Color/smoothing grids: pass a material-grid handle or CPU data; themes do not expose `webgl/Texture` in their public types.
- Geometry export: request CPU geometry and material samples asynchronously, rather than binding the resource as a framebuffer.

Define reference ownership where geometry is shared between views. Reject resources from another backend/device generation. A lost device invalidates its handles; recovery recreates them from CPU source or reruns their producer. Releasing a view must not destroy resources still owned by another view.

This is an operation/resource boundary, not a public replacement for every WebGL texture method.

### 4.4 Semantic capabilities and geometry services

Replace `webgl.extensions.*` and `isWebGL2` checks in representations with questions about the operation they need:

- Can this backend render sphere/cylinder impostors with the required picking/depth behavior?
- Can it create this volume field or isosurface for these dimensions and data types?
- Can it smooth this geometry's colors/materials?
- Which transparency and postprocessing modes are available?
- What capture dimensions and presentation modes are supported?

Report implemented capabilities and effective limits, not just hardware features. Support can depend on input size or format, so large operations need a suitability check as well as static capability metadata. Keep CPU/GPU selection heuristics near the algorithm; a backend translates physical limits into usable constraints.

Start `GeometryServices` with operations already crossing the boundary: volume preparation, Gaussian-density/surface generation, isosurface extraction, smoothing, and geometry/material readback. Pass numeric/geometry inputs such as positions, radii, grids, and transforms; avoid making the backend depend on plugin state or structure-provider classes.

Reuse existing task progress/cancellation for expensive operations. Preserve CPU alternatives where they exist. Return a structured unsupported result or an actionable error when no implementation exists; do not silently substitute a scientifically different result.

Extension-specific computation stays in its extension. For example, orbital density should have a typed operation with a WebGL implementation, not become a mandatory method on every renderer. The current generic GLSL grid-compute helper remains a WebGL extension facility until a separate backend implementation is written. A runtime shader-language translator is not part of this proposal.

### 4.5 Readback and initialization

WebGPU device acquisition and GPU-buffer mapping use promises; readback cannot be designed around immediate CPU access to a framebuffer. See the [official GPU-compute example](https://developer.chrome.com/docs/capabilities/web-apis/gpu-compute?hl=en).

Use asynchronous creation, capture, and geometry export. For interaction, retain a polling ticket similar to today's `asyncIdentify`: pending, completed hit/miss, failed, or cancelled. A WebGL ticket may resolve immediately. Avoid a promise allocation for every hover sample if the existing reusable ticket path suffices.

Associate picks with camera/viewport/scene versions, discard stale results, and cancel work on resize or disposal. Preserve exact object/instance/group identity and the established coordinate convention for loci lookup. Move input handling to the asynchronous path for portable operation; an explicitly WebGL-only synchronous shortcut can remain in interop while callers migrate. Do not promise synchronous fresh GPU picking on every backend.

Make plugin initialization actually await the backend factory. The existing async entry points and canvas-context injection provide a place to do this. The Viewer supplies the WebGL factory by default; backend-neutral canvas/representation modules do not import it. Explicit convenience entry points may choose WebGL.

### 4.6 Settings, coordinates, and optional presentation

Keep user-facing camera, lighting, material, clipping, and quality settings where their meaning is shared. Move defaults/types out of files that also implement GL passes so importing settings does not load a backend.

Separate requested settings from effective settings when capabilities require a fallback. Preserve existing snapshot fields where possible, and keep backend selection outside scientific scene state. Features tied to a particular algorithm can remain optional settings with support reporting; another backend need not implement all of WBOIT, DPOIT, illumination, or Hi-Z before drawing a scene.

Choose canonical camera and readback conventions. Adapters perform depth-range and viewport/image-orientation conversions; molecular geometry and loci code should not branch on backend type. Test clipping, depth reconstruction, winding/culling, transparency, premultiplied alpha, and screenshot orientation explicitly. Preserve current WebGL output as the compatibility reference.

Keep WebXR/session presentation, headless context creation, shader debugging, and detailed driver/resource statistics as optional adapter capabilities. Common diagnostics can report frame time, draw/object counts, memory estimates, and active features without requiring every backend to emulate GL program/VAO counts.

## 5. Blast radius by workstream

| Workstream | 6.0 abstraction work | Later WebGPU work | Risk |
| --- | --- | --- | --- |
| Device/context and resources | Factory, ownership, lifecycle, size/limits, headless adapter | Device/surface setup and backend resource implementation | Medium-high: shared contexts, loss, disposal |
| Render-object data and schemas | Separate shared values from GL factories; replace leaked handles | Binding layouts, upload/cache logic, per-kind pipelines | High: many representation producers |
| Scene, Renderer, passes | Extract backend-owned view, retaining WebGL implementation | Implement rendering, sorting, transparency, picking, effects | High: output and performance |
| Representation/visual context | Semantic capabilities and geometry services | Reuse callers; implement missing backend operations | Medium-high: constructors and update paths |
| GPU geometry and smoothing | Isolate existing acceleration; preserve CPU paths and opaque results | Native compute or deliberate CPU fallback per operation | High: ownership, numerics, memory |
| Camera/input/loci | Retain domain logic; adapt scheduling, coordinates, and async picking | Backend coordinate conversion and picking implementation | Medium |
| Plugin/UI/screenshots/export | Inject backend, await init, use capture/readback, report support | Mostly reuse neutral paths | Medium; direct WebGL consumers must migrate |
| XR, custom GLSL, debug integrations | Explicit WebGL adapters and capability gating | Separate implementations if required | Optional, potentially high |
| Parsers, molecular models, tasks, state | No new GPU dependency; retain existing contracts | Expected to remain unchanged | Low, beyond the existing v6 layer cleanup |

Most ordinary representations should retain their data selection, loci mapping, and CPU geometry algorithms. Shared visual helpers will still need edits because WebGL appears in their signatures. Specialized volume/GPU visuals and export paths require semantic changes, not only type substitutions.

The full backend migration is much larger than the interface extraction. The 91 shader modules and 21 pass modules explain why “WebGPU-ready” must mean a validated boundary, not implicit feature parity or a small future shader swap.

## 6. Minimal implementation sequence

Land this as an independent workstream before finalizing v6 public graphics/context contracts. Coordinate source relocations with the package-boundary audit; do not combine a new rendering algorithm, package moves, and interface extraction in one PR.

### A. Capture the baseline and separate data

- Record representative WebGL rendering, picking, memory, and performance results.
- Split render-object/value types from GL renderable factories and shader imports.
- Separate settings/capability types from concrete pass modules.
- Inventory `.webgl`, concrete ImagePass, and texture-handle consumers, including downstream examples.

Exit: the proposed shared data/contract entry points have no transitive WebGL implementation imports. Existing WebGL still runs through the old path.

### B. Introduce backend creation and views

- Wrap the existing context, scene, renderer, and passes behind the factory/view contract.
- Move the GPU-specific frame orchestration, helper rendering, and progressive-pass state into the WebGL view.
- Preserve Canvas3D's camera, input, representation membership, and scheduling responsibilities.
- Await initialization and preserve multiple views sharing one backend/context.

Exit: the default Viewer renders through the injected WebGL backend. Swapping a factory does not require editing the Canvas3D controller. Concrete pass/resource access no longer appears in its shared contract.

### C. Close the representation and readback leaks

- Move capability checks to semantic queries; replace raw context parameters in shared visual interfaces.
- Route GPU geometry, uploads, smoothing, and readback through services and opaque handles.
- Migrate screenshots, export, and interaction to the neutral asynchronous operations.
- Isolate custom GLSL/XR/debug integrations behind explicit WebGL adapters. Track every remaining exception with its owner and reason.

Exit: base plugin/UI, shared representation paths, and the portable rendering contract do not expose WebGL types. Existing WebGL accelerated features still work; unsupported features on a different backend are explicit.

### D. Validate the boundary and freeze the v6 API

- Supply a small recording/fake backend with delayed initialization and delayed pick/capture results; ensure controllers do not assume synchronous completion or access WebGL.
- Verify resource ownership and disposal, shared views, context loss/recovery, and cancellation.
- Build a neutral-consumer fixture with no WebGL implementation in its import graph.
- If practical, prototype one WebGPU mesh with depth and picking to challenge the contract. Keep this experiment out of the production feature promise.

Exit: the WebGL acceptance suite passes and the second/fake implementation uses the same contract without adding raw WebGL escape hatches to it. A fake verifies separation and scheduling only; it does not prove WebGPU rendering or performance. Without a real backend experiment, treat the API as an informed design awaiting that validation.

A factory wrapper alone is a useful intermediate PR, but does not satisfy the abstraction goal while shared representations still require WebGL textures/context.

## 7. Later WebGPU implementation

Port a deliberately small feature set first:

1. CPU-backed meshes, depth, transforms, and object/instance/group picking.
2. A useful molecular scene: ball-and-stick via existing mesh fallbacks, highlighting/selection, clipping, and capture.
3. Lines, points, text/images, and sphere/cylinder impostors, with explicit capability reporting during rollout.
4. Volume rendering, GPU isosurfaces/density, material smoothing, and geometry export.
5. Advanced transparency, postprocessing, illumination, XR, and extension-specific compute according to demand.

Retain WebGL as the production default until coverage and performance justify another policy. Backend choice/fallback happens during initialization; changing a live scene to another backend requires recreating backend-owned resources. Automatic fallback must report a changed backend or unsupported capability rather than silently losing requested features.

WGSL, pipeline caches, resource layouts, staging/readback, and pass implementations belong to this later work. Decide which algorithms to share after implementing representative workloads. Do not introduce a shader intermediate language or a universal render graph as a prerequisite for 6.0.

Other interactive targets can implement the same view contract with a different capability set. Offline renderers instead consume an extracted scene and return artifacts through an asynchronous job. Neither contract requires every target to emulate every Mol* effect.

### 7.1 Offline rendering with Blender

A future optional extension, for example `@molstar/blender-render-extension`, could render the current Mol* scene in Blender for final images and animations. Keep WebGL/WebGPU for interactive preview. Blender should not have to implement `RenderView`, hover picking, canvas presentation, or frame scheduling.

The shared foundation is **on-demand portable scene extraction**, built on render-object data and asynchronous geometry/material readback:

```text
Representations + camera/lighting settings
                  |
       Scene extraction + readback
                  |
       Portable scene snapshot
                  |
       Blender extension / render job
                  |
       Image, animation, or .blend file
```

A render scene snapshot contains resolved geometry, instance transforms, colors/materials, visibility, camera, lights, and relevant rendering settings. It is distinct from existing plugin-state snapshot JSON, which records application state and requires Mol* features to reconstruct representations. Snapshot assets must be portable data or resolvable references, with no WebGL/WebGPU handles, live `ValueCell`s, or backend objects.

Extract a consistent scene/camera version at a selected animation time. Retain stable object/instance identities and transforms where useful, and define units, axes, and camera conventions. CPU-backed geometry can be extracted directly; GPU-generated geometry and smoothed materials use readback from their owning backend. Pin the required resource versions until extraction completes, or detect changes and retry/fail explicitly. Do this only for export/render requests, without CPU copies or serialization in the normal frame loop.

The extension can start with this pipeline:

1. Extract the visible scene, resolving supported geometry and material data through the shared services.
2. Convert it to GLB plus Blender setup instructions for the camera, lights, materials, and output settings.
3. Submit a cancellable render job to a local Blender process or an explicitly configured rendering service.
4. Return output artifacts and a report of unsupported or approximated features. Animation repeats extraction at selected times, with asset reuse where practical.

The existing [GLB exporter](../src/extensions/geo-export/glb-exporter.ts) is a starting point for geometry/material conversion, not a complete Blender scene exporter. Camera/light transfer and render orchestration need new work. Its geometry recentering must also be applied to camera/light transforms, or disabled consistently. Blender supports background rendering through its [command-line interface](https://docs.blender.org/manual/en/dev/advanced/command_line/render.html).

Keep the offline contract small: scene input, output/render settings, support reporting, task progress/cancellation, and an asynchronous artifact result. Reuse existing task conventions. The extension owns Blender conversion, process/service transport, temporary assets, and cleanup; `@molstar/graphics` owns only the portable extraction/readback primitives it needs to share. A Node host can invoke an installed Blender; a browser needs a local companion or remote service. Blender and its runtime dependencies must remain optional.

Fidelity needs explicit translation rather than a promise of pixel-equivalent output:

| Mol* feature | Blender adaptation |
| --- | --- |
| Meshes and instances | Preserve geometry and transforms; translate materials |
| Sphere/cylinder impostors | Convert to meshes or Blender-native equivalents |
| GPU surfaces and smoothed colors | Read back geometry/material data; bake where necessary |
| Custom shaders, clipping, transparency, postprocessing | Translate supported semantics, bake, or report an approximation/unsupported feature |
| Direct volumes | Dedicated volume data/material path; not covered by mesh export |

Start with static meshes and ball-and-stick scenes, camera matching, and basic materials. Validate transforms, colors, instancing, output orientation, cancellation, and unsupported-feature reporting before adding volumes or animation. The v6 boundary checks should establish that export can obtain portable CPU data without raw WebGL access; they should not require a Blender installation or freeze a general offline-rendering framework. Implement and validate the snapshot/job format with the extension when that work begins.

## 8. Acceptance checks and migration

### 8.1 WebGL regression coverage

| Scenario | What it protects |
| --- | --- |
| SDF ball-and-stick and a large cartoon structure | Base geometry, camera, colors, materials, instancing |
| Marker/overpaint/transparency updates | ValueCell versions and partial uploads |
| Clipping, depth picking, and overlapping transparent objects | Coordinate/depth conventions and identity |
| CPU and GPU isosurfaces, direct volume, smoothed colors | Geometry services, resource handles, accelerated paths |
| WBOIT/DPOIT, postprocessing, progressive illumination | Existing feature output and scheduling |
| Hover/click with camera motion, resize, and delayed readback | Stale-pick rejection and interaction latency |
| Screenshots, MP4 frames, and geometry export | Pixel/geometry completion, orientation, materials |
| Shared contexts/views and headless rendering | Ownership, output isolation, sizing |
| Context loss, restore, disposal during pending work | Resource generations, cancellation, leaks |

Use appropriate tolerances for rendering comparisons and exact checks for IDs, geometry indices, and lifetime rules. Measure load/commit time, frame time, interaction latency, allocations, and CPU/GPU memory on representative datasets. No full-scene serialization or forced GPU readback belongs in the normal frame path.

### 8.2 Import and API migration

The v6 migration map should include replacements for:

- `canvas3d.webgl` and injected `Canvas3DContext.webgl/passes`.
- `RepresentationContext.webgl`, `VisualContext.webgl`, and visual callback parameters.
- `getImagePass()`, direct target binding/readPixels, and geometry-export texture reads.
- Concrete GL texture/program/renderable types in public extension APIs.
- Extension checks and backend-specific settings/statistics.

Mechanical imports can be migrated automatically. Resource ownership, synchronous readback, and custom shader integrations need a manual-work report. Existing WebGL-only integrations may use `backends/webgl/interop`, with runtime backend checks; shared code must not depend on that entry point.

Backend selection is separate from `PluginFeature` format/representation registration. Registering a representation does not guarantee backend support: its applicability check must account for capabilities and available fallbacks. Keep transformer identifiers, picking IDs, and scene semantics stable; include any unavoidable public API breaks in the migration guide.

### 8.3 Decisions to resolve during extraction

The following do not require a full WebGPU implementation, but should be settled before freezing the v6 contract:

- Exact render-object/value types that can remain shared, and the opaque handles replacing GPU textures.
- Ownership of shared generated geometry and its rebuild source after loss.
- Which semantic geometry operations need a shared service and which stay extension-owned.
- The canonical pixel/depth conventions and capture format.
- How existing synchronous `identify` callers migrate to asynchronous tickets or explicit WebGL interop.
- The supported initial capability set and behavior when a requested setting is unavailable.
- Which optional WebGL integrations remain in the interop exception list.

The recommended scope is the boundary, the working WebGL adapter, and its validation. Full WebGPU rendering and parity should have their own milestones after this foundation is proven.
