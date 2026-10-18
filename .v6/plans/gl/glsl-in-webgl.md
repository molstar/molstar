# Plan: GLSL sources only in the WebGL backend

## Status

- **Current**: done, steps 1–7. Step 5 (the WebGL-only compute in `gl/webgl/compute/`) was done by
  [backend-narrowing.md](backend-narrowing.md) step 3, step 7 (the bundle check) by
  [backend-selection.md](backend-selection.md) step 4 (`pnpm check:backend-bundles`).
- Steps 1–4 and 6 done: unique program names (`accumulate` split into `trace-accumulate` and
  `color-smoothing-accumulate`, test "program names are unique"); the GLSL assembly in `gl/webgl/shader-code.ts`;
  `ShaderCode` a descriptor without sources, GLSL looked up by name in `gl/webgl/programs.ts` (`registerGlslProgram` for
  grid3d and custom programs); sources moved to `gl/webgl/shader/`; the guard rejects nested `gl/webgl/...` imports
  (incl. `shader/`) outside the backend and imports between `gl/webgl/` and `gl/webgpu/`. Only grid3d (WebGL-only by
  construction, allowlisted) imports GLSL outside the backend. WebGL profiles byte-identical, WebGPU unchanged. The GLSL
  tests stayed in `gl/_test/` (they share the shader inventories with the WGSL tests). Step 5 (moving the WebGL-only
  compute) is deferred to the compute port ([webgpu-compute.md](webgpu-compute.md)): it gives these algorithms neutral
  entry points with per-backend implementations, so the WebGL ones move to `gl/webgl/compute/` then, changing the
  callers once. Step 7 waits for backend-selection step 4.
- Originally planned. Details the "program sources per backend" part of [backend-selection.md](backend-selection.md)
  step 2 (import hygiene). Related: [wgsl-strategy.md](wgsl-strategy.md) (WGSL registry by program name),
  [backend-split.md](backend-split.md), [core-gpu-compute.md](core-gpu-compute.md).

## Goal

`gl/shader/**` moves to `gl/webgl/shader/**` and nothing outside `gl/webgl/` imports GLSL. Neutral code (renderables,
passes, `canvas3d`, geometry, representations) references programs by a descriptor without sources. Each backend maps
descriptors to its own sources: GLSL in `gl/webgl/`, WGSL in `gl/webgpu/` (exists). A WebGPU-only build then contains no
GLSL, and a WebGL-only build contains no WGSL.

## Current state

- `gl/shader/` has 93 modules (chunks, graphics and pass shaders, compute shaders).
- `ShaderCode` (`gl/shader-code.ts`) carries the sources:
  `{ id, name, vert, frag, extensions, outTypes, ignoreDefine }`. Only the WebGL backend reads `vert`/`frag`
  (`webgl/program.ts`, `webgl/resources.ts`); the WebGPU render items only use `name` and look the WGSL up in
  `webgpu/shader/programs.ts`.
- 27 modules outside `gl/shader/` import GLSL:
  - `gl/shader-code.ts`: the 8 geometry programs (mesh, spheres, cylinders, points, lines, text, direct volume, image);
  - 16 passes: background, bloom, cas, dof, dpoit, fxaa, illumination, marking, multi-sample, outline, postprocessing,
    shadow, smaa, ssao, tracing, wboit;
  - GPU compute: `gl/compute/` (`util`, `grid3d`, `hi-z`, histogram pyramid `reduction` and `sum`, marching cubes
    `active-voxels` and `isosurface`), `geo/gaussian-density/gpu`, `geo/geometry/texture-mesh/color-smoothing`;
  - `gl/_test/shader-uniforms.test.ts`.
- `gl/shader-code.ts` mixes neutral and GLSL-only code:
  - **neutral:** the `ShaderCode` identity (id, name), the define types, `ignoreDefine` (filters defines out of program
    keys);
  - **GLSL-only:** source normalization, the GLSL 100/300 prefixes, `getContextDefineLookup`, the 300 fragment
    transform, `stripUniformDeclarations`/`stripSamplerDeclarations`/`stripAttributeDeclarations`, `addShaderDefines`.
  - The backend-neutral guard test allowlists it ("WebGLParameters, GLSL assembly").
- Program names are not unique:
  - `compose`: multi-sample, illumination `ComposeShaderCode` and illumination `MultiSampleComposeShaderCode` (three
    different sources);
  - `accumulate`: tracing and texture-mesh color smoothing.

  The WGSL registry is keyed by name, so these collide there as well.

- `grid3d` builds its GLSL at runtime from a spec with GLSL snippets; the alpha-orbitals extension uses it
  (`createGrid3dComputeRenderable`). No other code outside graphics creates a `ShaderCode`.

## Design

### Program descriptors

- `ShaderCode` becomes a neutral descriptor without sources: `{ id, name, extensions, ignoreDefine }`.
  - `extensions` stays as feature requirements (frag depth, draw buffers, texture LOD), which the WebGL backend turns
    into extension directives and both backends can check.
  - `outTypes` (WebGL2 fragment output types) moves into the GLSL entry; WGSL declares its outputs in the source.
- Descriptors stay where they are defined today (the passes, the geometry programs in `gl/shader-code.ts`), so neutral
  code keeps its imports.
- Program names are unique, and the name is the key into both registries (plain strings, matching the existing WGSL
  registry).

### Source registries per backend

- `gl/webgl/programs.ts` maps names to `{ vert, frag, outTypes }`, importing from `gl/webgl/shader/**`;
  `gl/webgpu/shader/programs.ts` is the WGSL counterpart (exists).
- The WebGL render items and resources look the GLSL up by name when they create a program, throwing on an unknown name.
- Only the backend contexts import their registries, so the registry of a backend not in a build drops out.
- Custom programs (extensions, apps) register their sources: `registerGlslProgram(name, source)` in `gl/webgl/`, and
  later `registerWgslProgram` for WebGPU.

### Runtime-generated GLSL

- `grid3d` generates its fragment shader from a spec with GLSL snippets, so it is WebGL-only by construction. It moves
  to `gl/webgl/compute/` and registers its generated program under the spec id.
- alpha-orbitals stays WebGL-only (it already narrows with `isWebGLContext`) and stays on the guard test allowlist until
  GPU compute becomes a backend service ([backend-selection.md](backend-selection.md) step 3).

### WebGL-only GPU compute

`gl/compute/**`, `geo/gaussian-density/gpu` and the texture-mesh GPU paths (color smoothing, GPU isosurfaces) only run
on WebGL; callers narrow with `isWebGLContext`. They move under `gl/webgl/compute/` (with their GLSL), and the callers
keep the narrowing. Turning them into backend compute services with CPU fallbacks is backend-selection step 3; this plan
only moves them.

## Steps

1. **Unique program names**: rename `compose` (`multi-sample-compose`, `illumination-compose`,
   `illumination-multi-sample-compose`) and `accumulate` (`trace-accumulate`, `color-smoothing-accumulate`); update the
   WGSL registry keys; add a test that all program names are unique. Renders unchanged.
2. **GLSL helpers into the backend**: move the GLSL-only half of `gl/shader-code.ts` to `gl/webgl/shader-code.ts`
   (prefixes, define lookup, stripping, `addShaderDefines`, normalization); the neutral file keeps descriptors and
   define types; drop the guard test allowlist entry.
3. **Registry and descriptors**: add `gl/webgl/programs.ts` with all graphics and pass programs; the WebGL backend
   resolves sources by name; remove `vert`/`frag`/`outTypes` from `ShaderCode`; pass and geometry modules stop importing
   GLSL.
4. **Move the sources**: `gl/shader/**` moves to `gl/webgl/shader/**`; GLSL-only tests (`shader-uniforms`,
   `minimum-limits`, the shader preprocessor inventory) move to `gl/webgl/_test/`.
5. **WebGL-only compute**: move `gl/compute/**`, `geo/gaussian-density/gpu` and the texture-mesh GPU paths to
   `gl/webgl/compute/`, with `grid3d` registering its generated programs; update the imports of callers and
   alpha-orbitals.
6. **Guard test**: `backend-neutral.test.ts` rejects imports of `gl/webgl/shader/`, `gl/webgl/programs` and
   `gl/webgl/compute/` from outside `gl/webgl/`, except the allowlisted WebGL-only callers (shrinking with
   backend-selection step 3), and imports between `gl/webgl/` and `gl/webgpu/`.
7. **Bundle check**: with backend-selection step 4, assert that a WebGPU-only bundle contains no `gl/webgl/shader/`
   modules and a WebGL-only bundle no `gl/webgpu/shader/` modules.

## Validation

- Every step: WebGL visual test profiles byte-identical, the `webgpu` profile unchanged, full type check, graphics unit
  tests.
- Steps 1 and 3: the program cache keys change (names), so check that program counts and the shader-manager warm-up
  (`linkPrograms`) are unchanged.

## Breaking changes (record in [breaking-v6-changes.md](../breaking-v6-changes.md))

- Import paths `@molstar/graphics/gl/shader/*` move to `@molstar/graphics/gl/webgl/shader/*`.
- `ShaderCode` no longer carries `vert`, `frag` and `outTypes`; custom programs register their GLSL with
  `registerGlslProgram`.
- `gl/compute/*` (including `grid3d`) moves to `gl/webgl/compute/*`.
- Renamed program names (`compose`, `accumulate`), visible in debug output and program statistics.

## Open questions

- Registry keys: name strings (proposed, like the WGSL registry) or the descriptor objects themselves (no name
  uniqueness needed, but registries must import the descriptors).
- Whether step 5 waits for the compute services of backend-selection step 3 instead of moving the modules twice.
