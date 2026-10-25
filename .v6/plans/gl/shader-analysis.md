# Plan: Shader analysis, GLSL and WGSL, for correctness and performance

## Status

- **Current**: steps 1–3 done but the WebGL runs' length (the variant matrix, the WGSL and GLSL compile walks). A
  systematic pass over every program on both backends: compile every reachable variant, review both languages side by
  side against a hazard list, then measure and analyse the passes that cost the most GPU time. Findings become fixes
  with tests, or measured changes recorded here.
- **Layout** (`tests/src/shader-variants/`, the `tests` package: Node and the browser share it):
  - `matrix.ts`: the sources (`buildSources`), one source compiled alone (`compileSource`: its precompile list, then its
    props as frames of the real scene, `prepareFrames`), all of them for coverage (`walkVariants`), the kernels
    (`prepareAllKernels`), and what a context compiles (`recordCompiles`: WebGPU, Tint's messages of every shader module
    and the pipelines; WebGL, every program asked for and the compile and link logs, queried late so that compiling
    stays parallel). Files through a reader (Node: the file system; the browser: the served repository);
  - `_test/`: Dawn (`webgpu-shader-variants.test.ts`) and native `gl` (`webgl-shader-variants.test.ts`), a test per
    source (`-t '<source>'` runs one), the kernels, the coverage. On demand only, never in a lane or CI (minutes per
    source): `SHADER_VARIANTS=1 pnpm native:run -- pnpm test:unit tests/src/shader-variants`;
  - the visual suite's on-demand cases (`tests/src/visual/cases/shader-variants.ts`, run when `--filter` or `--tags`
    select cases): `shaders/<source>`, `shaders/kernels`, `shaders/coverage`, on a context of the profile alone
    (`BackendProfile.createContext`: nothing compiled with it before). Coverage is required on profiles with every
    capability (WebGPU with `floatBlend`; WebGL2 with float targets, several targets and fragment depth).
- **The sources**: the presets (every applicable structure preset on a protein, a nucleic acid and a coarse model,
  every structure and volume representation type with its defaults), the MolViewSpec documents of the repository
  (snapshots made offline: downloads pointed at a missing file), the Mesoscale Explorer's spacefill per graphics mode
  and clip variant with its canvas props (rebuilt from `Spheres.LodPixelPresets`: apps can't be imported); the axes:
  each boolean define flipped and each string define through the values the sources reach for the kind (another option
  is a variant no object draws), with what a value needs (a palette, colors per group; clip objects, a count), in one
  transparency mode; each canvas option alone (nested groups too), plus dimming, illumination with multi-sampling and
  sharpening with each antialiasing. Every registered program and kernel is compiled (SMAA in the browser only: its
  lookup tables are `HTMLImageElement`s).
- **Findings**:
  - fixed: the precompile list of a MolViewSpec volume failed (`clip: undefined` in its params overwrote the default;
    an undefined param is the default now, `plugin/core/src/state/_test/precompile.test.ts`);
  - fixed: WebGL program items (precompile) were made with the schema as given, not as render items adapt it: no data
    textures (the accessors missing: any per-group color failed to compile), WebGL1 loop bounds as uniforms (GLSL ES
    1.00 needs constant bounds: clipped objects failed), no instance data of merged objects; precompiled programs
    differed from the render items' (`getProgramItemSchema`, `gl/_test/precompile.test.ts`);
  - fixed: GPU culling under the core limits: its instance copy binds 10 storage buffers, a stage may 8 there
    (`ComputeContext.isKernelSupported(code, schema)`, `canCullOnGpu`; `webgpu-core-limits.test.ts`);
  - no WGSL warning anywhere; the only WGSL errors were in define combinations no object draws (the axes then);
  - open: `Canvas3D.precompile` plans no Hi-Z pass (the Mesoscale Explorer turns it on), and misses `dim` and some pass
    variants `prepare` of the real scene makes (`postprocessing`, `fxaa`, `overlay`);
  - open: `prepare({ props })` plans Hi-Z with the canvas's props, not the given ones;
  - open: a headless canvas never syncs its viewport to Hi-Z (no resize): Hi-Z stays off until `handleResize`;
  - deferred, planned (WebGL runs, see WebGL runs below): the WebGL cases take long (about 95 ms a program on ANGLE
    D3D11, the define axes alone 110 s on WebGL2) and stall the visual runner (5 min a page) when every profile runs at
    once; no compile failure was found in the runs that finished;
  - to check: with SMAA unsupported (Node), sharpening (CAS) was not drawn either.
- **Hazard review** (step 6, WGSL against GLSL; every finding is in both languages, none is a port regression):
  - fixed, both languages:
    - SSAO blur: 0 / 0 at a blur depth bias of 0, the slider's minimum: no sample, not even the pixel itself, is within
      the bias, and SSAO disappeared. A bias of 0 now skips the blur passes (`postprocessing/occlusion/no-blur`), and a
      pixel with no sample within the bias is left unblurred.
    - Image: per-group transparency was read at group -1, for pixels without a group.
    - Color smoothing: a vertex of a cell with no grouped corner unpacks group 0 to -1, which wraps to `0xFFFFFFFF` in
      WGSL. Such a vertex now adds no color.
    - Path tracing, Russian roulette: `1 / p` with `p == 0` (a black hit, and `randomFloat` returning 0). The NaN stays
      in the accumulation, and the denoiser spreads it over 13×13 pixels.
    - Marching cubes: `normalize` of a zero gradient. A flat neighbourhood now gives no normal, which the callers
      already check for.
    - Cylinders: the two colors per group were mixed only in the color variants, so tracing and emissive got the first
      color only.
    - CAS: the limiter's ring left out the left neighbour `d` (FSR's RCAS uses `b, d, f, h`).
    - 0 / 0 that comes out right only where NaN is dropped (WGSL leaves it indeterminate), now guarded:
      - CAS at flat, black and white texels;
      - DOF at a blur spread of 0 when every sample is in focus (`postprocessing/dof/no-spread`);
      - shadows with no light.
  - open:
    - Spheres and cylinders: the cap depth bias `1 - 1e-7 / r` (and `2e-7 / r`) is lost to f32 rounding near 1
      (reversed Z), so for larger radii it is no bias at all.
    - Direct volume: group -1 from a packed group of 0 is used unguarded (reads and pick).
    - `pow(0, y)` (denoiser normal weight, x-ray falloff) is defined in GLSL and not in WGSL; `pow(0, 0)` is undefined
      in both.
  - latent, recorded: f32 indices past 2^24 (`packIntToRGB` too); the gaussian density's u32 index for huge grids;
    radius 0; an end-on cylinder or line, and a sphere quad with the camera inside; the cell cull's 1D dispatch above
    4.19M cells; SSAO 0 / 0 at an equal sample depth, and a tangent of a normal parallel to the kernel rotation; the
    shadow march with `w ≤ 0`; `normalize(normal + randomUnitVector)` at the exact opposite.
  - checked, no finding: reversed Z, Y orientation, uniformity (no derivatives or implicit LOD anywhere), edge reads
    (every `read_` clamps), WBOIT and DPOIT conventions, marking and pick packing, compose's premultiplied alpha.

## WebGL runs (planned)

The time of a WebGL run is the compilers', not the checks' (querying the status late changed nothing): ANGLE on D3D11
translates the GLSL in `compileShader` (fast) and compiles the generated HLSL in `linkProgram` (about 95 ms a program).

- **Compile only, once per source**: the per-source cases compile each distinct vertex and fragment source (after
  preprocessing; many define sets give the same text) once, with `COMPILE_STATUS` and the log, and link nothing. Every
  GLSL error found so far was of this stage (undeclared names, no matching overload, loop bounds of GLSL ES 1.00).
- **Link in the coverage case**: one program per distinct pair of sources, linked, the WebGL profiles one at a time.
  What only linking finds: mismatches between the stages (varyings, uniforms declared differently), the limits checked
  at link (uniform vectors, varyings, attributes: the minimal profiles), the D3D compiler's failures and timeouts on the
  generated HLSL, driver bugs.
- **Not `getError`**: compile and link failures are no GL errors (a status and a log); an error comes later, from using
  the program, with no variant or line.

## GPU time per pass (step 7)

`pnpm benchmark -- --filter 'passes/*' --gpu-timing` (`tests/src/benchmark/passes.ts`): 7qpd zoomed in to fill about
half of a 2560x1440 canvas, the geometry alone with the opaque pass only, and one pass at a time added to spacefill;
the direct-volume cases (`render/direct-volume/*`). GPUs by browser flags: `--browser-arg
--use-webgpu-power-preference=force-high-performance --browser-arg --force_high_performance_gpu` (or `force-low-power`
and `--force_low_power_gpu`); the adapter line of the output names the one used.

- **The AMD iGPU (Radeon, RDNA 2) decides the ranking**: its times are stable to a few percent across cases. The RTX
  4070 Ti's are not: the same opaque pass took between 0.48 and 1.17 ms in different cases on WebGPU. Other work
  loaded that GPU during the run, and light loads also let it clock down. To measure again on an idle GPU, with
  heavier scenes or clocks locked (`nvidia-smi -lgc`, as administrator).
- **Ranked by GPU time on the AMD iGPU, in ms** (WebGPU / WebGL2; NVIDIA in brackets, as an indication only):

  | Pass                                        | AMD           | NVIDIA      | Notes                                                    |
  | ------------------------------------------- | ------------- | ----------- | -------------------------------------------------------- |
  | illumination compose (denoiser)             | 153 / 156     | 4.9 / 2.4   | a 13×13 window with an `exp` per sample, at full size    |
  | illumination trace                          | 50 / 167      | 2.8 / 1.2   | WebGPU 3.4× faster on AMD, 2× slower on NVIDIA           |
  | depth of field                              | 57 / 56       | 2.3 / 1.0   | 81 taps of color and depth per pixel (`blurSize` 9)      |
  | direct volume (`wboit transparent`)         | 12–83 / 10–94 | 0.3–1.0     | WebGPU 20–50% slower on AMD, except the Gaussian volume  |
  | Gaussian surface, opaque                    | 43 / 32       | 0.6 / 0.6   | WebGPU 35% slower on AMD                                 |
  | spheres (impostors), opaque                 | 37 / 30       | 0.5 / 0.5   | WebGPU 26% slower on AMD                                 |
  | DPOIT transparent (two peels)               | 41 / 59       | 2.6 / 1.3   | WebGPU faster on AMD; WebGL2 also clears for 6 ms        |
  | SSAO                                        | 18 / 17       | 0.2 / 0.2   | multi-scale is cheaper (12 / 10): half-size levels       |
  | marking mask (the scene drawn again)        | 17 / 14       | 0.6 / 0.4   |                                                          |
  | WBOIT transparent (molecular surface)       | 12 / 12       | 0.6 / 0.4   |                                                          |
  | cylinders (impostors), opaque               | 8.7 / 8.6     | 0.1 / 0.2   |                                                          |
  | outlines                                    | 8.2 / 5.5     | 0.1 / 0.1   | WebGPU 48% slower on AMD                                 |
  | mesh (cartoon), opaque                      | 4.5 / 4.3     | 0.1 / 0.1   |                                                          |
  | volume composite                            | 3.3 / 2.1     | 0.1 / 0.04  | fixed per frame, WebGPU 55% slower                       |
  | SMAA weights                                | 1.3 / 3.5     | 0.05 / 0.1  | WebGL2 2.6× slower                                       |
  | FXAA, CAS, shadows, bloom, the compose pass | 0.5–2.2 each  | ≤ 0.1       |                                                          |

- **For step 8**, in this order:
  - the denoiser (by far the largest);
  - depth of field;
  - the gaps where WGSL is slower than GLSL on the same GPU: direct volume, Gaussian surface and sphere impostor
    fragments, outlines, volume composite;
  - illumination trace on NVIDIA, where WebGPU is 2× slower (to confirm on an idle GPU).

## Step 8: the hot passes

Measured on the AMD iGPU (`passes/illumination`, compose and trace in ms, WebGPU / WebGL2). The step 7 run was
disturbed by load on the other GPU: measured again, illumination trace is 37 / 70 (not 50 / 167) and compose 137 / 149.

- **The denoiser** (`smartDeNoise`, illumination compose): 169 samples per pixel, each reading an rgba32float color and a
  packed normal.
  - adopted: a background pixel returns at once. Its normal is zero, so every sample's normal weight is 0 and the
    result was the pixel itself after 169 samples. Compose 137 / 149 → 31.6 / 27.9, output unchanged.
  - adopted: the spatial weights from a table of 7, and the normal weight's 6th power by multiplication (no `pow`,
    whose `pow(0, y)` WGSL leaves undefined). 31.6 / 27.9 → 30.5 / 28.2, output unchanged.
  - not taken: skipping the color read of a sample with no normal weight. A branch per sample, slower (33.7 / 29.4).
  - measured, not taken: the circular kernel of the original (radius 6, 113 samples) changes the output. 30.5 → 26.7
    (WebGPU).
  - what is left: compose without the denoiser is 1.8 ms, so the denoiser is still about 29 ms. It is bound by texture
    reads: the same texels are read by up to 169 neighbours. Next candidate: a compute pass (WebGPU) that loads a tile
    and its 6-texel border into workgroup memory once.

## Today

- **Sources**:
  - GLSL: 36 programs and 35 chunks, about 6.3k lines (`gl/webgl/shader/`, `chunks/`, plus the `bloom`, `compute`,
    `histogram-pyramid`, `illumination`, `marching-cubes`, `marking`, `smaa` and `util` subfolders);
  - WGSL: 35 programs and 12 chunks, about 8.1k lines (`gl/webgpu/shader/`).

  Both go through the shared preprocessor ([wgsl-strategy.md](wgsl-strategy.md)). Variants come from defines (render
  variants, color, size and layer types, clipping, LOD, lights).
- **What checks them now**:
  - `shader-preprocessor.test.ts`: every shader in every render variant preprocesses (mesh quick, the rest medium);
  - `shader-uniforms.test.ts`, `wgsl-declarations.test.ts`, `data-accessors.test.ts`: identifiers, declarations,
    accessors;
  - `precompile.test.ts` and the Dawn precompile tests: the signatures scenes draw compile on WebGPU;
  - the visual suite: images on 7 profiles, with WebGPU compared against WebGL2.

  Nothing compiles every reachable variant on every backend. Nothing compares the two languages program by program.
  Nothing measures a shader's cost apart from its pass.
- **Measurements available**: GPU time per pass (`pnpm benchmark --gpu-timing`, [gpu-timer.md](gpu-timer.md)) on the
  benchmark scenes, on this machine's NVIDIA RTX (discrete) and AMD RDNA 2 iGPU
  (`--browser-arg --force_low_power_gpu`).

## Goals

- **Correctness**: every variant a scene can reach compiles on every backend and profile, with no WGSL warnings.
  GLSL and WGSL compute the same thing. Numerical hazards are found and fixed.
- **Performance**: the passes that cost the most GPU time are known per GPU, and their shaders are analysed (instruction
  and fetch counts, registers and occupancy, branches, overdraw and early-z). Changes are adopted only when measured on
  both GPUs.

## Phase 1: compile every reachable variant

1. **The variant matrix**: not the cross product of all defines, but the variants representations can reach. Take them
   from the precompile derivation (`getTreePrecompileList`, signatures from params) over the param trees of three
   sources: the default presets, the MolViewSpec examples and mesoscale scenes. Add each define axis varied alone from
   the mesh baseline. Output: program × variant × backend.
2. **WGSL**: compile each variant in Dawn (Tint) and collect `getCompilationInfo` errors and warnings, uniformity
   analysis included. One test that walks the matrix, on demand.
3. **GLSL**: compile and link each variant with the native `gl` package (WebGL1 and WebGL2 contexts, on demand), and
   in the browser with a visual-suite `check` per profile, so that `webgl1-minimal` and `webgl2-minimum` meet their
   limits (loop bounds, uniform and varying counts, dynamic indexing).
4. **Report**: a table of failures and warnings per backend, each turned into a fix with the test that found it.

## Phase 2: correctness review

1. **Parity, program by program**: preprocessed GLSL and WGSL of the same variant side by side (a script printing both
   for a program and its defines). What they compute should match:
   - the defines each handles;
   - constants;
   - depth conventions (reversed-Z, [depth-conventions.md](depth-conventions.md));
   - packing and unpacking (ids, depth);
   - the y orientation of texture reads;
   - branch structure.

   Differences are listed: intended (a backend's own path) or a bug.
2. **Hazards**, checked in every program in both languages:
   - division by zero, `normalize` of a zero vector, `sqrt`, `log`, `pow` and `acos` outside their domain, `atan`
     of (0, 0): NaN that reaches blending or depth;
   - indices and counts as floats beyond 2^24 (WebGL1's float math in the data-page accessors, group and instance
     ids), integer and float conversions, `packIntToRGB` precision;
   - out-of-bounds texel fetches and storage reads, clamping at grid and page edges;
   - derivatives and implicit-LOD samples in non-uniform control flow (an error in WGSL, undefined in GLSL), and
     after `discard`;
   - uninitialized variables, loop bounds;
   - writes to `frag_depth` that disagree with the depth test or the pick depth.
3. **Findings**: severity, the program and variant, how to reproduce. Each fix comes with a visual case, a Dawn test
   or a native GL test that fails without it.

## Phase 3: performance analysis

1. **Measure first**: GPU time per pass on the benchmark scenes, NVIDIA and AMD, WebGL2 and WebGPU. Ranked, this picks
   the passes worth analysing. Likely candidates:
   - sphere and cylinder impostor fragments;
   - mesh lighting;
   - SSAO, DOF and the postprocessing compose;
   - the OIT passes;
   - illumination and tracing.
2. **Cost of the hot variants on both vendors**, with tools run by hand, not dependencies:
   - the generated backend code: Dawn's shader dump (`dump_shaders`; HLSL, or SPIR-V where available) and ANGLE's
     translated HLSL for GLSL;
   - AMD: Radeon GPU Analyzer compiles that code offline: registers, occupancy and instruction counts per variant;
   - NVIDIA: Nsight Graphics captures a frame of Chrome (GPU Trace, the shader profiler) and measures on the hardware:
     registers, occupancy, warp stall reasons and their hot instructions, per pipeline (NVIDIA has no offline
     compiler).

   Per variant: ALU and texture instruction counts, registers, occupancy, and on NVIDIA where its warps stall.
3. **Checklist per hot shader**:
   - texture fetches per fragment, and dependent reads;
   - work per fragment that could be per vertex or per object;
   - repeated `normalize` and transforms;
   - uniform loads inside loops;
   - loops unrolled or not (the clip loop is no longer unrolled);
   - branches on uniforms against defines ([fewer-defines.md](fewer-defines.md) measured them);
   - `discard` and `frag_depth` writes that disable early depth testing;
   - overdraw (fragments shaded per pixel, from the GPU timing of a depth pre-pass);
   - precision (half precision paused in [shader-f16.md](shader-f16.md)).
4. **Changes**: each measured A/B with the benchmark (interleaved rounds, both GPUs, `--filter` on the affected cases),
   adopted only when faster without regressions elsewhere, and recorded here with its numbers, as are those not taken.

## Phase 4: guard rails

- The phase 1 compile walk stays as on-demand tests (`SHADER_VARIANTS=1`, never in a lane or CI): WGSL in Dawn, GLSL
  in native `gl`, and the visual suite's on-demand cases. WGSL warnings fail them.
- The parity script stays for reviews of shader changes: a change to one language prints the other's variant beside it.

## Steps

1. The variant matrix from the precompile derivation and single-axis variation. A test that it covers every program.
   Done (see Status).
2. The WGSL compile walk in the Dawn lane: errors and warnings. Done.
3. The GLSL compile walk: native `gl` (medium; done, runs where `gl` installs) and a visual-suite check per profile
   (done; the WebGL runs' length open).
4. Fixes for what 2 and 3 find.
5. The parity script, and the review of every program pair. Findings listed here.
6. The hazard review in both languages. Findings and fixes with tests.
7. GPU time per pass on both GPUs and backends, ranked. Done (GPU time per pass, above).
8. The top passes' hot variants analysed: shader dumps, Radeon GPU Analyzer (AMD), Nsight Graphics captures (NVIDIA).
   Candidates listed.
9. Candidates measured and adopted or recorded.

## Decisions

- **No external validator as a dev dependency** (glslang, naga): the backends' own compilers (ANGLE, Tint in Dawn)
  decide, and they compile what browsers run.
- **Mobile GPU behaviour is out of scope** (Mali, Adreno, Apple; `mediump` precision): no device here to check it on.
- **Reachable variants** come from three sources: the default presets, the MolViewSpec examples and mesoscale scenes.
  More sources are added where a finding shows a gap.
