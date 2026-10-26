# Plan: Shader analysis, GLSL and WGSL, for correctness and performance

## Status

- **Done**: every step done or decided (see Steps). A systematic pass over every program on both backends: compile
  every reachable variant, review both languages against a hazard list, then measure and analyse the passes that cost
  the most GPU time. Findings became fixes with tests, or measured changes recorded here. Left open: the WebGL runs'
  length (planned, WebGL runs), the cap depth bias and the latent hazards (Findings), the AMD iGPU's gaps (not
  pursued).
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
  - fixed: `Canvas3D.precompile` plans Hi-Z with each props entry's (also for the GPU's test where the list has objects
    culled on the GPU), and marking with a stand-in scene that has something marked, so the dim pass too; the marking
    pass takes the renderer props given (`colorMarker`, `dimStrength`) over the live renderer's. The coverage walk's
    programs first compiled by preparing or drawing frames: none (were Hi-Z, `dim`, `postprocessing`, `fxaa`,
    `overlay`);
  - fixed: `prepare({ props })` plans Hi-Z with the given props (`HiZPass.plan(b, camera, props, gpuCulling)`);
  - fixed: Canvas3D syncs its viewport to the renderer, picking and Hi-Z when created and restored, not only on a
    resize (a headless canvas never had Hi-Z levels; `plugin-headless/src/_test/webgpu-prepare-hiz.test.ts`);
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
  - fixed, both languages, output unchanged in the visual suite:
    - x-ray shading: `d^falloff` at the silhouette (`d = 0`) is its limit, 0 (1 for a falloff of 0, the parameter's
      minimum): WGSL's `pow` is undefined at 0, both languages' `pow(0, 0)`. The denoiser's `pow(·, 6)` is multiplied
      out (above).
    - direct volume: a cell without a group (packed 0 unpacks to -1: no atom within reach, density about 0) takes group
      0, as the clamped reads of index -1 did.
  - open:
    - spheres and cylinders: the cap depth bias `1 - 1e-7 / r` (and `2e-7 / r`) of the near-plane and clip-object caps.
      It was `0 + 1e-7 / r` before reversed Z ("spheres with a greater radius occlude smaller ones"); near 1, f32 steps
      are 6e-8, so for radii above about 1.7 it rounds away. Tried: a bias of float steps (`depth × (1 − 1.2e-7 ×
      (1 + 64 / (1 + r)))`), and a case of two overlapping solid spheres cut by one clip plane; old and new render
      alike, the smaller sphere's cap shows inside the larger's in both, so the caps there are not the tie the bias
      decides. Not changed without a scene that shows the problem.
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

- **The AMD iGPU (Radeon, RDNA 2) decides the ranking**: its times are stable to a few percent across cases. The first
  NVIDIA run was disturbed by other load on that GPU (the same opaque pass 0.48–1.17 ms); its column below is from
  that run, the idle run follows the table.
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

- **NVIDIA, idle** (RTX 4070 Ti, ms, WebGPU / WebGL2; stable: the spacefill opaque pass 0.49–0.53 in every case):
  every pass about equal or faster on WebGPU. Opaque spacefill 0.52 / 0.48, mesh 0.10 / 0.12, Gaussian surface
  0.56 / 0.58; SSAO 0.23 / 0.18; outlines 0.13 / 0.09; DOF 0.27 / 0.25 (separable on both: rows 0.17 / 0.12,
  columns 0.10 / 0.13); WBOIT transparent 0.25 / 0.38, DPOIT transparent 1.31 / 1.35 (WebGL2 also clears for
  0.21); illumination trace 0.76 / 1.19, compose 0.11 (with the denoise kernel 0.39) / 0.53; SMAA weights 0.04 /
  0.10, CAS 0.03 / 0.07; direct volumes within 10 %; volume composite 0.08 / 0.04. Frames carry about 0.1 ms more
  on WebGPU outside the timed passes.
- **The AMD iGPU's gaps are not pursued**: on NVIDIA none of them exists (above); the iGPU is a minor target.
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
  - compose without the denoiser is 1.8 ms, so the denoiser was still about 29 ms.
  - adopted (WebGPU): the denoiser as a kernel (`gl/compute/denoise.ts`): tiles of 16 × 8 pixels load their texels
    with a border of 6 into workgroup memory once, the normal decoded once per texel; the same weights in the same
    order, and the compose program reads its output (rgba16f). WebGL keeps the program's denoiser. Interleaved A/B,
    two rounds alike: NVIDIA 0.95 → 0.48 ms (kernel 0.37, compose 0.11: 2.0×), the AMD iGPU 30.2 → 26.0 ms (24.8 +
    1.2: 1.16×). The iGPU (2 compute units, about 0.5 TFLOPS) is bound by the arithmetic of 169 samples, not by
    reads; there only fewer samples help (the circular kernel above, 113).
- **NVIDIA (Nsight Graphics GPU Trace, by hand)**: `ngfx.exe --activity "GPU Trace Profiler"` launching Chrome on a
  served benchmark page (`?backend=webgpu&runs=5000&filter=<case>`, so that the case still draws when the trace
  starts). What it takes:
  - no `--platform` option: Qt takes it as its own (the window-system plugin) and ngfx fails to start;
  - Chrome with `--in-process-gpu` (the tool is injected into chrome.exe and sees the GPU process's D3D12 work only
    there), `--disable-gpu-compositing --disable-direct-composition` (else the window stays white and nothing is
    submitted), a fresh `--user-data-dir`, and the GPU flags of step 7;
  - `--no-timeout`, `--set-gpu-clocks base`, `--architecture Ada --metric-set-name "Throughput Metrics"`,
    `--real-time-shader-profiler`, `--auto-export`, and an existing `--output-dir`. The exported tables cover the
    whole trace: Dawn emits no PIX markers, so passes are not named; per-shader rows use Dawn's hashed entry points
    (also with `use_user_defined_labels_in_backend`). Shader names, registers and hot lines are in the `.ngfx-gputrace`
    report, read in the Nsight UI.
  - Chrome's `--enable-dawn-features=disable_symbol_renaming,dump_shaders` with `--enable-logging --log-file=...`
    writes each shader's WGSL, generated HLSL and DXIL with our names to the log: the code the hot shaders became.
    `emit_hlsl_debug_symbols` adds source lines for the UI's shader profiler; the DXIL is then unoptimized (the driver
    optimizes it again: the counters moved a few percent), so it is not for measuring.
- **`passes/illumination` on the RTX 4070 Ti** (WebGPU, base clocks, 500 ms): GPU active 90 %, SM throughput 33 %, not
  memory bound (L1 hit 84 %, L2 91 %, DRAM 10 %). Pixel shader warp launches stall on register allocation 27 % of the
  time, and 58 % of SM-active time has idle warp slots: occupancy is limited by registers. The busiest pipe is XU
  (transcendentals: `exp`, `rsqrt`, `sqrt`) at 25 %, above ALU 19 % and FMA 17 %. Two shaders take nearly all warps
  (11–13 % and 6 % of slots): the denoiser (compose, `smartDeNoise`) first, tracing (`trace`) second, as the UI's
  source view of each hash shows. The front end waits for idle 81 % of its cycles: passes are drained one after the
  other. The denoiser decodes a packed normal (octahedral, `normalize`: an `rsqrt`) and takes an `exp` for each of its
  169 samples, so a texel's normal is decoded by every neighbour: the case for the tiled compute denoiser on NVIDIA too
  (one decode and fetch per texel into workgroup memory).
- **Path tracing** (`trace`): on an idle NVIDIA GPU WebGPU is faster than WebGL2 (0.77 against 1.18 ms, two rounds;
  step 7's "2× slower" was the disturbed run), so there is no WGSL-specific cost. Measured, not taken: the ray march
  with an incremental projection (the clip position advanced by the projected step instead of a matrix product per
  step) and its invariants (bounds, the depth texture's size) hoisted out of the step loop: the AMD iGPU 37.5 → 36.7
  ms (2 %), NVIDIA unchanged; the compilers hoist them already. Its cost is the work itself (steps × renders per frame
  × bounces); the next look needs the shader profiler's hot lines (Nsight UI, `-Symbols`).
- **Sphere impostors** (`passes/spheres`, the opaque pass; two rounds each, alike): the AMD iGPU 37.0 (WebGPU) against
  29.4 ms (WebGL2), +26 %; meshes +5 % (4.5 / 4.3); NVIDIA equal (0.5 / 0.5). So AMD's D3D12 path against its D3D11
  path, for a fragment shader that writes depth and discards. Ruled out:
  - strict math (Dawn's default, `!dx.precise` on every operation): `--strict-math off` gives 37.0 → 35.1 ms, a
    quarter of the gap;
  - the demotion of `discard` (WGSL's semantics): Dawn's HLSL discards directly;
  - the compiler: FXC (`--disable-dawn-features=use_dxc`) instead of DXC, 37.1 / 37.3 ms;
  - the depth format: `depth32float-stencil8` in both (WebGL2's `DEPTH32F_STENCIL8`).

  - not the D3D12 path: Dawn's D3D11 backend (`--use-webgpu-adapter=d3d11`) is slower still, 42.7 ms (frame median;
    FXC at `OPTIMIZATION_LEVEL0` and IEEE strictness), where ANGLE takes 30.0 on the same driver. So it is the code
    our WGSL becomes against what ANGLE makes of our GLSL.
  - in the lighting: `passes/spheres-ignore-light` 16.3 / 14.4–15.7 ms, `passes/spheres-approximate` 13.1 / 11.0–11.3;
    lit 37.2 / 29.8. The bump branch (`perturbNormal`, never taken) is 0.7 ms of it.
  - the AMD ISA (Radeon GPU Analyzer 2.15, gfx1035; Dawn's HLSL and flags from `dump_shaders`, ANGLE's HLSL from
    `WEBGL_debug_shaders` wrapped as ANGLE links it): Dawn 38 VGPRs, 4420 bytes, 33 branches; ANGLE 46 VGPRs, 7852
    bytes (the bump noise unrolled), 14 branches. Dawn's arithmetic is unfused under IEEE strictness (247 `v_mul`, 141
    `v_add`, 7 `v_fma`; ANGLE 370 fused `v_fma*`); without `-Gis` Dawn's fuses too (142 `v_fma*`, 13 `v_add`). Yet
    strict math off (the sphere shader then compiles relaxed, checked in the dump) gains only 5 %: the remaining
    ~6 ms are not in what static counts show.

  Next: Radeon GPU Profiler's instruction timing of Dawn's shader on the iGPU (D3D12 only, so no ANGLE side), or a
  bisection inside the lit path with the same change on both backends. Capture tools (by hand, not in the
  repository): a Playwright script collecting both backends' shaders of a benchmark case, a splitter of Dawn's dump.
- **Sphere impostors on NVIDIA** (the PetWorld model, `merging/postsyn-gpu-hiz-zoom-*`, 3.68M atoms; two rounds
  alike): WebGPU is 16–27 % faster than WebGL2 in every case (opaque pass, ms: 640×360 10.4 / 14.2, rotating
  13.1 / 16.1, 2560×1440 24.0 / 28.8, LOD ultra 87 / 112). Not bound by fragments: 16× the pixels cost 2.3× the time,
  no lighting changes nothing (13.0), the approximate impostors are slower (15.4), the LOD preset moves it 20×. GPU
  Trace: GPU active 99 %, SM throughput 62 %; the vertex stage stalls on allocating its output memory (ISBE) 47 % of
  the time, the primitive engine's L2 traffic at 77 %: bound by the geometry of millions of quads.
  - measured, not taken: fewer vertex outputs (the unread model and view positions left out, radius and center in one
    flat output: from 7 slots to 4): no gain, NVIDIA 13.1 → 13.3, 23.9 → 24.6 ms; the driver prunes unread outputs.
    A triangle per sphere instead of a quad (3 vertices, more fragments) did not help much when tried before; fewer spheres (LOD) is what moves it.
  - open, latent (GLSL and WGSL): spheres, cylinders, lines, points and text assign their color and layer varyings
    before `vModelPosition`, so grid colors and grid layers sampled at the model position would read it unset (only
    meshes assign it first); today only meshes get those.
- **Depth of field**: a box blur weighted per sample, `blurSize`² taps of color and two depths each.
  - adopted: separable where half-float targets render. A horizontal pass writes a row's weighted color sums and
    weights (rgba16f), the DOF pass adds a column of them; `2 × blurSize` taps instead of `blurSize²`. DOF 57 / 56 →
    10.0 / 9.0 ms (6.7 + 3.3, 6.2 + 2.9). Other profiles keep the single pass (`dDofPass`).
  - adopted: taps centered (`i - (blurSize - 1) / 2`). With `i - blurSize / 2` every tap of an odd size sat on a texel
    edge, the kernel was half a texel off center, and rounding (which differs between shader variants) chose the
    texel: the split differed by up to 6/255 at edges. Centered, the split matches the single pass to 1/255.

## After the plan: large scenes and GPU surfaces

Profiled on the RTX 4070 Ti (and the AMD iGPU) after the steps above, as data for what comes next (the compute
rasterization of small spheres, the GPU surfaces).

- **Sphere impostors at scale** (the PetWorld model, `merging/postsyn-gpu-hiz-zoom-rotating`, 1280×720, 0.92 Mpx;
  GPU Trace, per frame from the trace's sums over its 32 frames: 19,753 draws at 620 a frame):

  | Stage                                                 | Per frame | Per pixel |
  | ----------------------------------------------------- | --------- | --------- |
  | sphere quads submitted (75 M triangles)               | 37 M      | 40        |
  | fragments rasterized                                  | 63 M      | 68        |
  | pixel shader threads launched (7.0 M warps × 32)      | 225 M     | 245       |
  | fragments reaching the depth test                     | 46 M      | 50        |
  | fragments passing it                                  | 3.25 M    | 3.5       |

  - the drawn spheres are pixel-sized: 1.7 fragments per sphere;
  - overshading 3.6× (2×2 quads on triangles of a pixel or two), and 93% of fragments fail the late depth test
    (`frag_depth` disables early depth): a visibility buffer would shade about 0.9 M pixels where 225 M threads run;
  - the front end carries it: 37 M quads, about 148 M vertex invocations, the vertex output memory (ISBE) stalled 47%.
- **64-bit atomics**: `atomic-vec2u-min-max` is offered by Chrome 157 (Beta) on both GPUs here, not by 156 (stable);
  both adapters also offer `subgroups`, `shader-f16`, `maxComputeWorkgroupStorageSize` 32 KiB.
- **GPU surfaces** (`mc/volume/*`, `mc/gaussian/*` with `--gpu-timing`, the kernels timed per pass,
  `withGpuPasses`), WebGPU, ms:

  | Case                                  | NVIDIA GPU | NVIDIA wall          | iGPU GPU | iGPU wall           |
  | ------------------------------------- | ---------- | -------------------- | -------- | ------------------- |
  | gyroid 256³ marching cubes (10.4 M △) | 9.1        | 11.2 (+21.8 upload)  | 40.6     | 45.3 (+21.6 upload) |
  | blobs 256³ marching cubes             | 0.9        | 3.7 (+13.7 upload)   | 36.6     | 40.9 (+21.9 upload) |
  | globule 100k atoms, 0.5 Å, density    | 4.8        | 20.4                 | 289      | 316                 |
  | same, marching cubes                  | 1.3        | 3.5                  | 64       | 67                  |
  | globule 1k atoms, 1 Å, density + mc   | 0.02 + 0.02| 4.1 + 3.2            | 0.6 + 0.2| 3.3 + 3.2           |

  - on NVIDIA the kernels are a small part: a fixed latency of about 3 ms a stage (waiting for the GPU, the count's
    readback), about 15 ms of CPU setup of the density for 100k atoms, the volume's upload (`writeBuffer`, about
    2.7 GB/s). The density benchmark keeps the GPU busy 12.5% of the time (GPU Trace); the kernel stalls on memory
    (L1 texture-path scoreboard, texture throttle), as a scatter with atomics does;
  - on the iGPU the kernels are the cost: density 289 ms for 100k atoms at 0.5 Å (60× NVIDIA, about the gap in
    throughput), marching cubes 30–40 ms at 256³.
  - **The stage latency is Chrome's**: a completion notice for real work (`onSubmittedWorkDone`, or a 4-byte copy and
    `mapAsync`) takes about 2.5–3 ms on both GPUs here (0.1 ms after an empty submit), the GPU process noticing the
    work finished. The app does not wait for it: GPU surfaces are drawn indirectly with the count on the GPU before
    anything is read back (`compute-isosurface.ts`); the benchmark's stages do (`StageTimer.stage`).
  - **The density's CPU setup** (`density cpu`, `StageTimer.cpuStage`, NVIDIA, 20 runs): 0.4 ms for 1k atoms, 1.4–1.6
    for 7qpd and 10k, 6.85 for 100k. The atoms are walked three times: collected (`getStructureConformationAndRadius`,
    per element the conformation's `x`, `y`, `z` and the size theme), repacked (`prepareGaussianDensityData`), and
    interleaved for the kernel (`computeGaussianDensity`).
    - adopted: the coordinates read directly from the unit's arrays where its operator is the identity (the mappings
      are three classes, their methods slow to call per element), and a range of indices indexed directly, not by
      `OrderedSet.getAt`; non-identity operators keep the mapping's `x`, `y`, `z` (exact). Interleaved A/B, three
      rounds alike: the density stage 16.8 / 16.5 / 17.2 → 15.9 / 15.5 / 16.2 ms for 100k atoms. A CPU profile
      overstated it (24 → 14 ms a run): its overhead inflates calls per element.
    - left: the three walks as one (the kernel's xyzr written while collecting), about 1–2 ms at 100k atoms; the size
      theme per element.
  - On the iGPU the density kernel's work per cell is what to look at.
- **Merged meshes** (reported: GTX 980, ANGLE on GL, 4V6X Gaussian surface 452 fps unmerged, 135 merged; benchmark
  `merging/{4v6x,3j3q,assembly}-gaussian-*`, CPU meshes on every backend, opaque pass in ms):

  | Case                             | NVIDIA WebGL2 | NVIDIA WebGPU | iGPU WebGL2 | iGPU WebGPU |
  | -------------------------------- | ------------- | ------------- | ----------- | ----------- |
  | 4V6X merged, before              | 1.10          | 1.28          | 59.4        | 43.2        |
  | 4V6X merged, after               | 0.32          | 0.42          | 24.2        | 16.9        |
  | 4V6X unmerged                    | 0.42          | 0.44          | 19.7        | 17.7        |
  | 3j3q (1356 chains) merged, after | 0.38          | 0.31          | 4.93        | 4.18        |
  | 3j3q unmerged                    | 4.41          | 0.61          | 4.94        | 4.38        |

  - the cost was the vertex work, the same at 320×180 as at 2560×1440: a merged item drew its solid interior caps every
    frame (`renderSolidInteriorMark` and `renderSolidInteriorFill`, all its vertices twice more). Its sphere spans the
    scene's, whose front the camera fits near to, and the near test took near nudged 1e-4 inward; for 3j3q the merged
    sphere (a union of the members') reached past the camera's by 8e-4. The test now takes the margin the other way and
    a merged item's member spheres (`Renderable.getSegmentSpheres`): the caps drawn where a member is cut, as unmerged.
  - WebGL2 without base instances (every Chrome: the extension is a draft) reads merged items' instance data from a
    texture (`dInstanceData`). Each use of `aInstance`, `aSegment`, `aTransform` (up to ten chunks) fetched again, with
    an integer division and modulo per texel: 23% on the iGPU for 4V6X, 17% for 1000 cartoon units. Now read once per
    vertex into globals (`readInstanceData`, first in `main`), rows of 256 instances addressed by shift and mask; iGPU
    opaque pass, before → after (the instanced attributes re-pointed per draw for comparison, 2–5 µs CPU a draw):

    | Case                | Before | After | Attributes | Unmerged |
    | ------------------- | ------ | ----- | ---------- | -------- |
    | 4V6X merged         | 24.2   | 19.5  | 20.5       | 19.7     |
    | 1000 cartoon units  | 9.80   | 7.81  | 8.19       |          |
    | 3j3q merged         | 4.93   | 4.26  | 4.36       | 4.94     |
    | 7qpd assembly       | 15.0   | 12.4  |            | 15.5     |
    | 1000 spacefill      | 2.76   | 2.70  | 2.77       |          |

  - open: a general analysis of vertex inputs and varyings (their count and size per geometry type and variant, e.g. a
    transform of 3 vec4s, packed group and instance indices), and a plan for reducing them.

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
3. The GLSL compile walk: native `gl` and a visual-suite check per profile, on demand. Done (the WebGL runs' length
   planned, see WebGL runs).
4. Fixes for what 2 and 3 find. Done.
5. The parity script, and the review of every program pair. Not done: the visual suite compares every WebGPU
   profile with WebGL2 per case, enough parity at the level that matters.
6. The hazard review in both languages. Done (findings above; the cap depth bias and the latent ones open).
7. GPU time per pass on both GPUs and backends, ranked. Done (GPU time per pass, above).
8. The top passes' hot variants analysed: shader dumps, Radeon GPU Analyzer (AMD), Nsight Graphics captures (NVIDIA).
   Done for the denoiser, depth of field, tracing and sphere impostors (step 8 above).
9. Candidates measured and adopted or recorded.

## Decisions

- **No external validator as a dev dependency** (glslang, naga): the backends' own compilers (ANGLE, Tint in Dawn)
  decide, and they compile what browsers run.
- **Mobile GPU behaviour is out of scope** (Mali, Adreno, Apple; `mediump` precision): no device here to check it on.
- **Reachable variants** come from three sources: the default presets, the MolViewSpec examples and mesoscale scenes.
  More sources are added where a finding shows a gap.
