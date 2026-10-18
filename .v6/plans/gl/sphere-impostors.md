# Plan: Sphere impostor performance

## Status

- **Current**: analysed (measurements below), nothing changed. The WGSL `discard` idea did not hold up (Tint emits a
  native `discard`); the WebGPU gap is in the shading code of the generated HLSL, half of it strict math (option 1).
  Conservative depth waits for measurements on other GPUs and a Chrome that compiles WGSL `fragment_depth` on D3D12.
- **Option 3 done: WebGPU has early-Z**, the earlier finding was wrong. A minimal WebGPU page (Chrome 155, D3D12,
  2560x1440, an expensive fragment shader drawn behind a near occluder, timestamp queries): hidden fragments cost
  nothing with every depth format (`depth24plus`, `depth32float`, with stencil), with and without depth writes, with
  `discard`, with stencil writes; only a `frag_depth` output disables it (hidden as costly as visible). In Mol\* (the
  globule as spacefill drawn once and twice, the second copy at equal depth, opaque pass, WebGL2 / WebGPU): approximate
  spheres (no depth output) 9.7 -> 19.4 / 11.0 -> 11.7 ms, so WebGPU rejects the hidden copy and WebGL (ANGLE D3D11,
  `discard`) does not; impostors (`frag_depth`) 13.5 -> 21.8 / 14.3 -> 28.5 ms, early-Z off on both. The earlier
  "WebGPU without the depth write gets no early-Z" kept the `frag_depth` member of the output struct (a depth output
  even when not assigned). WGSL conservative depth (`requires fragment_depth;`, `@builtin(frag_depth, greater)`) still
  fails pipeline creation on D3D12 in Chrome 155 (`DXC compile failed with error: E_FAIL`), shader validation passes:
  to report. Consequence for option 5: on WebGPU a shading pass without depth output keeps early-Z, also with
  `discard`, so a depth pre-pass can skip the shading of hidden impostors.
- **Option 1 measured for the report**: `pnpm benchmark --strict-math off` (Chrome with
  `--enable-webgpu-developer-features`, `strictMath: false` injected into every shader module) and the case
  `render/globule-300k-spacefill-large` (2560x1440). Opaque pass 14.3–14.4 ms strict, 12.4 ms relaxed (-13.5%, three
  alternating rounds of 9 frames each; frame 14.7–14.8 against 12.8 ms; WebGL2 10.4 ms). A minimal repro without Mol\*
  (one WGSL PBR fragment shader, full-screen draws at 2560x1440, timestamp queries) gives -19% in Chrome 155 and 156, 49
  of 14.7 M channel values off by 1, and 0% without the flag (the member is ignored). Repro and the text for the
  Chromium issue prepared outside the repository. Decided: not adopted, `strictMath` is a non-standard developer option;
  revisit if WebGPU gains a standard way to relax math (the report asks for one).

## Where things stand

- **Vertex stage** (`spheres.vert`, `spheres.wgsl`): six vertices per sphere (two triangles, no elements), each pulls
  the sphere from `tPositionGroup` and computes the projected bounding rectangle (`sphereProjection`: two square roots,
  four divisions) and `uInvProjection * gl_Position`; the quad lies at the sphere's center depth.
- **Fragment stage**: a ray-sphere intersection (`SphereImpostor`), `discard` where the ray misses (the quad corners,
  21% of the quad area for a circle), then the full material and lighting; the depth is written (`gl_FragDepth`,
  `frag_depth`). The shaders contain further runtime discards (opaque/transparent masks, LOD fade, clipping).
- **Measured** (`render/globule-300k-spacefill`, 300k spheres, RTX 4070 Ti, opaque pass, WebGL2 / WebGPU):

  | case                                                    | 640x360     | 1280x720  | 2560x1440   |
  | ------------------------------------------------------- | ----------- | --------- | ----------- |
  | impostors                                               | 1.3 / 1.4   | 4.2 / 4.5 | 10.8 / 14.8 |
  | impostors, `ignoreLight`                                |             | 3.1 / 3.1 | 10.4 / 10.6 |
  | approximate spheres (no depth write)                    |             | 3.0 / 3.3 | 9.7 / 10.9  |
  | approximate, trivial fragment shader                    |             | 1.45      | 4.3 / 4.4   |
  | approximate, trivial fragment shader, 128x72 (vertices) | 0.16 / 0.26 |           |             |

  Without lighting (`ignoreLight`) both backends take the same time: the WebGPU gap is the lighting code; at 2560x1440
  WebGL2 hides its lighting behind the rasterization and depth writes (10.4 ms with or without), at 1280x720 it costs
  1.1 ms there.

  Fragment-bound: about 0.2-0.4 ms vertices, 1.1 ms rasterizing the quads, 1.5 ms material and lighting, 1.1 ms the
  ray-sphere intersection and the depth write (1280x720, WebGL2).

- **Early depth testing does not happen**: with an expensive loop added to the shading, drawing the globule twice (the
  second copy fails the depth test everywhere) costs twice as much with impostors and still 70% more with approximate
  spheres without a sphere `discard`. The depth write disables early-Z and coarse depth culling; the remaining runtime
  discards together with depth writes force late depth tests (NVIDIA, D3D11/D3D12).
- **Conservative depth** (`EXT_conservative_depth`, `layout(depth_less) out float gl_FragDepth`, the quad moved to the
  sphere's front-most depth so that it bounds the fragment depths; depth is reversed): available in Chrome (ANGLE D3D11
  translates it, `GL_USES_FRAG_DEPTH_LESS`), but no gain here (4.12 ms; the expensive-loop test unchanged), also without
  the sphere `discard` (misses written at depth 0, which fails the test). WGSL `fragment_depth`
  (`requires fragment_depth;`, `@builtin(frag_depth, less)`, Chrome 156): every pipeline using it fails to compile on
  Windows D3D12 in Chrome Beta 156 (DXC and FXC `E_FAIL`, also a one-line shader): to report, to retry.
- **WebGPU's extra cost is in the material and lighting code** (2560x1440): the rasterization floor is equal (4.3 / 4.4
  ms); the intersection and depth write without material and lighting cost the same on both (10.4 / 10.5 ms); on WebGL2
  the shading is hidden behind them (10.4 ms with it), on WebGPU it adds 3.8 ms (14.3 ms). The WGSL lighting mirrors the
  GLSL (the same BRDF, the light loop bounded by a pipeline constant), so the difference is in the generated code: Dawn
  compiles with `D3DCOMPILE_IEEE_STRICTNESS` (no fast-math reassociation). With `strictMath: false` in
  `GPUShaderModuleDescriptor` (a developer option behind `--enable-webgpu-developer-features`,
  https://issues.chromium.org/issues/42241455) the pass takes 12.3 instead of 14.3-14.4 ms (-14%, alternating runs):
  about half of the gap; the other 1.9 ms remain. Not register pressure: the color variant's pixel shader compiled with
  FXC (O3) from Tint's HLSL and from ANGLE's HLSL needs 12 and 13 temporary registers, 868 and 818 instruction slots
  with IEEE strictness, 725 and 674 without (strictness adds about 20% on both). Tint's extra instructions are control
  flow (+10 `if`, +14 `endif`, +6 `else`; ANGLE uses `movc` and a `switch`) and moves (+19 `mov`); the arithmetic is the
  same (+7 each of `mad`, `add`, `mul`). WebGPU compiles with DXC (shader model 6, O3, strict; 335 `dx.op` calls), ANGLE
  with FXC: the rest of the gap is between these compilers and the D3D12 and D3D11 drivers. Not `discard`: Tint emits a
  native `discard` (dumped HLSL), returning after each `discard` of the entry points gained about 1% (14.2-14.3 ->
  14.0-14.2 ms, alternating runs); an earlier 4.47 -> 4.18 ms was run-to-run noise.
- **Without the depth write** WebGL2 gets early-Z (10.4 -> 5.3 ms at 2560x1440, the quad at the sphere center, so not a
  correct image); WebGPU seemed not to (14.3 -> 14.4 ms), but that measurement kept the `frag_depth` output (see Status,
  option 3): WebGPU has early-Z without a depth output, also with `discard`.

## Options

1. **`strictMath: false`** in the descriptors of the WebGPU shader modules (`createWebGPUProgram`, compute and the other
   modules): a developer option behind `--enable-webgpu-developer-features`
   (https://issues.chromium.org/issues/42241455), unknown descriptor members are ignored elsewhere, so it is harmless
   without the flag; -14% for the globule's opaque pass with the flag. Before adopting: the visual suite with the flag
   (fast-math reassociation may change impostor depths and lighting slightly), a `--webgpu-developer-features` switch
   for the benchmark and visual scripts. Report the measurement on the Chromium issue as a case for exposing the option
   without the flag.
2. **The rest of the WGSL shading cost** (1.9 ms at 2560x1440): not register pressure (above), not the branches: the
   short-circuit conditions of the fragment path (`||`/`&&` of the impostor miss, the transparency checks, the LOD fade,
   the bump flag, the order-independent transparency test) as non-short-circuit `|`/`&` measured the same (14.3 ms both,
   alternating runs); Dawn's FXC path is slower (16.3-16.4 ms, compiled at O0), so no lead there. What remains is DXC's
   code and the D3D12 driver, beyond Mol*'s reach except for `strictMath` (option 1).
3. **Early-Z on WebGPU**: done, it works (see Status); impostors lose it through their depth output, as on WebGL.
4. **Conservative depth** on both backends (front-most quad depth, `depth_less` / `frag_depth, less`): correct and
   cheap, no gain on NVIDIA D3D11 here; to measure on AMD, Intel and Apple (ANGLE Metal, `[[depth(less)]]`) before
   adopting, WGSL once Chrome compiles it.
5. **Early-Z through a depth pre-pass** for impostors: a depth-only pass (intersection and depth, no shading), then the
   shading pass without depth output at the front-most quad depth with `greater-equal` (coarse culling rejects hidden
   spheres) and an in-shader test against the pre-pass depth (WebGPU: read-only depth attachment; WebGL: a copy). Twice
   the vertex and rasterization work (about 1.5 ms here), saves shading of hidden fragments; pays off only with high
   overdraw, measured per scene. Planned and dropped in [sphere-depth-prepass.md](sphere-depth-prepass.md): the
   depth-only pass costs about as much as the impostors without lighting (10.4 of 14.2 ms for the globule), more than
   the lighting it saves.
6. **Fewer fragments**: an octagon instead of the quad (5% instead of 21% outside the circle), at the cost of more
   vertices; small spheres (a few pixels) are dominated by 2x2 quad rasterization anyway.
7. **Fewer vertex invocations**: four vertices with elements (the post-transform cache shares them) instead of six: at
   most a third of the 0.2-0.4 ms.

## Steps

1. Option 1: measured (`--strict-math off`), not adopted (a non-standard developer option); the measurement on the
   Chromium issue (prepared, to post).
2. Option 3: done (WebGPU has early-Z, the depth output disables it on both backends); option 2 explored without a
   lead.
3. Option 4 as a prop or automatic per GPU, once measured on other GPUs.
4. Options 5-7 only with measurements on scenes with high overdraw (assemblies, cellular scenes).
