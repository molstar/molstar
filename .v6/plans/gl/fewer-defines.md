# Plan: Fewer shader defines

## Status

- **Current**: paused, good enough for now: the color type, the instance granularity, the size type and (on WebGL) the
  clip object count converted (step 2, below). The shading switches (4) measured and not converted (below); open: the
  layers' presence (5) last if at all.
- Step 4 measured for cel shading and flip-sided, not adopted: as uniforms (`uCelShaded`, `uFlipSided`, one program
  each way, WebGL ignoring the defines) they cost every frame, cel shading off on WebGPU 2.1% of the fragment-bound
  globule (`render/globule-300k-spacefill-large`, 14.68 to 14.98 ms), both on WebGL2 0.3-0.9% and on WebGL1 1.0-1.7%
  (the globule as impostors and meshes, three interleaved rounds); with cel shading on within 1%. These switches rarely
  vary (one case each in the visual suite, none in the stories), so the programs saved do not pay for it. The other
  switches (x-ray, ignore-light, transparent backfaces, solid interior) gate more code in the same fragment paths, so
  not tried; `cel/*` benchmark cases (on demand) measure cel shading.
- Step 2 done for the clip object count on WebGL: the presence is the boolean define `dClipObjects`, the count the
  material uniform `uClipObjectCount` (`loopBound`: a define on WebGL1), WebGL ignores `dClipObjectCount`
  (`ignoreDefine`), so WebGL2 has one program for all counts; the GLSL clip loops are no longer unrolled. WebGPU keeps
  the count a pipeline constant (`oClipObjectCount` from the number define): a uniform bound cost it 15% with 1 and
  23% with 3 clip objects (`clip/globule-300k-spacefill-large-*`, on demand: inverted spheres that cut nothing, every
  fragment testing), a constant bound with a `break` still 6-10%, a guarded fixed loop 65%. Measured, three interleaved
  rounds (2560x1440): WebGL2 3 clip objects 33.2 to 29.0 ms, WebGL1 33.3 to 24.9 ms, 1 clip object and none unchanged;
  WebGPU unchanged. Checked: the visual suites (WebGL2, WebGL1, WebGPU, WebGL2 without uniform blocks), the headless
  WebGPU tests.
- Step 2 done for the size type, as the color type: `dSizeType` has `uniform` and `texture`, the object uniform
  `uSizeIndexType` selects the index (`SizeTypeIndex`: instance, group, groupInstance, vertex, vertexInstance; a
  `switch` in WGSL, an `if` chain in GLSL), the type itself is the value `sizeType` (`setSizeType`; read by merged
  objects, `getMaxSize`, the mesh exporter; the signature probe maps with `getSizeTypeDefine`). Frame cost, base and
  conversion alternated three times over the render cases (opaque pass): the large scenes unchanged on both backends
  (globule 300k as spacefill, large, as meshes; the 60-copy assemblies; 1000 units), the instanced spacefill assembly
  with SSAO +0.5% on WebGPU (2.44 → 2.45–2.46 ms), scenes under 0.5 ms within their spread. Checked: the visual suites
  on all profiles identical to the images before (`--save`/`--compare`), the headless WebGPU tests, a unit test of the
  size type values.
- Step 2 done for the instance granularity: the object uniform `uInstanceGranularity` (was the value
  `instanceGranularity`) selects the index of the marker and layer data, per instance or per group instance, at runtime
  (GLSL the macro `granularityIndex`, WGSL `select`; images and direct volumes in their fragment shaders). `dMarkerType`,
  `dClippingType` and `dWiggleType` are removed; the four layer type defines have `texture` (instance and
  groupInstance), `vertexInstance` and `volumeInstance` (`getLayerTypeDefine`), the types themselves are values
  (`markerType`, `overpaintType`, …, set with `setOverpaintType` etc.). The marker and clipping varyings are `flat` where
  groups vary, whatever the granularity (correct for instance data too). Signatures no longer resolve 'auto' instance
  granularity. Frame cost, measured before as a prototype: the layer cases (60-copy assembly and the 300k globule as
  spacefill, with overpaint and substance on everything), base and prototype alternated three times, opaque pass within
  noise on both backends (assembly WebGL2 0.94–1.00 → 0.95–0.97 ms, WebGPU 1.06–1.10 → 1.06–1.07 ms; globule ±1%; a first
  single run had suggested −9% and +18%). Fixed on the way: images and direct volumes read markers per group instance
  whatever the granularity, and direct volumes ignored overpaint with instance granularity. Checked: the visual suites on
  all profiles identical to the images before (`--save`/`--compare`, WebGPU too; only the cases that vary between runs,
  dpoit solid interiors, differ), the headless WebGPU tests (precompile derivation and snapshots).
- Step 2 done for the grid color types too: volume and volumeInstance share `grid`, `uColorIndexType` selects the
  position the grid is sampled at (the local position, or the model position); the grid sampling helpers of WGSL are
  emitted for `dColorType_grid`. The grid types stay apart from `texture`: they read `tColorGrid`, bound only for them
  (`usedIf`), and an always-bound second color texture would cost a vertex sampler in every program, scarce on the
  minimal WebGL profiles. Checked as before: the visual suites (WebGPU, WebGPU core), the WebGL profiles' images
  against those before the color type conversion unchanged (the same `blobs` difference), the headless tests.
- Step 2 done for the color type: `dColorType` has the values `uniform`, `texture`, `volume`, `volumeInstance` and
  `direct` (then; the grid types merged after, above); the five types read from `tColor` by index share `texture`, and
  the object uniform `uColorIndexType` selects the index at runtime (`ColorTypeIndex`: instance, group, groupInstance,
  vertex, vertexInstance; a `switch` in WGSL, an `if` chain in GLSL; direct volumes the same in their fragment shaders).
  The color type itself is the value `colorType`; `setColorType` keeps it, the define and the index in step (color data,
  mesh and texture mesh smoothing), and what read the type from the define reads it from `colorType` (merged objects,
  the units visual's updates, the signature probe, which derives define values with `getColorTypeDefine`). Dual colors
  (cylinders), only for the group types, are chosen at runtime too. Checked: the visual suites on all profiles, and the
  WebGL profiles' images against those before the change (`--save`/`--compare`): identical except `geometry/direct-
  volume/blobs` on `webgl1-minimal`, 4 pixels of 65536 (deterministic, likely the ray march's restructured branch); the
  headless precompile tests; a unit test of the color type values. Effect: the visual suite's drawn signatures 54 → 53,
  `kinase-story`'s drawn 10 → 9 (a scene rarely mixes texture types), its precompile list 24 → 18 signatures (the over-
  approximated color types: the multilayer theme's granularities), on WebGPU 120 → 69 programs and 151 → 89 pipelines,
  on WebGL 78 → 48 programs (1.9 → 1.3 s), playing still compiling nothing.

## Goal

Fewer programs per object: defines that select small branches become uniforms with uniform (coherent) branches, or
are merged into fewer variants. Every define less halves (or more) the program space, for drawing (programs compiled
when objects are created or change) as well as for precompiling (`Canvas3D.precompile`). The cost per frame stays
measured: impostors are fragment-bound ([sphere-impostors.md](sphere-impostors.md)) and more registers can lower
occupancy.

## Where things stand

A program is determined by its render object kind, its define values and the render variant (`ObjectSignature`,
`gl/precompile.ts`). On `kinase-story` the precompile list derives 71 signatures; together with the variants (color,
depth, pick, marking, emissive, tracing) and, on WebGPU, render states and targets, that is 211 programs and 664
pipelines (pipeline-readiness.md, 6f and 6g).

## How often they vary (step 1)

The signatures the visual suite's cases draw (WebGPU profile, 127 cases with objects; `Canvas3D.getPrecompileList` per
case): 54 distinct signatures, i.e. programs per variant (mesh 20, texture mesh 10, points 7, spheres 7, cylinders,
lines, direct volume, image, text). Ignoring one define at a time saves few, as each case varies one feature:
`dColorType` 5, `dClipObjectCount` 5, `dPointStyle` 3, `dXrayShaded` 2, and 1 each for `dApproximate`, `dCelShaded`,
`dClipVariant`, `dFlatShaded`, `dIgnoreLight`, `dPointSizeAttenuation`, `dSizeType`, `dSolidInterior`,
`dTransparentBackfaces`. Ignoring groups: the layers 54 → 51, the shading switches 54 → 46, the color and size types
54 → 48, all candidates 54 → 29.

A story varies differently: the 10 signatures `kinase-story` draws differ in `dColorType` (uniform, volume, group,
groupInstance), `dSizeType`, `dSolidInterior` and `dVaryingGroup`; no layer at all. Its precompile list
over-approximates the color types (the MVS multilayer theme's granularities). So the color type, the largest data
layout define, is also the one that multiplies programs most in practice, and the order below changes: the color type
(and size type) first, if its frame cost allows, then the shading switches; the layers last.

`dClipObjectCount`, a count (the clip objects' array size), varies too: a fixed maximum with a uniform count instead.

## The color type's frame cost (step 1)

Six color types read `tColor` or the uniform color and differ only in the index (uniform, instance, group,
groupInstance, vertex, vertexInstance); the grid types (volume, volumeInstance) read `tColorGrid`, whose dimension
changes the binding layout on WebGPU, so they stay a define. A prototype, GLSL and WGSL: for programs reading `tColor`,
the vertex shader selects among the six at runtime (WGSL `switch`, GLSL an `if` chain for GLSL ES 1.00) on a selector
the compiler can't fold, and the fragment shader chooses between the uniform color and the varying at runtime; the
programs stay split by the define, only the code is what the merged program would run. Visual profiles (WebGPU, WebGL2,
WebGL1, WebGL1 minimal) unchanged. GPU time of the opaque pass in the render benchmarks, before and after:

- large GPU-bound scenes within about ±1% on both backends (globule 300k as spacefill, large and as meshes; the
  60-copy cartoon assembly);
- the instanced impostor assemblies consistently +1.6–3.5% on both (60 copies of 7qpd as spacefill, rotating, with
  SSAO and multi-scale SSAO), a small real cost, likely registers of the impostors' vertex shader;
- scenes under 0.2 ms: noise either way (a Gaussian surface case first +13% on WebGPU, +3% repeated).

So merging all six costs up to about 3% in instanced impostor scenes, nothing measurable elsewhere.

Variant measured next: the uniform color kept as its own define (its objects need no color varying, the fragment code
stays as it is) and only the five texture index types merged (instance, group, groupInstance, vertex,
vertexInstance), against a baseline measured in the same session: the instanced spacefill assemblies +0.3–1.4%,
globule ±0.4%, cartoon assembly ±0.3% on both backends, i.e. within noise; scenes under 0.2 ms noise (`7qpd-spacefill`
first +50% on both, then 0.08 ms with and without, three runs each).

A third variant: the uniform color as a lookup too (a one-texel `tColor` at index 0), so that all six types share one
structure, and without a branch: the index `instance * a + group * b + vertex * c` with the factors per type (uniform
0, 0, 0; instance 1, 0, 0; group 0, 1, 0; groupInstance groupCount, 1, 0; vertex 0, 0, 1; vertexInstance vertexCount,
0, 1) in an object uniform; merged objects keep their segment-offset indices. Measured as before (the formula with
opaque factors, against a baseline of the same session; the render cases draw with group color types, `element-symbol`
and `chain-id`): within ±0.7% in the large scenes on both backends; the small cases repeated three times each way equal
(Gaussian surface 0.33–0.34 ms, 7qpd spacefill 0.08 ms). What uniform-colored objects would pay, a texel read and a
varying instead of `uColor`, measured as the difference of the group-colored baseline to one with uniform colors: the
globule's sphere impostors +2.3% (WebGPU, WebGL ±0), its meshes +4.4% (WebGPU) and +3.7% (WebGL), the assemblies about
0.

So: the texture types as one lookup cost nothing measurable (the index formula and the switch alike); the uniform color
as a lookup costs uniform-colored objects up to about 4% in large scenes, for the one structure and the fewest programs
(`dColorType`: a lookup or the two grid types). Decided: the texture types merged, the uniform color kept (no cost).

## Candidates

In order:

1. the color types: done (`texture`, `grid`, above);
2. **instance granularity**: done (above); `dMarkerType` and the layers' types (`dOverpaintType`, `dTransparencyType`,
   `dEmissiveType`, `dSubstanceType`, `dClippingType`, `dWiggleType`) where they are `instance` or `groupInstance`,
   which `instanceGranularity` decides for all of them together: one object uniform selects the instance or the
   groupInstance index at runtime (the varyings `flat` where groups vary, correct for instance values too), so that
   changing the granularity, or 'auto' resolving it by size, creates no programs (and the precompile lists need not
   resolve 'auto', pipeline-readiness.md 6d). The layers' grid types (`volumeInstance`, smoothed layers; `vertexInstance`
   without a context) stay their defines: they read the layers' grid textures, bound only for them;
3. the size type (`dSizeType`): done (above);
4. shading switches: `dXrayShaded`, `dCelShaded`, `dIgnoreLight`, `dTransparentBackfaces`, `dSolidInterior`,
   `dFlipSided`, `dDualColor`, `dUsePalette`;
5. the layers' presence (`dOverpaint`, `dTransparency`, `dEmissive`, `dSubstance`, `dClipping`, `dWiggle`): they gate
   texture bindings (`usedIf`) and blocks of code, as uniforms the layers' textures would always be bound (samplers),
   and transparency and emission also select passes; last if at all;
6. counts: `dClipObjectCount` (a fixed maximum, the count as a uniform).

## Steps

1. **Measure**: per candidate, how often it varies in the stories and the visual suite's scenes (programs it adds), and
   the frame time of the render benchmarks with it as a uniform branch, GLSL and WGSL (both backends, the impostors
   first).
2. **Convert** the candidates whose branch costs nothing measurable, one group at a time (in the order above): the
   define becomes a uniform (material or object uniforms), both shader languages; signatures and schemas follow.
3. **Check**: visual suite unchanged on all profiles, the benchmarks within noise, the precompile lists and programs per
   story counted again.
