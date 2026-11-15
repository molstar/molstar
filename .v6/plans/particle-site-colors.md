# Plan: A second coloring of particle target sites, blended in near the camera

## Status

- **Step 1 done**: the elements kept (`ParticleSites.elements`, `ParticleSitesQuantized.typeElements`, the pack's
  `siteElements`, version 5); the benchmark packs packed again (version 5, with elements).
- **Step 2 done**: the data, no shader reads it yet. `SiteColorSchema` on spheres (`dSiteColor` off, `siteType` or
  `group`; `tSiteColor`, alpha 0 a carbon kept the particle's; `uSiteColorRange`); the param `siteColor` of the
  structure targets' spacefill (`SiteColorParams`, on: element, carbons the particle's, 4-12 px), set for targets of
  sites (`setSitesTargetSiteColor`); a merged object's tables packed, a member's first entry its segment data's
  element 11 (`siteColorTableBases`, beside the size table's). Members with and without the second color do not merge
  (a define).
- **Step 3 done** (WebGPU checked, WebGL by CI's native suite: `gl` not installed here): the blend in the vertex stage
  (`site-color.wgsl.ts`, `assign-site-color.glsl.ts`) and in the sphere raster's resolve (Design 3 revised); a uniform
  color the object's (`object`) where the second color is on (the blend mixes into a texture's color); `tSiteColor` a
  WebGPU data layer (`DataLayers`, the sixteenth). Test: `plugin/headless/_test/particle-site-colors.test.ts` (a
  pack's oxygens red near, the particles' green far; impostors and raster).
- **Step 4 done**: no UI of its own: Prism's inspector shows a representation's params as a form
  (`particles.update`), `structure.siteColor` among them, and the assistant sets it by `particles.representation`'s
  params; no v6 docs page covers the target representation's params. The changelog's entry
  (`.v6/changelogs/particle-site-colors.md`).
- **Step 5 measured** (the yeast cell in Prism, WebGPU, Chrome headless, 1280×720): a forced frame 18.0 ms from 600 Å
  and from 80 Å, the second color on or off alike (the frames at the page's rate: no difference above it); the JS heap
  1.52 GB either way; its GPU data a table of 7 entries a member, 158 KB for the 5,639 targets. The GPU timer gave no
  results in headless Chrome: the shaders' cost not measured apart (a few operations a sphere in the vertex stage, a
  table read where blended).
- **Structure targets done** (Design 6): their spheres' groups are their elements, so the second color by group, a
  texel an element in the serial order from its `type_symbol` (`setStructureTargetSiteColor`, `site-color.ts`; a coarse
  unit's element none, the particle's color), `dSiteColor` 'group' in the precompile's signatures
  (`getStructureSignatureDefines`). The param and the shared parts in `site-color.ts`. Test: a RELION list's 1crn
  targets, red near and none far (`particle-site-colors.test.ts`).
- Decided: the blend by size on screen; the second theme element symbol with its carbons the particle's color; old
  packs (versions 2-4) without the option; on by default for lists with elements; the range 4 to 12 px.

## Why

- **Far and near want different colors**: a cell seen whole reads by its particles' theme (hierarchy, entity,
  compartment); zoomed in to a protein, by its atoms (element). Today a target representation has one color theme, its
  colors per particle (`createTargetColors`: instance granularity, or a constant the object's), the sites of a
  particle all one color.
- **The sites carry no element**: a pack's site type is a distinct radius (`pack.ts`: `siteRadii`, the type its
  index), not an element (S and P share 1.8 Å); the mmCIF reader reads `type_symbol` only for the radius
  (`formats/particles/mmcif.ts`, `VdwRadiusOfNumber`) and drops it.
- **The sizes already go by type**: quantized sites draw their radii from a table by type (`siteType` sizes,
  `createSiteTypeSize`), a merged object's members' tables packed and each segment's first entry in its segment data
  (`merged.ts` `siteType`, `CulledData` `SegmentSizeBase`). A color table by type is the same shape: no per-site color
  data for millions of sites.

## Design

1. **Elements of the sites** (model): `ParticleSitesQuantized.typeElements?: Uint8Array` (a type's atomic number) and
   `ParticleSites.elements?: Uint8Array` (a site's, for sites not quantized). The mmCIF reader keeps the atomic number
   it already looks up. A pack's types become the distinct (radius, element) pairs, the document's `siteElements` beside
   `siteRadii` (version 5; versions 2-4 read without: no elements). No cost a site in a pack (a type's table entry), a
   byte a site otherwise.
2. **The second color**: a target representation param, `siteColor`: off, or a theme of the sites (element symbol
   first) and the blend's range. Evaluated per type into a table (`tSiteColor`, a texel a type; merged: the members'
   tables packed as the size tables, the segment's first entry beside its size base in the segment data, vec4 2's `y`).
   Sites not quantized (WebGL1, a structure-less list's plain sites) by group (a texel a site, the group the site
   index).
3. **The blend**: per sphere in the vertex stage (flat, no fragment cost): `t = smoothstep(start, end, s)`, `s` the
   sphere's radius on screen in pixels (the projection's focal scale over the view depth, as the levels of detail
   (`LodView`) and the raster split (`MergedRasterParams.scale`) measure it, CSS pixels); `color = mix(primary,
   second, t * alpha)`. By size on screen, not distance: the same at any field of view, resolution and atom size. The
   sphere raster draws spheres up to 24 px (`sphereRasterRadius`), the whole blend's range: its resolve blends too (it
   builds the vertex stage's outputs, `resolveSphere`). The primary a texture's color (`dColorType` texture): a
   target's constant color the object's (`object`, a texel; merged, a member's by its segment), a uniform one made so
   where the second color is on.
4. **What stays**: picking and marking (ids, not colors); overpaint, transparency and substance layers over the
   blended color; illumination and tracing read the blended color (the vertex stage's).
5. **Backends**: WebGPU `spheres.wgsl` (culled: the segment's table base by `rasterSegment...`; not culled: the
   object's), WebGL `spheres` GLSL (`assign-color` chunk), both by a define (`dSiteColor`), set only where the option is
   on: no cost for objects without it.
6. **Structure targets** (`kind: 'structure'`): the same param, their elements by group
   (`createStructureElementSphereImpostor`, a group an element), done after the steps.

## Steps

1. **Elements kept**: `ParticleSites.elements`, `ParticleSitesQuantized.typeElements`; the mmCIF reader; the pack's
   `siteElements`, version 5, the older versions read without. Tests: the round trip, a version 4 pack read, S and P of
   one radius two types.
2. **The second color's data**: the param, the table by type (and by group for plain sites), the merged packing and
   segment base, the element symbol colors by atomic number. Tests: the table, a merged object's bases.
3. **The blend in the shaders**: WebGPU and WebGL spheres, the raster's resolve. Tests: the
   kernels and programs preprocess with the define; a headless render (medium), a sphere's color near and far.
4. **Prism**: the option in the target representation's controls; the docs.
5. **Measured**: the yeast cell near and far (frame time, memory: a table a type, nothing a site).

## Decisions

- **The measure**: a sphere's radius on screen in pixels (Design 3), not distance.
- **The second theme**: element symbol, its carbons the particle's color (`carbonColor`): most atoms keep the hue of
  the particle across the blend. Other site themes later where asked.
- **Old packs** (versions 2-4, no elements): the option unavailable (not shown), no guess from radii; packed again to
  have it.
- **Default**: on for lists whose sites have elements, off for the others.

## Open questions

- None. The range stays 4 to 12 px (decided by eye on the yeast cell).
