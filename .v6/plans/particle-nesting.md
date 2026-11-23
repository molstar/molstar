# Plan: Nested particle lists (instances of instances)

## Status

- **Done for grids** (steps 1 to 9; 5 on WebGPU and WebGL2, WebGL1 flattened; the sparse overrides within a budget):
  a particle may place a whole particle list (a template), stored once and expanded on the GPU only, each expanded
  instance addressable in each copy. Open: symmetric assemblies (a capsid, the NPC: no producer of nested packs from
  an assembly's operators, step 9 not measured for them). Parked: colors by level, deeper levels, the reader's form,
  templates and tiles ([Decisions open](#decisions-open)). A memory saving: a repeated part's placements stored once,
  not a copy each. Builds on [particle-slots.md](particle-slots.md) (step 6f, the page store and its cull that
  composes a pose a particle) and [gl/instance-poses.md](gl/instance-poses.md) (poses, not matrices).

## Question

A particle list stores every placement: about 48 B a particle on the CPU (target, entity, compartment, position,
rotation, radius) and its stored values on the GPU, whatever repeats. Mesoscale models repeat whole parts: the isg cell
copied on a grid (`isg_mature_grid-AxBxC`: 460k placements a cell, 29.4 M at the largest, 9.1 GB;
[mesoscale-scale.md](mesoscale-scale.md); made by a pack script gridding the input, every copy written flat, as
`scripts/particle-grid-pack.mjs` does the yeast's), symmetric assemblies (a capsid of 60 asymmetric units, the NPC's
eight spokes), mmCIF assemblies whose operators are a product (read flat today, `composeOperatorCombination`). Storing
such a part once and placing it by a pose a copy takes its memory from copies x placements to copies + placements.

What must hold: every expanded instance is a particle of its own to the user (picked, labelled, highlighted, selected,
overpainted alone, in one copy and not the others), and nothing a copy is expanded on the CPU (else the saving is lost
where it is needed most: the JS heap, [mesoscale-scale.md](mesoscale-scale.md)'s 4 GB a page).

## Design

- **A target of particles**: `ParticleTarget` kind `particles` (`{ particles: ParticleList, center? }`): a particle of
  that target places the whole list, its pose composed with each of the list's particles' (`parent ∘ child`). The
  template's own targets are what is drawn (its leaves). Two levels first (a top list and its templates, a template's
  targets leaves only); a template of templates later, where measurements ask for it (the same stage repeated, the
  identity's `local` decoded a level more), so the model keeps a depth (`MaxNestingDepth`, 1 now) and refuses cycles
  where a list is made. A template's particles are not
  in the list that places it: the list's `count` is its own particles, `expandedCount` all instances (by its targets').
- **An instance's identity**: its path, an index a level (`i0, i1, …, ik`). Kept as two words: the top particle `i0`
  and its index in the template's expansion (`local`: the template's particles' expansion counts as prefix sums, a
  level read by a binary search a level). A template's expansion within 2^32 (checked); the top list's count as now.
  The CPU decodes a path where it needs one (a label, a loci); the GPU carries the two words.
- **Loci by path sets**: a particle loci's elements are product sets, a path prefix and an index set at the next
  level (`OrderedSet`), a wildcard a level: one instance (`[i0, i1]`), a whole copy (`[i0, *]`), a template particle in
  every copy (`[*, i1]`), a flat list's as now (a prefix of none). Queries (by entity, by target, by compartment) give
  such sets without expanding (a template's matches times its copies). Labels name the path ("Particle 12 › 3,401").
- **Expanded on the GPU, per frame, the visible alone**: the page store holds the top list's pages and each template's
  pages once (not tiled). The cull gains an expansion stage before its count, scan and write: the top particles culled
  at their template's bounds (its radius scaled by the pose), the visible copies queued with their template's
  instance count, a scan, then a thread a queued copy and template particle: its pose composed, its leaf target's
  member and level, culled at its own radius, counted and written as a flat particle is today. A deeper level (later)
  repeats the stage (queue to queue). The composed pose record gains the identity words (top slot, local index: 4 B
  more).
- **The bound moves, it does not go away**: what is stored is copies + placements, but the composed poses of a frame
  are the visible expanded instances (27.6 M for 60 copies of the isg cell in view, at 36 B about 1 GB). So copies
  far away are drawn coarser (a template's coarse level, as the tiles' coarse nodes: a sample drawn larger; or the
  sphere raster), and the composed buffer has a budget (past it the farthest copies coarser, flagged, never a fault).
- **Per-instance state is sparse or by level, never dense**: an array an expanded instance would undo the saving.
  - Colors: by the leaf's own data (its entity, its target: the template's, shared by every copy) or by a level's
    (the copy's compartment: the top particle's); a theme says which level, resolved in the write kernel from that
    level's table.
  - Markers (decided, step 6): by copy, read in the vertex stage by the drawn instance's copy (its composed pose's
    word): per copy a marker for all of it, or a block of markers by its template's store indices while marked in part
    (made on its first such mark, let go once clear: a hover's copy at a time). Why not a sparse set searched in the
    write kernel: a mark then culls again, and a copy whole or every copy's particle (a selection) is a block or a
    word, not millions of pairs. A template particle in every copy (an entity's selection): the template's (the store's markers of its
    particles, or-ed with a copy's by every copy: no block), step 6b.
  - Layers (overpaint, transparency, emissive, substance, clipping, wiggle): the same, by level or sparse.
- **Picking by the composed slot**: the pick output keeps one instance word, the composed slot (the pose's index in the
  frame's composed buffer: within 2^32, the visible instances of a frame). Its identity is gathered on the GPU in the
  same submission as the pick pass: a small kernel reads the pick texels asked (a point, or a region's distinct slots)
  and copies their records' identity words into the readback buffer, read with the pick's ids (`mapAsync`, no second
  round trip). Not a read of the slot first and its record after: by then the next cull may have rewritten the buffer.
  Only nested objects pay it (their objects flagged in the pick's read); the pick targets and every other object
  unchanged. A second instance word in the pick output, the alternative, costs every pick draw of every object.
- **The CPU's readers by level**: bounds (a copy's: its pose and its template's sphere), counts (by multiplication),
  labels and loci (by path), export and stats (a template's once, times its copies), trajectories (a frame moving top
  particles moves copies; a template's own frame moves every copy, the point of sharing). A per-copy edit of a
  template's particles is not supported: it ends the sharing (the copy made a list of its own, by the user's tool).
- **WebGL: expanded in the vertex stage**: no compute, so no queue, scan or compaction; but WebGL2 decodes an
  instance in the vertex shader (`gl_InstanceID`, integer ops, `texelFetch`). The CPU culls the copies (few: the top
  list's particles, at their template's bounds) and writes the visible copies' poses to a small texture a frame; a draw
  a template and leaf target of `visible copies x that target's template particles` instances decodes its instance to
  a copy and a template particle (division by the count), composes the two poses and draws. The template's poses a
  texture stored once. No culling inside a copy (a copy's children all drawn, or a copy far away by the template's
  coarse level); markers and colors by level and the sparse set as the WebGPU path's (a binary search in a texture).
  Picking: the draw's instance id, decoded on the CPU by the visible-copy list it wrote (kept for the picked frame), no
  readback. Memory as WebGPU's (copies + placements); the vertex work the visible copies' children. WebGL's alone:
  on WebGPU it would drop the per-particle culling, the levels and the sphere raster a flat list has (step 5). While the compute expansion's kernels compile, a nested object waits for them
  (nothing drawn, no fallback), as every object culled on the GPU does ([gl/gpu-culled-views.md](gl/gpu-culled-views.md)).
- **The CPU-culled views on WebGPU**: by [gl/gpu-culled-views.md](gl/gpu-culled-views.md) (the thickness pass and a
  ray's pick culled on the GPU for every merged object); a nested object has no CPU culling of its children (as a page
  store's) and draws those views through it.

## Steps

1. **Measured first**: the cases and their numbers now: the isg grid, the yeast's grid (`particle-grid-pack.mjs`), a
   symmetric capsid, the NPC; for each the JS heap, the GPU bytes and the frames as flat lists, and what nesting
   stores (copies + placements). The grids are flat copies of one list made by a pack script (no product in the
   files): the grid script is the first producer. The isg grids measured (`merging/isg-pack-*`, WebGPU, Chrome, 720p,
   whole, levels by size on screen at 3 px; the list a view of its pack, its CPU bytes the pack's):

   | Grid  | Particles | Pack   | JS heap (scene) | GPU storage | GPU textures | Still   | Rotating |
   | ----- | --------- | ------ | --------------- | ----------- | ------------ | ------- | -------- |
   | 1x1x1 | 460k      | 16 MB  | 13 MB           | 61 MB       | 79 MB        | 0.53 ms | 0.84 ms  |
   | 2x2x2 | 3.68 M    | 65 MB  | 82 MB           | 435 MB      | 82 MB        | 0.73 ms | 1.80 ms  |
   | 3x3x3 | 12.4 M    | 198 MB | 275 MB          | 1,449 MB    | 90 MB        | 0.84 ms | 3.52 ms  |
   | 4x4x4 | 29.4 M    | 458 MB | 646 MB          | 3,426 MB    | 107 MB       | 1.29 ms | 7.84 ms  |

   Per particle, flat: about 16 B of pack, 22 B of JS heap (the merged object's copies of the instances, 19 B) and
   116 B of GPU storage; rotating, the frame bound by the GPU's cull of every particle (its wait 7.4 ms at 4x4x4).
   Nested, the 4x4x4 grid would store one cell (the 1x1x1 row) and 64 copies; what stays per expanded instance is the
   frame's composed poses of the visible (step 8's bound). Not yet measured: the yeast's grid, a capsid, the NPC.
2. **The model** (done; model `particle-nesting.ts`): target kind `particles` (`{ particles, center? }`, the center
   the template's boundary's by default), two levels (`MaxNestingDepth` 1, `ParticleNesting.check`: a template with
   templates, a list placing itself, a slot list as a template refused), the expansion (`getExpansion`: the count and
   the top particles' prefix sums, cached, none without templates), the path's encode and decode (`ParticlePath`,
   `getExpandedIndex`, `getPath`), `Particle.getBoundary` by the copies' spheres, `ParticleNesting.flatten` (the
   reference: the poses composed, the templates' ids moved past the top list's, an instance without an entity or
   compartment its top particle's, labels by path). Graphics: a template never a target geometry (its switches throw).
   Tests: `particle-nesting.test.ts`.
3. **A producer** (done; formats.md §3): the pack's templates (version 7 where a list has them, 6 otherwise: a
   target's `template` an index of `templates`, a document each whose blocks are in the same file, written once for
   every target placing it, read once into one list; not tile-ordered, no density grid of the top list's); the grid
   script (`particle-grid-pack.mjs`) nested by default (the input the template, a copy a top particle at its offset,
   about the origin, its radius the input's reach; `--nz` a third axis; `--flat` the copies written out, `--tiles` and
   `--levels` with it only). Until step 5 a nested list is drawn flattened by the representation (and by the plugin's
   transform, its themes and loci the flattened list's): the frame and GPU bytes as flat, the pack and the model's
   heap as nested. The isg grids nested: 17.2 MB a pack whatever the grid (the cell's 16.3 MB and its template's
   density grid; flat 65 to 458 MB); the 2x2x2's flatten 0.45 s, its 100 target and entity groups the flat pack's
   (counts, centroids, boxes, radii; within 3e-3 Å, the grids' quantization). Tests: `pack.test.ts` (a nested list
   written and read: its flattened instances the written one's, a template read once for the targets placing it, as
   views of the same bytes; refusals). Later, where a case needs it: mmCIF assemblies whose operator
   expression is a product, kept as a template and its copies (where the right side's placements are the same under
   each left operator).
4. **Loci, labels, queries by path** (done; model.md §5): `Particle.Loci`'s `partial` (the copies held in part, a
   template particle set each; `indices` the top particles held whole), canonical (`ParticleNesting.createLoci`);
   one instance, a copy, a template particle in every copy; `union`, `intersect`, `subtract`; queries by level
   (`selectTop`, `selectInTemplate`); size, bounds and labels by path; `toFlatLoci` and `fromFlatIndices` (by the
   indices' runs a copy). The particle and target representations draw the flattened list and give and mark the
   given list's loci (picked instances by path, a nested loci marked as its flattened one); the plugin's transform
   themes the flattened list. Theme locations by path moved to step 6 (what a theme reads where the copies are
   expanded on the GPU). Tests: `particle-nesting.test.ts` (every operation against the flattened list's indices of
   the same instances, deterministic subsets; canonical round trips; a query by a template's entity; bounds holding
   every instance), `nested.test.ts` (a representation's loci by path).
5. **The expansion** (WebGPU by the page store's cull, then WebGL by the vertex stage). Why not the vertex stage on
   both: on WebGPU it would drop what a flat list has, its particles culled each (frustum, Hi-Z), its levels and the
   sphere raster (the 4x4x4 flat grid still 1.29 ms, rotating 7.84 ms, bound by the cull), every particle of a copy in
   view drawn. The page store's cull (`paged-cull.ts`: a thread a particle of the active pages, decoded, culled at its
   target's radius, into its member's region a level, the raster's, the ghost's) already is that path; a copy becomes
   pages of its own, the template's particles stored once. Today a store page holds a slot list's particles in world
   coordinates, its poses composed once when written (`compose`: a pose slot a particle, 32 B) and the cull writing
   slots: copies' poses composed so would be a slot an expanded instance (0.94 GB at 4x4x4: less than flat's 116 B a
   particle, not copies + placements). So a copy's are composed a frame, of the visible alone. Every per-instance
   table (colors, sizes, markers, the layers) is read by `aInstance`, the pick by its instance word.
   - **5a. Templates in the store** (done; not yet drawn: the representation flattens till 5e). A template's
     particles a block of its leaf list (`ParticleNesting.getLeafList`: its targets, entities and compartments in the
     flattened list's ids, `getIdSpace`, shared with `flatten`: the themes, of the flattened list, and the members
     agree), its pages `TemplatePage` (the record's flags word): composed once in the template's frame, never an
     active page of their own (`writeActive`); its copies a table (`ParticlePageStore.setCopies`: the top particles and
     their poses, `t = p - R c` about the template's center, `q`: no particles), `composeCopies` the CPU's expansion
     (`t_c + R_c t_x`, `q_c q_x`); the top list's leaves a block as a flat list's (its copies' particles in no run).
     The paged object's `addNested` builds them (a template's block key below 0). The object's bounds the nested
     list's (the copies' spheres). Tests: `page-store.test.ts` (the template's pages once and flagged, the composed
     copies the flattened list's positions and rotations, a copy's targets its members'; refusals),
     `particle-nesting.test.ts` (the leaf list's values the flattened list's at every copy); the paged and slot
     WebGPU suites as before (native).
   - **5b, 5c. The cull by copy, the visible's poses** (done, the frame buffer's form). The active entries two words
     (a page, its copy or `NoCopy`): a page drawn as it is once, a template's once a copy (`collectActive`); the
     copies' poses a section of the tables (`uCopiesBase`). A group's box (its template's frame) placed by its copy
     against the frustum and the pyramid, so a copy out of view a test a group; a particle's position composed (the
     copy's rotation and offset) and culled as a flat one's (levels, raster split, ghost). The occlusion bits by entry
     (a template particle's shared by its copies). A visible copy instance's pose composed by the cull into the frame's
     slots after the persistent ones (`uFrameBase`, `uFrameRoom`: the copies' instances, at most 2^24, 512 MB, those
     past it not drawn), a slot a workgroup's (one global atomic a workgroup, its threads one entry's; the second
     phase's past the first's), the regions' slots the frame's: the draws, the raster and its resolve read them as any
     slot (no reader changed). The extras moved into the tables (a binding freed: the cull's eight with the poses read
     and written), their last section: grown in place (a copy, the store's new extras written; no other section moved
     or written again: a tile add as before), tested by a growth after the first sync. What it costs: tables and
     extras within one binding's limit (`maxStorageBufferBindingSize`, 128 MiB under the core limits), past it the
     buffer refused (the extras 8 B a particle of a page of many runs: millions of them). The members' instances a template's times its copies (`setCopies`; a template of no copies
     refused). The reader's form (no frame buffer, a pose more a vertex) not built: 5h measures the buffer's. Themes:
     a nested list's by `ParticleNesting.getThemeList` (a slot list of its top list's and its templates' leaves, the
     flattened list's tables: a theme colors a page by its own list, a template particle once for all copies; part of
     5e). Tests: `webgpu-paged-cull.test.ts` (every instance of three copies, one turned, and the leaves in its
     member's region, its pose the CPU's `composeCopies` within 0.05 Å; a narrow view the first copy's alone),
     `webgpu-particle-nested.test.ts` (a nested list drawn by the store's object: its frames, whole and zoomed on a
     copy, the flattened list's drawn flat, within 8 pixels; spheres of one depth not overlapping: another draw order
     colors their overlap otherwise); a quaternion product's order swapped fails them.
   - **5d. Identity** (done). The tables by store index (a template particle's: shared by the copies; step 6's levels
     and sparse set what tells them apart). A frame record's last word its segment (a member, 16 bits:
     `PagedSegmentMask`, a paged cull of more members refused) and the copy plus one above it (`PagedCopyShift`; past
     2^16 - 2 copies none, their instances not picked by copy), written by the cull as it composes. A page store's
     object (`dCulledPaged`) picks its composed pose's slot as its instance word, its group the geometry's (a site:
     kept, not repurposed; the raster's resolve its slot alike); the pick's identity gathered on the GPU in its own
     submission (`PickPass.resolvePaged`, `paged-pick.ts`: a dispatch a page store's object over the pick's texels, its
     record's store index into a resolved instance texture, its copy into a copy texture; the others' texels as
     picked), before the pick is read (the next cull rewrites a frame's slots); the readback five textures where
     resolved (`PickingId.copyId`); its kernel compiling, the pick incomplete (a ray's: none picked). The paged object's
     `getLoci`: a template page's store index and the copy (the copies in the store's order of its templates) its path
     (`ofInstance`). Tests: `webgpu-paged-cull.test.ts` (every composed instance's copy word its copy's, the leaves'
     none), `webgpu-particle-nested.test.ts` (picks through the canvas: every copy's instances by their copy, their
     group a site, the flattened list's instance of the path the nearest the hit; the leaf near, a line two pixels
     across smaller than the pick's texels whole); the cull not writing the copy fails both, the gather not dispatched
     fails the nested and the slot list's pick tests; the headless WebGPU suites (one timing test of render bundles
     failed under load, passed alone) and the visual sphere raster cases (WebGPU) as before.
   - **5e. The representation** (done). A list with templates drawn as it is by the store's
     object where its every leaf target is drawn as the store draws (spheres: `isNestedPagedDrawable`, again with the
     props), else flattened (`ParticleNesting.flatten`, a shape target's matrices, a mesh, WebGL until 5f): the
     flattened list made only there. Its params by its leaves' targets, its themes by `ParticleNesting.getThemeList`
     in the plugin's transform too (a slot list of its leaves by level: a theme colors a page by its own list, a
     template particle once for all copies; no copy expanded), the flattened draws themed by it as well (a theme of a
     slot list colors any list by its own). A slot list's spacing radius its largest page's (its own arrays empty).
     Its loci by path (`getAllLoci` its top particles whole); its leaves marked by their store indices, a copy's
     instances by step 6's copy markers; its instances' colors by an inherited entity or an index theme with step 6. Tests: `webgpu-particle-nested.test.ts` (never
     flattened where drawn so; a leaf marked; its targets as meshes flattened and its frame its flattened list's).
   - **5f. WebGL by the vertex stage** (done, WebGL2: the `instanceExpansion` capability; WebGL1, without an
     instance id, flattens). A list with templates drawn by its copies where not by a page store: its own leaves'
     visuals a flat list's (`getTopList`), a visual a leaf target of each template (`getLeafList`: its particles once,
     in the template's frame; `expanded.ts`), not merged. Its object's poses its template's copies
     (`instance-expansion.ts`: the affine rows `transform T(t) R(q) T(-c)`, the representation's transform folded in,
     and the template's instances' sphere placed by each); the CPU's cull (`SegmentedMdbList.expandedSegment`) tests a
     copy by its sphere (frustum, distance range, clip, occlusion) and draws a level of its instances times the copies
     whose sphere reaches the level's distances (`n x V`, the vertex stage drops an instance outside them), the copies
     seen written a level after the other into `tCopyPoses` (by the render item at draw time, `uCopyBase` a level's
     first); the instances read from the instance data texture (`dInstanceData`, as merged ones without base
     instances: `k = gl_InstanceID % n`, the copy `uCopyBase + gl_InstanceID / n`, its pose times the instance's),
     `aInstance` the template particle's (colors a template particle's, every copy's). Picked by the draw slot
     (`v n + aInstance`: the instance's own id, the grid reorders them), decoded by the copies the visible passes'
     last cull saw (`InstanceExpansion.seen`) into the
     path. Far copies at the template's coarse levels (`getTemplateLevels`, as the page store's): a level's sample a
     visual a target too (keyed past the target ids), a copy drawn by the objects of one level, its spacing's pixels at
     its center's distance (`getCopyBuckets`: below `ExpansionCoarsePixels`, 4, a level coarser each halving), the
     copies in the view by the template's sphere bucketed by level once a view for the template's objects; a template's
     poses one shared texture (`tCopyPoses`), an object's cull its copies' indices alone (`tCopyOrder`); a mark of a
     coarse level's by a template particle's index in it (an inverse a level). Markers by copy (step 6's `CopyMarkers`,
     the representation's templates' copies one after another, a template's particles by their index): a copy whole a
     word, in part a block; their buffer one texture of the representation's objects (`tCopyMarkers`, synced where
     drawn), the vertex stage or-ing an instance's copy's (its template particle's index in the instance data's spare
     word) into its own; a template's particles in every copy (the same in each) and every loci its objects' own
     markers. An object marked by a copy (one marked whole, or in part at its target) its uniform the texture's and its
     average a part's (its marking passes draw it), its own state kept apart (what `Visual.mark` reads). Its bounds the
     template's box's corners placed by each copy (as `Particle.getBoundary`'s), padded by an object's reach. Found on
     the way: WebGL2's asynchronous pick read checked five arrays of four started reads (the copy texture WebGPU's
     alone), every asynchronous pick failed; it fills what it read, as WebGPU's reader. The CPU's cull was an object's
     pass over its level's copies (its distance ranges by its own sphere: targets times copies, at 16x16x16 400 objects
     a frame over 1,024 copies each, 13 ms): a level's copies now ordered by their centers' distance once a view
     (`getCopyBuckets`), an object's levels a slice each of them (`getCopyRange`: the copies within its reach of a
     level's distances, `getCopyReach`, the most any copy places its instances' sphere from the template's center; a
     superset, the vertex stage drops an instance outside), where neither clipped nor occluded (else the pass as
     before). Rotating 16x16x16 22.4 to 9.99 ms, 12x12x12 13.7 to 10.9, 8x8x8 8.55 to 7.67, the smaller and still
     frames as before (the GPU's). The copy order an object's own texture uploaded a frame each (400 rows of 16 KB at
     16x16x16, 13% of its frame): an unmoved object's levels now draw from their slice of the template's order
     (`getCopyOrderImage`, `uCopyBase` the slice's start), one image and shared texture a template's view, synced
     where drawn; kept for four views at once (a pick's, another canvas'), an object whose view's image was since
     rewritten draws nothing and is culled again by the next frame (`loadPoses`); a pick's pose read from the image.
     Rotating 16x16x16 10.5 to 9.04 ms, 12x12x12 11.0 to 9.07, 8x8x8 7.54 to 6.36, 4x4x4 7.12 to 5.74 (the submit
     8.58 to 5.47 at 16x16x16), still frames as before. A moved object (a drag's preview transform, a
     representation's objects alike) had its own buckets sorted and its own order uploaded each cull: the buckets and
     the template's orders now kept by the view and the transform, its objects one sort and one image a frame
     (`merging/isg-grid-nested-<n>-drag`, another transform a frame): 16x16x16 237 to 9.08 ms, 12x12x12 107 to
     10.5, 8x8x8 35.6 to 6.60, 4x4x4 12.5 to 5.97, as rotating.
   - **5g. Tests** (done): the visual suite's `particles/nested/{spheres,mesh}` (three copies, two turned, and a
     leaf) and `particles/nested/far` (copies of a cube of 512 receding, the far ones at its coarse level: copies seen
     at a coarse level and the whole, a check), `particles/nested/marked-{spheres,mesh}` (an instance of one copy
     selected, a copy highlighted whole, a template particle selected in every copy, the leaf highlighted); their
     reference WebGL2's expansion, WebGPU's page store and its flattened mesh the same image, the camera given;
     WebGL1's and the minimum WebGL2's flattened alone), each pick on a grid one instance by its path, its
     flattened position at the pick, every top particle picked; the reference drawn by its copies, not flattened (a
     check). Headless WebGPU as before (`webgpu-particle-nested.test.ts`); the headless GL suite is WebGL1 (no
     expansion).
     WebGL2 (Chrome, ANGLE D3D11, RTX 4070 Ti, 720p, `merging/isg-grid-nested-*`), flattened before (every copy
     expanded on the CPU) and by its copies:

     | Grid     | Form      | Add and commit | JS heap (scene) | GPU bytes | Still   | Rotating | Hover   |
     | -------- | --------- | -------------- | --------------- | --------- | ------- | -------- | ------- |
     | 2x2x2    | flattened | 843 ms         | 697 MB          | 593 MB    | 4.69 ms | 4.63 ms  | 4.99 ms |
     |          | copies    | 259 ms         | 136 MB          | 185 MB    | 4.82 ms | 5.74 ms  | 0.84 ms |
     | 3x3x3    | flattened | 2,706 ms       | 2,166 MB        | 1,773 MB  | 13.0 ms | 12.9 ms  | 17.3 ms |
     |          | copies    | 201 ms         | 118 MB          | 213 MB    | 6.68 ms | 5.86 ms  | 0.93 ms |
     | 4x4x4    | flattened | 6,453 ms       | 5,350 MB        | 4,238 MB  | 28.5 ms | 28.4 ms  | 43.4 ms |
     |          | copies    | 211 ms         | 117 MB          | 213 MB    | 5.40 ms | 6.44 ms  | 0.67 ms |
     | 8x8x8    | copies    | 221 ms         | 128 MB          | 213 MB    | 6.27 ms | 7.44 ms  | 0.68 ms |
     | 12x12x12 | copies    | 212 ms         | 68 MB           | 213 MB    | 10.2 ms | 10.5 ms  | 0.73 ms |
     | 16x16x16 | copies    | 230 ms         | 95 MB           | 213 MB    | 8.83 ms | 16.4 ms  | 0.74 ms |

     8x8x8 and larger flattened refused (an array buffer's allocation failed). The memory and the setup the template's
     and its levels' (the same from 2x2x2 to 16x16x16; the JS heap between runs the collector's; WebGPU's page store's
     GPU storage 926 MB at 4x4x4, step 5h, the frame's slots most of it). Every copy its template whole (no coarse
     level) the frame the flattened list's: 4x4x4 29.1 ms, 16x16x16 1.8 s (its hover's device lost: D3D11 removed it, a
     draw past Windows' two seconds); with them still frames 5 to 10 ms at every size, rotating its CPU's cull
     (16x16x16: 37.6 ms with every pose written an object a frame, 23 ms with a template's poses one texture, 16.4 with
     the frustum tested once a template). WebGPU's page store (step 9) at 16x16x16 rotating 3.17 ms. With copy markers
     the hover 0.96, 1.27 and 0.88 ms (2x2x2, 4x4x4, 16x16x16: a hovered copy's block one object's marking draw, not
     every object's as when any marked copy drew all of them, about 8 ms).
   - **5h. Measured** (done, the frame buffer's form; WebGPU, Chrome, RTX 4070 Ti, 720p, whole, levels by size on screen
     at 3 px; `merging/isg-grid-*`, a case a page, the camera of every form the nested list's sphere): the isg grids
     written by the grid script flat (the merged culler), flat tile-ordered (`--tiles 16`: the page store's cull, as
     the nested one's without copies) and nested.

     | Grid  | Form   | Pack   | Add and commit | JS heap (scene) | GPU storage | Still   | Rotating |
     | ----- | ------ | ------ | -------------- | --------------- | ----------- | ------- | -------- |
     | 2x2x2 | flat   | 69 MB  | 58 ms          | 13 MB           | 435 MB      | 0.80 ms | 1.68 ms  |
     |       | tiled  | 70 MB  | 100 ms         | 100 MB          | 634 MB      | 0.81 ms | 3.20 ms  |
     |       | nested | 17 MB  | 24 ms          | 22 MB           | 140 MB      | 0.67 ms | 1.52 ms  |
     | 3x3x3 | flat   | 209 MB | 204 ms         | 70 MB           | 1,449 MB    | 1.04 ms | 2.80 ms  |
     |       | tiled  | 210 MB | 191 ms         | 173 MB          | 1,513 MB    | 1.02 ms | 4.73 ms  |
     |       | nested | 17 MB  | 24 ms          | 25 MB           | 407 MB      | 0.91 ms | 3.75 ms  |
     | 4x4x4 | flat   | 481 MB | 291 ms         | 180 MB          | 3,426 MB    | 1.38 ms | 5.27 ms  |
     |       | tiled  | 483 MB | 297 ms         | 196 MB          | 2,216 MB    | 1.38 ms | 7.63 ms  |
     |       | nested | 17 MB  | 26 ms          | 2-25 MB         | 926 MB      | 1.29 ms | 7.99 ms  |

     The pack the cell's once at every size, its add a tenth of the flat's at 4x4x4, its JS heap the cell's (no copy's;
     between runs 2 to 25 MB, the collector's). Still frames at most the flat's. Rotating: the tiled list's about the
     nested one's (the same cull: 7.63 and 7.99 ms at 4x4x4), slower than the merged culler's (5.27): the page store's
     cull and its second phase (at 4x4x4 rotating, nested: cull 1.99 and 3.09 ms, raster 0.98 and 1.22; tiled: 2.67
     and 2.92, 0.97 and 0.22; flat: 1.84 and 1.40, 0.93 and 0.36). The nested second phase's raster five times the
     tiled one's: no occluder grid (a flat list's the pack's density grid, a tiled one's its tiles' opacities prime
     the first phase's pyramid; with them off theirs as the nested one's, tiled's 1.20 ms, flat's 1.84). Fixed:
     `ParticleNesting.getDensityGrid` (its templates' grids placed by their copies, into the blocks a copy covers whole:
     a block on a copy's edge none, so the gaps between copies stay open; first placed with every block the seams
     were lost, more than the flat grid's own: a block, about 320 Å, coarser than a gap, about 250 Å) primes a nested
     list drawn whole: at 4x4x4 rotating 6.84 ms (the second phase's raster 0.25), 3x3x3 3.42; its image the
     occluder-less one's (0 pixels apart; the flat list's grid closes some seams of its own: a separate question); the
     representation's setup 28 to about 50 ms at 4x4x4 (once a list: a list's block a copy covers the average of the
     template blocks it spans, from a table of the template grid's sums; the first form, every template block placed
     once a copy, 335 ms). GPU storage the frame's slots most of it (every instance visible
     here: 29.4 M x 32 B, 897 MiB of 926 at 4x4x4): what step 8's budget bounds, or 5c's reader form removes (a pose
     more a vertex instead). Found and fixed on the way: a store's format the first list's (a top list without
     rotations refused its templates'), now of all its lists (`getStoreFormatOf`); the frame's slots past their room
     dropped (at 2^24 the 4x4x4's 29.4 M: the still frame 6.5 ms), now 2^25 (dropped silently still: step 8 flags
     them); a nested list's bounds by its copies' boxes' corners, not their spheres (the camera's fit a seventh
     farther); a page store's `aInstance` none (the shared serial of its instances, 4 B each: 118 MB at 4x4x4); the
     benchmark's byte notes made a lazy `aSegment` (470 MB) they counted.
6. **Markers and colors by level and sparse**: a highlight of one instance in one copy, of a copy, of a template
   particle in every copy; colors by the leaf's and by a level's data, theme locations by path. Tests: the composed flags and colors as the
   flattened list's.
   - **6a. Markers by copy** (done, WebGPU). `copy-markers.ts`: two words a copy (its marker, or `CopyMarkerBlock` and
     its block's first byte less its template's first store index), then the blocks' bytes (a block its template's
     store pages' span; let go once nothing in it is marked, reused by the template's other copies), uploaded by the
     culler where changed (`sCopyMarkers`, `dCopyMarkers` where the store has copies). The vertex stage and the sphere
     raster's resolve read an instance's copy from its pose's word (`PagedCopyShift`, as the pick) and its marker
     there, not its store index's. The object's marker status per instance and its average above 0 while a copy is
     marked (a mark moves its uniform: the object updated, the culler syncing). A loci all of a nested list by its
     instances' count, not its top particles'. Tests: `copy-markers.test.ts`; `webgpu-particle-nested.test.ts` (one
     instance by its path its copy's pixels alone, the same particle of another copy others, a copy whole; a copy's
     raster-drawn spheres tinted without the marking pass).
   - **6b. A template's marks** (done). The same particles of every copy of a template (`ParticleNesting.ofTemplateParticles`,
     a selection by entity: every copy's locals the same) marked in the store's markers at the template's particles'
     store indices (its coarse levels' by their samples), read by every copy: the vertex stage and the resolve or the
     store's marker with the copy's. An action taking away from a part of them (a copy's, one instance) makes the
     template's marks each copy's first (`CopyMarkers.include`), then applies; taking away from all of them takes away
     from a copy's own too; a toggle by copy. Tests: `copy-markers.test.ts` (include); `webgpu-particle-nested.test.ts`
     (a particle of every copy in each copy's pixels, no copy's word set; one instance of it taken away its copy's alone).
7. **The layers** by level and sparse (overpaint, transparency, clipping first).
   - **7a. A target's visuals' layers** (done). The target representation took none (its state's overpaint,
     transparency, emissive, substance, clipping and wiggle dropped, a flat list's too): its visuals' now
     (`TargetVisual.setLayers`, by the loci of its particles as its marks, `createTargetLociApply`), again where its
     objects are made (before they merge); a list with templates flattened by its flat loci; drawn by its copies, an
     object's the particles of its template every copy holds (`getTemplateWideLocals`, its coarse level's by its
     inverse). Found with it: a layer of four instances or fewer made from an object's empty layer wrote into that
     array, shared by every object (`createDataArray` reused it), so one object's cleared another's (overpaint,
     transparency, clipping, emissive, substance, wiggle, markers). Tests: `layer-data.test.ts`; visual
     `particles/layers` (overpaint, transparency, a clipping layer exempting from the clip object).
   - **7b. By copy on WebGL2** (done). `copy-layers.ts`: a copy's own layers (overpaint, transparency, clipping,
     emissive, substance, wiggle) by the layers in order (`getCopyLayerStates`): a copy has one once a layer reaches it
     by copy, whole a value, in part a block by template particle (from its value, or from the template's layers so
     far); a layer on the same particles of every copy the template's (its objects' own) and those copies' with one.
     Their buffer one texture of the representation's objects (`tCopyLayers`, by the markers' copies: six words a copy,
     none, `CopyLayerValue` or `CopyLayerBlock` and its first byte, then the bytes), set again with the state, read
     where a copy has one (`dCopyLayers`): the vertex stage's layer its copy's where it has one (`copyLayerOr`), else
     its object's. An object keeps its layer where it reaches none of its own (its define on), its transparency,
     emission and wiggle averages at least a half where its copies have them (the passes that draw them, a wiggle's
     frames drawn anew). Tests: `copy-layers.test.ts`; visual `particles/nested/layers-mesh` (a copy overpainted whole,
     a template particle in every copy, the leaf transparent, two instances of a clipped copy exempt): WebGL2's copies
     as WebGPU's flattened list.
   - **7c. The page store's** (WebGPU, done). Its layers by store index (`TargetSlotsObject.setLayers`: overpaint,
     transparency, clipping, emission, substance, wiggle; a flat list's slots, a template's particles in every copy at
     each of its levels), the merged object's own cells read where their define is on (a source's layer wired while off,
     `SourceOptionalLayers`: turned on later by its define alone); its copies' own as 7b's (`getCopyLayerStates`), by
     the copy markers' copies (`sCopyLayers`: two words a copy's kind, a value or a block by store index over its
     block's span, then their bytes), read in the vertex stage and the sphere raster's resolve by the pose's copy and
     its template particle's store index (`copyLayerOr`, generated with the culled input reader: a program's storage
     bindings now by what that reader references too). Kept with the state, applied again where the store updates: all
     where its texels are laid out again or its templates or their copies changed, else the pages written alone (their
     texels cleared, written by the layers and uploaded; the copies' own kept), none where none was written (a theme);
     a layer's loci by an earlier version of the slot list (a tile view's step: its slots the same particles) the
     store's. Tests: `webgpu-particle-slots.test.ts` (a slot list overpainted by slot, its view's tiles let go and
     added again still overpainted, cleared again); visual
     `particles/nested/layers-spheres` (WebGPU's page store as WebGL2's copies). The layers' cases wiggle too (a copy
     whole, two particles of the flat list; animated, at time 0: their pixels alone moved).
   - **7d. Measured** (isg nested, WebGPU and WebGL2, 720p; `merging/isg-grid-nested-<n>-layers`: a copy and a thousand
     template particles in every copy overpainted, set once; `-layers-update`: another copy overpainted a frame, the
     state set and the canvas told, as the plugin does: committed by the frame). Still frames with them: WebGPU 4x4x4
     0.79 ms (0.74 without), 16x16x16 0.84 (0.83); WebGL2 6.74 (5.71), 9.90 (8.85), the GPU's wait the more. A change's
     update was led by the precompile recorder signing the representation's objects again on the canvas' update
     (`holdForUpdated`: `getObjectSignature`'s defines by every schema key, `getSignatureKey`), any state set paying it:
     an object's signature not made again where the cells it reads (its defines', the values it takes) are the same
     cells at the same versions (`PrecompileRecorder`, the defines' keys a schema once), WebGL2's update 8.9 to 3.6 ms
     at 4x4x4. A frame with a change: WebGL2 7.47 ms at 4x4x4, 10.0 at 16x16x16; WebGPU 1.13 and 2.80 (its update 2.07
     ms at 16x16x16). (Committed at once by the update, `commit(true)`, it was 13.4 ms on WebGL2: the wait for the GPU
     the longer, not the frame's passes; what a change uploads small, the uniform buffers' flushes 29 KB a frame, 16 KB
     without, the copies' layers' texture 32 KB, written in place no faster.)
8. **The frame's bound**: far copies by the template's coarse level, the composed buffer's budget.
   - **8a. The frame's slots by what a cull needs** (done). Measured first: a frame's copies' instances written
     (the cursor read back) at 4x4x4 0.56 M a still frame, 0.9 M a rotating one on average and 1.2 M at most, of
     29.4 M (the pyramid's occlusion hiding the rest; 3x3x3 at most 0.61 M of 12.4 M). So the room at first
     `PagedFrameInitialSlots` (2^21, 64 MB), the cursor of each cull read back (one read at a time), past the room: the
     view culled again (the culler's `frameVersion`, in the merged culler's data version) and a draw asked for (a
     settled promise tracked as pending work), the room grown to half again what it needed; not grown where a pyramid
     is to come (`PagedCullView.pyramidToCome`: the view culled in two phases, `Renderable.prepareCull`'s `twoPhase`,
     none built yet: an app's first frames need every instance in view, 29.4 M; grown so, the budget would be none);
     the first cull to fit after one that did not culled once more (its pyramid of a frame with instances missing:
     a still view kept its lists, at 4x4x4 still 1.72 ms against 1.40). The room at 4x4x4 4.8 M (154 MB, from 897
     MiB), at 2x2x2 and 3x3x3 2 M; frames and images as before (0 pixels apart; 4x4x4 still 1.40 ms, rotating 6.84;
     3x3x3 1.05 and 3.38; 2x2x2 0.78 and 1.60); the first frames' overflow not drawn whole (a frame or two). Not
     shrunk once grown. Within what the device's limit leaves after the persistent slots (a zoom into the yeast 2x2's
     copies: its template's 36 M slots, 1.16 GB, and a frame's 32 M, 1.07 GB, past 2 GiB, the poses' buffer not made);
     a room at its cap overflowed asks for another cull while each needs fewer than the one before, none past that
     (a still view's every frame culled else; a zoom's first cull, of the last view's pyramid, the isg 16x16x16's 87.6
     M of 33.5 M: the rest not drawn till the camera moved). That cull's slots went to whichever copies' workgroups came
     first (a third of the copies in view drawn, far ones too: its first frame 369 K pixels of background); the select
     now orders its entries near first (`PagedCullOrderBuckets`: 256 buckets of a copy's center's depth over the copies'
     range, a counting sort in two passes), the dropped copies the far ones (about 220 px along the cells' seams), and
     an eighth of a room below the copies' instances (a template with coarse levels) is the fallback's
     (`PagedFallbackShare`): a copy group a phase dropped drawn after the last a level coarser, without the pyramid (one
     primed from the occluder grid hides the dropped copies' own instances: 11.1 M of 14.3 M), none left. Without the
     order a third was at its edge (10.9 M of 11.2 M wanted, holes where over), the coarsest level fit (1.4 M) but drew
     the whole frame in molecule-sized spheres; a continuous zoom's frames as before (median 23.7 ms against 23.2).
     Rocking near (0.35 of the distance at 1280x1130, `merging/isg-grid-nested-16x16x16-near-rotating`): 45 ms frames in
     Chrome, the GPU busy (Nsight GPU Trace: 99 % active, compute warps 87 % of the slots, 36 % of their issue stalls at
     barriers), each phase testing about 350 M template particles, 98 % of those tested one by one hidden. The pack
     written before templates were tile-ordered (a group a target across the cell: none culled whole); written again
     39.9 → 29.2 ms (16 tiles an axis 28.8: the default kept). A copy the pyramid hides whole left out by the select
     (marked; the second phase's select appends those the frame's pyramid shows after the first's entries, their groups
     tested whole) 31.9 → 28.9 ms. The groups culled first, a thread a group (`paged-cull-groups`: tested whole, those
     not culled whole the survivors, by their entry's order bucket, the cull a workgroup a survivor: 630 k of 1.05 M
     workgroups a phase culled whole no more launched) 28.9 → 23.3 ms, the cull 8.35 + 9.9 → 5.4 + 6.0 ms and the group
     passes 0.45 + 0.77; appended at one cursor, the order the workgroups ran, the first frames of a pyramid needed 17 M
     frame slots (not 9 M) and grew the rooms by them (poses 665 MB, not 410). The survivors' room the selected entries'
     groups (read back, grown), not every active entry's (2.7 M: 340 MB). Tests: `webgpu-paged-cull.test.ts` (a room of
     100: the copies' instances past it not drawn, a cull again asked for, not grown while a pyramid is to come, grown
     where none is, then every instance; not growing fails it; a limit below the copies' instances: the room capped, a
     cull again asked for once, and where one needs fewer, none for the same again; the copies selected near first from
     either side; a leaf block made resident past the room: its copies dropped drawn at their coarse level; a copy the
     pyramid hides whole left out, added by the second phase; the survivors' room grown). 5c decided: the frame buffer's
     form (budgeted, its readers unchanged; the reader's form not built).
   - **8b. What else holds a copy's every instance**, GPU storage after the runs (the 5h table's at the first frame:
     the lazy buffers not made yet; `--byte-notes` gives both now), 4x4x4 rotating: nested 1,534 MB (the paged cull's
     indices 449 MB: a block of regions a level, every expanded instance's room at each of four levels; the sphere
     raster's chunks 899 MB, by its instances; the poses 167 MB), flat 4,346 MB, tiled 4,498 MB (each case after
     the ones before it: what they left counted; alone, nested 1,533 MB). Done: the regions' rooms by what a cull
     needs (as the frame's slots): a region's room at first a sixteenth of its member's instances, at least 4,096
     (all where fewer; the ghost's too: its every visible instance, at 4x4x4 118 MB where full), packed one after another (no block a level), a slot past it not written
     (the kernel's check, the draw records cut at it), the counts read back with the frame's cursor (the second
     phase's with the first's), grown to half again, the view culled again; not grown while a pyramid is to come.
     The raster's first chunk capacity the raster regions' rooms (it was the store's instances). Nested 4x4x4 alone
     591 MB (from 1,533: the cull's buffers 104 MB from 461, the raster's chunks 313 MB from 899, its readback
     growing them), frames as before (still 1.41 ms, rotating 6.85), images whole. Tests: `webgpu-paged-cull.test.ts`
     (20,000 instances of one member, two levels and the ghost: rooms of 4,096, what was written within them, not grown while a pyramid is to
     come, then grown, every instance; the tests reading exact regions with full rooms). The first rooms within
     `PagedRegionInitialSlots` (2^24, 64 MB) in all, the shares scaled down: a sixteenth of every instance outgrew the
     device's buffer limit at 16x16x16 (1.9 G instances), the regions past it not written, nothing drawn; now drawn,
     GPU storage 16x16x16 4.5 GB to 0.78, 12x12x12 4.0 GB to 0.67. Next: far copies by the
     template's coarse level; the raster's chunks' cap (every raster's, not nesting's: 8f).
   - **8c. Far copies by the template's coarse level** (done). `ParticleNesting.getTemplateLevels`: per template up to
     four levels, `Particle.sampleByTarget` (each target's every `8^L`-th particle in Morton order over its box, at
     least one, while a level keeps at most half the one before), drawn larger by the square root of seen over kept
     (the area covered kept: the cube root, the volume kept as the tiles' coarse nodes, left the isg cell's thin outer
     layer holed, its interior seen through on a grid of copies). A level a template block of its own with the same
     copies (`StoreCopies.level`, `levels`, the template's `center` and `spacing`: its bounding sphere's volume a
     particle, a cube's side); an entry's page word carries its level and the template's levels; the kernel picks a
     copy's level by the template's spacing in pixels at its center's distance (`copyLevel`: a level a halving below
     `PagedCoarsePixels`, 4), a workgroup of another level's page none. Colors, picks (a level's particle its
     template's: `coarse` in `paged.ts`) and copy markers (a top particle's copy at each level, a path's particle in a
     level's sample) by the template's. The sample order and the targets in it once for a template's levels (typed
     arrays: the isg cell's levels 68 ms from 144). Measured (isg nested, rotating; still about as before or
     better): 4x4x4 7.57 to 2.89 ms, 8x8x8 53.5 to 10.4 ms (GPU storage 3.9 to 2.7 GB); 12x12x12 30.3 ms, 16x16x16
     69.2 ms, their cull (both phases 28 and 66 ms) by every copy's entries of every level, those of another level
     launched and leaving. Setup (representations) 130 to 180 ms (38 before: the levels and their blocks). Tests:
     `particle-nesting.test.ts` (the sample, the levels), `webgpu-particle-nested-far.test.ts` (near whole as
     flattened; far a sample of about its color and cover, picked by a sampled particle's path, marked).
   - **8d. The cull by the selected entries** (done). `paged-cull-select` (a thread an active entry): a page drawn as
     it is, a template's page at its copy's level (`copyLevel`, no longer per workgroup in the cull) where the copy's
     sphere (its template's radius about its center, `StoreCopies.radius`, padded by the largest reach) is in the
     frustum (not with the ghost), compacted into the tables' selected section (`uSelectedBase`: a header vec4, their
     count, then an entry's index a word) by an atomic count; `paged-cull-select-args` (one thread) the count into the
     header and the cull's indirect dispatch (rows of `PagedCullMaxRow` workgroups). Both phases dispatched
     indirectly over the same selection, a workgroup's entry the selected one, the occlusion bits by selected index.
     The cull's eight storage bindings unchanged (the selection in `tables`, its count in a buffer of the select's).
     Measured (isg nested, rotating): 4x4x4 2.89 to 1.57 ms, 8x8x8 10.4 to 1.65, 12x12x12 29.9 to 2.77, 16x16x16
     67.6 to 3.91 (cull both phases 64 to 2.0 ms, the select 0.02); still frames and images as before.
   - **8e. Setup and the bounds** (done). A level from the one before (`sampleByTarget`'s `finer`: every eighth of a
     target's kept), the far copies' spacing by the template's own cached sphere (not its leaf list's: a pass), setup
     (representations) 130-180 to 96-141 ms. The padding by the largest reach times the largest scale of all (a
     template's fourth coarse level's 63) made every group's, copy's and the object's sphere loose: a group's now by the
     largest reach times its own largest scale (the compose's, `groupBoxes`' w), a copy's by its block's (the copies'
     fourth vec4), the object's by the blocks drawn as they are (not the coarse levels': a far copy's coarse particles
     at the scene's edge past it by a little): isg nested 4x4x4 50,726 to 21,574 Å (the copies' spheres 19,931),
     16x16x16 109,761 to 80,609. A level kept where it keeps eight particles a target on average at least
     (`MinLevelKeptPerTarget`): the isg cell's fourth (two a target, drawn up to 63 times larger) none, its third
     (up to 22.6) the coarsest; images at 8x8x8 to 16x16x16 as before (their views at the third or finer). Tried and not taken: a target's scale by its box-counting dimension (its cells at the sample grid's
     resolutions; the square root for a layer, the cube root for a volume): its particles drawn smaller (the third
     level's median scale 6.8, from 16.9), the isg cell's interior seen through again as with the cube root.
   - **8f. The cull's memory** (done). GPU storage at 16x16x16 rotating 778 MB: the paged cull's 399 MB (its occlusion
     bits 285 of them, a bit for every active entry's particle: every copy's every level's pages), the sphere raster's
     chunks 284 MB, the poses 88 MB. The occlusion bits by the selected entries the cull runs over
     (`PagedOcclusionInitialEntries`, 4,096 at first: 4 MB), the select's count read back where they hold fewer than
     the active entries, grown to half again, the view culled again; an entry past them not tested against the
     pyramid (drawn in the first phase, none in the second: culled less, none hidden). 16x16x16 rotating 778 to
     480 MB, 3.84 to 3.17 ms; 4x4x4 to 12x12x12 458-513 MB, frames as before. The raster's chunks left as they are
     (280 MB at 4x4x4 and 16x16x16 alike): their first capacity the raster regions' rooms, bounded by the regions'
     budget and the device's limit; grown by a frame's need, never shrunk (a shrink a path of the shared scratch
     buffer's: every sphere raster's, not nesting's). Test: `webgpu-paged-cull.test.ts` (bits for one entry: the others
     drawn in the first phase, grown, then as with bits for all).
9. **Measured** (done for the isg grids and the yeast's 2x2; a capsid and the NPC no packs). WebGPU, Chrome, 720p,
   whole, levels by size on screen at 3 px (`merging/isg-grid-*`, `merging/yeast-grid-*`; still, rotating, `-hover`:
   per frame another instance highlighted, the one before not, redrawn over the last image). GPU storage after the
   runs, each case alone; the JS heap retained by the scene; load the representations, add and commit, first frame:

   | Grid    | Form   | Pack   | JS heap | GPU storage | Load (ms)    | Still   | Rotating | Hover   |
   | ------- | ------ | ------ | ------- | ----------- | ------------ | ------- | -------- | ------- |
   | 2x2x2   | flat   | 69 MB  | 14 MB   | 550 MB      | 28+44+803    | 0.82 ms | 1.70 ms  | 5.66 ms |
   | 2x2x2   | tiled  | 70 MB  | 84 MB   | 728 MB      | 102+189+371  | 0.80 ms | 3.31 ms  | 13.8 ms |
   | 2x2x2   | nested | 17 MB  | 57 MB   | 194 MB      | 99+15+99     | 0.77 ms | 1.64 ms  | 0.66 ms |
   | 3x3x3   | flat   | 209 MB | 83 MB   | 1,837 MB    | 35+113+90    | 1.16 ms | 2.85 ms  | 20.8 ms |
   | 3x3x3   | tiled  | 210 MB | 110 MB  | 1,804 MB    | 103+280+83   | 1.03 ms | 4.90 ms  | 27.5 ms |
   | 3x3x3   | nested | 17 MB  | 48 MB   | 319 MB      | 130+12+55    | 0.93 ms | 1.47 ms  | 1.06 ms |
   | 4x4x4   | flat   | 481 MB | 128 MB  | 4,346 MB    | 50+240+239   | 1.35 ms | 5.25 ms  | 40.0 ms |
   | 4x4x4   | tiled  | 483 MB | 76 MB   | 3,109 MB    | 98+358+63    | 1.37 ms | 7.84 ms  | 27.0 ms |
   | 4x4x4   | nested | 17 MB  | 47 MB   | 493 MB      | 81+12+49     | 1.03 ms | 1.76 ms  | 0.56 ms |
   | 8x8x8   | nested | 17 MB  | 58 MB   | 466 MB      | 65+26+27     | 1.33 ms | 2.83 ms  | 0.45 ms |
   | 12x12x12| nested | 17 MB  | 46 MB   | 516 MB      | 72+14+56     | 1.57 ms | 2.36 ms  | 0.71 ms |
   | 16x16x16| nested | 17 MB  | 27 MB   | 515 MB      | 86+22+62     | 1.37 ms | 3.18 ms  | 0.72 ms |
   | yeast   | flat 1 | 714 MB | 378 MB  | 4,655 MB    | 495+399+797  | 1.79 ms | -        | -       |
   | yeast 2x2 | nested | 714 MB | 550 MB | 4,252 MB   | 6203+811+383 | 0.87 ms | 1.65 ms  | 5.97 ms |

   (flat 1: one yeast cell, tile-ordered with levels, `merging/yeast-pack`; the flat 2x2, 2.8 GB, past what the
   page fetches.) Nested, the isg grids' memory flat from 4x4x4 on (the template once, the cull's and the raster's
   buffers by what a frame draws), their frames within 2x of the 2x2x2's to 16x16x16 (1.9 G instances); the flat and
   tiled forms' hover redraws by every instance (5.7-40 ms). Found and fixed here: a hover's frame at 16x16x16 7.2 ms,
   the paged cull's sync collecting every copy's active entries and poses again where nothing it reads changed (a
   marker's change: the object updated), now none (2.0 ms). The yeast nested (a template of 31.5 M particles): a
   hover 73 ms, every mark scanning the template's markers whole (`getMarkersAverage` 5.7 s of a profile's 20, the
   copy markers' `countMarked` 2.7, `applyMarkerAction` over a copy's block 2.4) and the copy markers' every block
   uploaded (a template's size and its levels' a mark); now the store's markers their own aggregates (`storeCells`: the object's by them and the
   copies', a clear store's clear none), a copy's block its marked positions (an action on the copy whole by them,
   past a sixteenth of the block its every byte), the bytes changed alone uploaded (`CopyMarkers.takeChanged`): 73 to
   6.0 ms (the isg grids' hovers 1.6-2.0 to 0.45-1.06 ms, the table's). Open: its 4.3 ms waited past the GPU's timed
   passes (0.75 ms). Not the draws' count: the marking redraw's 3,552 indirect draws are the records with instances
   (compacting them as a still view's next frame does changes none), the still frame replays the same from a recorded
   list (its 19 the stats' count of that). Its wait the GPU process's CPU, Dawn's handling of the redraw's draws
   encoded directly (29 chunks a frame drawn directly, 3,559 draws): by a trace (`scripts/benchmark.mjs --trace`), a
   hover frame's Dawn commands 5.8 ms on the GPU process's main thread (its submit, the D3D12 lists recorded, 3.1 of
   them; the rest the commands as they arrive), the still frame's, the same records replayed, 0.76 ms (0.41); the
   GPU's passes 0.59 ms. The store's chunk (one row) was drawn directly as a chunk of few rows changing every frame,
   its row 3,559 draws; and its flush failed, its recording's globals (a frame's mask pass, at another ring offset than
   a marking redraw's) in the digest it compares. Now a chunk of few rows drawn directly only where it drew few draws
   (64), the digests without the globals (a replay compares them), and a chunk recorded again for its globals alone
   not taken for draws changing: the hover flushed and replayed, 6.0-8.0 to 2.1 ms (its wait 0.60), the other hovers
   and recolors as before. Its setup 6.2 s (the template's levels and pages, the isg cell's 0.1), 7.4 s in a later
   session: of it the template's levels 6.5 s (profiled), the sample order's radix sort 3.8 (its four passes each
   reading a key out of order, through the index) and the positions decoded into a copy. Now the sort by bucket of
   the top digit, each bucket's passes with its keys in cache (`radixSortIndices`, the same order), the positions
   decoded where read: 7.3-7.6 to 3.9-4.2 s. The levels' sample picks where it reads a particle's target by the
   order, into runs by target id put in order after, not a target a particle written in that order and read again
   (`pickByTarget`, the same sample): 1.39 to 1.06 s, the setup 3.6 s.
   Left: that pass's read out of order (0.7 s) and the template's sphere 0.9 s (`getBoundary`'s pass over its
   particles). Not taken for it: the pack's target bounds (a target an ingredient spread over the cell: the sphere 38%
   larger for the yeast, 61% for MG, 12% for the isg cell) nor the box's center and the largest distance from it
   (0.19 s, but 10% larger for the yeast and the isg cell, 12% for MG): a camera's fit and the copies' spacing by it.
   Then, from a profile of its load (the index theme): the template's colors filled for it and each of its three
   levels (a gather from one fill now, a word a color), and a template written by target in its levels' sample order
   with its boundary (the pack's `sampleOrdered` and `boundary`, formats.md: no sort, no pass, the picks the same): the
   benchmark's representations (uniform colors) 3.68 s to 3.35 by the fill, to 1.11-1.15 by the pack (the pack 200 B
   larger); its frames as before (still 0.89 and 0.88 ms, rotating 1.73 and 1.68). The first level's pick of such a
   template by its runs (a target's every stride-th, not a pass over its 31.5 M targets: 0.4 s in a profile): 1.17-1.25
   to 0.71-0.80 s. Not taken: the first sync's page writes merged by run (its positions and rotations, about 576 MB in
   about 9,000 writes): the add and commit 0.49-0.51 to 0.36-0.38 s, the first frame 0.35-0.39 to 0.47-0.49 (the same
   bytes copied later: by the bytes, not the writes).

## Decisions open

- **Decided**: two levels first, deeper by measurements; WebGL by the vertex stage, not flattened; the CPU-culled
  views by [gl/gpu-culled-views.md](gl/gpu-culled-views.md); picking by the composed slot gathered on the GPU;
  WebGPU by the page store's cull, not the vertex stage (its particles culled each, its levels and the sphere raster
  kept: step 5).
- **The visible's poses** (decided, step 8a): composed into the frame's buffer, its room by what a cull needs (at
  4x4x4 154 MB, not 897 MiB), its readers unchanged; the reader's form not built. A template's particles have no
  composed poses of their own ([particle-residency.md](particle-residency.md) step 2: the cull composes a copy's visible
  instance from the values as stored; the yeast's 1 GB of template poses none).
- **The extras in the tables** (decided, kept; step 5b, reviewed at the end): a binding freed for the cull's poses, at
  the cost of one binding's limit for tables and extras. The reader's form not built (8a: the frame buffer's), so the
  poses stay the cull's eighth binding; extras of their own its ninth, past the core limits' eight (the device
  requests the adapter's maxima: this desktop's 16 a stage, 2 GiB a binding; a core-limits fallback needed for it).
  And the binding's limit not the first a store reaches: its positions one buffer under the same limit (8 B a
  particle), the extras at most as many bytes (8 B a particle of a page of many runs), the tables' other sections a
  page's, a run's, a member's and a copy's. Steps 6 and 7's buffers (copy markers and layers) the draw stage's, not
  the cull's.
- **Sparse overrides' scale** (decided: (a), after the measurements): markers are by copy (6a): a block a copy marked
  in part, a byte a store index of its template's pages (the isg cell's 57 pages 466 KB, the yeast's 3,841 pages
  31.5 MB); the same particles of every copy the template's (6b: no block). The copies' own layers alike (7b: a block a
  copy layered in part, its template's particles times the kind's bytes). Neither has a cap: the copy markers' and
  layers' buffers grow by half again to what the blocks need (a copy layer's offset 30 bits, 1 GiB), past the device's
  binding limit not made (what the store draws then not measured). What makes many blocks: a selection differing by
  copy (a sphere about a point across a grid of copies, a query by a copy's own property); a hover none (a copy at a
  time, let go). Plan:
  - **Measured first** (done, on the CPU: the packs' templates and copies, the blocks by `createCopyMarkers`' layout,
    its bytes checked against it; the levels' blocks where a level has a selected particle). Now, (b), (c) (a word a
    template page a copy in part, the pages touched whole):

    | Selection                                         | In part | Pages touched | Now      | (b) 2 bits | (c) by page |
    | ------------------------------------------------- | ------- | ------------- | -------- | ---------- | ----------- |
    | isg 16x16x16, a sphere about the center, r 0.05   | 8       | 54.5/57       | 4.4 MB   | 1.1 MB     | 4.2 MB      |
    | the same, r 0.15 (8 copies whole)                 | 80      | 37.4          | 43.9 MB  | 11.0 MB    | 29.2 MB     |
    | the same, r 0.3 (208 whole)                       | 344     | 40.7          | 187 MB   | 46.8 MB    | 136 MB      |
    | isg, an entity (2,757 particles) every 16th copy  | 256     | 2             | 141 MB   | 35.1 MB    | 10.6 MB     |
    | the same, every 2nd copy                          | 2,048   | 2             | 1,124 MB | 281 MB     | 84.3 MB     |
    | yeast 2x2, a sphere in copy 0, r 0.1 of its reach | 1       | 1,212/3,841   | 36.0 MB  | 9.0 MB     | 12.7 MB     |
    | the same, r 0.33                                  | 1       | 1,976         | 36.0 MB  | 9.0 MB     | 19.4 MB     |
    | yeast, an entity (5,216 particles) in 2 copies    | 2       | 2             | 71.9 MB  | 18.0 MB    | 0.1 MB      |

    (r a fraction of the grid's side; the levels' blocks 15-18% of the leaves'.) A sphere's particles on most of the
    isg cell's pages (a template's order by target), (c) saving 5-33% there, 46-65% in the yeast copy; an entity's on
    two pages, (c) a thirteenth of now (isg) to under a hundredth (yeast). An entity in half the isg 16's copies past
    the core limit (128 MiB) by 8x, within this machine's 2 GiB (its buffer grown half again: 1.69 GB). Past the
    limit: the poses' buffer showed it (a zoom into the yeast 2x2, fixed in 8a): the buffer not made, the copy into
    it failing validation; the copy markers' not run past it.
  - **(a) A budget, a coarser granularity past it** (chosen): the blocks' bytes within a budget (a share of the binding
    limit); past it a copy marked in part drawn marked whole (its word, no block), its layer the value of its first
    part, a warning; the CPU's state exact (loci, picks, the next marks), what is drawn coarser. Bounded, nothing
    refused; a selection of thousands of copies in part shown by copy.
    - **Markers** (done; `copy-markers.ts`): the budget half the device's copy marker binding (`getCopyBlockBudget`:
      WebGPU its storage buffer limit, 1 GiB at 2 GiB; WebGL the copy texture's rows, 4,096 texels of four bytes
      by its largest texture size, 128 MiB at 16,384; 64 MiB without a device), read at each layout (a
      representation's graphics set later); the buffer's growth by half again capped by it. A copy marked in part with
      no block within it (none of its template's let go to reuse) marked whole by its word, a removal clearing it
      whole; a template's marks or-ed in one value; counted (`coarseCount`), a warning once a layout. Tests:
      `copy-markers.test.ts` (copies past the budget marked whole, a block let go reused, a template's marks into one
      value, the budget read again at a layout; the device's budgets).
    - **Layers** (done; `copy-layers.ts`): the same budget for the copies' own layers (their binding's), shared by a
      layout's kinds and templates (`CopyLayerRoom`), the kinds not given again holding theirs; a block's cost the
      expanded objects' a template's particles times the kind's bytes, the paged store's its pages' spans (its
      levels' too) times them; checked as a block is made (`getCopyLayerStates`: none past it on the CPU either). A
      copy with no room its value, that layer and the ones after it on the copy applied whole (in order: the last
      wins); counted, a warning once a layout. Tests: `copy-layers.test.ts` (a copy past the room its value, a
      block's cost given; a layout's budget across kinds, the image's words a value past it). The paged store's not
      tested past its budget (the same states, its cost alone its own).
  - **(b) Markers two bits a particle** (highlight, select: four a byte): a block a quarter (the isg cell's 117 KB),
    read by a shift in the vertex stage; layers not (their bytes a kind's value).
  - **(c) Blocks by page**: a copy's block only of the store pages it marks (a page's offset a word, none for a clear
    page): bounded by the pages a selection touches, not the template's size; small where a selection keeps to few
    pages (a target's run: the template's pages by its order, by target), not for one by space (its particles across
    most pages).
  - (b) and (c) where (a)'s budget is reached by selections in use, by the measurements; (a) stays the bound.
- **Colors by level** (parked: no data needs it). The nested packs drawn (the isg grids', the yeast's) have their
  template particles' own entities and compartments (the isg cell's top particles none): drawn as flattened. Taken up
  where a list's template particles have none and its top particles do (a copy's entity its top particle's), checked
  as such a pack comes (the flattened list's colors against the nested one's, `particle-nesting.test.ts`'s form). Plan
  then, as sketched:
  - **The gap**: a template particle without an entity or compartment of its own takes its copy's top particle's where
    flattened; drawn nested the default color by its template.
  - **By groups of copies**: the copies grouped by what their top particles give, a variant of the template colored
    by the theme as it is a group (`getLeafList` with those values), its colors a block after the store's; a copy's
    word (a third in the copy markers' buffer) its group's block, read as a copy's markers are; nothing for data
    without (one group).
  - **Not covered**: themes by instance (a particle's index) not by level: every copy its template's colors, the
    flattened list's each instance's (dense: not done).
- **Templates and tiles**: templates stored whole (not tiled, not streamed); a template too large for that tiled as
  a list of its own, or not nested. Measured (the yeast nested 2x2, [particle-residency.md](particle-residency.md)): a
  view needs 1.3% of its 4,391 pages whole (its copies on coarse levels), 46-71% near after the Hi-Z pyramid (a page a
  target's particles spread over the cell), 13-55% were its leaves by space. Taken up there: the coarse levels
  resident, the leaf pages by what the cull draws, a template tile-ordered measured.
- **Copies past 2^16 - 2 picked** (parked: 65k copies enough for now). A frame record's words are full (the pose's
  translation and scale; the rotation's two words, the store index, the member and the copy's 16 bits). The option:
  the rotation in its smallest three components (16 bits each, the largest's index in 2 bits, rebuilt from them: it is
  at least a half), 14 bits freed for the copy's high bits (2^30 copies); by its quantization (estimated, not
  measured) a rotation about 1e-4 rad off at worst (0.2 Å at 2,000 Å from the center; about twice the four components'
  now). Not its w rebuilt from x, y, z: near w = 0 off by up to 0.015 rad (29 Å at 2,000 Å).
