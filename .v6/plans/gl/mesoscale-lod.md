# Plan: levels of detail for mesoscale scenes

## Status

- **Current**: proposed, nothing built. Decisions open (see Open).
- What the meso\* viewer adapter needs from graphics to draw a whole scene (its render IR: a list a placement set,
  instanced targets, transient looks a frame) at every zoom, from a cell in view to a protein's atoms, scalable,
  holding each datum once, without popping, within a frame budget.
- Builds on: [screen-space-lod.md](screen-space-lod.md) (sphere levels by pixel size),
  [../particle-nesting.md](../particle-nesting.md) (templates, copies' coarse levels),
  [../particle-tiles.md](../particle-tiles.md) and [../particle-slots.md](../particle-slots.md) (tiles, slots, coarse
  tiles as Morton ranges), [../particle-residency.md](../particle-residency.md) (pages by what the cull wants, the
  fallback), [gaussian-splats-lod.md](gaussian-splats-lod.md) (streamed splats), [../streaming.md](../streaming.md)
  (one streamer for every kind: sources, scheduling, budgets, upload pacing, fallback and fades by unit), [two-phase-culling.md](two-phase-culling.md),
  [sphere-compute-raster.md](sphere-compute-raster.md).

## Question

A mesoscale scene spans four orders of magnitude of scale on one screen: a cell 10⁴ Å across in view, its proteins
3-100 Å, their atoms 1.7 Å. What reaches a pixel is one of:

| Size on screen (a target's radius) | What shows                         | What should draw it                                    |
| ---------------------------------- | ---------------------------------- | ------------------------------------------------------ |
| > ~30 px                           | atoms, chains, ligands             | the target's atoms (impostors), or a structure's full representation (cartoon, ball and stick) for a few |
| ~3-30 px                           | a molecule's shape                 | a surface (blob, Gaussian), a sphere level (strides)   |
| ~0.5-3 px                          | a blob of the target's color       | one primitive a particle (sphere, ellipsoid, splat)    |
| < ~0.5 px                          | a region's color and density       | an aggregate a region (a coarse tile, a copy's coarse level, a density) |

Each row exists in part today, each as its own mechanism, chosen by its own metric, with no transition between rows
and, between representations, a copy of the instances each.

## What exists

- **Sphere levels within a geometry** (`lodPixelSize`, `LodLevelsValue`): nested strides 1, 10, 100, 1000, a level
  every stride-th item scaled up, chosen per sphere by size on screen (CSS px, `LodView`), culled per cell and drawn per
  level (`uLod`, `uLodLevel`, `aSegmentLod`). No fade (decided in [screen-space-lod.md](screen-space-lod.md): "no level
  chain needed one").
- **An object's distance range** (`lod`: min, max, overlap in Å): the whole object, a cell or an instance culled
  outside it (CPU `cullSegment`, the merged culler's `uLod`), faded in its overlap by a 4×4 ordered dither
  (`fade-lod`, `dLodFade`). The only way to switch representations today: two representations with adjacent ranges,
  each with its own instances, transforms and culler (the transforms were 64 B an instance a copy,
  [../mesoscale-scale.md](../mesoscale-scale.md)), in Å, not pixels; the two sides' dithers not exact complements.
- **Copies' coarse levels** (nested, `PagedCullView.coarse`, `PagedCoarsePixels` 4): a template's samples drawn
  larger where its spacing is under 4 px, a level each halving; the fallback for leaf pages not resident. Its own
  metric (spacing), no fade.
- **A particle budget for copies** (`maxParticles`, `PagedCullView.maxInstances`, WebGPU): the copies' coarse pixels
  doubled a step of a bias while a cull read back needs more, begun at the coarsest, as many steps finer a read back as
  their seen growths keep within ([../particle-residency.md](../particle-residency.md) 4d.5). A read back about a
  second where the GPU is busy: a few to settle.
- **Coarse tiles** (Morton ranges, `ParticleTileVisibility`): which tiles load; a coarse tile the union of its tiles'
  runs. Nothing draws a tile's aggregate.
- **Geometry**: spheres (impostors, the compute raster), mesh, texture mesh (GPU surfaces), cylinders, lines, points,
  text, image, direct volume, Gaussian splats (sorted, streamed SOG planned). No mesh levels; the target surfaces
  (blob, Gaussian, isosurface) at one resolution.
- **One instance source a target representation**: the page store (positions, rotations a store index), composed on
  the GPU into the frame's poses (a 32 B pose and a 4 B index a visible instance), read by index (`dCulledPaged`).
  One geometry kind a target representation; markers and layers by store index and by copy (sparse, budgeted).

## Goals and how they are measured

- **Scalable**: the frame's cost by what is on screen, not by the scene: drawn primitives a pixel bounded (the
  postsyn zoomed drew 48 spheres a pixel at `quality`; 0.45 ns a sphere), the yeast cell and the isg 16³ nested within
  the frame budget at every zoom.
- **Memory**: an instance's data once, whatever represents it; a target's geometry once a level; a level's geometry
  only while some instance draws it (resident by use, as pages are). Accounted: bytes a particle and a target by band.
- **Smooth**: no pop moving the camera, none as data streams in. Measured: a slow dolly and a pan, the image change
  between frames past what the motion explains (the frame against the previous one reprojected; a pop is a block of
  pixels changing at once), and by eye (images at the default level: a metric alone misjudges dithers and blurs).
- **Fast**: a transition costs no second draw of the whole (only the instances inside a transition draw twice), no
  per-frame CPU work by instance, no stall while a level or a page loads.
- **Live**: a scene that moves (motion tracks, physics, dynamics) costs by what moves and how, not by re-uploading
  every pose a frame; what a motion can be written down as is evaluated where it is drawn.

## Fit with meso\*

What the meso\* architecture (spec, IRs, runtimes, the viewer adapter over the render IR) asks of this plan, and
what this plan asks of it. A "copy" in meso\* is one placed instance of a population (a particle here); a copy here is
a placed template.

- **What fits**: levels, bands and fades are derived from the view by the renderer, never state in the spec or an IR
  ("derived data is never edited"); overrides (6) are the player's transient transforms and looks over committed
  lists, a freeze the list's next version; bands share the one instance source a population's layer edits act on
  (a recolour a layer edit at every band); a flat list's pages are views of its arrays, so the compiler's cached
  placement set can be the store's CPU data; a `primitive` model is a target whose one band is the proxy (5); the
  streaming runtime's tile-source port is [../particle-tiles.md](../particle-tiles.md) and
  [../particle-residency.md](../particle-residency.md).
- **The render IR must keep what repeats nested** (asked of meso\*): `resolveDefinables` expands a definition at
  its uses and the scene pass resolves a population to flat positions and rotations. An assembly, a definition used
  many times or a lift is a template and its copies ([../particle-nesting.md](../particle-nesting.md)); expanded
  before graphics sees it, the saving (copies + placements, not copies × placements) is lost. Required of meso\*
  (Asked of meso\* 1): the render IR keeps both forms, the expanded placement the runtimes read and, per instance, its
  definition and use, from which the adapter builds a template's list once and its copies' poses.
- **Bands are the spec's, evaluated by the renderer**: meso\*'s view language has bands and context rules
  ("colours that follow the zoom") evaluated on the CPU into a keyed diff of looks. A representation chosen by size
  on screen is lowered into this plan's band list (2), never switched by a context rule (that would be a commit, a
  rebuild and a pop); context rules drive looks alone (color, opacity, emissive). The two meanings of "band" are one
  concept: meso\*'s written, this plan's drawn. The spec's bands give the representations and their order, the renderer
  each boundary by error (1); a boundary the spec pins is kept (Asked of meso\* 2).
- **Looks that follow the zoom cost nothing a frame by instance**: a context rule evaluated on the CPU at each view
  change is a diff of looks every frame of a dolly; by population or member they are uniforms, by instance a layer's
  upload a frame. So zoom-following looks are by population or member (or by depth on the GPU, as a band's `t`), with
  a margin between switching on and off (a camera resting at a threshold flickers otherwise) (Asked of meso\* 3).
- **Bands by target, not by representation**: a model is one list with a layer per population, its targets of every
  kind and size; the bands (and their boundaries) are a table by target (member), part of the list's version.
- **Many lists**: the render IR has one list a placement set, a scene hundreds; each a merged object culled on its
  own ([many-culled-objects.md](many-culled-objects.md): passes by object were the cost) and budgeted on its own (11).
  Open: lists merged into one object a scene where they share targets, or a frame's budget across objects.
- **Motion is meso\*'s runtimes'**: animation engines (per copy: explode, spiral, tween, fade, pairs by source),
  physics and dynamics move every copy of a population each frame. What graphics gives them: live poses (7) and
  motion evaluated on the GPU (8), behind the engine port (a GPU engine one implementation, passing the engine's
  suite). meso-runtime imports nothing of Mol\*: the adapter lends it a buffer handle, never the reverse. A buffer is
  its device's: a GPU engine or world writes the renderer's poses only on the page, on the renderer's device the
  adapter lends it; in a worker, on a device of its own, its poses are read back and uploaded (8).
- **One writer a value, overrides included**: meso\* decides for a run who writes each value (physics leaves out what a
  track moves). An override (6) is a writer too: a copy it holds leaves the physics world or the track while it does,
  and rejoins when the override is dropped or committed (Asked of meso\* 4).
- **A new version keeps what did not change**: an edit of the copies is a list's new version, updated in place by the
  adapter. Today a flat list's new version resets the page store (every page written again, its residency and fade
  stamps lost: an edit or a freeze reloads and pops). A version changed keeps the pages whose particles did not change
  (their residency, their stamps), by a diff of the pages (as `ParticleSlots.getChanges` gives a slot list's), from
  the changed ranges the render IR says where it can (Asked of meso\* 5).
- **A frame at rest**: meso\* requests no frame at rest; this plan's readbacks (residency's wants, promotion's
  candidates) arrive a frame late, fades and band geometry over several. After the camera stops, frames until they
  settle (no wants unmet, no fade running, no band building), bounded, then none (4). The player's 0 B a frame
  holds for override writes, readbacks and motion uniforms too (no array or promise a frame): work, not a property
  today (the paged cull's read backs allocate their promises and arrays).
- **Identity**: loci by path sets (`[i0, i1]`, `[*, i1]`) round-trip with MesoQL locators and addresses (the
  adapter's map: a population's list and its template paths). Promotion needs a `Structure` from the molecular IR's
  handles on demand; a promoted representation is session state outside the spec, rebuilt from the view after a
  restore or an undo (the adapter keeps no record of what it drew).

## Design

1. **One metric: error on screen in CSS pixels.** Every level of every kind declares its geometric error in Å (a
   stride's spacing and scale, a mesh level's deviation, a copy's sample spacing, a tile aggregate's cell, a proxy's
   deviation from its target's shape); a level is drawn where its error is under `D` pixels, `D` the one quality knob
   (`lodPixelSize`, the graphics mode's). `LodView` converts error to distance ranges once a frame, as now. The
   object's `lod` range gains a pixel form (a target radius in px) and the Å form is kept for what is not instanced.
   Why: the four mechanisms then agree where they meet (a copy's coarse level and a sphere level switch at the same
   pixel size), one knob trades quality for speed everywhere, and an adaptive `D` (later) steers all of it.

2. **Bands: several representations of one target over one instance source.** A target representation takes an
   ordered list of bands by target (a table: a porin and a ribosome differ), each a geometry and its props and where
   it ends: e.g. atoms as impostors, then a Gaussian surface, then a proxy (the list and its order the spec's; a
   boundary it pins kept). A band ends where the next one's error is under `D` pixels (1), not at a fixed radius on
   screen: two targets of one radius can differ 100× in atoms, a
   surface worth drawing earlier for the larger. The cull decides an instance's band (by its depth), and writes it into
   that band's regions: one dispatch, one pose buffer, one set of
   markers, layers and pick ids by store index; each band's renderable (a geometry kind's, merged by target as now)
   reads the poses by index (`dCulledPaged`). Within a band its own levels (strides, mesh levels) as today. What a band
   adds: its geometry a target (built where an instance first draws it: the bands' geometry resident by use) and its
   draw records. Why one source: a representation per band copies the instances (64 B and a culler each), culls them
   again, and lets markers, layers and picks diverge between the bands. The CPU-culled path (WebGL2, few instances)
   decides a band per cell (the cell's nearest depth), the cell drawn by one band, faded by fragment.
   - **A band not ready** (its geometry building, its pipelines compiling) is drawn by the nearest ready coarser band,
     the proxy at least (built with the list: a few floats a target), never by nothing; faded in once ready (4).

3. **Exact crossfades by a shared dither.** An instance within a transition (its pixel size within a band's fade
   width, `overlap` in px) is written to both bands; each fragment keeps or discards by one threshold `ξ` (screen
   position, the frame's jitter index) against the instance's `t`: the finer side keeps `ξ < t`, the coarser
   `ξ ≥ t`, so every pixel is drawn by exactly one side (no hole, no double coverage), opaque, depth-tested, picked as
   drawn (the pick pass the same `ξ`: of the full-resolution pixel a pick pixel stands for, the pick buffer smaller,
   `PickScale`; else a pick lands on the side discarded). `t` from the instance's depth in the vertex stage (as
   `assign_lod_factor`), not stored. With temporal accumulation (the multi-sample offset, TAA) the dither resolves
   to a blend. Why stochastic: no sort, no second pass, both backends, and only the instances in a transition pay
   twice.
   - **Transparency** (meso\*'s opacity tracks, looks' transparency): an instance whose layer is transparent fades by
     alpha (the finer side's `a t`, the coarser's `a (1 - t) / (1 - a t)`: composited over each other exactly the
     layer's `a`, no dip mid-transition as `a t` and `a (1 - t)` would make, nor twice it; weighted blending
     approximates it), drawn in the transparent pass by both bands; opaque instances by the dither. An instance switches
     between the passes where its alpha reaches 1 (an opacity track's end): faded over the frames of (4), not at once.
   - **What reads dithered pixels**: SSAO, outlines, the illumination's input, depth of field and the marking's
     edges see the pattern of a still frame without accumulation; XR draws each eye with one `ξ` (the views' same
     pixels the same side). The selection outline's shape differs by band: it changes where the band does.
   - **Nested sphere levels fade without a second draw**: level `j` is a subset of level `j - 1`, so within the
     transition the finer level is drawn with the items not in `j` dithered out by `t` and every item's radius
     interpolated from `j - 1`'s scale to `j`'s (a geomorph). Reopens [screen-space-lod.md](screen-space-lod.md)'s "no
     fade" (Open).
   - **Copies' coarse levels and the residency fallback** by the same rule (a copy's level by its error, faded
     across the halving).

4. **Transitions in time for discrete events.** A page or tile becoming resident, a band's geometry built, a promotion
   (6) are events, not distances: what replaces a coarser draw fades in over a short time (~200 ms) by the same dither,
   `t` from the time since the page's or copy's first draw (a stamp a page and a copy in the cull's tables, a frame's
   clock uniform). Hysteresis where a decision is discrete (a band per cell on the CPU path, a promotion): a margin
   in pixels between switching up and down.
   - **Settling**: a frame is pending while a fade runs, a readback is outstanding (the residency's wants, the
     promotion's candidates) or a band builds; the canvas requests frames while any is, at most a bound (e.g. 2 s)
     past the last input, then none. A still view's last frame is its settled one (what an image export waits for, as
     [resource-readiness.md](resource-readiness.md)'s pending work). The bound only stops asking for frames: what waits
     for a final image waits on the pending work itself (remote reads and read backs take seconds: the particle
     budget's settled in four read backs, 3.7 s).

5. **Proxies: one primitive a particle.** The coarsest band of every target: a sphere (its radius the target's) or
   an ellipsoid impostor (its axes the target's gyration, oriented by the particle's rotation: a fiber, an antibody, a
   membrane protein keep their outline) or a splat (opaque, no sort: its footprint an ellipse). Its parameters a
   target (a few floats), its color by member (the target's, entity's or compartment's: no per-particle layer). Why:
   below a few pixels a target's atoms (hundreds to thousands of spheres) draw into one or two pixels; the proxy draws
   one quad. Replaces the particles' own spacefill for that purpose. Its data by regime: a flat list's pages, a
   copy's stored coarse levels (not its leaf pages: residency keeps a far copy's out), a tile set's coarse tiles.

6. **Instance overrides: a few instances posed, pinned or taken out.** A sparse table a representation, its entries
   keyed by an instance's address at a level, as the loci's path sets are: a store index (a flat list's), one instance
   in one copy (`[i0, i1]`), a template particle in every copy (`[*, i1]`); a whole copy (`[i0, *]`) is its pose in
   the store's copy table written (`setCopies`), no entry. An entry:
   - **a pose**: replacing the stored one, or composed after it (a delta: a wiggle, a nudge); a copy's applied where
     the cull reads the copy's pose (its template culled at it), a template particle's in every copy's expansion;
   - **flags**: `hidden` (no band draws it: the app draws something bespoke in its place), `pinned` (a band and level
     fixed: always atoms, the finest stride, whatever the distance), `noFade`.

   The cull skips an overridden instance where it is stored and draws it from the table: a bit a cull group (256
   particles, a page's group test already) says the group has entries, a binary search over the sorted entries only
   inside such groups; a kernel over the entries composes each (its override, or its stored pose and the delta),
   culls it by its own sphere and writes it into its band's (or its pinned band's) regions, or nowhere where hidden.
   Why from the table and not in place: a moved instance leaves its page's, group's and cell's bounds, so it would be
   culled with them; and a pose written into its page rewrites 8,192 particles and their bounds for one.
   - **Cost**: a frame's change writes the entries (`writeBuffer` of `K` × ~48 B), no page, no CPU work by instance
     beyond `K`; the object's bounds grown by the entries' spheres (the near and far planes, the whole object's cull).
     Identity unchanged: picking, markers, layers and loci by the store index as stored. Sized for up to about 10⁴-10⁵
     entries; past it a frame of poses (posed pages, a list's next version) is cheaper, and copies are few (a copy's
     override moves thousands of instances for one entry: an animation by copy, the player's per-copy motion).
   - **One writer**: a copy an override holds is out of the physics world and the tracks that move it while it does
     (meso\*'s one writer a value), drawn by the override alone; rejoined when it is dropped or committed.
   - **Transient, then committed**: an override is a look (a drag, a motion track's value, a bespoke placement),
     never written into the list; to keep it, the list's next version with the page posed again, and the entry
     dropped in the same commit (no frame drawn twice or not at all). As meso\*: transient values over committed
     lists. A whole list moved is [object-transforms.md](object-transforms.md)'s, not an entry each.
   - **Opting out of levels** by instance is `pinned` (the instance still culled by the frustum and occlusion) or
     `hidden`; by population a band list of one band and no levels (the representation's props, as now).
   - **WebGL2**: the entries a texture, culled on the CPU (they are few), drawn as an extra instance range of each
     band; the hidden by the group bit and the search in the vertex stage (collapsed, as a culled instance).

   **Promotion** builds on it: near the camera the few largest instances (by pixel size, at most `K`, e.g. 16-64)
   drawn by a full structure representation (cartoon, ball and stick, labels, interactions: what needs a `Structure`,
   made from the target on demand, [../mesoscale-scale.md](../mesoscale-scale.md)'s open part) placed by the
   instance's pose. Graphics provides:
   - **the candidates**: the cull's readback of the instances over a pixel threshold, nearest first, capped (as the
     residency's wants: one read at a time, a frame late), with their address and composed pose;
   - **the swap**: an entry `hidden`, the promoted representation fading in while the instance fades out (4), then
     the entry holds. Which instances to promote, the representation built and its budget are the app's (the
     adapter's).
   The same entries serve "show these instances differently" in general (a focus, a selection drawn as a structure,
   a particle dragged out of its packing), without a second instance list.

7. **Moving instances and spatial order.** The store's culling rests on space: pages tile-ordered, cull groups of
   256 particles with a bounding sphere each, cells, the residency's pages by their bounds. An instance that moves
   leaves them, and is culled wrongly with its group. Three regimes, by how far and how it moves:
   - **Bounded motion** (a wiggle, Brownian jitter about a site, a tumble, a membrane protein diffusing a few Å): its
     reach `ρ` known when the motion is (an amplitude, a step's bound over the run); the groups' (and cells', pages')
     spheres grown by `ρ` while it runs, the order unchanged. A uniform a motion, nothing stored, no bounds pass.
   - **Motion by a known law** (8: explode, spiral, a tween between two placements, a copy along a curve): the
     group's bounds by the same law where it is affine in the group (an explode about `c` by `s(t)`: the sphere's
     center moved by `s(t)` from `c`, its radius scaled), else the sweep's bounds over the track's range (a tween's:
     both ends' spheres and the path between), computed when the motion is set, not a frame.
   - **Free motion** (physics, dynamics, live input, a drag of many): the poses written each frame (a live pose
     buffer: the runtime's, lent to it by the adapter, or its upload of what moved), the groups' spheres recomputed on
     the GPU from them (a reduction a group: one pass over the moving poses, ~16 B an instance read), and the spatial
     order kept as an index, not as storage: a list of store indices sorted by Morton key (4 B an instance), the cull's
     groups runs of it; re-sorted on the GPU (the splat sort's counting sort, by 16-bit keys of a coarse grid) when the
     groups' spheres have grown past a share of their sorted size, not a frame. Poses, markers, layers and picks stay
     by store index: nothing that identifies an instance moves. Why an index: storage reordered would renumber the
     store indices every identity is keyed by (loci, marks, layers, overrides), a remap of each at every sort.
   - **Residency under motion**: pages that move freely are resident while in the motion (their bounds a frame's);
     a tile's particles moving out of the tile stay its (in the tile of their first frame,
     [../particle-tiles.md](../particle-tiles.md) 3), its reach grown by the motion's.
   - **What freely moving costs**: a frame's upload of the moving poses (28 B an instance as float32, 16 B as a
     quantized delta of its group: 400k moving 6-11 MB a frame), the groups' reduction, a sort now and then. Measured
     in step 9 on a packing relaxed live and a cellPAINT scene's physics.

8. **Motion on the GPU.** What a motion is written as (meso\*'s animation language: primitives and combinators run
   per copy), evaluated where the pose is composed (the cull's kernel; the vertex stage on WebGL2) at the player's
   time: a frame uploads `t` and a motion's few parameters, not a pose an instance. The precedent: the wiggle and tumble
   animations evaluated by `uTime` in the vertex stage (`animation.wgsl.ts`), nothing a frame from the CPU.
   - **What it evaluates from**: the stored pose (the copy's placement), the instance's index (a seed: a phase, a
     delay, a direction by a hash; nothing stored), the population's or group's parameters (a center, an axis, a
     curve's control points), and, where the motion needs it, a second placement (a tween's other end, pairs by
     source: a second pose an instance, stored once, uploaded once). Order: stored, then the motion at `t`; a copy an
     override holds drawn by the override alone (6: one writer).
   - **Primitives in the shader, combinators by a table**: each primitive (explode, spiral, tween, orbit, along a
     curve, wiggle, fade as `t` of (3)) a function in one shader chunk; a motion's combinators (sequence, blend, a
     delay by copy, a selection by entity or compartment) a small program in a storage buffer, interpreted. Why not a
     shader generated a motion: a pipeline a motion, compiled as it changes (a spec edit a stall,
     [pipeline-readiness.md](pipeline-readiness.md)); an interpreter's few ops cost less than the compile. Generated
     where the interpreter measures too slow.
   - **Keyed frames** (a simulation's trajectory, a physics run recorded): two frames resident, the pose between
     them interpolated in the kernel (position linearly, rotation by nlerp), the next frame uploaded once an interval
     and not a render frame: a 10 fps trajectory at 60 fps a sixth of the traffic, quantized (a delta of the group:
     16 B an instance) a further half. A frame is a store's pages posed ([../particle-slots.md](../particle-slots.md)
     5), the two a page each.
   - **The CPU agrees**: the same motion evaluated on the CPU for what reads one instance (a label's position, a
     focus, a pick's 3D point, the bounds of a promoted instance), and the engine port's suite run against the GPU
     engine (headless WebGPU, within a tolerance); a freeze reads the CPU's (or the GPU's poses back once).
   - **Bounds**: by the motion's law or its reach (7); a motion whose bounds cannot be written is free motion.
   - **Not GPU**: physics and dynamics stepped on the CPU (their worlds' state is the runtime's); their poses are
     free motion (7). A GPU world (the packing's relax, already WGSL kernels) writes the lent pose buffer itself, no
     upload at all, where it runs on the page on the renderer's device (lent by the adapter: a buffer is its device's);
     in a worker, on a device of its own, its poses are read back and uploaded as free motion.

9. **Mesh levels by index ranges.** A target surface's levels in one vertex and index buffer, a level an index range
   (its first index and count): the paged and merged cullers already write a level's draw count and first index from
   the member's `lodLevels`, so a mesh level is a level as a sphere's is. Made by kind:
   - density surfaces (blob, Gaussian, a volume's isosurface): marching cubes at coarser grids (the resolution halved
     a level: a level a quarter of the triangles; all levels 1.33× the finest), on the GPU where the surface is;
   - other meshes (shapes, imported meshes): simplified (edge collapse by quadric error) in `@molstar/core`
     (`math/geometry/`, with tests; a new API, asked first) in a worker, the error a level's metric (1).
   A mesh level's error its grid spacing or its quadric error. Why index ranges: one buffer a target, the existing
   draw records, no new renderable.

10. **Aggregates a region, far away.** Below the proxy (sub-pixel particles), what shows is a region's color and
   density: a tile's or a copy's aggregate, drawn in place of its particles. Candidates: the copy's coarsest template
   level as proxies (exists, extended to proxies), a tile's splats (a splat a block of its density grid, the tiles'
   opacity blocks already computed for visibility), a compartment's surface. One band of the region, not of a target;
   fades as (3). Its data by tile, made when the tile set is written (the pack) or once loaded.

11. **Budgets.** The composed slots and the copy blocks are bounded today; bands add: a band's geometry resident by
    use (let go after a cooldown, as pages), the promotions' `K`, the moving poses' upload a frame, and a frame budget
    of drawn primitives, past it the farthest drawn a band coarser (flagged, never a fault): a multiplier of `D` (1),
    the copies' particle budget its first form (What exists). The budget is an
    object's: a scene of many lists (Fit with meso\*) trades none against another (Open). An adaptive `D` by frame
    time is later (not planned yet).

## Memory

What a band, a level or a fade adds, by design: no instance data (the store and its markers and layers shared); a
band's geometry a target (the proxy a few floats; a surface's levels 1.33× its finest; the atoms as now); a fade
nothing stored (`t` from depth or a stamp a page); overrides `K` entries and a bit a cull group (the yeast's 31.5 M
particles 15 KB); a tile's aggregate a block of its density grid. Motion: by law nothing an instance (a seed is its
index), a tween a second placement, keyed frames two pages a page, free motion its spatial index (4 B an instance)
and the moving poses. What the scene holds stays copies + placements ([../particle-nesting.md](../particle-nesting.md)) and
pages resident by view ([../particle-residency.md](../particle-residency.md)).

## Steps

Each measured on the PetWorld postsyn (flat), the isg 16³ (nested), the yeast cell (tiled and nested 2×2): frame time
zoomed, whole and in a dolly, memory kept and peak, the popping metric, images by eye.

1. **The baseline**: the dolly and pan sequences (a benchmark case, headless and in Chrome), the popping metric
   (frame against the previous reprojected), the drawn primitives a pixel by kind; for every scene at `D` 3. Says
   where popping is worst (expected: sphere strides, copies' coarse levels, page loads) and what each later step buys.
2. **One metric**: a level's error (`LodLevelsValue`), the copies' coarse levels by it (not their own 4 px), the
   object's `lod` in pixels; images unchanged at the default (the same boundaries, re-expressed).
3. **The shared dither**: one chunk on both backends (GLSL, WGSL), exact complements, the pick pass's; the object's
   `lod` fade moved to it (two representations' overlap without holes); then the nested sphere levels' fade
   (geomorph and dither, no second draw), measured: its cost in the transition and the popping metric; the
   transparent instances' alpha fade; SSAO, outlines and XR over dithered pixels by eye.
4. **Proxies**: the sphere and ellipsoid proxy (a target's gyration), color by member; the coarsest band of the
   particle representations' default, before bands (what draws a band not ready); measured against the strides'
   coarsest level (primitives, frame, images).
5. **Bands in the cull**: the band table by target; the paged cull's band per instance and regions per band;
   renderables of several geometry kinds over one `CulledData`; markers, layers, picking across bands; the CPU path
   by cell; a band not ready drawn by the next ready coarser. First bands: atoms (impostors) and the existing
   surfaces; then the overlap (3).
6. **Time fades and settling**: the stamps a page and a copy, the clock uniform; residency loads and the coarse
   fallback faded; the pending state (fades, readbacks, builds) and the bounded frames after input.
7. **Mesh levels**: index-range levels in the merged and paged draw records; density surfaces at coarser grids; then
   the simplifier in core (asked first).
8. **Instance overrides**: (a) the table by address level (store index, template particle; a copy its pose), the group
   bits, the override kernel into the bands' regions, `hidden` and `pinned`, the bounds grown; tests: a moved instance
   drawn outside its page's bounds, a hidden one in no band, picked and marked by its store index; measured: a frame of
   1, 10³ and 10⁵ entries moving against a page posed again. (b) The commit (the list's version and the entry dropped
   together). (c) WebGL2. (d) Promotion: the candidates' readback, the swap's fade; an example in the plugin (the
   nearest particles as cartoons, a particle dragged), the adapter's part described for meso\*.
9. **Moving instances** (7): (a) bounded motion's reach growing the groups', cells' and pages' spheres (the
   existing wiggle and tumble first: today their bounds are not grown); (b) the live pose buffer, its upload of what
   moved, the groups' spheres reduced on the GPU; (c) the spatial index (Morton-sorted store indices, the groups its
   runs) and its re-sort by the groups' growth; tests: an instance moved across the scene drawn and picked, never
   culled with its old group; measured: a packing relaxed live, a cellPAINT scene's physics (the upload, the
   reduction, the sort's frequency and cost, the cull's efficiency against the order's age).
10. **Motion on the GPU** (8): (a) the primitives' chunk and the combinators' interpreter in the cull's compose (and
    the vertex stage on WebGL2), `t` from the player; (b) bounds by law and by sweep; (c) the CPU's evaluation and the
    engine port's suite against the GPU engine; (d) keyed frames, two resident, interpolated; measured: the explode
    and spiral of the isg grid's copies and a 10 fps trajectory at 60 fps against posing the pages each frame
    (upload bytes a frame, frame time). The engine port's GPU implementation described for meso\* with the adapter.
11. **Aggregates**: a tile's or copy's far band (proxies of the coarsest level, then tile splats), the yeast's whole
    view against today's.
12. **Budgets**: band geometry residency and cooldown, the primitives budget, the moving poses' upload; adaptive `D`
    left to an adaptive-quality plan.

## Asked of meso\*

What this plan needs meso\* to change or decide, required where the design rests on it:

1. **The render IR keeps what repeats nested** (required): `resolveDefinables` stays the one form the runtimes read
   (packing, physics, the freeze); beside it, each instance's definition and use (its path in them), so the adapter
   builds a definition used many times as a template once and its copies' poses. Without it every meso\* scene draws
   copies × placements (the yeast cell 2x2: a pack of 714 MB nested, 2.8 GB flat).
2. **Bands give representations and order, not boundaries**: the spec's band list says what represents a target and
   in what order; the renderer sets each boundary by error (Design 1); a boundary the spec pins is kept.
   `meso-convert`'s default bands written without boundaries.
3. **Zoom-following looks by population or member**: a context rule's looks that follow the zoom are per population
   or member (uniforms), not per copy (a layer's upload every frame of a dolly), or written as a band's look the GPU
   evaluates by depth; their thresholds with a margin between on and off.
4. **One writer includes overrides**: the engine port and the worlds release a copy an override holds (a held set the
   adapter gives them) and take it back when it is dropped or committed, as they leave out what a track moves.
5. **A version says what changed**: a placement set's new version with the ranges of copies it changed (or a stable
   order of the unchanged copies), so the store keeps the pages that did not change (their residency and fades)
   instead of a reset at every edit or freeze.
6. **A GPU world on the renderer's device**: the engine port able to take a lent `GPUDevice` (a world on the page
   writing the renderer's pose buffer: no upload); a world in a worker hands its poses over (read back, uploaded).
7. **A frame at rest after settling**: meso\* requests no frame at rest; the renderer asks for frames after input
   until its fades, read backs and builds settle (bounded, Design 4), and an export waits on that pending work.

## Open

- **Fades between sphere levels**: [screen-space-lod.md](screen-space-lod.md) decided none (no chain needed one).
  Smoothness is now a goal; step 1's popping metric decides whether its strides pop enough to need (3)'s geomorph.
- **The band's boundary** (proposed: by error, 2): each band's error needs every geometry's (a surface's grid, a
  proxy's deviation from its target); the radius on screen is simpler but blind to a target's cost. Whether a band
  should also weigh cost (atoms a pixel), not error alone.
- **Repetition the spec does not say** (the nested form itself: Asked of meso\* 1): who detects it in a flat
  placement (an assembly's operators, a grid), the compiler or the adapter.
- **Many lists**: one merged object a scene's lists of shared targets, or a frame budget across objects
  ([many-culled-objects.md](many-culled-objects.md)'s pool, measured and removed for its coupling); by the number of
  lists meso\* scenes have.
- **The spatial index's cost**: the indirection read a particle a cull (4 B, a gather where the store was read in
  order) against a cull's efficiency kept under motion; whether a static list keeps its stored order (no index) and
  only lists that move freely get one.
- **The motion interpreter**: its op set against meso\*'s animation language (every primitive and combinator, or a
  subset on the GPU and the rest on the CPU as free motion); its cost a particle a cull against a generated shader.
- **The live pose buffer's owner**: a GPU world lent the store's buffer writes the renderer's memory; its layout
  (float32 or the store's quantized form) the store's, so a runtime behind the port depends on it.
- **The dither**: the 4×4 ordered matrix (exists, visible as a pattern in a still frame) or blue noise and the frame's
  jitter (needs the accumulation to resolve it; a still frame without it noisier). SSAO and outlines over dithered
  pixels: measured in step 3.
- **Bands as renderables or one renderable**: proposed separate renderables (a geometry kind each) over one cull
  and one pose buffer; one renderable of mixed geometry would need a program per band anyway.
- **The proxy's shape**: sphere, ellipsoid or opaque splat; an ellipsoid's cost against a sphere's in the raster.
- **Overrides' lookup**: the group bits and a search (proposed), a bit a store index (dense: 4 MB for the yeast,
  no search), or a hash; by step 8's measurement of the cull's cost with no entries (must be none) and with many.
- **Overrides' reach**: an entry a template particle in every copy (`[*, i1]`) needs the group bits in the
  template's store, read in every copy's expansion; whether that is wanted beyond hiding (a pose for every copy of
  one particle) or left to editing the template.
- **Promotion's representation**: built by the adapter as a structure representation (one per promoted instance) or
  a merged object of the promoted (a few at a time: the former, simpler).
- **Aggregates' form**: tile splats, a density volume (direct volume a tile), or the coarsest copy level only.
- **Fibers and curves** (copies along curves, ropes: `fibers.ts`): bands and levels for them, or drawn as now.
- **WebGL1**: bands by object distance (no per-instance band), no time fades; or not supported for bands.
