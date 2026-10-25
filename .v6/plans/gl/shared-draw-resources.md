# Plan: Shared draw resources on WebGPU

## Status

- Steps 1 (geometry in shared GPU pages), 3.1–3.2 (theme layers in storage pages, updated by range) and 4 (shared
  object groups) done; step 2 (the geometry store) deferred: measured, it gains little (below). Layer B of
  [webgpu-draw-submission.md](webgpu-draw-submission.md), ordered by what the measurements below say it gains, with a
  store that geometry builders, theme layer producers and the GPU share (no extra copies, grouped by representation and
  material, updated by range); [draw-table.md](draw-table.md)'s row buffer (step 6) has no consumer since GPU culling
  went per merged object.

## Why

Replaying a render bundle costs Dawn (the GPU process: commands translated into the D3D12 command list) about 0.4–0.5
µs per draw for each buffer or bind group the draw binds of its own, but little per command. Measured with
`draw-overhead/*` (raw WebGPU, 10k tiny cubes in one bundle, NVIDIA, D3D12), per frame:

| How the draws bind                                                         | ms   |
| -------------------------------------------------------------------------- | ---- |
| as Mol\* does: two bind groups and six buffers (5 vertex, index) their own | 37   |
| only the buffers their own                                                 | 24   |
| only the bind groups their own                                             | 10   |
| shared buffers, bound per draw at an offset; shared bind groups            | 2.2  |
| shared bind group with a dynamic offset per draw; shared buffers           | 1.1  |
| shared bindings set once, objects by `firstInstance`                       | 1.0  |
| one instanced draw                                                         | 0.2  |

In Mol\*, `draw-table/boxes-10k-drawlists` waits about 45–60 ms per frame for the replay of its 10k draws, the CPU
submitting in 0.15 ms (still) to 2.3 ms (rotating). Per draw a render item binds the global group (shared), its
material group (shared per material without material textures), its object group (its own: the object block's slot,
the draw block's slot, its data textures and storage buffers) and its vertex and index buffers (its own).

## Design

### Geometry in shared buffers (first)

- Vertex attributes and indices of render items are sub-allocated in shared GPU buffers of the context (pages): an
  item binds its page with its range's offset and size (`setVertexBuffer(slot, page, offset, size)`,
  `setIndexBuffer(page, format, offset, size)`): no change to the shaders, the index data or the draw calls. A page is
  one buffer with vertex, index (and later storage) usage: offsets aligned to 4 bytes (256 for storage bindings).
- Updates write into the range (`writeBuffer`); a grown attribute gets a new range and frees the old one. A recording
  binding the old range differs in what it binds (`GraphicsContext.captureDraws` records buffers with offsets) and
  records again.
- Instanced attributes (transforms, instance ids) the same. Buffers produced on the GPU (texture meshes in the chunked
  vertex arena, indirect draws, storage-pulled geometry) stay as they are.
- Expected: of the 37 ms above, the 24 of the buffers go to about 1.

### A geometry store: builders and GPU sharing memory

Today a mesh reaches the GPU through up to three copies: the builder writes into chunked arrays sized from an estimate
(`MeshBuilder.createState`), `ChunkedArray.compactTrimmed` copies them into one array unless a single chunk is mostly
full, the WebGPU render item interleaves attributes only beyond `maxVertexBuffers` (rare), and uploads each attribute
into a buffer of its own. The geometry stays on the CPU as well (bounding spheres, picking positions, export, CPU
smoothing, themes reading positions).

- **Pages in CPU and GPU memory alike**: the store holds pages, each an `ArrayBuffer` with a GPU buffer of the same
  size (WebGPU), allocated per owner; ranges are typed-array views into a page (`subarray`), the GPU buffer bound at the
  same offsets. The CPU side is the geometry: what code reads today keeps working on views.
- **Builders write into the store**: a builder allocates its ranges (vertices, normals, groups, indices) from the store
  for its estimate and writes in place; past the estimate the range grows (a new range, the written part copied once,
  as compaction copies today); at the end the unused tail is freed instead of copied away. One write by the builder,
  one upload of the written range, no compaction copy.
- **Owners: the representation and its material**: a representation's visuals allocate from the pages of their
  representation (the material id groups them today: draws are sorted by program and material), so that its geometry
  is contiguous (few, large uploads; consecutive draws bind the same page) and freed together when it goes. A page of
  an owner is reused by its own updates (trajectory frames, theme or parameter changes rebuild into freed ranges).
- **Uploads by range**: writes mark their range dirty; before drawing, the store uploads each page's dirty ranges
  (merged when adjacent) with `writeBuffer` from the page's `ArrayBuffer`: the only copy, as today's per-item upload.
  Updates in place (positions of a trajectory frame) upload just that range.
- **Lifetime**: ranges are freed explicitly with their geometry (visuals destroying their render objects, geometry
  replaced by an update); a `FinalizationRegistry` frees what was missed and counts it (a leak to fix). A representation
  removed frees its pages.
- **Context loss**: the CPU pages stay; a new context creates the GPU pages and uploads them whole.
- **WebGL**: the CPU side is the same (views into pages); it uploads per item from the views as today, or binds shared
  buffers with attribute offsets (`vertexAttribPointer`) later. Nothing it needs is WebGPU-only.
- **Merged renderables** (`renderMerged`) copy their members' geometry into merged values today; members contiguous in
  one page could be drawn as ranges of it instead (multi-draw lists), without that copy: later, measured.
- **Theme layers** live in the same store (below).

### Theme layers in the store

The theme layers (color, size, marker, overpaint, transparency, emissive, substance, clipping, wiggle, palette, the
sphere positions) are data textures of each render object today: their producers (`createColors`, `createOverpaint`,
`createMarkers`, ...) allocate a `TextureImage` sized to a 2D texture layout, and every change uploads the whole
texture (`ValueCell.update`), also a hover changing a few markers of a large structure.

- **`DataSpec` values as store ranges**: [storage-buffers.md](storage-buffers.md) makes them flat typed arrays with an
  element count (`DataSpec`, read through accessors); here their arrays are views into the owner's pages, next to its
  geometry. The producers write into them in place, as they write into the texture arrays today.
- **Updates by range**: changes mark the ranges they touch (`lociApply` hands out intervals: a marker or overpaint of
  some groups marks those), and the store uploads only them. A hover over a large structure uploads the marked bytes,
  not the layer. Marks coalesce: ranges closer than a gap merge, and past a few ranges or a low fill (every other group
  changed) the layer uploads whole in one write, never many short ones.
- **Layers present or not**: a layer that no visual of an object uses (no overpaint) takes no range; a cleared layer
  (`KnownClearLayers` today) can free its range and bind a shared zero range.
- **WebGPU**: the page is bound as one storage buffer, the object block holds the sections' offsets (u32 words, packed
  as [storage-buffers.md](storage-buffers.md) says); the object group is shared per program and page (below).
- **WebGL**: the data arrays in ranges of data textures shared per context and format, as WebGPU's pages
  (`webgl/data-pages.ts`, see [storage-buffers.md](storage-buffers.md)), their offsets object uniforms, the changed
  elements written from views of the arrays (`texSubImage2D` by row); draws of items in the same pages skip binding them
  (WebGL2 1000 unmerged units 3.7× faster per frame).
- **Large layers**: an item's layers together beyond the binding limit (default 128 MiB) are spread over up to 4
  bound pages, each within it ([storage-buffers.md](storage-buffers.md) step 5); a single layer beyond it is not drawn
  (warned).
- **Merged render objects** write their members' layers by part into the merged item's range, without a merged copy
  ([merged-render-objects.md](merged-render-objects.md) step 4); members have no GPU layers of their own, so the part
  writes are the only upload, nothing left to share. Measured: in the merged recolor of `assembly-60-spacefill` the
  writes (composing words across parts) are 3% of the CPU, the color theme most of the rest.

### Shared object groups (then)

- The object group was an item's own because of its data textures, storage buffers and uniform arena (a buffer per
  item). With the theme layers in the store's pages (above) and the object and draw blocks in shared uniform pages, the
  object group is shared per layout and bound buffers, a draw sets its dynamic offsets (the object and draw blocks'
  slots).
- Shaders read the data arrays from storage with the section offsets instead of from data textures (the
  [storage-buffers.md](storage-buffers.md) design). Textures that stay textures (grids, images, fonts, GPU-produced
  data) keep groups of their items for those programs only.
- Expected: the remaining 10 ms of bind groups to about 1.

### Not needed for submission

- Object records indexed by `firstInstance` (the original layer B): dynamic offsets per draw cost about the same (1.1
  against 1.0 ms). Records come back if GPU culling (draw-table step 7) wants them.

## Steps

1. **Geometry in shared buffers** (WebGPU, done): the render items' vertex and index data in ranges of shared pages
   (`WebGPUDrawContext.geometryPages`, `buffer-pages.ts`: pages of 32 MiB with vertex and index usage, ranges of 16-byte
   alignment first fit, larger than a quarter page a page of their own; a freed range reusable after the next
   submission, empty pages but one dropped), bound at offsets; the pass state skips rebinding the same buffer and
   offset. Measured back to back (per frame): `draw-table/boxes-10k-drawlists` 45 to 17 ms (still and rotating, the
   replay's wait), `boxes-10k-rotating` direct 76 to 68 ms, `units-1356-gaussian-drawlists-rotating` 2.8 to 1.1 ms,
   `render/boxes-1000` 6.3 to 5.1 ms. Byte counts are the ranges' (attribute, elements); the pages' free space is not
   counted yet.
2. **Geometry store** (deferred): measured (WebGPU, `render/*` and `merging/*`), the builders' copies are small: of
   300k spacefill meshes' 378 ms build the compaction is at most 3 ms (the estimates fit), the upload 133 ms (stays);
   chunks allocated against final bytes 625 / 625 MB (300k spacefill meshes), 19.0 / 18.7 (7qpd cartoon), 60.5 / 54.8
   (1000 cartoon units, 22 MB copied by compaction, transient). The CPU keeps the geometry the GPU has (625 MB each),
   and the store would keep it too. The one large waste, merged objects' copy of their members' geometry (+52 MB at
   1000 cartoon units), went without a store ([merged-render-objects.md](merged-render-objects.md) step 3). Taken up
   again where a reader of the CPU geometry can go (then GPU-only data, `mappedAtCreation`), or for the fragmentation
   of the pages. The plan was: CPU pages behind the GPU pages, owners by representation and material, ranges as typed-array
   views; the mesh builder writing into it first (`MeshBuilder`: estimate, growth, trimming), the others after
   (spheres, cylinders, lines, points, text, image); explicit frees with their geometry. Measured: copies and peak
   memory of building (`render/*`, the structure representation cases of `pnpm benchmark`), upload sizes.
3. **Theme layers in the store**, in parts, WebGPU first:
   1. (done) **Layers in storage pages**: the WebGPU backend reads `tColor`, `tSize`, `tMarker`, `tOverpaint`,
      `tTransparency`, `tEmissive`, `tSubstance`, `tClipping` and `tWiggle` from shared storage pages
      (`WebGPUDrawContext.dataPages`, `data-layers.ts`) instead of data textures: an item's layers in one range, their
      word offsets object uniforms of a schema the backend extends (`getDataLayerSchema`), WGSL accessors
      (`read_tColor(i)`: `rgba8` one word per element, `r8` four) generated for the referenced layers of that schema.
      The layers are `DataSpec`s now ([storage-buffers.md](storage-buffers.md) steps 1 and 2: flat arrays, accessors
      in GLSL too, WebGL laying them out); direct volumes read the pages as well; the palette stays a texture
      (sampled). Measured: `boxes-10k-drawlists` 17.6 to 13.4 ms
      per frame, `units-1356-gaussian-drawlists-rotating` 1.25 to 0.94 ms. Layers together beyond the binding limit
      read from up to 4 pages ([storage-buffers.md](storage-buffers.md) step 5).
   2. (done) **Updates by range**: producers log the element ranges an update changed together with the update
      (`TextureChanges.update`, a log of 8 updates per image by the version each produced; an update without it, a
      replaced image or a cleared layer upload all); consumers merge the ranges since their last upload (`get`: gaps
      under 256 elements merged, more than 8 ranges or half the layer upload all, so every other group changed is one
      write), WebGPU writing them (whole words) into the layer's range, WebGL the rows they touch (`Texture.loadData`).
      Logged by the marker (hover) and by overpaint, transparency, emissive, substance and clipping applied without
      clearing. A change of the image without an update of its cell (the marker while a status uniform stands for it)
      is pending (`defer`), part of the next update's ranges: else a mark cleared to an all-clear status stayed in the
      uploaded layer and showed again with the next partial mark (`mark-uploads`). Checked on Dawn: single-atom hovers uploaded by range draw what a whole upload draws (the test fails
      with the range writes left out).
   3. **The store's CPU pages**: producers writing into store views, on the `DataSpec`s of
      [storage-buffers.md](storage-buffers.md) (done), with step 2.
4. **Shared object groups** (done): the object and draw blocks of render items in ranges of shared uniform pages
   (`WebGPUDrawContext.uniformPages`, 8 MiB, the dynamic offset alignment; freed with the item), and object groups
   without textures shared by the items binding the same buffers (`getSharedGroup` in `render-item.ts`: a trie of
   weak maps by layout, uniform, draw and storage buffers). Items with object textures (grids, images, text, the
   palette, GPU-produced data) keep groups of their own. Measured back to back (per frame):
   `draw-table/boxes-10k-drawlists` 11.1 to 4.75 ms, `-drawlists-sweep` 18.3 to 8.9 ms, the `units-1356-gaussian`
   cases about equal (1356 objects).
5. Then [draw-table.md](draw-table.md) steps 6 and 7 (rows on the GPU, GPU culling and indirect draws): superseded,
   GPU culling per merged object ([merged-render-objects.md](merged-render-objects.md) step 6).

## Open questions

- Page size: large pages bind once for many items but grow in big steps; small pages split items across more bind
  points (a draw binding another page costs like a buffer of its own).
- Fragmentation of the pages where geometry is updated often (trajectories, surfaces recomputed): compaction, or
  pages per size class.
- Whether sorting draws by page (as by program and material today) is needed once pages are few.
- Growth past the estimate: a new range and one copy (as compaction today), or chained ranges drawn as several ranges
  of one item (multi-draw lists, no copy, indices relative to each range).
- How a geometry learns its owner: passed by the visual (its representation's store and material) to the builder;
  geometry built outside representations (helpers, shapes, tests) gets a default owner.
- `mappedAtCreation` for GPU-only data (never read on the CPU): building straight into a mapped buffer saves the CPU
  page; not for geometry the CPU reads.
