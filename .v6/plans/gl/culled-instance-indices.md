# Plan: GPU-culled regions as instance indices

## Status

- **Built** (steps 1 to 4), measured (see Measured): the cost of the reads per vertex is small, the opaque pass the same
  or less but where impostors fill a moving view. Follows `cull-workgroup-limit.md` (its next limit).

## Question

The merged culler copied each drawn instance's attributes into its regions: a region per level, the raster's regions,
the ghost's, then the identity region, each sized for every instance. Per instance and region 68 bytes in three buffers
(`aTransform` 48, `aInstance` 4, `aSegment` 16), so `aTransform` alone was `(regions + 1) × instances × 48` bytes and
passed WebGPU's 2 GiB a buffer at about 4.5 M instances with ten regions (the isg model 3×3×2: 2.98 GB, refused).
Mesoscale models are wanted at this scale and past it.

## Design

Regions hold a 4-byte index into the instances' data, stored once, instead of copies of the attributes.

- **Memory.** Per region and instance 4 bytes instead of 68 (12× for the transforms' buffer, the one that overflowed;
  17× in all). The instances' data once: 56 bytes an instance (68 before the segments as indices,
  `mesoscale-scale.md`). Ten regions, 6.2 M instances: 4.2 GB in three buffers (2.98 GB the largest) -> 0.35 GB of data
  and 0.25 GB of indices. The data reaches 2 GiB at about 38 M instances, the indices at 50 M with ten regions.
- **The data** (`CulledData`, `culled-data.ts`): one buffer (one binding: the vertex stage's storage buffers are few, the
  sphere resolve has three and the theme layers up to four pages) in blocks: the transforms' affine rows (three vec4s
  each, from 0), `aInstance` (four a vec4), each instance's segment index (four a vec4), each segment's groups (a vec4:
  `aSegment` composed from those two, `culledSegment`). The transforms first: the sphere raster and its resolve read
  them as they read an unmerged object's (stride three), unchanged.
- **The culled program variant** (`dCulled`, `CulledSchema` in `merged.ts`, WebGPU): the merged schema's instanced
  attributes are used `!dCulled`; the culled programs read the index from storage at `instance_index` (which includes
  the record's first instance; a `u32` vertex input would need an attribute kind of its own, WebGL's too), then the
  instance's data (`sCulledData`, `uCulledInstanceCount` for the blocks), into the program's `VertexInput` by a
  generated reader; a generated `vertexMain` calls the program's own, renamed (`webgpu/culled.ts`). The program sources
  stay as they are. As the dropped version did (`3f3081d93`), with three rows (not seven vec4s) in one buffer. A culled
  item binds no instanced vertex buffers.
- **The culler** writes the data where it wrote the identity region (`syncInstances`), the identity region's indices
  `0..n-1` with the layout; the copy kernel writes `indices[slot] = i` where it copied the attributes, reading each
  instance's sphere from the data; the records are unchanged. The copy kernel binds 6 storage buffers (8 before).
- **The sphere raster**: the regions kernel reads a slot's index and appends the chunk with it (the chunk's slot is
  the instance's index into the data), so the raster kernel and the resolve read the transforms as before. The ghost's
  region as the levels'.
- **Precompile**: `dCulled` is not a signature define: a program item sets it from `gpuCulled` and whether the context
  culls on the GPU (`canCullOnGpu`). `gpuCulled` is part of merged signatures' key (culled and unculled objects draw with
  different programs; merging them lost the unculled one).

## Measured

WebGPU, NVIDIA (Lovelace), rotating unless said; frame median and the opaque pass's GPU time, copies -> indices:

| Case                                         | Frame             | Opaque pass       |
| -------------------------------------------- | ----------------- | ----------------- |
| `merging/instances-1000-spacefill-merged-gpu` | 0.85 -> 0.82 ms   | 0.02 -> 0.02 ms   |
| `merging/instances-1000-cartoon-merged-gpu`   | 1.23 -> 1.08 ms   | 0.87 -> 0.73 ms   |
| `mesoscale-split/postsyn-whole-noraster`      | 3.42 -> 3.40 ms   | 1.21 -> 1.20 ms   |
| `mesoscale-split/postsyn-whole-auto`          | 4.76 -> 4.34 ms   | 0.86 -> 0.56 ms   |
| `mesoscale-split/postsyn-zoom0.3-noraster`    | 8.47 -> 8.36 ms   | 2.72 -> 3.08 ms   |
| the same, still                               | 5.70 -> 5.80 ms   | 2.47 -> 2.47 ms   |

The per-vertex reads cost where impostors fill the view moving (zoomed PetWorld without the raster: the opaque pass
+13%, the frame CPU-bound and the same); elsewhere the same or less (meshes 16% less: no instanced vertex buffers, the
copy kernel's bandwidth gone; with the raster its culling 0.25 -> 0.15 ms). The merged spacefill's storage 32.2 ->
20.4 MB at 100k instances.

## Cost

The vertex stage reads its instance's index and data from storage per vertex. Measured when culling was built
(`merged-render-objects.md`, `merging/instances-1000-spacefill-*`, 100k instances, WebGPU): the opaque pass 0.65 ms
read through indices against 0.55 ms with copied attributes (and 0.58 ms culled on the CPU). That version read seven
vec4s; this one reads the index and five values. The copy kernel writes 4 bytes an instance and region instead of 68.

## Steps

1. Done: measured against the copies (Measured).
2. Done: the culled variant (schema, reader, wrapper, precompile).
3. Done: the culler's data and index buffers, the copy kernel writing indices, the identity region's indices.
4. Done: the sphere raster through the indices.
5. Tests pass: `webgpu-gpu-culling` (the GPU-culled image as the CPU's), the raster's native tests, precompile, the
   native WebGPU suite. Open: the WebGPU visual cases of merging, the sphere raster and the ghost; the isg 3×3×2 model
   in the viewer.
