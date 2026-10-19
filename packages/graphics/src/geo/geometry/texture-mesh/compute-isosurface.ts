/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Sphere3D } from '@molstar/core/math/geometry';
import type { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import type { ComputeContext } from '@molstar/graphics/gl/compute';
import {
  type ChunkedIsosurfaceResult,
  type ChunkedSurface,
  computeChunkedIsosurface,
  getChunkByteCount,
  getVertexArena,
  growVertexArena,
} from '@molstar/graphics/gl/compute/marching-cubes/chunks';
import type { IsosurfaceInput } from '@molstar/graphics/gl/compute/marching-cubes/compute';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { TextureMesh } from './texture-mesh.js';
import { isTimingMode } from '@molstar/core/util/debug';

/**
 * Where an extraction writes: pages of the vertex arena, `surface` reused; `estimate`, the expected vertex count, grows
 * the arena beforehand when needed
 */
export type IsosurfaceTarget = { readonly surface?: ChunkedSurface; readonly estimate?: number };

/**
 * Extracts an isosurface into `target` (see `extractIsosurface`). Called again when the surface did not get all its
 * pages, so it must not rely on inputs that may have changed in between (e.g. shared density textures).
 */
export type IsosurfaceExtraction = (target: IsosurfaceTarget) => Promise<ChunkedIsosurfaceResult>;

/** `computeChunkedIsosurface` into `target` */
export function extractIsosurface(
  graphics: GraphicsContext & { compute: ComputeContext },
  input: IsosurfaceInput,
  gridDim: Vec3,
  isoValue: number,
  transform: Mat4,
  invert: boolean,
  target: IsosurfaceTarget,
): Promise<ChunkedIsosurfaceResult> {
  const { compute } = graphics;
  return computeChunkedIsosurface(
    compute,
    input,
    gridDim,
    isoValue,
    transform,
    invert,
    target.surface,
    target.estimate,
    isTimingMode ? graphics.timer : undefined,
  );
}

/** Headroom of the expected count over the last vertex count, so small changes fit */
const Headroom = 1.25;
const MinCapacity = 3 * 1024;

/** The last vertex count of a texture mesh, read back after its extraction */
const lastCounts = new WeakMap<TextureMesh, number>();
/** The latest extraction into a texture mesh, so that the read back of an earlier one isn't handled */
const latestExtractions = new WeakMap<TextureMesh, number>();
let nextExtraction = 0;

/** Whole triangles */
const toCapacity = (vertexCount: number) => Math.max(MinCapacity, Math.ceil(vertexCount / 3) * 3);

/**
 * The isosurface in pages of the vertex arena, drawn indirectly with the count on the GPU: nothing is read back before
 * drawing. The count, read back afterwards, is drawn directly from then on (`TextureMesh.setDrawnVertexCount`), and keeps
 * the upper bound of the texture mesh (`vertexCount`) and the position iterator's data right; a surface that did not get all its pages (the arena's free pages ran out) is extracted again
 * after the arena grew. The re-extraction handlers of the mesh (`TextureMesh.setReextractionHandler`) run after either
 * changed the vertices data derived from them saw. The expected count, which grows the arena beforehand, is the last
 * count of `textureMesh` (with headroom), else `estimate` (a function: called only then, e.g. a pass over the grid,
 * `estimateIsosurfaceVertexCount`).
 */
export async function computeTextureMeshIsosurface(
  graphics: GraphicsContext & { compute: ComputeContext },
  extract: IsosurfaceExtraction,
  groupCount: number,
  boundingSphere: Sphere3D,
  textureMesh?: TextureMesh,
  estimate?: number | (() => number),
): Promise<TextureMesh> {
  const { compute } = graphics;
  const surface = textureMesh?.doubleBuffer.getChunks(graphics);
  const last = textureMesh && lastCounts.get(textureMesh);
  const expected = last !== undefined ? last * Headroom : typeof estimate === 'function' ? estimate() : estimate;
  const result = await extract({ surface, estimate: expected });
  const mesh = TextureMesh.createChunked(
    toCapacity(expected ?? 0),
    groupCount,
    result.surface,
    compute,
    boundingSphere,
    textureMesh,
  );
  mesh.meta.readVertexCount = undefined;
  const extraction = nextExtraction++;
  latestExtractions.set(mesh, extraction);

  // a newer extraction, or the mesh destroyed (which returns its pages)
  const isCurrent = () =>
    latestExtractions.get(mesh) === extraction && mesh.doubleBuffer.getChunks(graphics) === result.surface;
  const checked = result.count.then(async (state) => {
    if (!isCurrent()) return;
    let { count, drawn, pages } = state;
    let reextracted = false;
    if (state.missingPages > 0) {
      if (isTimingMode) {
        console.log(`Isosurface extracted again: ${state.missingPages} pages missing, the vertex arena grown first`);
      }
      growVertexArena(compute, getVertexArena(compute), state.missingPages);
      const again = await extract({ surface: result.surface });
      // the surface's table grows with the arena
      if (isCurrent()) TextureMesh.bindChunkedTable(mesh, again.surface);
      ({ count, drawn, pages } = await again.count);
      if (!isCurrent()) return;
      reextracted = true;
    }
    lastCounts.set(mesh, count);
    mesh.meta.readVertexCount = drawn;
    mesh.meta.chunkByteCount = getChunkByteCount(pages) + result.surface.table.byteLength;
    // the upper bound, for draws without the GPU count and the samples of data derived from the vertices (e.g.
    // smoothing, which then left out the vertices beyond it)
    const beyondBound = drawn > mesh.vertexCount;
    if (beyondBound) TextureMesh.setVertexCount(mesh, toCapacity(drawn));
    // the count known: drawn directly, not with the count on the GPU
    TextureMesh.setDrawnVertexCount(mesh, drawn);
    // data derived from incomplete vertices
    if (reextracted || beyondBound) for (const handler of mesh.meta.reextracted?.values() ?? []) handler();
  });
  graphics.resources.pending.track(checked, 'isosurface vertex count');
  mesh.meta.extracted = checked.catch(() => {});
  return mesh;
}
