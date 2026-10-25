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
  getIsosurfacePages,
  type IsosurfaceEstimate,
  getVertexArenaStats,
  isVertexArenaExhausted,
  releaseVertexArenaPages,
  reserveVertexArenaPages,
  type VertexArenaReservation,
} from '@molstar/graphics/gl/compute/marching-cubes/chunks';
import type { IsosurfaceInput } from '@molstar/graphics/gl/compute/marching-cubes/compute';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { TextureMesh } from './texture-mesh.js';
import { isTimingMode } from '@molstar/core/util/debug';

/**
 * Where an extraction writes: pages of a vertex arena, `surface` reused; `estimate`, the expected size, sizes
 * the page table; `reservation`, the arena and the pages reserved for it (`reserveVertexArenaPages`)
 */
export type IsosurfaceTarget = {
  readonly surface?: ChunkedSurface;
  readonly estimate?: IsosurfaceEstimate;
  readonly reservation?: VertexArenaReservation;
};

/**
 * Extracts an isosurface into `target` (see `extractIsosurface`). Called after awaiting the reservation of its pages
 * and again when the surface did not get all its pages, so it must compute inputs others may change in between (e.g.
 * shared density textures) itself, and record its commands without awaiting after.
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
    target.reservation,
  );
}

/** Headroom of the expected count over the last vertex count, so small changes fit */
const Headroom = 1.25;
const MinCapacity = 3 * 1024;

/** The last count (3 per triangle) and vertices of a texture mesh, read back after its extraction */
const lastCounts = new WeakMap<TextureMesh, { readonly count: number; readonly vertices: number }>();
/** The latest extraction into a texture mesh, so that the read back of an earlier one isn't handled */
const latestExtractions = new WeakMap<TextureMesh, number>();
let nextExtraction = 0;

/** Whole triangles */
const toCapacity = (vertexCount: number) => Math.max(MinCapacity, Math.ceil(vertexCount / 3) * 3);

/** Extractions after the first, for a surface that did not get all its pages */
const MaxReextractions = 2;

/** How long after the last incomplete surface they are reported together */
const ReportDelayMs = 500;
const reports = new WeakMap<
  GraphicsContext,
  { surfaces: number; pages: number; tooLarge: number; timer?: ReturnType<typeof setTimeout> }
>();

const formatBytes = (bytes: number) =>
  bytes >= 2 ** 30 ? `${(bytes / 2 ** 30).toFixed(1)} GiB` : `${Math.ceil(bytes / 2 ** 20)} MiB`;

/**
 * Reports surfaces drawn incomplete (`GraphicsResources.outOfMemory`), together once no more came for a while;
 * `tooLarge`: a surface of more vertices than a storage buffer holds (else the memory ran out)
 */
function reportIncompleteSurface(
  graphics: GraphicsContext & { compute: ComputeContext },
  missingPages: number,
  tooLarge: boolean,
) {
  const report = reports.get(graphics) ?? { surfaces: 0, pages: 0, tooLarge: 0 };
  report.surfaces += 1;
  report.pages += missingPages;
  if (tooLarge) report.tooLarge += 1;
  clearTimeout(report.timer);
  report.timer = setTimeout(() => {
    reports.delete(graphics);
    const stats = getVertexArenaStats(graphics.compute);
    const missing = formatBytes(getChunkByteCount(report.pages));
    const surfaces = report.surfaces === 1 ? '1 surface' : `${report.surfaces} surfaces`;
    const reasons: string[] = [];
    if (report.tooLarge < report.surfaces) {
      reasons.push(`the GPU has no memory for more than the ${formatBytes(stats?.byteCount ?? 0)} of surface vertices`);
    }
    if (report.tooLarge > 0) {
      reasons.push(
        `${report.tooLarge} of them larger than a GPU storage buffer (${formatBytes(stats?.maxByteCount ?? 0)})`,
      );
    }
    graphics.resources.outOfMemory.next(
      `Not enough GPU memory: ${surfaces} drawn incomplete, about ${missing} of vertices missing; ` +
        `${reasons.join(', ')}. Use a coarser resolution or a lower quality.`,
    );
  }, ReportDelayMs);
  // not keeping Node (headless rendering, tests) alive: its timers are objects, browsers' numbers
  if (typeof report.timer === 'object') report.timer.unref();
  reports.set(graphics, report);
}

/** `extract` into `target`, its reservation returned when it fails (else when its count is read back) */
async function extractReserved(extract: IsosurfaceExtraction, target: IsosurfaceTarget) {
  try {
    return await extract(target);
  } catch (e) {
    if (target.reservation) releaseVertexArenaPages(target.reservation);
    throw e;
  }
}

/**
 * The isosurface in pages of a vertex arena, drawn indirectly with the count on the GPU: nothing is read back before
 * drawing. The count, read back afterwards, is drawn directly from then on (`TextureMesh.setDrawnVertexCount`), and keeps
 * the upper bound of the texture mesh (`vertexCount`) and the position iterator's data right. The pages of the
 * expected count are reserved before (`reserveVertexArenaPages`: an arena grown or added, the memory checked); a
 * surface that did not get all its pages (its table too small, or its arena's free pages ran out) is extracted again
 * with its count, its pages reserved first, which sizes the table and may move it to another arena; one the arenas
 * can't hold (out of memory, or larger than a storage buffer) is drawn incomplete and reported
 * (`GraphicsResources.outOfMemory`). The re-extraction handlers of the mesh (`TextureMesh.setReextractionHandler`)
 * run after either changed the vertices data derived from them saw. The expected count is the last count of
 * `textureMesh` (with headroom), else `estimate` (a function: called only then, e.g. a pass over the grid,
 * `estimateIsosurfaceVertexCount`). `boundingSphere` as a function is called after the first extraction (e.g. of what
 * the extraction computed).
 */
export async function computeTextureMeshIsosurface(
  graphics: GraphicsContext & { compute: ComputeContext },
  extract: IsosurfaceExtraction,
  groupCount: number,
  boundingSphere: Sphere3D | (() => Sphere3D),
  textureMesh?: TextureMesh,
  estimate?: number | (() => number),
): Promise<TextureMesh> {
  const { compute } = graphics;
  const surface = textureMesh?.doubleBuffer.getChunks(graphics);
  const last = textureMesh && lastCounts.get(textureMesh);
  const expected: IsosurfaceEstimate | undefined =
    last !== undefined
      ? { count: last.count * Headroom, vertices: Math.ceil(last.vertices * Headroom) }
      : typeof estimate === 'function'
        ? estimate()
        : estimate;
  const expectedCount = typeof expected === 'object' ? expected.count : expected;
  const reservation =
    expected !== undefined
      ? await reserveVertexArenaPages(compute, getIsosurfacePages(expected).pages, surface)
      : undefined;
  const result = await extractReserved(extract, { surface, estimate: expected, reservation });
  const mesh = TextureMesh.createChunked(
    toCapacity(expectedCount ?? 0),
    groupCount,
    result.surface,
    compute,
    typeof boundingSphere === 'function' ? boundingSphere() : boundingSphere,
    textureMesh,
  );
  mesh.meta.readVertexCount = undefined;
  const extraction = nextExtraction++;
  latestExtractions.set(mesh, extraction);

  // a newer extraction, or the mesh destroyed (which returns its pages)
  const isCurrent = () =>
    latestExtractions.get(mesh) === extraction && mesh.doubleBuffer.getChunks(graphics) === result.surface;
  const checked = result.count.then(async (first) => {
    if (!isCurrent()) return;
    let state = first;
    let reextracted = false;
    let tooLarge = false;
    for (let i = 0; i < MaxReextractions && state.missingPages > 0; ++i) {
      const size = { count: state.count, vertices: state.vertices };
      const again = await reserveVertexArenaPages(compute, getIsosurfacePages(size).pages, result.surface);
      // pointless when no arena can take the pages (until surfaces are destroyed)
      if (!isCurrent() || (state.arenaFull && !again.fits)) {
        releaseVertexArenaPages(again);
        if (!isCurrent()) return;
        tooLarge = !isVertexArenaExhausted(compute);
        break;
      }
      if (isTimingMode) {
        const reason = state.arenaFull ? 'the pages reserved first' : 'the page table sized to the count';
        console.log(`Isosurface extracted again: ${state.missingPages} pages missing, ${reason}`);
      }
      const extracted = await extractReserved(extract, {
        surface: result.surface,
        estimate: size,
        reservation: again,
      });
      // the table grown, the arena changed
      if (isCurrent()) TextureMesh.bindChunkedTable(mesh, extracted.surface);
      state = await extracted.count;
      if (!isCurrent()) return;
      reextracted = true;
    }
    if (state.missingPages > 0) reportIncompleteSurface(graphics, state.missingPages, tooLarge);
    const { count, drawn, pages } = state;
    lastCounts.set(mesh, { count, vertices: state.vertices });
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
