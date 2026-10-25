/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Box3D } from '@molstar/core/math/geometry';
import { Mat3, type Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util';
import type { ComputeContext, DrawIndirect, StorageBuffer } from '../../compute.js';
import { getScratchBuffer } from '../scratch.js';
import { IsosurfaceChunkSize, IsosurfaceTableHeader, IsosurfaceVertexBytes } from './chunk-size.js';
import {
  getInputValues,
  getIsosurfaceKernels,
  getShared,
  type IsosurfaceInput,
  McChunkAllocateKernel,
  McChunkAllocateSchema,
  McChunkCountKernel,
  McChunkCountOutput,
  McChunkEmitKernel,
  McChunkEmitOutput,
  McChunkFreeKernel,
  McChunkFreeSchema,
  McChunkPushKernel,
  McChunkPushSchema,
  WorkgroupSize,
} from './compute.js';
import type { GraphicsTimer } from '../../timer.js';

/**
 * Marching cubes into chunks of vertex arenas, without reading the vertex count back first: an arena is a pool of
 * pages (`IsosurfaceChunkSize` vertices, or as many triangles), a context has a list of them shared by all chunked
 * surfaces; a surface is indexed (a vertex per crossed edge, shared by the triangles around it: about a sixth of the
 * vertices drawn), has a page table (its counts, its box and the ids of its vertex and triangle pages) and its pages in
 * one arena, the GPU takes the pages its counts need from the arena's free pages (`mc-chunk-allocate`) between
 * counting (`mc-chunk-count`, one atomic reservation per workgroup) and emitting (`mc-chunk-emit`). The vertex and
 * triangle order varies between runs (workgroups and cells reserve in the order they run). The count and the arena's free pages are read back afterwards, never awaited before drawing.
 * The pages an extraction expects are reserved beforehand (`reserveVertexArenaPages`: an arena grown, or one added,
 * the device's memory checked before any command uses it); a surface that did not get all its pages is extracted again
 * with its count.
 */

const TableHeader = IsosurfaceTableHeader;
/** Most steps of a quantized position (16 bits per axis), see `iso_vertex_wgsl` */
const QuantizedSteps = 65535;
/** Pages of the first arena (256 K vertices, 3 MB) */
const InitialPages = 64;
/**
 * Bytes of an arena (the first grows to it, others are created with it): a growth copies, the old and the new buffers
 * at once, more arenas copy nothing; larger for a surface that needs more, up to the storage buffer limit
 */
const ArenaBytes = 256 * 1024 * 1024;
/** Free pages kept beyond the estimate of an extraction */
const HeadroomPages = 16;
/**
 * Page ids of a table beyond the estimate (at 4 bytes a page, a table of the arena's pages was 128 KiB per surface at
 * a 2 GiB arena): a surface exceeding them is extracted again with its count (`ChunkedSurfaceCount.missingPages`)
 */
const TableHeadroom = 2;
const TableSlackPages = 16;
const GrowthFactor = 1.5;
/**
 * Vertices per triangle of an estimate: about 1/2 for a closed surface (Euler), more where the grid's faces cut it and
 * for small surfaces, see `getIsosurfacePages`
 */
const VerticesPerTriangle = 0.6;
/**
 * How long an arena stays without surfaces before it is released: longer than a frame of a trajectory played back,
 * whose surfaces may be created again each frame (the pages of the old ones freed, taken again by the new ones), so
 * that its arenas are kept while playing and released when it stops
 */
const ArenaIdleMs = 3000;

/** Pages of a context's surfaces in one buffer */
export interface VertexArena {
  /**
   * 3 `u32` per vertex (`IsosurfaceVertexBytes`): the position quantized in the surface's box, the normal octahedral
   * and the packed group, see `iso_vertex_wgsl`
   */
  readonly vertices: ValueCell<StorageBuffer>;
  /** the free page stack: its size, then the page ids */
  free: StorageBuffer;
  pageCount: number;
  /** the most pages it grows to */
  readonly maxPages: number;
  /** free pages at the last readback */
  freeRead: number;
  /** pages extractions since that readback may take */
  reserved: number;
  /** surfaces with pages here (or about to take them): an arena without any is released, see `ArenaIdleMs` */
  surfaces: number;
  releaseTimer?: ReturnType<typeof setTimeout>;
  readonly owner: VertexArenas;
}

interface VertexArenas {
  readonly compute: ComputeContext;
  readonly list: VertexArena[];
  /** an arena could not be grown or created for want of memory, until a surface returns its pages */
  exhausted: boolean;
  /** reservations in turn: one may grow or create an arena, which the next then sees */
  lock: Promise<unknown>;
  /** called after an arena grew (its buffers replaced), see `onVertexArenaGrow` */
  readonly listeners: Set<() => void>;
}

const contexts = new WeakMap<ComputeContext, VertexArenas>();

/** Pages of a buffer at the storage buffer limit */
function getMaxPages(compute: ComputeContext) {
  return Math.floor(compute.maxStorageBufferSize / getChunkByteCount(1));
}

function getArenaPages(compute: ComputeContext) {
  return Math.min(getMaxPages(compute), Math.floor(ArenaBytes / getChunkByteCount(1)));
}

function getArenas(compute: ComputeContext): VertexArenas {
  let arenas = contexts.get(compute);
  if (!arenas) {
    arenas = { compute, list: [], exhausted: false, lock: Promise.resolve(), listeners: new Set() };
    contexts.set(compute, arenas);
  }
  if (arenas.list.length === 0) {
    // small, created without awaiting the device (as before the first extraction, or after all were released)
    const pages = Math.min(InitialPages, getMaxPages(compute));
    const vertices = compute.createBuffer(getChunkByteCount(pages), 'isosurface vertices');
    const free = compute.createBuffer((pages + 1) * 4, 'isosurface free pages');
    pushPages(compute, free, 0, pages);
    arenas.list.push({
      vertices: ValueCell.create(vertices),
      free,
      pageCount: pages,
      maxPages: getArenaPages(compute),
      freeRead: pages,
      reserved: 0,
      surfaces: 0,
      owner: arenas,
    });
  }
  return arenas;
}

/**
 * Releases `arena` once it has had no surfaces and no reservations for `ArenaIdleMs`, or keeps it when it has: called
 * whenever either changed
 */
function updateIdle(arena: VertexArena) {
  if (arena.surfaces > 0 || arena.reserved > 0) {
    clearTimeout(arena.releaseTimer);
    arena.releaseTimer = undefined;
    return;
  }
  if (arena.releaseTimer !== undefined) return;
  arena.releaseTimer = setTimeout(() => {
    arena.releaseTimer = undefined;
    if (arena.surfaces === 0 && arena.reserved === 0) releaseArena(arena);
  }, ArenaIdleMs);
  // not keeping Node (headless rendering, tests) alive: its timers are objects, browsers' numbers
  if (typeof arena.releaseTimer === 'object') arena.releaseTimer.unref();
}

/** Destroys the buffers of `arena`, which no surface uses (none draws them), and drops it from its context's arenas */
function releaseArena(arena: VertexArena) {
  const { list } = arena.owner;
  const index = list.indexOf(arena);
  if (index < 0) return;
  list.splice(index, 1);
  // destroyed after the commands recorded with them
  arena.vertices.ref.value.destroy();
  arena.free.destroy();
  arena.owner.exhausted = false;
}

function addReserved(arena: VertexArena, pages: number) {
  arena.reserved += pages;
  updateIdle(arena);
}

function removeReserved(arena: VertexArena, pages: number) {
  arena.reserved -= pages;
  updateIdle(arena);
}

/** `surface` (new, or moving) into `arena`, out of its previous one */
function setArena(surface: { arena?: VertexArena }, arena: VertexArena) {
  if (surface.arena === arena) return;
  if (surface.arena) {
    surface.arena.surfaces -= 1;
    updateIdle(surface.arena);
  }
  surface.arena = arena;
  arena.surfaces += 1;
  updateIdle(arena);
}

/** The vertex arenas of `compute`, the first created if there is none */
export function getVertexArenas(compute: ComputeContext): readonly VertexArena[] {
  return getArenas(compute).list;
}

/** Free pages to expect: at the last readback, minus those reserved since */
export function getAvailablePages(arena: VertexArena) {
  return arena.freeRead - arena.reserved;
}

function pushPages(compute: ComputeContext, free: StorageBuffer, first: number, count: number) {
  const kernel = compute.createKernel(McChunkPushKernel, McChunkPushSchema, {
    uFirst: ValueCell.create(first),
    uCount: ValueCell.create(count),
    free: ValueCell.create(free),
  });
  compute.beginComputePass({ label: 'isosurface arena pages' });
  kernel.dispatch(1);
  compute.endComputePass();
  kernel.destroy();
}

/** Both buffers of an arena of `pages`, or none when the device ran out of memory */
async function createArenaBuffers(compute: ComputeContext, pages: number) {
  const vertices = await compute.tryCreateBuffer(getChunkByteCount(pages), 'isosurface vertices');
  const free = vertices && (await compute.tryCreateBuffer((pages + 1) * 4, 'isosurface free pages'));
  if (!vertices || !free) {
    vertices?.destroy();
    return;
  }
  return { vertices, free };
}

/**
 * Grows `arena` by at least `pages` (by `GrowthFactor` at least, up to its `maxPages`): new buffers with the contents
 * of the old ones (page ids stay valid), the new pages pushed onto the free stack; returns whether it grew (not at its
 * most pages, the memory there)
 */
export async function growVertexArena(compute: ComputeContext, arena: VertexArena, pages: number) {
  const old = arena.pageCount;
  const count = Math.min(arena.maxPages, Math.max(Math.ceil(old * GrowthFactor), old + pages));
  if (count <= old) return false;
  const buffers = await createArenaBuffers(compute, count);
  if (!buffers) return false;
  // grown meanwhile (callers reserve in turn, see `reserveVertexArenaPages`): kept as it is
  if (arena.pageCount !== old) {
    buffers.vertices.destroy();
    buffers.free.destroy();
    return true;
  }
  const { vertices, free } = buffers;
  arena.vertices.ref.value.copyTo(vertices);
  arena.free.copyTo(free);
  pushPages(compute, free, old, count - old);
  // destroyed after the commands recorded with them
  arena.vertices.ref.value.destroy();
  arena.free.destroy();
  ValueCell.update(arena.vertices, vertices);
  arena.free = free;
  arena.pageCount = count;
  arena.freeRead += count - old;
  for (const listener of getArenas(compute).listeners) listener();
  return true;
}

/** An arena of `pages` added to the arenas of `compute`, if the device has the memory */
async function addVertexArena(compute: ComputeContext, pages: number) {
  const buffers = await createArenaBuffers(compute, pages);
  if (!buffers) return;
  pushPages(compute, buffers.free, 0, pages);
  const arenas = getArenas(compute);
  const arena: VertexArena = {
    vertices: ValueCell.create(buffers.vertices),
    free: buffers.free,
    pageCount: pages,
    maxPages: pages,
    freeRead: pages,
    reserved: 0,
    surfaces: 0,
    owner: arenas,
  };
  arenas.list.push(arena);
  // reserved by the caller right after, see `reserveVertexArenaPages`
  return arena;
}

/** Pages of an arena for a surface, reserved (`VertexArena.reserved`) until the extraction's count is read back */
export interface VertexArenaReservation {
  readonly arena: VertexArena;
  readonly pages: number;
  /** whether `arena` is expected to have them; else the most it could give (the memory, or a surface of more pages than a storage buffer holds) */
  readonly fits: boolean;
}

/**
 * The arena for a surface of `pages` pages (`surface`, if given, moves only when its arena can't take them): one
 * with the free pages, else grown, else a new one (created with the device's memory checked, so no command uses a
 * buffer that failed); pages beyond the surface's current ones reserved, see `computeChunkedIsosurface`
 */
export function reserveVertexArenaPages(compute: ComputeContext, pages: number, surface?: ChunkedSurface) {
  const arenas = getArenas(compute);
  const run = arenas.lock.then(() => reserve(compute, arenas, pages, surface));
  arenas.lock = run.catch(() => {});
  return run;
}

async function reserve(
  compute: ComputeContext,
  arenas: VertexArenas,
  pages: number,
  surface?: ChunkedSurface,
): Promise<VertexArenaReservation> {
  const needed = (arena: VertexArena) => Math.max(0, pages - (surface?.arena === arena ? surface.pages : 0));
  const fits = (arena: VertexArena) => getAvailablePages(arena) >= needed(arena) + HeadroomPages;
  const take = (arena: VertexArena, fit: boolean): VertexArenaReservation => {
    const reserved = needed(arena);
    addReserved(arena, reserved);
    return { arena, pages: reserved, fits: fit };
  };
  // its own arena, any with the free pages, one grown, a new one
  const ordered = surface ? [surface.arena, ...arenas.list.filter((a) => a !== surface.arena)] : arenas.list;
  for (const arena of ordered) if (fits(arena)) return take(arena, true);
  for (const arena of ordered) {
    if (arena.pageCount >= arena.maxPages || arena.maxPages < needed(arena) + HeadroomPages) continue;
    const missing = needed(arena) + HeadroomPages - getAvailablePages(arena);
    if (await growVertexArena(compute, arena, missing)) {
      if (fits(arena)) return take(arena, true);
    } else {
      arenas.exhausted = true;
    }
  }
  const maxPages = getMaxPages(compute);
  if (!arenas.exhausted && pages + HeadroomPages <= maxPages) {
    const added = await addVertexArena(compute, Math.max(getArenaPages(compute), pages + HeadroomPages));
    if (added) return take(added, true);
    arenas.exhausted = true;
  }
  // the most the arenas give: its own arena when it has the most
  let best = ordered[0];
  for (const arena of ordered) if (getAvailablePages(arena) > getAvailablePages(best)) best = arena;
  return take(best, false);
}

/** Returns the pages of `reservation` (an extraction that did not take place) */
export function releaseVertexArenaPages(reservation: VertexArenaReservation) {
  removeReserved(reservation.arena, reservation.pages);
}

/** Whether the arenas of `compute` ran out of memory (an arena not grown or created, see `reserveVertexArenaPages`) */
export function isVertexArenaExhausted(compute: ComputeContext) {
  return !!contexts.get(compute)?.exhausted;
}

/** Calls `listener` whenever an arena of `compute` grew (e.g. to bind its new buffers); returns the unsubscribe */
export function onVertexArenaGrow(compute: ComputeContext, listener: () => void) {
  const { listeners } = getArenas(compute);
  listeners.add(listener);
  return () => {
    listeners.delete(listener);
  };
}

/** Without a reservation: the surface's arena, else the first with the free pages, else the one with the most */
function pickArena(arenas: VertexArenas, pages: number, surface?: ChunkedSurface) {
  if (surface) return surface.arena;
  let best = arenas.list[0];
  for (const arena of arenas.list) {
    if (getAvailablePages(arena) >= pages + HeadroomPages) return arena;
    if (getAvailablePages(arena) > getAvailablePages(best)) best = arena;
  }
  return best;
}

/** A surface in an arena: its page table and draw arguments */
export interface ChunkedSurface {
  /** the arena of its pages; changed by an extraction with a reservation of another (its pages moved there) */
  arena: VertexArena;
  /** the counts and the box of the positions (`TableHeader` words), then the ids of the vertex and triangle pages */
  table: StorageBuffer;
  /** vertex pages of the table, the triangle pages after them */
  vertexCapacity: number;
  /**
   * tables the surface grew out of: what draws the surface may still bind one, released when it binds the current
   * table (`releaseRetiredTables`)
   */
  readonly retired: StorageBuffer[];
  readonly drawArgs: StorageBuffer;
  /** pages of the surface at its last readback, of vertices and triangles */
  pages: number;
}

/**
 * An expected size of an isosurface: its count (3 per triangle, the vertices drawn), with the vertices of the
 * indexed surface when known (e.g. read back), else about `VerticesPerTriangle` a triangle
 */
export type IsosurfaceEstimate = number | { readonly count: number; readonly vertices: number };

/** The vertex and triangle pages of a surface of `estimate` */
export function getIsosurfacePages(estimate: IsosurfaceEstimate) {
  const count = typeof estimate === 'number' ? estimate : estimate.count;
  const triangles = Math.ceil(count / 3);
  const vertices = typeof estimate === 'number' ? Math.ceil(triangles * VerticesPerTriangle) : estimate.vertices;
  const vertexPages = getChunkPages(vertices);
  const trianglePages = getChunkPages(triangles);
  return { vertexPages, trianglePages, pages: vertexPages + trianglePages };
}

/** Destroys the tables `surface` grew out of, once what draws it binds the current one */
export function releaseRetiredTables(surface: ChunkedSurface) {
  for (const table of surface.retired) table.destroy();
  surface.retired.length = 0;
}

/** The state of a surface after an extraction, read back */
export interface ChunkedSurfaceCount {
  /** the count of the surface, 3 per triangle (more than drawn when it did not get the pages it needed) */
  readonly count: number;
  /** the vertices of the indexed surface */
  readonly vertices: number;
  readonly drawn: number;
  /** of vertices and triangles */
  readonly pages: number;
  /** pages missing, 0 if the surface got all it needed */
  readonly missingPages: number;
  /** pages were missing because the arena's free pages ran out (else the table was too small for the count) */
  readonly arenaFull: boolean;
}

export interface ChunkedIsosurfaceResult {
  readonly surface: ChunkedSurface;
  readonly drawIndirect: DrawIndirect;
  /** the state after the extraction, read back (not needed to draw it) */
  readonly count: Promise<ChunkedSurfaceCount>;
}

function createTable(compute: ComputeContext, vertexPages: number, trianglePages: number) {
  return compute.createBuffer((TableHeader + vertexPages + trianglePages) * 4, 'isosurface page table');
}

/**
 * Extracts the isosurface of `input` into an arena (`surface` reuses its pages and draw arguments): into the arena of
 * `reservation` (`reserveVertexArenaPages`, awaited before, so that this records its commands without awaiting),
 * else the surface's own or one with the free pages; `estimate`, the expected size, sizes the page table (without it,
 * the table spans the arena's pages for either) and the pages reserved without a reservation
 */
export async function computeChunkedIsosurface(
  compute: ComputeContext,
  input: IsosurfaceInput,
  gridDim: Vec3,
  isoValue: number,
  transform: Mat4,
  invert: boolean,
  surface?: ChunkedSurface,
  estimate?: IsosurfaceEstimate,
  /** times the extraction's passes (timing mode) */
  timer?: GraphicsTimer,
  reservation?: VertexArenaReservation,
): Promise<ChunkedIsosurfaceResult> {
  // compiled in the background, not blocking the queue at the first dispatch; awaited only when not ready, the
  // extraction then records its commands without awaiting
  const kernels = getIsosurfaceKernels(input);
  if (!compute.hasKernels(kernels)) await compute.prepareKernels(kernels);

  const arenas = getArenas(compute);
  const expected = estimate !== undefined ? getIsosurfacePages(estimate) : undefined;
  const estimatePages = expected?.pages ?? 0;
  // a reservation of no pages doesn't hold its arena, which may have been released
  const reserved = reservation && arenas.list.includes(reservation.arena) ? reservation : undefined;
  const arena = reserved?.arena ?? pickArena(arenas, estimatePages, surface);
  let reserve = reserved?.pages ?? 0;
  if (!reserved) {
    reserve = Math.max(0, estimatePages - (surface?.arena === arena ? surface.pages : 0));
    addReserved(arena, reserve);
  }
  if (surface && surface.arena !== arena) {
    // its pages returned to its arena, taken from the new one by `mc-chunk-allocate`
    freePages(compute, surface);
    setArena(surface, arena);
  }

  const capacity = (pages: number) => Math.min(arena.maxPages, pages * TableHeadroom + TableSlackPages);
  const vertexCapacity = expected ? capacity(expected.vertexPages) : arena.pageCount;
  const triangleCapacity = expected ? capacity(expected.trianglePages) : arena.pageCount;
  if (!surface) {
    surface = {
      arena,
      table: createTable(compute, vertexCapacity, triangleCapacity),
      vertexCapacity,
      retired: [],
      drawArgs: compute.createBuffer(16, 'isosurface draw arguments'),
      pages: 0,
    };
    arena.surfaces += 1;
    updateIdle(arena);
  } else {
    const triangles = surface.table.byteLength / 4 - TableHeader - surface.vertexCapacity;
    if (surface.vertexCapacity < vertexCapacity || triangles < triangleCapacity) {
      // the page ids at other offsets of a new table: the pages returned (in queue order, before they are taken
      // again), not copied
      freePages(compute, surface);
      // what draws the surface binds the old table until it rebinds (e.g. a re-extraction after the arena grew)
      surface.retired.push(surface.table);
      surface.table = createTable(compute, vertexCapacity, triangleCapacity);
      surface.vertexCapacity = vertexCapacity;
    }
  }
  const { table, drawArgs } = surface;

  const [dx, dy, dz] = gridDim;
  // per grid point (its edges), the cells by their lower point
  const groups = (n: number, size: number) => Math.max(1, Math.ceil(n / size));
  const dispatch = [groups(dx, WorkgroupSize[0]), groups(dy, WorkgroupSize[1]), groups(dz, WorkgroupSize[2])] as const;
  const workgroupCount = dispatch[0] * dispatch[1] * dispatch[2];
  const { triCounts, triTable } = getShared(compute);
  const { schema: inputSchema, values: inputValues } = getInputValues(input, gridDim, isoValue);
  // the positions quantized in the box of the grid's vertices (as `TextureMesh.getGridVertexBox`)
  const box = Box3D.create(Vec3.create(0, 0, 0), Vec3.create(dx - 1, dy - 1, dz - 1));
  Box3D.transform(box, box, transform);
  const boxSize = Box3D.size(Vec3(), box);
  const boxStep = Vec3.scale(Vec3(), boxSize, 1 / QuantizedSteps);
  const toSteps = (size: number) => (size > 0 ? QuantizedSteps / size : 0);
  const boxScale = Vec3.create(toSteps(boxSize[0]), toSteps(boxSize[1]), toSteps(boxSize[2]));
  // kept between calls (see `getScratchBuffer`): nothing is awaited from here on
  const bases = getScratchBuffer(compute, 'mc chunk bases', workgroupCount * 4);
  const pointBases = getScratchBuffer(compute, 'mc chunk point bases', dx * dy * dz * 4);

  const countKernel = compute.createKernel(McChunkCountKernel, { ...inputSchema, ...McChunkCountOutput }, {
    ...inputValues,
    triCounts: ValueCell.create(triCounts),
    table: ValueCell.create(table),
    bases: ValueCell.create(bases),
    pointBases: ValueCell.create(pointBases),
  } as never);
  const allocateKernel = compute.createKernel(McChunkAllocateKernel, McChunkAllocateSchema, {
    uBoxMin: ValueCell.create(box.min),
    uBoxStep: ValueCell.create(boxStep),
    uVertexCapacity: ValueCell.create(surface.vertexCapacity),
    table: ValueCell.create(table),
    free: ValueCell.create(arena.free),
    drawArgs: ValueCell.create(drawArgs),
  });
  const emitKernel = compute.createKernel(McChunkEmitKernel, { ...inputSchema, ...McChunkEmitOutput }, {
    ...inputValues,
    uInvert: ValueCell.create(invert),
    uTransform: ValueCell.create(transform),
    // the normal matrix, as on WebGL
    uTransformAdjoint: ValueCell.create(Mat3.adjointFromMat4(Mat3(), transform)),
    uBoxMin: ValueCell.create(box.min),
    uBoxScale: ValueCell.create(boxScale),
    triCounts: ValueCell.create(triCounts),
    triTable: ValueCell.create(triTable),
    bases: ValueCell.create(bases),
    pointBases: ValueCell.create(pointBases),
    table: ValueCell.create(table),
    vertices: ValueCell.create(arena.vertices.ref.value),
  } as never);

  timer?.group('computeChunkedIsosurface');
  // the vertex and triangle counts from 0, the pages of the previous extraction kept (reused, or freed by
  // `mc-chunk-allocate`)
  table.clear(0, 4);
  table.clear(10 * 4, 4);
  compute.beginComputePass({ label: 'marching cubes chunk count' });
  countKernel.dispatch(...dispatch);
  compute.endComputePass();
  compute.beginComputePass({ label: 'marching cubes chunk allocate' });
  allocateKernel.dispatch(1);
  compute.endComputePass();
  compute.beginComputePass({ label: 'marching cubes chunk emit' });
  emitKernel.dispatch(...dispatch);
  compute.endComputePass();
  timer?.groupEnd('computeChunkedIsosurface');
  countKernel.destroy();
  allocateKernel.destroy();
  emitKernel.destroy();

  const tableRead = table.readAsync(0, TableHeader * 4);
  const freeRead = arena.free.readAsync(0, 4);
  const current = surface;
  const vertexPageCapacity = surface.vertexCapacity;
  const trianglePageCapacity = table.byteLength / 4 - TableHeader - vertexPageCapacity;
  const count = Promise.all([tableRead, freeRead]).then(([t, f]) => {
    const header = new Uint32Array(t);
    const [vertices, vertexPages, status, drawn] = header;
    const triangles = header[10];
    const trianglePages = header[11];
    arena.freeRead = new Uint32Array(f)[0];
    removeReserved(arena, reserve);
    const pages = vertexPages + trianglePages;
    current.pages = pages;
    const neededVertices = getChunkPages(vertices);
    const neededTriangles = getChunkPages(triangles);
    const missing = neededVertices - vertexPages + neededTriangles - trianglePages;
    return {
      count: triangles * 3,
      vertices,
      drawn,
      pages,
      missingPages: status ? Math.max(1, missing) : 0,
      arenaFull:
        !!status &&
        (vertexPages < Math.min(neededVertices, vertexPageCapacity) ||
          trianglePages < Math.min(neededTriangles, trianglePageCapacity)),
    };
  });
  return { surface, drawIndirect: { buffer: drawArgs, byteOffset: 0 }, count };
}

/** Returns the pages of `surface` to its arena (in queue order), its counts reset */
function freePages(compute: ComputeContext, surface: ChunkedSurface) {
  const { arena } = surface;
  const kernel = compute.createKernel(McChunkFreeKernel, McChunkFreeSchema, {
    table: ValueCell.create(surface.table),
    free: ValueCell.create(arena.free),
  });
  compute.beginComputePass({ label: 'isosurface free pages' });
  kernel.dispatch(1);
  compute.endComputePass();
  kernel.destroy();
  // free from now on (in queue order); a readback of an earlier extraction may still report fewer, which only
  // grows the arena sooner
  arena.freeRead += surface.pages;
  surface.pages = 0;
  const arenas = contexts.get(compute);
  if (arenas) arenas.exhausted = false;
}

/** Returns the pages of `surface` to its arena and destroys its buffers */
export function destroyChunkedSurface(compute: ComputeContext, surface: ChunkedSurface) {
  // the arenas of another context (before a restore or backend switch) are gone
  if (contexts.get(compute)?.list.includes(surface.arena)) {
    freePages(compute, surface);
    surface.arena.surfaces -= 1;
    updateIdle(surface.arena);
  }
  surface.table.destroy();
  releaseRetiredTables(surface);
  surface.drawArgs.destroy();
}

/**
 * The first `count` vertices drawn of the surface of `table` (3 per triangle, unindexed), gathered from its pages and
 * read back: positions (`w` 1), normals and groups (packed rgba bytes) per vertex, as the texture readback of a
 * texture mesh gives them
 */
export async function readChunkedVertices(
  compute: ComputeContext,
  chunks: { readonly table: StorageBuffer; readonly vertices: StorageBuffer },
  count: number,
) {
  const { table, vertices: arenaVertices } = chunks;
  const header = new Uint32Array(await table.readAsync());
  const triangleCount = Math.min(Math.floor(count / 3), Math.floor(header[3] / 3));
  const vertexPages = header[1];
  const trianglePages = Math.min(header[11], getChunkPages(triangleCount));
  const vertexIds = header.subarray(TableHeader, TableHeader + vertexPages);
  const triangleIds = header.subarray(TableHeader + header[12], TableHeader + header[12] + trianglePages);
  const box = new Float32Array(header.buffer, 4 * 4, 6);
  const pageBytes = getChunkByteCount(1);
  const pages = vertexPages + trianglePages;
  const gathered = compute.createBuffer(Math.max(1, pages) * pageBytes, 'isosurface read vertices');
  try {
    const ids = [...vertexIds, ...triangleIds];
    for (let i = 0; i < pages; ++i) {
      arenaVertices.copyTo(gathered, ids[i] * pageBytes, i * pageBytes, pageBytes);
    }
    const words = new Uint32Array(await gathered.readAsync());
    const indexOffset = vertexPages * IsosurfaceChunkSize * 3;
    const vertexCount = triangleCount * 3;
    const vertices = new Float32Array(vertexCount * 4);
    const normals = new Float32Array(vertexCount * 4);
    const groups = new Uint8Array(vertexCount * 4);
    const normal = Vec3();
    for (let i = 0; i < vertexCount; ++i) {
      const v = words[indexOffset + i];
      const w0 = words[v * 3];
      const w1 = words[v * 3 + 1];
      const w2 = words[v * 3 + 2];
      vertices[i * 4] = box[0] + (w0 & 0xffff) * box[3];
      vertices[i * 4 + 1] = box[1] + (w0 >>> 16) * box[4];
      vertices[i * 4 + 2] = box[2] + (w1 & 0xffff) * box[5];
      vertices[i * 4 + 3] = 1;
      decodeOctahedral(normal, (w1 >>> 16) | ((w2 >>> 24) << 16));
      normals[i * 4] = normal[0];
      normals[i * 4 + 1] = normal[1];
      normals[i * 4 + 2] = normal[2];
      // the bytes of the packed rgb, opaque as in group textures
      groups[i * 4] = (w2 >>> 16) & 0xff;
      groups[i * 4 + 1] = (w2 >>> 8) & 0xff;
      groups[i * 4 + 2] = w2 & 0xff;
      groups[i * 4 + 3] = 255;
    }
    return { vertices, normals, groups, vertexCount };
  } finally {
    gathered.destroy();
  }
}

/** The unit vector of the 24-bit octahedral `bits` (12 per axis, as `iso_vertex_wgsl`) */
export function decodeOctahedral(out: Vec3, bits: number) {
  const x = ((bits & 0xfff) / 4095) * 2 - 1;
  const y = ((bits >>> 12) / 4095) * 2 - 1;
  const z = 1 - Math.abs(x) - Math.abs(y);
  const t = Math.max(-z, 0);
  return Vec3.normalize(out, Vec3.set(out, x + (x >= 0 ? -t : t), y + (y >= 0 ? -t : t), z));
}

/** Bytes of `pages` pages of the arena */
export function getChunkByteCount(pages: number) {
  return pages * IsosurfaceChunkSize * IsosurfaceVertexBytes;
}

/** The arenas of `compute`, if created: their pages, the free ones at the last readback, and their bytes */
export function getVertexArenaStats(compute: ComputeContext) {
  const arenas = contexts.get(compute);
  if (!arenas) return;
  let pages = 0;
  let used = 0;
  for (const arena of arenas.list) {
    pages += arena.pageCount;
    used += arena.pageCount - arena.freeRead;
  }
  return {
    arenas: arenas.list.length,
    pages,
    usedPages: used,
    byteCount: getChunkByteCount(pages),
    usedByteCount: getChunkByteCount(used),
    /** the bytes of the most pages a storage buffer allows */
    maxByteCount: getChunkByteCount(getMaxPages(compute)),
    exhausted: arenas.exhausted,
  };
}

/** Pages of `count` vertices, or of as many triangles */
export function getChunkPages(count: number) {
  return Math.ceil(count / IsosurfaceChunkSize);
}
