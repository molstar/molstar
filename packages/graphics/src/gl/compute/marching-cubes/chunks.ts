/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat3, type Mat4, type Vec3 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util';
import type { ComputeContext, DrawIndirect, StorageBuffer } from '../../compute.js';
import { getScratchBuffer } from '../scratch.js';
import { IsosurfaceChunkSize } from './chunk-size.js';
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
 * Marching cubes into chunks of a vertex arena, without reading the vertex count back first: the arena is a pool of
 * pages (`IsosurfaceChunkSize` vertices each) per context, shared by all chunked surfaces; a surface has a page table
 * (its counts and page ids), the GPU takes the pages its count needs from the arena's free pages (`mc-chunk-allocate`)
 * between counting (`mc-chunk-count`, one atomic reservation per workgroup) and emitting (`mc-chunk-emit`). The vertex
 * order varies between runs (workgroups and cells reserve in the order they run). The count and the arena's free pages
 * are read back afterwards, never awaited before drawing: they keep the arena large enough (grown before an extraction
 * whose estimate the free pages may not cover) and tell when a surface did not get all its pages (it is extracted
 * again after the arena grew).
 */

/** Words of a page table before its page ids: vertex count, page count, status, drawn vertex count */
const TableHeader = 4;
/** Bytes per vertex in the arena's buffers */
const PositionBytes = 16;
const NormalBytes = 8;
/** Pages of a new arena (256 K vertices, 6 MB) */
const InitialPages = 64;
/** Free pages kept beyond the estimate of an extraction */
const HeadroomPages = 16;
const GrowthFactor = 1.5;

/** The pages of the chunked surfaces of a context */
export interface VertexArena {
  /** `vec4<f32>` per vertex: the position, the group's packed rgba bits in `w` */
  readonly positions: ValueCell<StorageBuffer>;
  /** `vec2<u32>` per vertex: the normal as snorm16 (`pack2x16snorm` of xy and of z) */
  readonly normals: ValueCell<StorageBuffer>;
  /** the free page stack: its size, then the page ids */
  free: StorageBuffer;
  pageCount: number;
  /** free pages at the last readback */
  freeRead: number;
  /** pages extractions since that readback may take */
  reserved: number;
  /** called after the arena grew (its buffers replaced), see `onVertexArenaGrow` */
  readonly listeners: Set<() => void>;
}

const arenas = new WeakMap<ComputeContext, VertexArena>();

function getMaxPages(compute: ComputeContext) {
  return Math.floor(compute.maxStorageBufferSize / (IsosurfaceChunkSize * PositionBytes));
}

/** The vertex arena of `compute`, created with `InitialPages` */
export function getVertexArena(compute: ComputeContext): VertexArena {
  let arena = arenas.get(compute);
  if (!arena) {
    const pages = Math.min(InitialPages, getMaxPages(compute));
    arena = {
      positions: ValueCell.create(
        compute.createBuffer(pages * IsosurfaceChunkSize * PositionBytes, 'isosurface vertex positions'),
      ),
      normals: ValueCell.create(
        compute.createBuffer(pages * IsosurfaceChunkSize * NormalBytes, 'isosurface vertex normals'),
      ),
      free: compute.createBuffer((pages + 1) * 4, 'isosurface free pages'),
      pageCount: 0,
      freeRead: 0,
      reserved: 0,
      listeners: new Set(),
    };
    arenas.set(compute, arena);
    pushPages(compute, arena.free, 0, pages);
    arena.pageCount = pages;
    arena.freeRead = pages;
  }
  return arena;
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

/**
 * Grows the arena by at least `pages` (by `GrowthFactor` at least, up to the storage buffer limit): new buffers with
 * the contents of the old ones (page ids stay valid), the new pages pushed onto the free stack; returns whether it grew
 */
export function growVertexArena(compute: ComputeContext, arena: VertexArena, pages: number) {
  const maxPages = getMaxPages(compute);
  const old = arena.pageCount;
  const count = Math.min(maxPages, Math.max(Math.ceil(old * GrowthFactor), old + pages));
  if (count <= old) return false;
  const positions = compute.createBuffer(count * IsosurfaceChunkSize * PositionBytes, 'isosurface vertex positions');
  const normals = compute.createBuffer(count * IsosurfaceChunkSize * NormalBytes, 'isosurface vertex normals');
  const free = compute.createBuffer((count + 1) * 4, 'isosurface free pages');
  arena.positions.ref.value.copyTo(positions);
  arena.normals.ref.value.copyTo(normals);
  arena.free.copyTo(free);
  pushPages(compute, free, old, count - old);
  // destroyed after the commands recorded with them
  arena.positions.ref.value.destroy();
  arena.normals.ref.value.destroy();
  arena.free.destroy();
  ValueCell.update(arena.positions, positions);
  ValueCell.update(arena.normals, normals);
  arena.free = free;
  arena.pageCount = count;
  arena.freeRead += count - old;
  for (const listener of arena.listeners) listener();
  return true;
}

/** Calls `listener` whenever the arena of `compute` grew (e.g. to bind its new buffers); returns the unsubscribe */
export function onVertexArenaGrow(compute: ComputeContext, listener: () => void) {
  const { listeners } = getVertexArena(compute);
  listeners.add(listener);
  return () => {
    listeners.delete(listener);
  };
}

/** A surface in the arena: its page table and draw arguments */
export interface ChunkedSurface {
  /** the counts (`TableHeader` words), then the page ids */
  table: StorageBuffer;
  /**
   * tables the surface grew out of: what draws the surface may still bind one, released when it binds the current
   * table (`releaseRetiredTables`)
   */
  readonly retired: StorageBuffer[];
  readonly drawArgs: StorageBuffer;
  /** pages of the surface at its last readback */
  pages: number;
}

/** Destroys the tables `surface` grew out of, once what draws it binds the current one */
export function releaseRetiredTables(surface: ChunkedSurface) {
  for (const table of surface.retired) table.destroy();
  surface.retired.length = 0;
}

/** The state of a surface after an extraction, read back */
export interface ChunkedSurfaceCount {
  /** the vertex count of the surface (more than drawn when it did not get the pages it needed) */
  readonly count: number;
  readonly drawn: number;
  readonly pages: number;
  /** pages missing, 0 if the surface got all it needed */
  readonly missingPages: number;
}

export interface ChunkedIsosurfaceResult {
  readonly surface: ChunkedSurface;
  readonly drawIndirect: DrawIndirect;
  /** the state after the extraction, read back (not needed to draw it) */
  readonly count: Promise<ChunkedSurfaceCount>;
}

function createTable(compute: ComputeContext, pages: number) {
  return compute.createBuffer((TableHeader + pages) * 4, 'isosurface page table');
}

/**
 * Extracts the isosurface of `input` into the arena (`surface` reuses its pages and draw arguments); `estimate`, the
 * expected vertex count, grows the arena beforehand when its free pages may not cover it
 */
export async function computeChunkedIsosurface(
  compute: ComputeContext,
  input: IsosurfaceInput,
  gridDim: Vec3,
  isoValue: number,
  transform: Mat4,
  invert: boolean,
  surface?: ChunkedSurface,
  estimate?: number,
  /** times the extraction's passes (timing mode) */
  timer?: GraphicsTimer,
): Promise<ChunkedIsosurfaceResult> {
  // compiled in the background, not blocking the queue at the first dispatch; awaited only when not ready, the
  // extraction then records its commands without awaiting
  const kernels = getIsosurfaceKernels(input);
  if (!compute.hasKernels(kernels)) await compute.prepareKernels(kernels);

  const arena = getVertexArena(compute);
  const estimatePages = estimate !== undefined ? Math.ceil(estimate / IsosurfaceChunkSize) : 0;
  const reserve = Math.max(0, estimatePages - (surface?.pages ?? 0));
  const available = getAvailablePages(arena);
  if (available < reserve + HeadroomPages) growVertexArena(compute, arena, reserve + HeadroomPages - available);
  arena.reserved += reserve;

  if (!surface) {
    surface = {
      table: createTable(compute, arena.pageCount),
      retired: [],
      drawArgs: compute.createBuffer(16, 'isosurface draw arguments'),
      pages: 0,
    };
  } else if (surface.table.byteLength < (TableHeader + arena.pageCount) * 4) {
    // a surface can take any page of the arena: its table grows with it
    const table = createTable(compute, arena.pageCount);
    surface.table.copyTo(table);
    // what draws the surface binds the old table until it rebinds (e.g. a re-extraction after the arena grew)
    surface.retired.push(surface.table);
    surface.table = table;
  }
  const { table, drawArgs } = surface;

  const [dx, dy, dz] = gridDim;
  const groups = (n: number, size: number) => Math.max(1, Math.ceil((n - 1) / size));
  const dispatch = [groups(dx, WorkgroupSize[0]), groups(dy, WorkgroupSize[1]), groups(dz, WorkgroupSize[2])] as const;
  const workgroupCount = dispatch[0] * dispatch[1] * dispatch[2];
  const { triCounts, triTable } = getShared(compute);
  const { schema: inputSchema, values: inputValues } = getInputValues(input, gridDim, isoValue);
  // kept between calls (see `getScratchBuffer`): nothing is awaited from here on
  const bases = getScratchBuffer(compute, 'mc chunk bases', workgroupCount * 4);

  const countKernel = compute.createKernel(McChunkCountKernel, { ...inputSchema, ...McChunkCountOutput }, {
    ...inputValues,
    triCounts: ValueCell.create(triCounts),
    table: ValueCell.create(table),
    bases: ValueCell.create(bases),
  } as never);
  const allocateKernel = compute.createKernel(McChunkAllocateKernel, McChunkAllocateSchema, {
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
    triCounts: ValueCell.create(triCounts),
    triTable: ValueCell.create(triTable),
    bases: ValueCell.create(bases),
    table: ValueCell.create(table),
    positions: ValueCell.create(arena.positions.ref.value),
    normals: ValueCell.create(arena.normals.ref.value),
  } as never);

  timer?.group('computeChunkedIsosurface');
  // the vertex count from 0, the pages of the previous extraction kept (reused, or freed by `mc-chunk-allocate`)
  table.clear(0, 4);
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
  const count = Promise.all([tableRead, freeRead]).then(([t, f]) => {
    const [vertices, pages, status, drawn] = new Uint32Array(t);
    arena.freeRead = new Uint32Array(f)[0];
    arena.reserved -= reserve;
    current.pages = pages;
    const needed = Math.ceil(vertices / IsosurfaceChunkSize);
    return { count: vertices, drawn, pages, missingPages: status ? Math.max(1, needed - pages) : 0 };
  });
  return { surface, drawIndirect: { buffer: drawArgs, byteOffset: 0 }, count };
}

/** Returns the pages of `surface` to the arena and destroys its buffers */
export function destroyChunkedSurface(compute: ComputeContext, surface: ChunkedSurface) {
  const arena = arenas.get(compute);
  if (arena) {
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
  }
  surface.table.destroy();
  releaseRetiredTables(surface);
  surface.drawArgs.destroy();
}

/**
 * The first `count` vertices of the surface of `table`, gathered from its pages and read back: positions (`w` 1),
 * normals and groups (packed rgba bytes) per vertex, as the texture readback of a texture mesh gives them
 */
export async function readChunkedVertices(compute: ComputeContext, table: StorageBuffer, count: number) {
  const arena = getVertexArena(compute);
  const header = new Uint32Array(await table.readAsync());
  const pages = Math.min(header[1], Math.ceil(count / IsosurfaceChunkSize));
  const ids = header.subarray(TableHeader, TableHeader + pages);
  const positionPage = IsosurfaceChunkSize * PositionBytes;
  const normalPage = IsosurfaceChunkSize * NormalBytes;
  const positions = compute.createBuffer(Math.max(1, pages) * positionPage, 'isosurface read positions');
  const normals = compute.createBuffer(Math.max(1, pages) * normalPage, 'isosurface read normals');
  try {
    for (let i = 0; i < pages; ++i) {
      arena.positions.ref.value.copyTo(positions, ids[i] * positionPage, i * positionPage, positionPage);
      arena.normals.ref.value.copyTo(normals, ids[i] * normalPage, i * normalPage, normalPage);
    }
    const [p, n] = await Promise.all([positions.readAsync(), normals.readAsync()]);
    const position = new Float32Array(p);
    const groupBits = new Uint8Array(p);
    const normal = new Uint32Array(n);
    const vertexCount = Math.min(count, pages * IsosurfaceChunkSize);
    const vertices = new Float32Array(vertexCount * 4);
    const outNormals = new Float32Array(vertexCount * 4);
    const groups = new Uint8Array(vertexCount * 4);
    const snorm = (v: number) => Math.max(-1, ((v << 16) >> 16) / 32767);
    for (let i = 0; i < vertexCount; ++i) {
      vertices[i * 4] = position[i * 4];
      vertices[i * 4 + 1] = position[i * 4 + 1];
      vertices[i * 4 + 2] = position[i * 4 + 2];
      vertices[i * 4 + 3] = 1;
      for (let k = 0; k < 4; ++k) groups[i * 4 + k] = groupBits[i * 16 + 12 + k];
      const xy = normal[i * 2];
      const z = normal[i * 2 + 1];
      outNormals[i * 4] = snorm(xy & 0xffff);
      outNormals[i * 4 + 1] = snorm(xy >>> 16);
      outNormals[i * 4 + 2] = snorm(z & 0xffff);
    }
    return { vertices, normals: outNormals, groups, vertexCount };
  } finally {
    positions.destroy();
    normals.destroy();
  }
}

/** Bytes of `pages` pages of the arena */
export function getChunkByteCount(pages: number) {
  return pages * IsosurfaceChunkSize * (PositionBytes + NormalBytes);
}

/** The arena of `compute`, if created: its pages, the free ones at the last readback, and their bytes */
export function getVertexArenaStats(compute: ComputeContext) {
  const arena = arenas.get(compute);
  if (!arena) return;
  const used = arena.pageCount - arena.freeRead;
  return {
    pages: arena.pageCount,
    usedPages: used,
    byteCount: getChunkByteCount(arena.pageCount),
    usedByteCount: getChunkByteCount(used),
  };
}

/** Pages a surface of `vertexCount` vertices takes */
export function getChunkPages(vertexCount: number) {
  return Math.ceil(vertexCount / IsosurfaceChunkSize);
}
