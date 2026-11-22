/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Box3D } from '@molstar/core/math/geometry/primitives/box3d';
import { GridPositions, type Positions } from '@molstar/core/math/geometry/grid-positions';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { Task } from '@molstar/core/task';
import { CustomPropertyDescriptor } from '../custom-property.js';
import { type ParticleList, ParticleValues } from './particle-list.js';
import {
  ParticleTileCatalog,
  ParticleTileGrid,
  type ParticleTile,
  type ParticleTileOpacity,
  type ParticleTileSet,
} from './particle-tiles.js';

// How opaque tiles are (`ParticleTile.opacity`), where the format does not say: offline, a pack writer's. By the
// content's density: molecules filling a volume fraction φ, of mean chord ℓ (4V/S), let a line of length L pass with
// probability exp(-φ L / ℓ) (randomly placed objects, the mean over lines: a line's own count of molecules is noise);
// per block of the cell its optical depth a length, φ / ℓ: its molecules' cross-sections S/4 over its volume. A
// molecule a solid sphere of its volume (its mass's, `ParticleTileCatalog.getTargetVolume`): fewer, larger objects
// than its atoms let more through, so the estimate errs towards passing. What is finer than a block (a hole, a
// membrane) is in its mean only. Rays gather the blocks' optical depths (`ParticleTileVisibility`): tiles that each
// let much through hide together what is behind them.

/** A molecule of a tile's content: a sphere at x, y, z of `radius` and `crossSection` (S/4: π r² for a sphere) */
export type ParticleTileMoleculeSink = (x: number, y: number, z: number, radius: number, crossSection: number) => void;

/** Optical densities (1/Å) on a grid of `dimensions` cells over `min`-`max`, x fastest */
export interface ParticleDensityGrid {
  readonly min: Vec3;
  readonly max: Vec3;
  readonly dimensions: readonly [number, number, number];
  readonly density: Float32Array;
  /** this grid's state of `density` (`updateDensityGrid` writes the same array in place, a new version each time) */
  readonly version?: number;
  /**
   * what an update changed (`updateDensityGrid`): the version it was made from and the boxes of cells it wrote,
   * `[x0, y0, z0, x1, y1, z1)` six a box. Why: a copy of the grid (an upload) kept up to date by those cells alone, not
   * all of it (a view's step: a few tiles of thousands)
   */
  readonly update?: { readonly from: number; readonly boxes: Int32Array };
  /** the tiles whose opacities it holds (`updateDensityGrid`: none left, no grid; kept by the changes alone) */
  readonly tileCount?: number;
}

let densityGridVersion = 0;

export interface ParticleTileOpaqueOptions {
  /** the blocks an axis: 4 by default */
  readonly blocks?: number;
}

const _p = Vec3();

/** of the particles' positions: a derived list's (a move's, a subset's) computed again */
const DensityDescriptor = CustomPropertyDescriptor({ name: 'particle-density-grid' });

export namespace ParticleDensityGrid {
  /** `list`'s optical densities over its box, known (a pack's): what hides what is behind as its whole is drawn */
  export function set(list: ParticleList, grid: ParticleDensityGrid) {
    list.customProperties.add(DensityDescriptor);
    list._propertyData[DensityDescriptor.name] = grid;
  }

  export function get(list: ParticleList): ParticleDensityGrid | undefined {
    return list._propertyData?.[DensityDescriptor.name];
  }

  /**
   * The optical densities of `list`'s particles (`ParticleTileOpaque.computeDensity`: a particle a sphere of its
   * target's volume, `volume`, the catalog's by default) over its particles' centers' box, `longest` blocks (64 by
   * default) along its longest side; undefined where none has a volume
   */
  export function compute(list: ParticleList, options: ComputeOptions = {}): ParticleDensityGrid | undefined {
    if (list.count === 0) return undefined;
    const box = ParticleTileOpaque.getCentersBox(list.coordinates, list.count);
    const { dimensions, volume } = prepare(list, box, options);
    const density = ParticleTileOpaque.computeDensity(box, dimensions, (sink) =>
      ParticleTileOpaque.forEachMolecule(list, volume, sink),
    );
    return density && { min: Vec3.clone(box.min), max: Vec3.clone(box.max), dimensions, density };
  }

  /** `compute` as a task: its particles in chunks between the task's updates (the main thread not held) */
  export function computeTask(list: ParticleList, options: ComputeOptions = {}): Task<ParticleDensityGrid | undefined> {
    return Task.create('Particle density grid', async (ctx) => {
      const { count, coordinates } = list;
      if (count === 0) return undefined;
      // the centers' box, then the molecules: a pass each, in chunks
      const box = Box3D.setEmpty(Box3D());
      const p = Vec3();
      for (let start = 0; start < count; start += ChunkSize) {
        if (ctx.shouldUpdate) await ctx.update({ message: 'Density grid', current: start, max: 2 * count });
        const end = Math.min(count, start + ChunkSize);
        for (let i = start; i < end; ++i) Box3D.add(box, GridPositions.get(coordinates, i, p));
      }
      const { dimensions, volume } = prepare(list, box, options);
      const density = ParticleTileOpaque.createDensity(box, dimensions);
      const targetOf = (i: number) => ParticleValues.get(list.targets, i);
      for (let start = 0; start < count; start += ChunkSize) {
        if (ctx.shouldUpdate) await ctx.update({ message: 'Density grid', current: count + start, max: 2 * count });
        const end = Math.min(count, start + ChunkSize);
        ParticleTileOpaque.forEachMoleculeOf(start, end, targetOf, coordinates, volume, density.sink);
      }
      const values = density.finish();
      return values && { min: Vec3.clone(box.min), max: Vec3.clone(box.max), dimensions, density: values };
    });
  }

  export interface ComputeOptions {
    /** the blocks along the box's longest side: 64 by default */
    readonly longest?: number;
    /** a target's volume: the catalog's by default */
    readonly volume?: (id: number) => number | undefined;
  }

  /** particles a chunk of `computeTask` (about 35 ms of molecules in Chrome) */
  const ChunkSize = 1 << 18;

  function prepare(list: ParticleList, box: Box3D, options: ComputeOptions) {
    const dimensions = ParticleTileOpaque.getGridDimensions(box, options.longest);
    let volume = options.volume;
    if (!volume) {
      const targets = ParticleTileCatalog.fromList(list).targets;
      volume = (id) => targets.get(id)?.volume;
    }
    return { dimensions, volume };
  }
}

export namespace ParticleTileOpaque {
  /**
   * The molecules of `list` (a tile's): a particle a sphere of its target's volume (`volume`) at its position; a
   * target without one nothing
   */
  export function forEachMolecule(
    list: ParticleList,
    volume: (id: number) => number | undefined,
    sink: ParticleTileMoleculeSink,
  ) {
    forEachMoleculeOf(0, list.count, (i) => ParticleValues.get(list.targets, i), list.coordinates, volume, sink);
  }

  /**
   * The molecules of particles `[start, end)` of `positions` (`targetOf` a particle's target): a particle a sphere of
   * its target's volume (`volume`) at its position (a pack's tile as written: `writeParticlePack`)
   */
  export function forEachMoleculeOf(
    start: number,
    end: number,
    targetOf: (i: number) => number,
    positions: Positions,
    volume: (id: number) => number | undefined,
    sink: ParticleTileMoleculeSink,
  ) {
    let id = -1;
    let radius = 0;
    let crossSection = 0;
    for (let i = start; i < end; ++i) {
      const t = targetOf(i);
      if (t !== id) {
        id = t;
        const v = volume(t) ?? 0;
        radius = Math.cbrt((3 * v) / (4 * Math.PI));
        crossSection = Math.PI * radius * radius;
      }
      if (!(radius > 0)) continue;
      GridPositions.get(positions, i, _p);
      sink(_p[0], _p[1], _p[2], radius, crossSection);
    }
  }

  /**
   * The opacity of `cell`'s content `forEach` gives (`forEachMolecule`): per block the cross-sections of the molecules
   * in it over its volume (1/Å): a molecule's shared by the blocks whose centers it covers, else its center's block's.
   * Undefined where every block is empty
   */
  export function computeOpacity(
    cell: Box3D,
    forEach: (sink: ParticleTileMoleculeSink) => void,
    options: ParticleTileOpaqueOptions = {},
  ): ParticleTileOpacity | undefined {
    const b = Math.max(1, options.blocks ?? 4);
    const density = computeDensity(cell, [b, b, b], forEach);
    return density && { dimensions: [b, b, b], density };
  }

  /**
   * The optical densities (1/Å) of what `forEach` gives on a grid of `dimensions` blocks over `box` (x fastest, as
   * `computeOpacity`'s); undefined where every block is empty
   */
  export function computeDensity(
    box: Box3D,
    dimensions: readonly [number, number, number],
    forEach: (sink: ParticleTileMoleculeSink) => void,
  ): Float32Array | undefined {
    const density = createDensity(box, dimensions);
    forEach(density.sink);
    return density.finish();
  }

  /**
   * `computeDensity` in parts: the molecules given to `sink` (as many calls as wanted, e.g. a list's in chunks between
   * a task's updates), the densities by `finish`
   */
  export function createDensity(box: Box3D, dimensions: readonly [number, number, number]) {
    const [bx, by, bz] = dimensions;
    const o = box.min;
    const s = [(box.max[0] - o[0]) / bx, (box.max[1] - o[1]) / by, (box.max[2] - o[2]) / bz];
    const sums = new Float64Array(bx * by * bz);
    const bxy = bx * by;
    let any = false;
    const sink: ParticleTileMoleculeSink = (x, y, z, radius, crossSection) => {
      // spread over the blocks whose centers it covers (a molecule larger than a block)
      const x0 = Math.max(0, Math.ceil((x - radius - o[0]) / s[0] - 0.5));
      const x1 = Math.min(bx - 1, Math.floor((x + radius - o[0]) / s[0] - 0.5));
      const y0 = Math.max(0, Math.ceil((y - radius - o[1]) / s[1] - 0.5));
      const y1 = Math.min(by - 1, Math.floor((y + radius - o[1]) / s[1] - 0.5));
      const z0 = Math.max(0, Math.ceil((z - radius - o[2]) / s[2] - 0.5));
      const z1 = Math.min(bz - 1, Math.floor((z + radius - o[2]) / s[2] - 0.5));
      const r2 = radius * radius;
      let covered = 0;
      for (let pass = 0; pass < 2; ++pass) {
        const share = pass === 1 ? crossSection / covered : 0;
        for (let bk = z0; bk <= z1; ++bk) {
          const dz = o[2] + (bk + 0.5) * s[2] - z;
          for (let bj = y0; bj <= y1; ++bj) {
            const dy = o[1] + (bj + 0.5) * s[1] - y;
            for (let bi = x0; bi <= x1; ++bi) {
              const dx = o[0] + (bi + 0.5) * s[0] - x;
              if (dx * dx + dy * dy + dz * dz > r2) continue;
              if (pass === 0) covered += 1;
              else sums[bi + bj * bx + bk * bxy] += share;
            }
          }
        }
        if (covered === 0) break;
      }
      if (covered > 0) {
        any = true;
        return;
      }
      // smaller than a block's reach: all in the block of its center
      const i = Math.floor((x - o[0]) / s[0]);
      const j = Math.floor((y - o[1]) / s[1]);
      const k = Math.floor((z - o[2]) / s[2]);
      // a center on the box's far face in its last block; outside the box (a trajectory's frame) not counted
      const ic = i === bx && x <= box.max[0] ? bx - 1 : i;
      const jc = j === by && y <= box.max[1] ? by - 1 : j;
      const kc = k === bz && z <= box.max[2] ? bz - 1 : k;
      if (ic < 0 || ic >= bx || jc < 0 || jc >= by || kc < 0 || kc >= bz) return;
      sums[ic + jc * bx + kc * bxy] += crossSection;
      any = true;
    };
    const finish = (): Float32Array | undefined => {
      if (!any) return undefined;
      const volume = s[0] * s[1] * s[2];
      const density = new Float32Array(sums.length);
      for (let i = 0; i < sums.length; ++i) density[i] = sums[i] / volume;
      return density;
    };
    return { sink, finish };
  }

  /** The box of the centers of particles `[0, count)` of `positions` (a density grid's: `ParticleDensityGrid`) */
  export function getCentersBox(positions: Positions, count: number): Box3D {
    const box = Box3D.setEmpty(Box3D());
    for (let i = 0; i < count; ++i) Box3D.add(box, GridPositions.get(positions, i, _p));
    return box;
  }

  /**
   * A grid's blocks an axis over `box`: `longest` along its longest side, the others as many of about that size
   * (near-cubic blocks, at least one)
   */
  export function getGridDimensions(box: Box3D, longest = 64): [number, number, number] {
    const e = Vec3.sub(Vec3(), box.max, box.min);
    const side = Math.max(e[0], e[1], e[2]) / longest;
    if (!(side > 0)) return [1, 1, 1];
    return [0, 1, 2].map((k) => Math.max(1, Math.min(longest, Math.round(e[k] / side)))) as [number, number, number];
  }

  /**
   * The optical densities (1/Å) of `tiles` of `set` as one grid over the set's grid: a tile's opacity blocks at its
   * cell's place (`ParticleTileOpacity`, of the same blocks an axis all), none (0) where no tile of them is; undefined
   * where none has one. What the GPU's occlusion takes as what hides what is behind (the depth pyramid primed from it)
   */
  export function getDensityGrid(set: ParticleTileSet, tiles: Iterable<number>): ParticleDensityGrid | undefined {
    let b = 0;
    const withOpacity: ParticleTile[] = [];
    for (const id of tiles) {
      const tile = set.tiles.get(id);
      const opacity = tile?.opacity;
      if (!tile || !opacity) continue;
      const [x, y, z] = opacity.dimensions;
      if (x !== y || y !== z) return undefined;
      if (b === 0) b = x;
      else if (x !== b) return undefined;
      withOpacity.push(tile);
    }
    if (b === 0) return undefined;
    const [gx, gy, gz] = set.grid.dimensions;
    const nx = gx * b;
    const ny = gy * b;
    const nz = gz * b;
    const density = new Float32Array(nx * ny * nz);
    for (const tile of withOpacity) writeTileDensity(set, density, b, tile, false);
    const { origin, cellSize } = set.grid;
    return {
      min: Vec3.clone(origin),
      max: Vec3.create(origin[0] + cellSize[0] * gx, origin[1] + cellSize[1] * gy, origin[2] + cellSize[2] * gz),
      dimensions: [nx, ny, nz],
      density,
      version: ++densityGridVersion,
      tileCount: withOpacity.length,
    };
  }

  /** `tile`'s opacity blocks (`b` an axis) at its cell's place in `density` (a grid over `set`'s), or 0 (`clear`) */
  function writeTileDensity(
    set: ParticleTileSet,
    density: Float32Array,
    b: number,
    tile: ParticleTile,
    clear: boolean,
    /** the box of cells written appended (six a box) */
    boxes?: number[],
  ) {
    const [gx, gy] = set.grid.dimensions;
    const nx = gx * b;
    const ny = gy * b;
    const x = tile.id % gx;
    const y = Math.floor(tile.id / gx) % gy;
    const z = Math.floor(tile.id / (gx * gy));
    boxes?.push(x * b, y * b, z * b, (x + 1) * b, (y + 1) * b, (z + 1) * b);
    const d = tile.opacity!.density;
    for (let bk = 0; bk < b; ++bk) {
      for (let bj = 0; bj < b; ++bj) {
        const to = x * b + (y * b + bj) * nx + (z * b + bk) * nx * ny;
        if (clear) {
          density.fill(0, to, to + b);
          continue;
        }
        const from = bj * b + bk * b * b;
        for (let bi = 0; bi < b; ++bi) density[to + bi] = d[from + bi];
      }
    }
  }

  /**
   * `getDensityGrid` of `tiles` from `last` (the grid of `lastTiles` of `set`): the tiles let go cleared and those added
   * written in place, not a grid made again (a following view's steps: a few tiles each of a grid of millions of
   * cells); a new object (its users tell a change by identity), its array the last's. Made again where `last` is none
   * or of other blocks
   */
  export function updateDensityGrid(
    set: ParticleTileSet,
    tiles: ArrayLike<number>,
    last?: { readonly grid: ParticleDensityGrid; readonly tiles: ArrayLike<number> },
  ): ParticleDensityGrid | undefined {
    if (!last) return getDensityGrid(set, Array.from(tiles));
    const [gx, gy, gz] = set.grid.dimensions;
    const [nx, ny, nz] = last.grid.dimensions;
    const b = nx / gx;
    if (!Number.isInteger(b) || ny !== gy * b || nz !== gz * b) return getDensityGrid(set, Array.from(tiles));
    const next = new Set<number>();
    for (let i = 0; i < tiles.length; ++i) next.add(tiles[i]);
    const before = new Set<number>();
    for (let i = 0; i < last.tiles.length; ++i) before.add(last.tiles[i]);
    const withOpacity = (id: number) => {
      const opacity = set.tiles.get(id)?.opacity;
      return opacity && opacity.dimensions.every((d) => d === b) ? set.tiles.get(id)! : undefined;
    };
    // the tiles added and let go alone (a step: a few of thousands shown), the count of those with opacities kept
    const added: number[] = [];
    for (const id of next) {
      if (before.has(id)) continue;
      const tile = set.tiles.get(id);
      if (tile?.opacity && !withOpacity(id)) return getDensityGrid(set, next);
      added.push(id);
    }
    const { density } = last.grid;
    const boxes: number[] = [];
    let count = last.grid.tileCount;
    for (const id of before) {
      if (next.has(id)) continue;
      const tile = withOpacity(id);
      if (!tile) continue;
      writeTileDensity(set, density, b, tile, true, boxes);
      if (count !== undefined) count -= 1;
    }
    for (const id of added) {
      const tile = withOpacity(id);
      if (!tile) continue;
      writeTileDensity(set, density, b, tile, false, boxes);
      if (count !== undefined) count += 1;
    }
    // a grid of no count (made elsewhere): its tiles counted
    if (count === undefined) {
      count = 0;
      for (const id of next) if (withOpacity(id)) count += 1;
    }
    if (count === 0) return undefined;
    // a grid of no version (made elsewhere): no update to go by, a copy of it made again whole
    const from = last.grid.version;
    return {
      min: last.grid.min,
      max: last.grid.max,
      dimensions: last.grid.dimensions,
      density,
      version: ++densityGridVersion,
      update: from !== undefined ? { from, boxes: Int32Array.from(boxes) } : undefined,
      tileCount: count,
    };
  }

  /**
   * `set` with its tiles' opacities computed (`computeOpacity` of each tile's cell, a trajectory's each block's least
   * over its frames), its tiles read from its provider (a task run), the targets' volumes from its catalog
   */
  export function computeOpacities(
    set: ParticleTileSet,
    options: ParticleTileOpaqueOptions = {},
  ): Task<ParticleTileSet> {
    return Task.create('Particle tile opacities', async (ctx) => {
      const resolve = async <T>(v: T | Task<T>) => (Task.is(v) ? await v.runInContext(ctx) : v);
      const volume = (id: number) => set.catalog.targets.get(id)?.volume;
      const tiles = new Map<number, ParticleTile>();
      const cell = Box3D();
      let done = 0;
      for (const tile of set.tiles.values()) {
        if (ctx.shouldUpdate) await ctx.update({ message: 'Opacities', current: done, max: set.tiles.size });
        ParticleTileGrid.getBox(set.grid, tile.id, cell);
        let opacity: ParticleTileOpacity | undefined;
        for (let f = 0; f < set.frameCount; ++f) {
          const list = await resolve(set.provider.getTile(tile.id, f));
          const next = computeOpacity(cell, (sink) => forEachMolecule(list, volume, sink), options);
          if (!next) {
            opacity = undefined;
            break;
          }
          if (opacity)
            for (let i = 0; i < next.density.length; ++i)
              next.density[i] = Math.min(next.density[i], opacity.density[i]);
          opacity = next;
        }
        tiles.set(tile.id, { ...tile, opacity });
        done += 1;
      }
      return { ...set, tiles };
    });
  }
}
