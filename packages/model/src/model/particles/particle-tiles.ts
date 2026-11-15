/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Box3D } from '@molstar/core/math/geometry/primitives/box3d';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { BoundaryHelper } from '@molstar/core/math/geometry/boundary-helper';
import { GridPositions } from '@molstar/core/math/geometry/grid-positions';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { Task } from '@molstar/core/task';
import { getSharedSerial } from '@molstar/core/util/array';
import { CustomPropertyDescriptor } from '../custom-property.js';
import type { Cell } from '@molstar/core/math/geometry/spacegroup/cell';
import {
  Particle,
  type ParticleCompartmentInfo,
  type ParticleEntityInfo,
  type ParticleList,
  ParticleSites,
  ParticleTarget,
  type ParticleTargetInfo,
  ParticleValues,
} from './particle-list.js';
import type { ParticleTrajectory } from './particle-trajectory.js';
import { Grid } from '../volume/grid.js';
import { DefaultAtomWeight, ElementAtomWeights } from '../structure/model/properties/atomic/measures.js';

// Tiles of a particle list: a mesoscale model subdivided on a grid, its tiles and its targets'
// data loaded on demand; what tiles share (the targets, their info, entities and compartments) known up front and one
// object each for all tiles, so whatever is keyed by a target (its visual, its geometry, its GPU data) is one.

/**
 * What is known of a target before its data is loaded (`ParticleTileCatalog`): its kind, and what decides the programs
 * that draw it (sites quantized, their elements known: the precompile's signatures); its radius and its data's bytes
 * where known
 */
export interface ParticleTileTargetInfo {
  readonly id: number;
  readonly kind: ParticleTarget['kind'];
  readonly sites?: { readonly quantized: boolean; readonly elements: boolean };
  /**
   * how far its data reaches from the point a particle places (`ParticleTarget.getCenter`), in Å: a tile's reach past
   * its box and its size on screen (`ParticleTileVisibility`); none for a shape (placed by its geometry)
   */
  readonly radius?: number;
  /**
   * its molecule's volume, Å³ (`ParticleTileCatalog.getTargetVolume`): how much of what it covers it fills, so how
   * opaque its particles are together (`ParticleTileOpaque`); none for a shape or a volume
   */
  readonly volume?: number;
  readonly byteCount?: number;
}

/**
 * What the tiles of a tile set share, known up front and light (kilobytes for thousands of targets): per target its
 * kind (`targets`) and info (`targetInfo`, as a list's), the entities and compartments, the cell; the targets' data
 * on demand (`ParticleTileProvider.getTarget`)
 */
export interface ParticleTileCatalog {
  readonly targets: ReadonlyMap<number, ParticleTileTargetInfo>;
  readonly targetInfo: ReadonlyMap<number, ParticleTargetInfo>;
  readonly entityInfo?: ReadonlyMap<number, ParticleEntityInfo>;
  readonly compartmentInfo?: ReadonlyMap<number, ParticleCompartmentInfo>;
  readonly cell?: Cell;
}

/** A target's reach and volume where known before (`ParticleTileTargetInfo.radius`, `volume`): a pack's */
export interface ParticleTargetMeasures {
  readonly radius?: number;
  readonly volume?: number;
}

export namespace ParticleTileCatalog {
  /** of the targets' data: a derived list's (a move's, a subset's) computed again */
  const MeasuresDescriptor = CustomPropertyDescriptor({ name: 'particle-target-measures' });

  /** `list`'s targets' reach and volume, known (a pack's): the catalog's without a pass over their data */
  export function setTargetMeasures(list: ParticleList, measures: ReadonlyMap<number, ParticleTargetMeasures>) {
    list.customProperties.add(MeasuresDescriptor);
    list._propertyData[MeasuresDescriptor.name] = measures;
  }

  /** The catalog of `list`'s targets: their info, entities and compartments the list's own objects */
  export function fromList(list: ParticleList): ParticleTileCatalog {
    const targets = new Map<number, ParticleTileTargetInfo>();
    const measures = list._propertyData?.[MeasuresDescriptor.name] as
      | ReadonlyMap<number, ParticleTargetMeasures>
      | undefined;
    for (const id of list.targetInfo.keys()) {
      const target = list.targetMapping?.get(id);
      if (target) targets.set(id, getTargetInfo(id, target, measures?.get(id)));
    }
    return {
      targets,
      targetInfo: list.targetInfo,
      entityInfo: list.entityInfo,
      compartmentInfo: getCompartmentEntities(list),
      cell: list.cell,
    };
  }

  /**
   * `list`'s compartments' info with each one's entities (in `entityInfo`'s order) and whether some of its particles
   * have none (`ParticleCompartmentInfo.entities`): a view of some tiles colors by them as the whole; the info as it is
   * without compartments
   */
  export function getCompartmentEntities(list: ParticleList): ReadonlyMap<number, ParticleCompartmentInfo> | undefined {
    const { compartmentInfo, compartments, entities, entityInfo, count } = list;
    if (!compartmentInfo || !compartments) return compartmentInfo;
    const byCompartment = new Map<number, Set<number>>();
    const withoutEntity = new Set<number>();
    // a pair seen once a flag of a table where the ids are small (a pack's: its runs 2.29 M in the yeast tile-ordered,
    // its pairs a few hundred), the sets' adds a pair, not a run
    let maxC = -1;
    let maxE = -1;
    for (const c of compartmentInfo.keys()) if (c > maxC) maxC = c;
    if (entityInfo) for (const e of entityInfo.keys()) if (e > maxE) maxE = e;
    const stride = maxE + 2;
    const seen = maxC >= 0 && (maxC + 1) * stride <= 1 << 22 ? new Uint8Array((maxC + 1) * stride) : undefined;
    const add = (c: number, e: number) => {
      if (c < 0) return;
      if (seen && c <= maxC && e <= maxE) {
        const k = c * stride + (e < 0 ? 0 : e + 1);
        if (seen[k]) return;
        seen[k] = 1;
      }
      if (e < 0) {
        withoutEntity.add(c);
        return;
      }
      let set = byCompartment.get(c);
      if (!set) byCompartment.set(c, (set = new Set()));
      set.add(e);
    };
    if (entities) ParticleValues.forEachRun2(compartments, entities, count, (_s, _e, c, e) => add(c, e));
    else ParticleValues.forEachRun(compartments, count, (_s, _e, c) => add(c, -1));
    // an entity's place in the info's order, those it does not hold after them
    const rank = new Map<number, number>();
    entityInfo?.forEach((_, e) => rank.set(e, rank.size));
    const rankOf = (e: number) => rank.get(e) ?? rank.size + e;
    const out = new Map<number, ParticleCompartmentInfo>();
    for (const [c, info] of compartmentInfo) {
      const ids = [...(byCompartment.get(c) ?? [])].sort((a, b) => rankOf(a) - rankOf(b));
      out.set(c, { ...info, entities: ids, withoutEntity: withoutEntity.has(c) });
    }
    return out;
  }

  /** What the catalog knows of `target`: its reach and volume `known`'s where given, else of its data */
  export function getTargetInfo(
    id: number,
    target: ParticleTarget,
    known?: ParticleTargetMeasures,
  ): ParticleTileTargetInfo {
    const radius = known?.radius ?? getTargetRadius(target);
    const volume = known?.volume ?? getTargetVolume(target);
    if (target.kind !== 'sites') return { id, kind: target.kind, radius, volume };
    const { sites } = target;
    const q = sites.quantized;
    return {
      id,
      kind: 'sites',
      sites: { quantized: !!q, elements: ParticleSites.hasElements(sites) },
      radius,
      volume,
      // 8 bytes a quantized site, else 16 (a position and a radius)
      byteCount: q ? q.data.byteLength : sites.count * 16,
    };
  }

  const _p = Vec3();
  /**
   * How far `target`'s data reaches from the point a particle places: its sites' farthest sphere, a structure's bounding
   * sphere around it, a volume's box; none for a shape
   */
  export function getTargetRadius(target: ParticleTarget): number | undefined {
    if (target.kind === 'shape') return undefined;
    const center = ParticleTarget.getCenter(target)!;
    if (target.kind === 'sites') {
      const { sites } = target;
      let radius = 0;
      for (let i = 0; i < sites.count; ++i) {
        const d = Vec3.distance(ParticleSites.getPosition(sites, i, _p), center) + ParticleSites.getRadius(sites, i);
        if (d > radius) radius = d;
      }
      return radius;
    }
    const sphere =
      target.kind === 'structure' ? target.structure.boundary.sphere : Grid.getBoundingSphere(target.volume.grid);
    return Vec3.distance(sphere.center, center) + sphere.radius;
  }

  /** Å³ a dalton of protein: its partial specific volume, 0.73 cm³/g (nucleic acids' 0.55: theirs overestimated) */
  export const VolumePerDalton = 1.21;
  /** The part of its sites' spheres a molecule of sites without elements fills (coarse beads overlap) */
  export const SitesPacking = 0.7;

  /**
   * The volume of `target`'s molecule, Å³: by its mass where its atoms' elements are known (their weights, a hydrogen a
   * heavy atom where it has none, times `VolumePerDalton`: a molecule's atoms overlap, so not their spheres' sum), else
   * its sites' spheres' volume times `SitesPacking`; a structure's coarse elements not counted; none for a shape or
   * a volume
   */
  export function getTargetVolume(target: ParticleTarget): number | undefined {
    if (target.kind === 'shape' || target.kind === 'volume') return undefined;
    let mass = 0;
    let heavy = 0;
    let hydrogens = 0;
    const add = (z: number) => {
      mass += ElementAtomWeights[z] ?? DefaultAtomWeight;
      if (z === 1) hydrogens += 1;
      else heavy += 1;
    };
    if (target.kind === 'sites') {
      const { sites } = target;
      if (!ParticleSites.hasElements(sites)) {
        let v = 0;
        for (let i = 0; i < sites.count; ++i) v += (4 / 3) * Math.PI * ParticleSites.getRadius(sites, i) ** 3;
        return v * SitesPacking;
      }
      for (let i = 0; i < sites.count; ++i) add(ParticleSites.getElement(sites, i) || 6);
    } else {
      for (const unit of target.structure.units) {
        // atomic units only (Unit.Kind.Atomic)
        if (unit.kind !== 0) continue;
        const { atomicNumber } = unit.model.atomicHierarchy.derived.atom;
        const { elements } = unit;
        for (let i = 0, il = elements.length; i < il; ++i) add(atomicNumber[elements[i]] || 6);
      }
    }
    if (hydrogens === 0) mass += heavy * (ElementAtomWeights[1] ?? 1);
    return mass * VolumePerDalton;
  }
}

/**
 * A regular grid of tiles: `dimensions` cells of `cellSize` from `origin`; a tile a cell, its id the cell's linear
 * index (`x + y nx + z nx ny`)
 */
export interface ParticleTileGrid {
  readonly origin: Vec3;
  readonly cellSize: Vec3;
  readonly dimensions: readonly [number, number, number];
}

export namespace ParticleTileGrid {
  /** The grid of `dimensions` cells over the box `min`-`max` (an axis of no extent a cell of size 1) */
  export function create(min: Vec3, max: Vec3, dimensions: readonly [number, number, number]): ParticleTileGrid {
    const cellSize = Vec3();
    for (let a = 0; a < 3; ++a) {
      const extent = max[a] - min[a];
      cellSize[a] = extent > 0 ? extent / dimensions[a] : 1;
    }
    return { origin: Vec3.clone(min), cellSize, dimensions: [dimensions[0], dimensions[1], dimensions[2]] };
  }

  export function cellCount(grid: ParticleTileGrid) {
    const [nx, ny, nz] = grid.dimensions;
    return nx * ny * nz;
  }

  /** The id of cell `x`, `y`, `z` */
  export function getId(grid: ParticleTileGrid, x: number, y: number, z: number) {
    const [nx, ny] = grid.dimensions;
    return x + y * nx + z * nx * ny;
  }

  /** The id of the cell `position` is in, clamped to the grid (a position on its far face in the last cell) */
  export function getIdOf(grid: ParticleTileGrid, x: number, y: number, z: number) {
    const { origin, cellSize, dimensions: d } = grid;
    const cx = Math.min(d[0] - 1, Math.max(0, Math.floor((x - origin[0]) / cellSize[0])));
    const cy = Math.min(d[1] - 1, Math.max(0, Math.floor((y - origin[1]) / cellSize[1])));
    const cz = Math.min(d[2] - 1, Math.max(0, Math.floor((z - origin[2]) / cellSize[2])));
    return getId(grid, cx, cy, cz);
  }

  /** The box of cell `id` into `out` */
  export function getBox(grid: ParticleTileGrid, id: number, out: Box3D): Box3D {
    const [nx, ny] = grid.dimensions;
    const x = id % nx;
    const y = Math.floor(id / nx) % ny;
    const z = Math.floor(id / (nx * ny));
    const { origin, cellSize } = grid;
    Vec3.set(out.min, origin[0] + x * cellSize[0], origin[1] + y * cellSize[1], origin[2] + z * cellSize[2]);
    Vec3.set(out.max, out.min[0] + cellSize[0], out.min[1] + cellSize[1], out.min[2] + cellSize[2]);
    return out;
  }
}

/**
 * A tile's opacity (`ParticleTileOpaque.computeOpacity`): its grid cell in `dimensions` blocks (x fastest), per block
 * the optical depth a length (1/Å) every axis-aligned line through it reaches at least: a ray through it passes with
 * probability about exp(-density × its path in it)
 */
export interface ParticleTileOpacity {
  readonly dimensions: readonly [number, number, number];
  readonly density: Float32Array;
}

/** A tile of a tile set: a cell of the grid with particles */
export interface ParticleTile {
  readonly id: number;
  /** where its particles are: its cell's box, a trajectory's where they go in any frame (neighbours' may overlap) */
  readonly box: Box3D;
  readonly count: number;
  /**
   * the index in the whole of its first particle (the counts of the tiles before it: of smaller ids, a tile-ordered
   * list's in its order, `ParticleTileDirectory`): its particles' keys
   */
  readonly base: number;
  /** the ids of the targets its particles are of, ascending */
  readonly targets: Int32Array;
  /**
   * how opaque its content is, by blocks of its cell (a trajectory's for all frames): what hides the tiles behind it
   * together with others' (`ParticleTileVisibility`)
   */
  readonly opacity?: ParticleTileOpacity;
  /** an estimate of a frame's bytes */
  readonly byteCount?: number;
}

/**
 * What serves a tile set's data, sync where in memory, a task where fetched or decoded: a tile at a frame (its
 * particles, their targets by the catalog's ids, no targets' data: `targetMapping` none; its frames the same particles
 * in the same order, other poses, of one line of `Particle.withPoses`), and a target's data
 */
export interface ParticleTileProvider {
  getTile(id: number, frame: number): ParticleList | Task<ParticleList>;
  getTarget(id: number): ParticleTarget | Task<ParticleTarget>;
}

/**
 * The coarse levels of a tile set (an octree over its grid, a pack's `tiles.levels`): level `L` (1 to `count`) a node
 * a cell of the grid of `dimensions >> L` cells an axis, a tile of the set of its own: every `8^L`-th of its `2^L`³
 * tiles' particles (copies of theirs, spread over it), drawn larger by `2^L` (`ParticleList.scale`: an eighth of them
 * a level fill about the space all do). A view shows a node in place of its tiles where they are far
 */
export interface ParticleTileLevels {
  /** the finest tiles alone (the grid's cells): what a camera sees (`ParticleTileVisibility`) */
  readonly base: ParticleTileSet;
  readonly count: number;
  /** the tile id of level `level`'s node at `cell` of its grid (`x + y d + z d²`, `d = dimensions >> level`), none where it has none */
  getNode(level: number, cell: number): number | undefined;
}

/**
 * A particle list subdivided on a grid: the catalog, the grid, the tiles with particles (by id,
 * ascending), the frames (1 without a trajectory) and what serves them; a particle in the tile of its center (of the
 * first frame), in one tile only
 */
export interface ParticleTileSet {
  readonly label?: string;
  readonly catalog: ParticleTileCatalog;
  readonly grid: ParticleTileGrid;
  readonly tiles: ReadonlyMap<number, ParticleTile>;
  /** the particles of all tiles */
  readonly count: number;
  /**
   * the sphere its particles are all within (of the first frame; its tiles' bounds' union): a view's, whatever tiles it
   * shows (its object's bounds, the camera's reset to the whole)
   */
  readonly boundingSphere?: Sphere3D;
  readonly frameCount: number;
  readonly provider: ParticleTileProvider;
  /** its tiles levels of another's (`ParticleTiles.withLevels`): their ids not the grid's cells */
  readonly levels?: ParticleTileLevels;
}

/** a box of nothing, its min above its max */
const EmptyBox = [Infinity, Infinity, Infinity, -Infinity, -Infinity, -Infinity];

/**
 * The tiles a list holds its particles in (a tile-ordered pack's, `ParticleTiles.setDirectory`): its particles by tile,
 * the tiles in Morton order of their cells, a tile's by target (its targets a few runs of the list's); a tile a cell of
 * a grid of `dimensions` cells an axis (a power of two) over `min`-`max` (its particles' box). A tile set of that grid,
 * or of a coarser power of two (a coarse tile's cells consecutive in Morton order), made from it without a pass over
 * the particles, a tile a range of the list
 */
export interface ParticleTileDirectory {
  readonly min: Vec3;
  readonly max: Vec3;
  readonly dimensions: number;
  /** per tile with particles, in the list's order: its cell's id (`ParticleTileGrid.getId`), first particle, count */
  readonly ranges: Int32Array;
  /** per tile its particles' centers' box: min x, y, z, max x, y, z */
  readonly boxes: Float32Array;
  /**
   * per tile, in `ranges`' order, its opacity where known (a pack's, computed when written): `blocks`³ optical
   * densities a tile, x fastest (`ParticleTileOpacity`); a coarser set's tile its cells' blocks side by side
   */
  readonly opacity?: { readonly blocks: number; readonly density: Float32Array };
}

/** Values a tile takes in `ParticleTileDirectory.ranges` and `boxes` */
export const ParticleTileDirectoryRange = 3;
export const ParticleTileDirectoryBox = 6;

export namespace ParticleTiles {
  export interface FromListOptions {
    /** the tile set's label, the list's by default */
    readonly label?: string;
    /** a tile and a target served after this many ms (a task), as a fetch would be: sync by default */
    readonly delay?: number;
  }

  /** by position: not kept by a move (`Particle.move`), its tiles computed again */
  const DirectoryDescriptor = CustomPropertyDescriptor({ name: 'particle-tile-directory' });

  /** `list`'s particles ordered by `directory`'s tiles (a tile-ordered pack's) */
  export function setDirectory(list: ParticleList, directory: ParticleTileDirectory) {
    list.customProperties.add(DirectoryDescriptor);
    list._propertyData[DirectoryDescriptor.name] = directory;
  }

  export function getDirectory(list: ParticleList): ParticleTileDirectory | undefined {
    return list._propertyData?.[DirectoryDescriptor.name];
  }

  /**
   * The tile set of `list` on a grid of `dimensions` cells over its particles' box; its tiles and targets served from
   * the list, a tile's particles a view of a range of it: the list's own where it is tile-ordered for that grid
   * (`getDirectory`: its tiles' or a power of two coarser), else the list ordered by tile once
   */
  export function fromList(
    list: ParticleList,
    dimensions: readonly [number, number, number],
    options: FromListOptions = {},
  ): ParticleTileSet {
    const directory = getDirectory(list);
    const set = directory ? fromDirectory(list, directory, dimensions, options) : undefined;
    return set ?? fromFrames(list, 1, () => list, dimensions, options);
  }

  /**
   * The tile set of `trajectory` on a grid of `dimensions` cells over its first frame's particles' box, the particles
   * in the tiles of their centers there; a tile's box where its particles go in the frames given sync, its frames
   * served from the trajectory's
   */
  /** particles a tile the automatic grid aims at: a page of a page store (`StorePageSize`), what a tile's load writes */
  export const AutoTileParticles = 8192;

  /**
   * The cells along each axis of a tile set of `list` that draw it best: a tile-ordered list's written grid (its tiles
   * views of its ranges: another grid orders it again, a pass and a copy of every particle), else the power of two
   * whose cells hold about `AutoTileParticles` particles each (1-64)
   */
  export function getAutoCells(list: ParticleList) {
    const directory = getDirectory(list);
    if (directory) return directory.dimensions;
    const n = Math.cbrt(Math.max(1, list.count) / AutoTileParticles);
    return Math.min(64, Math.max(1, 2 ** Math.round(Math.log2(Math.max(1, n)))));
  }

  export function fromTrajectory(
    trajectory: ParticleTrajectory,
    dimensions: readonly [number, number, number],
    options: FromListOptions = {},
  ): ParticleTileSet {
    const first = trajectory.getFrameAtIndex(0);
    const list = Task.is(first) ? trajectory.representative : first;
    return fromFrames(list, trajectory.frameCount, (i) => trajectory.getFrameAtIndex(i), dimensions, options);
  }

  /**
   * The tile set of tile-ordered `list` from its directory where `dimensions` are its grid's or a power of two coarser
   * on every axis (a coarse tile its cells' ranges, consecutive), else undefined: no pass over the particles, no copy
   */
  function fromDirectory(
    list: ParticleList,
    directory: ParticleTileDirectory,
    dimensions: readonly [number, number, number],
    options: FromListOptions,
  ): ParticleTileSet | undefined {
    const n = directory.dimensions;
    const d = dimensions[0];
    if (dimensions[1] !== d || dimensions[2] !== d || !(d >= 1) || d > n || n % d !== 0) return;
    const ratio = n / d;
    if ((ratio & (ratio - 1)) !== 0) return;
    const shift = Math.log2(ratio);
    const grid = ParticleTileGrid.create(directory.min, directory.max, [d, d, d]);
    const cells = ParticleTileGrid.cellCount(grid);
    const firsts = new Int32Array(cells).fill(-1);
    const counts = new Int32Array(cells);
    const centers = new Float64Array(cells * 6);
    for (let t = 0; t < cells; ++t) centers.set(EmptyBox, t * 6);
    const { ranges, boxes } = directory;
    const tileCount = ranges.length / ParticleTileDirectoryRange;
    // a coarse tile's opacity its cells' blocks side by side: (ratio × blocks)³
    const opacity = directory.opacity;
    const opacities = opacity ? new Map<number, ParticleTileOpacity>() : undefined;
    const b = opacity?.blocks ?? 0;
    const bb = b * ratio;
    let end = 0;
    for (let k = 0; k < tileCount; ++k) {
      const fine = ranges[k * 3];
      const first = ranges[k * 3 + 1];
      const count = ranges[k * 3 + 2];
      const x = (fine % n) >> shift;
      const y = (Math.floor(fine / n) % n) >> shift;
      const z = Math.floor(fine / (n * n)) >> shift;
      const id = ParticleTileGrid.getId(grid, x, y, z);
      // a coarse tile's cells one after another, its range theirs (else not of this grid: ordered again)
      if (firsts[id] < 0) firsts[id] = first;
      else if (firsts[id] + counts[id] !== first) return;
      if (first !== end) return;
      end = first + count;
      counts[id] += count;
      const c = id * 6;
      const o = k * ParticleTileDirectoryBox;
      for (let a = 0; a < 3; ++a) {
        if (boxes[o + a] < centers[c + a]) centers[c + a] = boxes[o + a];
        if (boxes[o + 3 + a] > centers[c + 3 + a]) centers[c + 3 + a] = boxes[o + 3 + a];
      }
      if (opacity && opacities) {
        let tile = opacities.get(id);
        if (!tile) opacities.set(id, (tile = { dimensions: [bb, bb, bb], density: new Float32Array(bb * bb * bb) }));
        // the cell's place in its coarse tile, its blocks there
        const lx = (fine % n) - x * ratio;
        const ly = (Math.floor(fine / n) % n) - y * ratio;
        const lz = Math.floor(fine / (n * n)) - z * ratio;
        const from = k * b * b * b;
        for (let bk = 0; bk < b; ++bk) {
          for (let bj = 0; bj < b; ++bj) {
            for (let bi = 0; bi < b; ++bi) {
              const to = lx * b + bi + (ly * b + bj) * bb + (lz * b + bk) * bb * bb;
              tile.density[to] = opacity.density[from + bi + bj * b + bk * b * b];
            }
          }
        }
      }
    }
    if (end !== list.count) return;
    return createSet(list, list, 1, () => list, grid, firsts, counts, centers, undefined, options, opacities);
  }

  function fromFrames(
    list: ParticleList,
    frameCount: number,
    getFrame: (i: number) => ParticleList | Task<ParticleList>,
    dimensions: readonly [number, number, number],
    options: FromListOptions,
  ): ParticleTileSet {
    const n = list.count;
    const bounds = new Float64Array(7);
    Particle.computeRangeBounds(list.coordinates, 0, n, bounds);
    const grid = ParticleTileGrid.create(
      Vec3.create(bounds[0], bounds[1], bounds[2]),
      Vec3.create(bounds[3], bounds[4], bounds[5]),
      dimensions,
    );

    // each particle's tile by its center, the particles of a tile in their order (a counting sort)
    const p = Vec3();
    const tileOf = new Int32Array(n);
    const counts = new Int32Array(ParticleTileGrid.cellCount(grid));
    // per tile its particles' centers' box (min, max: six a tile), in the same pass: its bounds without a pass of its own
    const centers = new Float64Array(counts.length * 6);
    for (let t = 0; t < counts.length; ++t) centers.set(EmptyBox, t * 6);
    for (let i = 0; i < n; ++i) {
      GridPositions.get(list.coordinates, i, p);
      const id = ParticleTileGrid.getIdOf(grid, p[0], p[1], p[2]);
      tileOf[i] = id;
      counts[id] += 1;
      const b = id * 6;
      if (p[0] < centers[b]) centers[b] = p[0];
      if (p[1] < centers[b + 1]) centers[b + 1] = p[1];
      if (p[2] < centers[b + 2]) centers[b + 2] = p[2];
      if (p[0] > centers[b + 3]) centers[b + 3] = p[0];
      if (p[1] > centers[b + 4]) centers[b + 4] = p[1];
      if (p[2] > centers[b + 5]) centers[b + 5] = p[2];
    }
    const firsts = new Int32Array(counts.length);
    for (let t = 0, first = 0; t < counts.length; ++t) {
      firsts[t] = first;
      first += counts[t];
    }
    const order = new Int32Array(n);
    const next = firsts.slice();
    for (let i = 0; i < n; ++i) order[next[tileOf[i]]++] = i;
    const orders = new Map<number, Int32Array>();
    for (let id = 0; id < counts.length; ++id) {
      if (counts[id] > 0) orders.set(id, order.subarray(firsts[id], firsts[id] + counts[id]));
    }

    // the list in tile order (a stable sort: by target within a tile where the list is by target, a pack's), gathered
    // once: a tile's particles a range of it, its list a view (`Particle.range`), not a gather a tile
    const ordered = Particle.subset(list, order, { label: options.label ?? list.label });
    return createSet(list, ordered, frameCount, getFrame, grid, firsts, counts, centers, orders, options);
  }

  /**
   * the distinct targets of the particles `[base, base + count)`, ascending: by their runs (a pack tile's a run a
   * target, by target: ascending already, no set or sort; the yeast's 2.29 M runs), sorted and made distinct otherwise
   */
  function getRangeTargets(targets: ParticleValues<Int32Array>, base: number, count: number) {
    let out = new Int32Array(16);
    let k = 0;
    let ascending = true;
    ParticleValues.forEachRun(ParticleValues.slice(targets, base, count), count, (_, __, t) => {
      if (k > 0 && t <= out[k - 1]) ascending = false;
      if (k === out.length) {
        const grown = new Int32Array(k * 2);
        grown.set(out);
        out = grown;
      }
      out[k++] = t;
    });
    if (ascending) return out.slice(0, k);
    const sorted = out.subarray(0, k).sort();
    let d = 0;
    for (let i = 0; i < k; ++i) if (d === 0 || sorted[i] !== sorted[d - 1]) sorted[d++] = sorted[i];
    return sorted.slice(0, d);
  }

  /**
   * The tile set of `list` whose particles are `ordered`'s (tile-ordered: tile `id`'s `[firsts[id], + counts[id])`,
   * its centers' box `centers`' six); `orders`: a tile's particles' indices in `list` (its frames'), none where
   * `ordered` is `list`
   */
  function createSet(
    list: ParticleList,
    ordered: ParticleList,
    frameCount: number,
    getFrame: (i: number) => ParticleList | Task<ParticleList>,
    grid: ParticleTileGrid,
    firsts: Int32Array,
    counts: Int32Array,
    centers: Float64Array,
    orders: Map<number, Int32Array> | undefined,
    options: FromListOptions,
    opacities?: ReadonlyMap<number, ParticleTileOpacity>,
  ): ParticleTileSet {
    const n = list.count;
    const p = Vec3();
    let boundingSphere = getTilesSphere(list, centers, counts);
    const bytesPerParticle =
      (GridPositions.is(list.coordinates) ? 8 : 12) +
      (list.rotations ? list.rotations.BYTES_PER_ELEMENT * 4 : 0) +
      (list.attributes ? list.attributes.size * 4 : 0);
    const label = options.label ?? list.label;
    const tiles = new Map<number, ParticleTile>();
    for (let id = 0; id < counts.length; ++id) {
      const count = counts[id];
      if (count === 0) continue;
      const base = firsts[id];
      const box = ParticleTileGrid.getBox(grid, id, Box3D());
      tiles.set(id, {
        id,
        box,
        count,
        base,
        targets: getRangeTargets(ordered.targets, base, count),
        opacity: opacities?.get(id),
        byteCount: count * bytesPerParticle,
      });
    }
    // a trajectory's tile where its particles go: in every frame given sync
    if (frameCount > 1 && orders) {
      for (let f = 0; f < frameCount; ++f) {
        const frame = getFrame(f);
        if (Task.is(frame)) continue;
        for (const [id, tile] of tiles) {
          const tileOrder = orders.get(id)!;
          for (let i = 0; i < tileOrder.length; ++i) {
            GridPositions.get(frame.coordinates, tileOrder[i], p);
            Box3D.add(tile.box, p);
          }
        }
      }
      // the set's sphere of its tiles' boxes (where their particles go in every frame), not of the first frame's centers
      const boxes = new Float64Array(counts.length * 6);
      for (const [id, tile] of tiles) {
        boxes.set(tile.box.min, id * 6);
        boxes.set(tile.box.max, id * 6 + 3);
      }
      boundingSphere = getTilesSphere(list, boxes, counts);
    }

    // a tile's list a view of the ordered list's range (kept: the same object a request), its keys its particles'
    // indices in tile order (`base + i`: a view of the shared serial array); its frames the same particles posed
    let serial: ReturnType<typeof getSharedSerial> | undefined;
    const firstLists = new Map<number, ParticleList>();
    const getFirst = (id: number) => {
      let first = firstLists.get(id);
      if (!first) {
        const tile = tiles.get(id);
        if (!tile) throw new Error(`particle tiles: no tile ${id}`);
        serial ??= getSharedSerial(n);
        const keys = new Int32Array(serial.buffer, serial.byteOffset + tile.base * 4, tile.count);
        const view = Particle.range(ordered, tile.base, tile.count, {
          label: `${label ?? 'Particles'} tile ${id}`,
          keys,
        });
        first = { ...view, targetMapping: undefined };
        // its bounds: its centers' box padded by its largest radius (by its radii's runs), not a pass over it
        let radius = 0;
        if (first.radii) {
          ParticleValues.forEachRun(first.radii, first.count, (_, __, r) => {
            radius = Math.max(radius, r);
          });
        }
        const b = id * 6;
        const box = Box3D.create(
          Vec3.create(centers[b] - radius, centers[b + 1] - radius, centers[b + 2] - radius),
          Vec3.create(centers[b + 3] + radius, centers[b + 4] + radius, centers[b + 5] + radius),
        );
        const sphere = Sphere3D.create(
          Vec3.create(
            (centers[b] + centers[b + 3]) / 2,
            (centers[b + 1] + centers[b + 4]) / 2,
            (centers[b + 2] + centers[b + 5]) / 2,
          ),
          Math.hypot(centers[b + 3] - centers[b], centers[b + 4] - centers[b + 1], centers[b + 5] - centers[b + 2]) /
            2 +
            radius,
        );
        Particle.setBoundary(first, { box, sphere });
        firstLists.set(id, first);
      }
      return first;
    };
    const tileAt = (id: number, frame: ParticleList) => {
      const first = getFirst(id);
      if (frame === list) return first;
      const posed = Particle.subset(frame, orders!.get(id)!);
      return Particle.withPoses(first, posed.coordinates, posed.rotations);
    };
    const later = <T>(name: string, value: () => T | Task<T>): T | Task<T> => {
      const { delay } = options;
      if (!delay) return value();
      return Task.create(name, async (ctx) => {
        await new Promise((resolve) => setTimeout(resolve, delay));
        const v = value();
        return Task.is(v) ? await v.runInContext(ctx) : v;
      });
    };

    const provider: ParticleTileProvider = {
      getTile: (id, f) =>
        later(`Particle tile ${id}`, () => {
          if (f < 0 || f >= frameCount) throw new Error(`particle tiles: no frame ${f}`);
          const frame = getFrame(f);
          if (!Task.is(frame)) return tileAt(id, frame);
          return Task.create(`Particle tile ${id}`, async (ctx) => tileAt(id, await frame.runInContext(ctx)));
        }),
      getTarget: (id) =>
        later(`Particle target ${id}`, () => {
          const target = list.targetMapping?.get(id);
          if (!target) throw new Error(`particle tiles: no target ${id}`);
          return target;
        }),
    };

    return {
      label,
      catalog: ParticleTileCatalog.fromList(list),
      grid,
      tiles,
      count: n,
      boundingSphere,
      frameCount,
      provider,
    };
  }

  /** a node's scale at level `level` (`ParticleTileLevels`): its particles drawn larger by `2^level` */
  export function getLevelScale(level: number) {
    return 2 ** level;
  }

  /** the cell of base tile `id` (of a grid of `dimensions` cells an axis) in level `level`'s grid */
  export function getLevelCell(dimensions: number, id: number, level: number) {
    const d = dimensions >> level;
    const x = (id % dimensions) >> level;
    const y = (Math.floor(id / dimensions) % dimensions) >> level;
    const z = Math.floor(id / (dimensions * dimensions)) >> level;
    return x + y * d + z * d * d;
  }
}

/**
 * The sphere of the tiles' spheres (each its particles' centers' box's, `centers`: six a tile, padded by the largest
 * radius of `list`'s particles, by its radii's runs): not a pass over the particles
 */
function getTilesSphere(list: ParticleList, centers: Float64Array, counts: Int32Array) {
  let radius = 0;
  if (list.radii) {
    ParticleValues.forEachRun(list.radii, list.count, (_, __, r) => {
      radius = Math.max(radius, r);
    });
  }
  const helper = new BoundaryHelper('98');
  const sphere = Sphere3D();
  const tileSphere = (t: number) => {
    const b = t * 6;
    Vec3.set(
      sphere.center,
      (centers[b] + centers[b + 3]) / 2,
      (centers[b + 1] + centers[b + 4]) / 2,
      (centers[b + 2] + centers[b + 5]) / 2,
    );
    sphere.radius =
      Math.hypot(centers[b + 3] - centers[b], centers[b + 4] - centers[b + 1], centers[b + 5] - centers[b + 2]) / 2 +
      radius;
    return sphere;
  };
  helper.reset();
  for (let t = 0; t < counts.length; ++t) if (counts[t] > 0) helper.includeSphere(tileSphere(t));
  helper.finishedIncludeStep();
  for (let t = 0; t < counts.length; ++t) if (counts[t] > 0) helper.radiusSphere(tileSphere(t));
  return helper.getSphere(Sphere3D());
}
