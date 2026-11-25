/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Box3D } from '@molstar/core/math/geometry/primitives/box3d';
import { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { Clip } from '@molstar/core/util/clip';
import { type ParticleTile, ParticleTileGrid, type ParticleTileSet } from './particle-tiles.js';

// Which tiles of a tile set would show something to a camera, so that only those are loaded: in view, not sub-pixel,
// not hidden behind the tiles in front of them (their content's optical depth gathered along the rays). A tile wrongly
// kept rather than one wrongly left out: the densities a lower bound, a pixel opaque only where all its corners' rays
// are, at the deepest back of what made them so.

/** A camera as numbers (the model is below graphics) */
export interface ParticleTileCamera {
  /** the eye */
  readonly position: Vec3;
  /** where it looks, unit */
  readonly direction: Vec3;
  /** world to clip space */
  readonly viewProjection: Mat4;
  readonly orthographic?: boolean;
  /**
   * the pixels an Å makes at a distance of 1 (the projection's focal scale times half the viewport's height), or
   * anywhere where orthographic
   */
  readonly pixelScale: number;
  /** the near plane's distance: nothing nearer is drawn, so nothing nearer hides, 0 by default */
  readonly near?: number;
}

export interface ParticleTileVisibilityOptions {
  /** the smallest size on screen worth drawing, in pixels (a tile's largest target's radius): 1 by default */
  readonly minPixels?: number;
  /** whether tiles hide those behind them, true by default */
  readonly occlusion?: boolean;
  /** the optical depth a ray gathers to be stopped: 3 by default (under 5% passes) */
  readonly opticalDepth?: number;
  /** the depth buffer's pixels a side, 64 by default */
  readonly resolution?: number;
  /**
   * the times a density block a clip surface cuts is split in eight (`ClipSplitDepth` by default): fewer, a quicker
   * pick hiding less near the surface (a cut being moved)
   */
  readonly clipSplitDepth?: number;
  /**
   * what draws the tiles (a representation each): none, one drawing them whole and opaque. A tile shows something
   * where any keeps some of its reach; what hides is what every opaque one keeps (`ParticleTileDrawing`)
   */
  readonly drawings?: readonly ParticleTileDrawing[];
  /** whether a tile is drawn opaque (no entity of it hidden, no transparency): only those hide; all by default */
  readonly isOpaque?: (tile: ParticleTile) => boolean;
}

/**
 * What draws a set's tiles: its clip objects (`Clip.Objects`, any shapes; a point clipped by any, none: drawn whole),
 * whether it draws them opaque (only then do they hide). Its kept part of a density block hides: cut exactly by the
 * planes, by the other shapes a block kept whole, left out whole, or split in eight and each tested again
 * (`ClipSplitDepth` times; a part cut by a surface still then hides nothing: a tile wrongly kept, not one left out)
 */
export interface ParticleTileDrawing {
  readonly clip?: Clip.Objects;
  readonly opaque: boolean;
}

/** a block cut by a clip surface split in eight at most this many times (a part a 512th of the block at 3) */
const ClipSplitDepth = 3;

/** a pick's candidates by slot (`ParticleTileVisibility.get`): arrays grown to the largest set seen, reused */
interface CandidateScratch {
  readonly index: Int32Array;
  readonly size: Float64Array;
  readonly distance: Float64Array;
  readonly scale: Float64Array;
  readonly exposure: Float64Array;
  readonly depth: Float64Array;
  readonly kept: Uint32Array;
  /** slots `0..count-1` in order, in an array of its own of the scratch (two at once: the candidates, their order) */
  slots(count: number): Uint32Array;
}

let candidateScratch: (CandidateScratch & { readonly capacity: number }) | undefined;
function getCandidateScratch(n: number): CandidateScratch {
  if (candidateScratch && candidateScratch.capacity >= n) return candidateScratch;
  const capacity = Math.max(n, 1024);
  const lists = [new Uint32Array(capacity), new Uint32Array(capacity)];
  let next = 0;
  candidateScratch = {
    capacity,
    index: new Int32Array(capacity),
    size: new Float64Array(capacity),
    distance: new Float64Array(capacity),
    scale: new Float64Array(capacity),
    exposure: new Float64Array(capacity),
    depth: new Float64Array(capacity),
    kept: new Uint32Array(capacity),
    slots(count: number) {
      const out = lists[next].subarray(0, count);
      next = 1 - next;
      for (let i = 0; i < count; ++i) out[i] = i;
      return out;
    },
  };
  return candidateScratch;
}

/** A tile that would show something */
export interface ParticleTileVisible {
  readonly id: number;
  /** its largest target's radius on screen, in pixels (Infinity: the camera within its reach, or a radius unknown) */
  readonly size: number;
  /** the camera's distance to its reach (0 inside) */
  readonly distance: number;
  /**
   * the pixels an Å makes at its box (the camera's distance to its cell, not its reach: the same for tiles alike far
   * whatever they hold), Infinity within it
   */
  readonly scale: number;
  /**
   * how much of it is seen: the light through the tiles in front of it to its pixels, the mean of exp(-their optical
   * depth), 1 with nothing in front (or without occlusion)
   */
  readonly exposure: number;
}

export interface ParticleTileVisibility {
  /**
   * by priority: the least occluded first (`exposure`, by tenths: what fills the screen before what shows through
   * gaps), then larger on screen, then nearer
   */
  readonly tiles: ParticleTileVisible[];
  /** the tiles left out, by why */
  readonly culled: {
    readonly frustum: number;
    readonly clipped: number;
    readonly size: number;
    readonly occluded: number;
  };
}

const _block = Box3D();
const _sphere = Sphere3D();
const _cell = Box3D();
const _reach = Box3D();
const NoDrawings: readonly ParticleTileDrawing[] = [];
const NoBits: readonly number[] = [];
const _clip = new Float64Array(4);
const _invVP = Mat4();
const _ndc = Vec3();
const _frustum = Frustum3D();

export namespace ParticleTileVisibility {
  /** The tiles of `set` that would show something to `camera` (its frame's boxes: a trajectory's cover all frames) */
  export function get(
    set: ParticleTileSet,
    camera: ParticleTileCamera,
    options: ParticleTileVisibilityOptions = {},
  ): ParticleTileVisibility {
    const minPixels = options.minPixels ?? 1;
    const drawings = options.drawings ?? [WholeOpaque];
    const occluders = drawings.filter((d) => d.opaque);
    const occlusion = (options.occlusion ?? true) && occluders.length > 0;
    // what hides: the parts every opaque drawing keeps, their planes cut exactly, the other shapes by classification
    const planes = occluders.flatMap((d) => (d.clip ? getClipPlanes(d.clip) : []));
    const shaped = occluders.filter((d) => d.clip && hasClipShapes(d.clip));
    const shapedPlaneBits = shaped.map((d) => getClipPlaneBits(d.clip!));
    const splitDepth = options.clipSplitDepth ?? ClipSplitDepth;
    const { position: eye, direction: dir, pixelScale } = camera;
    Frustum3D.fromProjectionMatrix(_frustum, camera.viewProjection);
    const culled = { frustum: 0, clipped: 0, size: 0, occluded: 0 };

    // in view and on screen, by the reach: the box grown by the largest radius of the tile's targets. The candidates by
    // their index in the set's tiles, their values in arrays reused across calls, their reach read again where tested:
    // no object nor box made a tile (thousands in view, a pick every few frames while moving)
    const { tiles: tileList, reaches, radii: tileRadius } = getTileReaches(set);
    const scratch = getCandidateScratch(tileList.length);
    const { index, size: sizes, distance: distances, scale: scales, exposure: exposures, depth: depths } = scratch;
    const loadReach = (k: number, out: Box3D) => {
      const b = k * 6;
      Vec3.set(out.min, reaches[b], reaches[b + 1], reaches[b + 2]);
      Vec3.set(out.max, reaches[b + 3], reaches[b + 4], reaches[b + 5]);
      return out;
    };
    // a reach clipped whole only where every drawing clips
    const allClip = drawings.every((d) => !!d.clip);
    let count = 0;
    for (let k = 0; k < tileList.length; ++k) {
      const tile = tileList[k];
      const radius = tileRadius[k];
      const reach = loadReach(k, _reach);
      if (!Frustum3D.intersectsBox3D(_frustum, reach)) {
        culled.frustum += 1;
        continue;
      }
      // all of its reach clipped by every drawing: nothing of it shows
      if (allClip && drawings.every((d) => clipsWhole(d.clip!, reach, 0))) {
        culled.clipped += 1;
        continue;
      }
      const distance = distanceToBox(reach, eye);
      let size = Infinity;
      if (!Number.isNaN(radius)) {
        if (camera.orthographic) size = radius * pixelScale;
        else if (distance > 0) size = (radius * pixelScale) / distance;
      }
      if (size < minPixels) {
        culled.size += 1;
        continue;
      }
      const boxDistance = camera.orthographic ? 1 : distanceToBox(tile.box, eye);
      index[count] = k;
      sizes[count] = size;
      distances[count] = distance;
      scales[count] = camera.orthographic ? pixelScale : boxDistance > 0 ? pixelScale / boxDistance : Infinity;
      exposures[count] = 1;
      count += 1;
    }

    // the visible candidates (their slots), in their order
    let visible: ArrayLike<number> = scratch.slots(count);
    let opaqueSome = false;
    for (let c = 0; c < count && !opaqueSome; ++c) if (tileList[index[c]].opacity) opaqueSome = true;
    if (occlusion && opaqueSome) {
      const buffer = OpticalBuffer.get(options.resolution ?? 64, camera, options.opticalDepth ?? 3);
      const blocks = getTileBlocks(set);
      // front to back: each tile tested against what is in front of it, then its blocks added where it is kept
      const order = scratch.slots(count);
      for (let c = 0; c < count; ++c) depths[c] = nearestDepth(loadReach(index[c], _reach), eye, dir);
      order.sort((x, y) => depths[x] - depths[y] || distances[x] - distances[y] || x - y);
      const kept = scratch.kept;
      let keptCount = 0;
      for (let o = 0; o < count; ++o) {
        const c = order[o];
        const tile = tileList[index[c]];
        const reach = loadReach(index[c], _reach);
        // its pixels once for its tests (the buffer unchanged between them)
        const r = buffer.rect(reach);
        if (buffer.hidesRect(r, depths[c])) {
          culled.occluded += 1;
          continue;
        }
        exposures[c] = buffer.exposureRect(r);
        kept[keptCount++] = c;
        const tileBlocks = blocks.get(tile.id);
        if (!tileBlocks || (options.isOpaque && !options.isOpaque(tile))) continue;
        // nothing to add where every ray through it is opaque already
        if (r && buffer.isOpaqueRect(r)) continue;
        // its cell against the shapes first: kept whole or clipped away whole by them, its blocks not tested one by one
        // (a shape about a whole cell, thousands of tiles in view)
        const tileClass = classifyBox(shaped, shapedPlaneBits, ParticleTileGrid.getBox(set.grid, tile.id, _cell));
        if (tileClass === Clip.SphereClass.All) continue;
        // no clip object cutting it: its blocks by the rays' walk through its cell (a block's own of a uniform one)
        const opacity = tile.opacity;
        if (planes.length === 0 && tileClass === Clip.SphereClass.None && tileBlocks.length > 1 && opacity) {
          buffer.addGrid(_cell, opacity.dimensions, opacity.density);
          continue;
        }
        const byBlock = tileClass === Clip.SphereClass.Some ? shaped : NoDrawings;
        const byBlockBits = tileClass === Clip.SphereClass.Some ? shapedPlaneBits : NoBits;
        for (const { box, density } of tileBlocks) {
          if (clipBox(box, planes, _block)) addKept(buffer, _block, density, byBlock, byBlockBits, splitDepth, 0);
        }
      }
      visible = kept.subarray(0, keptCount);
    }

    const tiles: ParticleTileVisible[] = new Array(visible.length);
    for (let v = 0; v < visible.length; ++v) {
      const c = visible[v];
      tiles[v] = {
        id: tileList[index[c]].id,
        size: sizes[c],
        distance: distances[c],
        scale: scales[c],
        exposure: exposures[c],
      };
    }
    const tenths = (e: number) => Math.round(e * 10);
    tiles.sort(
      (a, b) => tenths(b.exposure) - tenths(a.exposure) || b.size - a.size || a.distance - b.distance || a.id - b.id,
    );
    return { tiles, culled };
  }

  /**
   * How `box` lies against the shapes (not the planes) of `drawings`: `All` where one clips all of it, `None` where
   * none clips any of it, `Some` else (by its bounding sphere: `Clip.classifySphere`)
   */
  function classifyBox(drawings: readonly ParticleTileDrawing[], planeBits: readonly number[], box: Box3D) {
    let some = false;
    for (let i = 0; i < drawings.length; ++i) {
      const c = Clip.classifySphere(drawings[i].clip!, boxSphere(box, _sphere), planeBits[i]);
      if (c === Clip.SphereClass.All) return c;
      if (c === Clip.SphereClass.Some) some = true;
    }
    return some ? Clip.SphereClass.Some : Clip.SphereClass.None;
  }

  /**
   * `box` (of `density`, cut by the planes already) added where every one of `drawings` keeps it from its shapes other
   * than planes: whole where none clips any of it, not where one clips all, else split in eight and each tested again
   * (`ClipSplitDepth` times, then left out: a part a surface cuts hides nothing)
   */
  function addKept(
    buffer: OpticalBuffer,
    box: Box3D,
    density: number,
    drawings: readonly ParticleTileDrawing[],
    planeBits: readonly number[],
    maxDepth: number,
    depth: number,
  ) {
    // none of the rays through it (off screen), or all opaque already: nothing to add, nor to test
    const r = buffer.rect(box);
    if (r && (r[0] > r[2] || r[1] > r[3] || buffer.isOpaqueRect(r))) return;
    // its size on screen (not the pixels it touches: a part across a pixel's border touches two however small)
    const extent = r ? buffer.extent : Infinity;
    let whole = true;
    for (let i = 0; i < drawings.length; ++i) {
      const c = Clip.classifySphere(drawings[i].clip!, boxSphere(box, _sphere), planeBits[i]);
      if (c === Clip.SphereClass.All) return;
      if (c === Clip.SphereClass.Some) whole = false;
    }
    if (whole) {
      buffer.add(box, density, r);
      return;
    }
    // split no finer than the rays: a part across a ray at most changes none it would not
    if (depth >= maxDepth || extent <= 1) return;
    const parts = getOctants(depth);
    for (let k = 0; k < 8; ++k) {
      addKept(buffer, getOctant(box, k, parts[k]), density, drawings, planeBits, maxDepth, depth + 1);
    }
  }

  /**
   * `box` cut by `planes` into `out`: per plane shrunk along its normal's largest axis until all of it is on the kept
   * side (for a given coordinate on that axis, the other two at their worst); false where nothing is left
   */
  export function clipBox(box: Box3D, planes: readonly Plane3D[], out: Box3D): boolean {
    Box3D.copy(out, box);
    for (const plane of planes) {
      const n = plane.normal;
      let a = 0;
      if (Math.abs(n[1]) > Math.abs(n[a])) a = 1;
      if (Math.abs(n[2]) > Math.abs(n[a])) a = 2;
      if (n[a] === 0) continue;
      // the plane's distance at the other axes' worst corner, less the axis' own term
      let rest = plane.constant;
      for (let b = 0; b < 3; ++b) {
        if (b !== a) rest += Math.min(n[b] * out.min[b], n[b] * out.max[b]);
      }
      // n[a] x + rest >= 0
      const bound = -rest / n[a];
      if (n[a] > 0) out.min[a] = Math.max(out.min[a], bound);
      else out.max[a] = Math.min(out.max[a], bound);
      if (out.min[a] >= out.max[a]) return false;
    }
    return true;
  }
}

/**
 * Whether `clip` clips all of `box`: a plane with all of it on its clipped side (its farthest corner), else by
 * classification of its bounding sphere, where some of it is clipped split in eight and each tested again
 * (`ClipSplitDepth` times; still some there: not all)
 */
function clipsWhole(clip: Clip.Objects, box: Box3D, depth: number): boolean {
  for (const plane of getClipPlanes(clip)) {
    const n = plane.normal;
    const corner = Vec3.create(
      n[0] > 0 ? box.max[0] : box.min[0],
      n[1] > 0 ? box.max[1] : box.min[1],
      n[2] > 0 ? box.max[2] : box.min[2],
    );
    if (Plane3D.distanceToPoint(plane, corner) < 0) return true;
  }
  const c = Clip.classifySphere(clip, boxSphere(box, Sphere3D()));
  if (c !== Clip.SphereClass.Some) return c === Clip.SphereClass.All;
  if (depth >= ClipSplitDepth) return false;
  for (let k = 0; k < 8; ++k) if (!clipsWhole(clip, getOctant(box, k, Box3D()), depth + 1)) return false;
  return true;
}

/** per split depth its eight octants' boxes (`addKept`'s, a depth's in use until its parts are done) */
const _octants: Box3D[][] = [];
function getOctants(depth: number) {
  while (_octants.length <= depth) _octants.push(Array.from({ length: 8 }, () => Box3D()));
  return _octants[depth];
}

/** octant `k` of `box` (bit a of `k` its upper half along axis a) into `out` */
function getOctant(box: Box3D, k: number, out: Box3D) {
  for (let a = 0; a < 3; ++a) {
    const mid = (box.min[a] + box.max[a]) / 2;
    const high = (k >> a) & 1;
    out.min[a] = high ? mid : box.min[a];
    out.max[a] = high ? box.max[a] : mid;
  }
  return out;
}

/** a set's tiles drawn whole and opaque (no drawing given) */
const WholeOpaque: ParticleTileDrawing = { opaque: true };

/** `box`'s bounding sphere into `out` */
function boxSphere(box: Box3D, out: Sphere3D) {
  Vec3.set(out.center, (box.min[0] + box.max[0]) / 2, (box.min[1] + box.max[1]) / 2, (box.min[2] + box.max[2]) / 2);
  out.radius = Vec3.distance(box.min, box.max) / 2;
  return out;
}

const isPlane = (clip: Clip.Objects, i: number) => clip.type[i] === Clip.Type.plane;

/** the bits of `clip`'s planes (`Clip.classifySphere`'s `flag`: ignored there, cut exactly by `getClipPlanes`) */
function getClipPlaneBits(clip: Clip.Objects) {
  let bits = 0;
  for (let i = 0; i < clip.count; ++i) if (isPlane(clip, i)) bits |= 1 << i;
  return bits;
}

/** whether `clip` has shapes other than planes (none switched off) */
function hasClipShapes(clip: Clip.Objects) {
  for (let i = 0; i < clip.count; ++i) if (!isPlane(clip, i) && clip.type[i] !== Clip.Type.none) return true;
  return false;
}

/**
 * `clip`'s planes in the world, what is kept where `Plane3D.distanceToPoint` is at least 0: a plane clips where its
 * frame's y is at least 0 (`Clip.getSignedDistance`: -y), inverted where below
 */
function getClipPlanes(clip: Clip.Objects): Plane3D[] {
  const out: Plane3D[] = [];
  for (let i = 0; i < clip.count; ++i) {
    if (!isPlane(clip, i)) continue;
    const f = i * 16;
    // the frame's y of a point: the matrix's second row
    const s = clip.invert[i] ? 1 : -1;
    const n = Vec3.create(s * clip.frame[f + 1], s * clip.frame[f + 5], s * clip.frame[f + 9]);
    const length = Vec3.magnitude(n);
    if (!(length > 0)) continue;
    out.push(Plane3D.create(Vec3.scale(n, n, 1 / length), (s * clip.frame[f + 13]) / length));
  }
  return out;
}

/**
 * the largest radius of `set`'s targets (those known; undefined where none is): what a set's detail is by, the same
 * for all its tiles (a tile's own largest would make tiles alike far differ in coarseness by what they hold)
 */
export function getParticleTileSetRadius(set: ParticleTileSet): number | undefined {
  let radius: number | undefined;
  for (const r of getTileRadii(set).values()) {
    if (r !== undefined && (radius === undefined || r > radius)) radius = r;
  }
  return radius;
}

/** per tile set (its tiles and catalog) its tiles' largest target radius, undefined where one is unknown */

const tileRadii = new WeakMap<ParticleTileSet['tiles'], Map<number, number | undefined>>();
function getTileRadii(set: ParticleTileSet) {
  let radii = tileRadii.get(set.tiles);
  if (radii) return radii;
  radii = new Map();
  for (const tile of set.tiles.values()) {
    let radius: number | undefined = 0;
    for (const t of tile.targets) {
      const r = set.catalog.targets.get(t)?.radius;
      if (r === undefined) {
        radius = undefined;
        break;
      }
      if (r > radius) radius = r;
    }
    radii.set(tile.id, radius);
  }
  tileRadii.set(set.tiles, radii);
  return radii;
}

/**
 * per tile set its tiles in an array, their reach boxes (a tile's box grown by its targets' largest radius, six a
 * tile) and radii (`getTileRadii`, NaN: none known)
 */
const tileReaches = new WeakMap<
  ParticleTileSet['tiles'],
  { tiles: ParticleTile[]; reaches: Float64Array; radii: Float64Array }
>();
function getTileReaches(set: ParticleTileSet) {
  let r = tileReaches.get(set.tiles);
  if (r) return r;
  const radii = getTileRadii(set);
  const tiles = [...set.tiles.values()];
  const reaches = new Float64Array(tiles.length * 6);
  const radius = new Float64Array(tiles.length);
  for (let k = 0; k < tiles.length; ++k) {
    const { id, box } = tiles[k];
    const value = radii.get(id);
    radius[k] = value ?? Number.NaN;
    const grow = value ?? 0;
    for (let a = 0; a < 3; ++a) {
      reaches[k * 6 + a] = box.min[a] - grow;
      reaches[k * 6 + 3 + a] = box.max[a] + grow;
    }
  }
  r = { tiles, reaches, radii: radius };
  tileReaches.set(set.tiles, r);
  return r;
}

interface Block {
  readonly box: Box3D;
  readonly density: number;
}

/** the factor within which a tile's blocks' densities are one, its least (the spread of a cytoplasm's: about 1.2) */
const UniformFactor = 1.5;

/** per tile set its tiles' blocks of some density (their boxes in the world), from their opacity and cells */
const tileBlocks = new WeakMap<ParticleTileSet['tiles'], Map<number, Block[]>>();
function getTileBlocks(set: ParticleTileSet) {
  let blocks = tileBlocks.get(set.tiles);
  if (blocks) return blocks;
  blocks = new Map();
  const cell = Box3D();
  for (const tile of set.tiles.values()) {
    const opacity = tile.opacity;
    if (!opacity) continue;
    ParticleTileGrid.getBox(set.grid, tile.id, cell);
    const [bx, by, bz] = opacity.dimensions;
    const size = [(cell.max[0] - cell.min[0]) / bx, (cell.max[1] - cell.min[1]) / by, (cell.max[2] - cell.min[2]) / bz];
    const list: Block[] = [];
    for (let k = 0; k < bz; ++k) {
      for (let j = 0; j < by; ++j) {
        for (let i = 0; i < bx; ++i) {
          const density = opacity.density[i + j * bx + k * bx * by];
          if (!(density > 0)) continue;
          const v = [i, j, k];
          const box = Box3D();
          for (let a = 0; a < 3; ++a) {
            box.min[a] = cell.min[a] + v[a] * size[a];
            // the cell's own face at its far end: neighbours' blocks meet without a gap
            box.max[a] = v[a] + 1 === opacity.dimensions[a] ? cell.max[a] : cell.min[a] + (v[a] + 1) * size[a];
          }
          list.push({ box, density });
        }
      }
    }
    if (!list.length) continue;
    // all of a density within a factor (a cytoplasm's): one box of the least, the cell's
    let least = Infinity;
    let most = 0;
    for (const b of list) {
      least = Math.min(least, b.density);
      most = Math.max(most, b.density);
    }
    const uniform = list.length === opacity.density.length && most <= UniformFactor * least;
    blocks.set(tile.id, uniform ? [{ box: Box3D.clone(cell), density: least }] : list);
  }
  tileBlocks.set(set.tiles, blocks);
  return blocks;
}

function distanceToBox(box: Box3D, p: Vec3) {
  let d = 0;
  for (let a = 0; a < 3; ++a) {
    const v = p[a] < box.min[a] ? box.min[a] - p[a] : p[a] > box.max[a] ? p[a] - box.max[a] : 0;
    d += v * v;
  }
  return Math.sqrt(d);
}

/** the box's nearest depth along `dir` from `eye` (a linear function's minimum: at a corner) */
function nearestDepth(box: Box3D, eye: Vec3, dir: Vec3) {
  let d = -eye[0] * dir[0] - eye[1] * dir[1] - eye[2] * dir[2];
  for (let a = 0; a < 3; ++a) d += Math.min(box.min[a] * dir[a], box.max[a] * dir[a]);
  return d;
}

/** a ray's floats: origin, direction, the origin's depth, the depth a unit of t, the length a unit of t */
const RayStride = 9;

/**
 * Rays through a small grid of the camera's (pixel corners over normalized device x, y in [-1, 1]), each gathering the
 * optical depth of the blocks added (their density times its path in them, past the near plane), opaque from the
 * deepest back of the blocks it crossed when it reaches `opticalDepth`: the blocks added front to back by tile, so
 * about in its order. A pixel opaque from the deepest of its corners' depths (where all are): a gap narrower than the
 * rays' spacing is not seen. A box tested by its projection's pixels at its nearest depth
 */
class OpticalBuffer {
  /** per pixel corner (n + 1)², x fastest: a ray (`RayStride`) */
  private readonly rays: Float64Array;
  /** per ray the optical depth gathered, the deepest back of the blocks crossed, and where it is opaque from */
  private readonly sum: Float64Array;
  private readonly back: Float64Array;
  private readonly opaque: Float64Array;

  /** the last buffer of a resolution: its arrays reset for the next pick (`get`), not made again a pick */
  private static last: OpticalBuffer | undefined;

  /** a buffer of `n` pixels a side for `camera`: the last one's arrays where of `n`, reset */
  static get(n: number, camera: ParticleTileCamera, opticalDepth: number) {
    const last = OpticalBuffer.last;
    const buffer = new OpticalBuffer(n, camera, opticalDepth, last && last.n === n ? last : undefined);
    OpticalBuffer.last = buffer;
    return buffer;
  }

  constructor(
    readonly n: number,
    readonly camera: ParticleTileCamera,
    readonly opticalDepth: number,
    reuse?: OpticalBuffer,
  ) {
    const m = n + 1;
    this.rays = reuse ? reuse.rays : new Float64Array(m * m * RayStride);
    this.sum = reuse ? reuse.sum.fill(0) : new Float64Array(m * m);
    this.back = (reuse ? reuse.back : new Float64Array(m * m)).fill(-Infinity);
    this.opaque = (reuse ? reuse.opaque : new Float64Array(m * m)).fill(Infinity);
    const { position: eye, direction: dir, orthographic } = camera;
    Mat4.invert(_invVP, camera.viewProjection);
    const r = this.rays;
    for (let j = 0; j < m; ++j) {
      for (let i = 0; i < m; ++i) {
        // a point on the ray inside the frustum (depth 0.5: within both clip depth conventions)
        Vec3.set(_ndc, (i / n) * 2 - 1, (j / n) * 2 - 1, 0.5);
        Vec3.transformMat4(_ndc, _ndc, _invVP);
        const o = (j * m + i) * RayStride;
        if (orthographic) {
          // its origin on the camera's plane, along the view direction
          const t = (_ndc[0] - eye[0]) * dir[0] + (_ndc[1] - eye[1]) * dir[1] + (_ndc[2] - eye[2]) * dir[2];
          r[o] = _ndc[0] - t * dir[0];
          r[o + 1] = _ndc[1] - t * dir[1];
          r[o + 2] = _ndc[2] - t * dir[2];
          r[o + 3] = dir[0];
          r[o + 4] = dir[1];
          r[o + 5] = dir[2];
          r[o + 6] = 0;
          r[o + 7] = 1;
          r[o + 8] = 1;
        } else {
          r[o] = eye[0];
          r[o + 1] = eye[1];
          r[o + 2] = eye[2];
          r[o + 3] = _ndc[0] - eye[0];
          r[o + 4] = _ndc[1] - eye[1];
          r[o + 5] = _ndc[2] - eye[2];
          r[o + 6] = 0;
          r[o + 7] = r[o + 3] * dir[0] + r[o + 4] * dir[1] + r[o + 5] * dir[2];
          r[o + 8] = Math.hypot(r[o + 3], r[o + 4], r[o + 5]);
        }
      }
    }
  }

  /** `rect`'s, its result: read before the next call (no array a box: millions of parts a pick with clip shapes) */
  private readonly _rect = new Int32Array(4);
  /** the last `rect`'s box on screen, across its larger side, in pixels (not rounded to the pixels it touches) */
  extent = 0;

  rect(box: Box3D): Int32Array | undefined {
    const { n } = this;
    const m = this.camera.viewProjection;
    let x0 = Infinity,
      y0 = Infinity,
      x1 = -Infinity,
      y1 = -Infinity;
    for (let c = 0; c < 8; ++c) {
      const x = c & 1 ? box.max[0] : box.min[0];
      const y = c & 2 ? box.max[1] : box.min[1];
      const z = c & 4 ? box.max[2] : box.min[2];
      _clip[0] = m[0] * x + m[4] * y + m[8] * z + m[12];
      _clip[1] = m[1] * x + m[5] * y + m[9] * z + m[13];
      _clip[3] = m[3] * x + m[7] * y + m[11] * z + m[15];
      if (_clip[3] <= 0) return undefined;
      const px = ((_clip[0] / _clip[3]) * 0.5 + 0.5) * n;
      const py = ((_clip[1] / _clip[3]) * 0.5 + 0.5) * n;
      if (px < x0) x0 = px;
      if (px > x1) x1 = px;
      if (py < y0) y0 = py;
      if (py > y1) y1 = py;
    }
    this.extent = Math.max(x1 - x0, y1 - y0);
    const out = this._rect;
    out[0] = Math.max(0, Math.floor(x0));
    out[1] = Math.max(0, Math.floor(y0));
    out[2] = Math.min(n - 1, Math.ceil(x1) - 1);
    out[3] = Math.min(n - 1, Math.ceil(y1) - 1);
    return out;
  }

  /** Whether every pixel `box` touches is opaque nearer than its nearest depth (none of it on screen: hidden) */
  hides(box: Box3D) {
    const { position: eye, direction: dir } = this.camera;
    const nearest = nearestDepth(box, eye, dir);
    if (nearest <= (this.camera.near ?? 0)) return false;
    return this.hidesRect(this.rect(box), nearest);
  }

  /** `hides` of a box of pixels `r` (its `rect`) and nearest depth `nearest` */
  hidesRect(r: Int32Array | undefined, nearest: number) {
    if (nearest <= (this.camera.near ?? 0)) return false;
    if (!r) return false;
    const [x0, y0, x1, y1] = r;
    const { n, opaque } = this;
    const m = n + 1;
    for (let y = y0; y <= y1; ++y) {
      for (let x = x0; x <= x1; ++x) {
        const c = y * m + x;
        const depth = Math.max(opaque[c], opaque[c + 1], opaque[c + m], opaque[c + m + 1]);
        if (!(depth < nearest)) return false;
      }
    }
    return true;
  }

  /**
   * How much light reaches `box` through what was added (the blocks in front: added front to back): the mean over the
   * rays of its pixels of exp(-their optical depth so far); 1 where none of it is on screen
   */
  exposure(box: Box3D) {
    return this.exposureRect(this.rect(box));
  }

  /** `exposure` of a box of pixels `r` (its `rect`) */
  exposureRect(r: Int32Array | undefined) {
    if (!r) return 1;
    const [x0, y0, x1, y1] = r;
    if (x0 > x1 || y0 > y1) return 1;
    const m = this.n + 1;
    let light = 0;
    let rays = 0;
    for (let j = y0; j <= y1 + 1; ++j) {
      for (let i = x0; i <= x1 + 1; ++i) {
        light += Math.exp(-this.sum[j * m + i]);
        rays += 1;
      }
    }
    return light / rays;
  }

  /** Whether every ray through `box`'s pixels is opaque (at any depth) */
  isOpaque(box: Box3D) {
    const r = this.rect(box);
    return !!r && this.isOpaqueRect(r);
  }

  /** whether every ray of pixels `r` (`rect`'s) is opaque already */
  isOpaqueRect(r: ArrayLike<number>) {
    const [x0, y0, x1, y1] = [r[0], r[1], r[2], r[3]];
    const m = this.n + 1;
    for (let j = y0; j <= y1 + 1; ++j) {
      for (let i = x0; i <= x1 + 1; ++i) {
        if (this.opaque[j * m + i] === Infinity) return false;
      }
    }
    return true;
  }

  /**
   * A tile's blocks (`dims` an axis over `cell`, `density` x fastest, 0: none) added to the rays through its cell: per
   * ray one slab test of the cell, then its blocks front to back (a 3D DDA), not a block's rectangle and slabs each. A
   * ray opaque from the deepest back of the blocks of some density it crossed, where its sum reaches the threshold after
   * the tile's (as `add` per block, never nearer)
   */
  addGrid(cell: Box3D, dims: readonly [number, number, number], density: Float32Array) {
    const { n, rays, sum, back, opaque, opticalDepth } = this;
    const m = n + 1;
    const r = this.rect(cell);
    const [x0, y0, x1, y1] = r ? [r[0], r[1], r[2], r[3]] : [0, 0, n - 1, n - 1];
    if (x0 > x1 || y0 > y1) return;
    const near = this.camera.near ?? 0;
    const [bx, by, bz] = dims;
    const sx = (cell.max[0] - cell.min[0]) / bx;
    const sy = (cell.max[1] - cell.min[1]) / by;
    const sz = (cell.max[2] - cell.min[2]) / bz;
    for (let j = y0; j <= y1 + 1; ++j) {
      for (let i = x0; i <= x1 + 1; ++i) {
        const c = j * m + i;
        if (opaque[c] !== Infinity) continue;
        const o = c * RayStride;
        let tn = Math.max(0, (near - rays[o + 6]) / rays[o + 7]);
        let tf = Infinity;
        for (let a = 0; a < 3; ++a) {
          const origin = rays[o + a];
          const d = rays[o + 3 + a];
          if (d === 0) {
            if (origin < cell.min[a] || origin > cell.max[a]) tf = -Infinity;
            continue;
          }
          const t1 = (cell.min[a] - origin) / d;
          const t2 = (cell.max[a] - origin) / d;
          tn = Math.max(tn, Math.min(t1, t2));
          tf = Math.min(tf, Math.max(t1, t2));
        }
        if (!(tf > tn)) continue;
        const ox = rays[o];
        const oy = rays[o + 1];
        const oz = rays[o + 2];
        const dx = rays[o + 3];
        const dy = rays[o + 4];
        const dz = rays[o + 5];
        // the block of the entry (its middle's, clear of the faces), the parameters of the next faces per axis
        const tm = tn + (tf - tn) * 1e-6;
        let ix = Math.min(bx - 1, Math.max(0, Math.floor((ox + dx * tm - cell.min[0]) / sx)));
        let iy = Math.min(by - 1, Math.max(0, Math.floor((oy + dy * tm - cell.min[1]) / sy)));
        let iz = Math.min(bz - 1, Math.max(0, Math.floor((oz + dz * tm - cell.min[2]) / sz)));
        const stepX = dx > 0 ? 1 : -1;
        const stepY = dy > 0 ? 1 : -1;
        const stepZ = dz > 0 ? 1 : -1;
        const nextX = dx === 0 ? Infinity : (cell.min[0] + (ix + (dx > 0 ? 1 : 0)) * sx - ox) / dx;
        const nextY = dy === 0 ? Infinity : (cell.min[1] + (iy + (dy > 0 ? 1 : 0)) * sy - oy) / dy;
        const nextZ = dz === 0 ? Infinity : (cell.min[2] + (iz + (dz > 0 ? 1 : 0)) * sz - oz) / dz;
        const deltaX = dx === 0 ? Infinity : sx / Math.abs(dx);
        const deltaY = dy === 0 ? Infinity : sy / Math.abs(dy);
        const deltaZ = dz === 0 ? Infinity : sz / Math.abs(dz);
        let tx = nextX;
        let ty = nextY;
        let tz = nextZ;
        let t = tn;
        let added = 0;
        let exit = -Infinity;
        while (t < tf) {
          const end = Math.min(tx, ty, tz, tf);
          const value = density[ix + iy * bx + iz * bx * by];
          if (value > 0 && end > t) {
            added += value * (end - t);
            exit = end;
          }
          if (end >= tf) break;
          t = end;
          if (tx <= ty && tx <= tz) {
            ix += stepX;
            tx += deltaX;
            if (ix < 0 || ix >= bx) break;
          } else if (ty <= tz) {
            iy += stepY;
            ty += deltaY;
            if (iy < 0 || iy >= by) break;
          } else {
            iz += stepZ;
            tz += deltaZ;
            if (iz < 0 || iz >= bz) break;
          }
        }
        if (added === 0) continue;
        sum[c] += added * rays[o + 8];
        back[c] = Math.max(back[c], rays[o + 6] + exit * rays[o + 7]);
        if (sum[c] >= opticalDepth) opaque[c] = back[c];
      }
    }
  }

  /** `box` of `density` (1/Å) added to the rays through it (`rect` its pixels where known, `rect(box)`'s) */
  add(box: Box3D, density: number, rect?: ArrayLike<number>) {
    const { n, rays, sum, back, opaque, opticalDepth } = this;
    const m = n + 1;
    const r = rect ?? this.rect(box);
    const [x0, y0, x1, y1] = r ? [r[0], r[1], r[2], r[3]] : [0, 0, n - 1, n - 1];
    if (x0 > x1 || y0 > y1) return;
    const near = this.camera.near ?? 0;
    for (let j = y0; j <= y1 + 1; ++j) {
      for (let i = x0; i <= x1 + 1; ++i) {
        const c = j * m + i;
        // opaque already: what is added does not move where from (the blocks come about front to back)
        if (opaque[c] !== Infinity) continue;
        const o = c * RayStride;
        // the slabs' latest entry and earliest exit, the entry no nearer than the near plane
        let tn = Math.max(0, (near - rays[o + 6]) / rays[o + 7]);
        let tf = Infinity;
        for (let a = 0; a < 3; ++a) {
          const origin = rays[o + a];
          const d = rays[o + 3 + a];
          if (d === 0) {
            if (origin < box.min[a] || origin > box.max[a]) tf = -Infinity;
            continue;
          }
          const t1 = (box.min[a] - origin) / d;
          const t2 = (box.max[a] - origin) / d;
          tn = Math.max(tn, Math.min(t1, t2));
          tf = Math.min(tf, Math.max(t1, t2));
        }
        if (!(tf > tn)) continue;
        sum[c] += density * (tf - tn) * rays[o + 8];
        back[c] = Math.max(back[c], rays[o + 6] + tf * rays[o + 7]);
        if (opaque[c] === Infinity && sum[c] >= opticalDepth) opaque[c] = back[c];
      }
    }
  }
}
