/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util/value-cell';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import type { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import type { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import type { LodView } from './renderable.js';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { BoundaryHelper } from '@molstar/core/math/geometry/boundary-helper';
import { AffineRowsItemSize } from './affine-rows.js';
import { DefineSpec, TextureSpec, UniformSpec, ValueSpec, type RenderableSchema } from './renderable/schema.js';
import type { TextureImage } from './renderable/util.js';

// Instances composed in the vertex stage (the `instanceExpansion` capability, WebGL2): an object's instances (a
// template's particles, in its frame) drawn once per pose of a list (a template's copies), the pose of each culled on
// the CPU by its sphere (the object's instances placed by it); a draw `n x V` instances (`n` the object's, `V` the
// poses seen at the draw's level), the vertex stage composing pose `v` with instance `k` (`gl_InstanceID`, divided by
// `n`). The poses seen written each cull, by level, into a texture (`tCopyPoses`); an instance picked by its draw slot
// (`v n + k`), the poses seen kept until the next cull to tell its pose

/**
 * Copies a row of `tCopyOrder` (the copies a cull saw, their indices, four a texel: 1,024 texels wide), indices as
 * floats (exact below 2^24)
 */
export const CopyOrderPerRow = 4096;

/** Texels a row of `tCopyMarkers` and `tCopyLayers` (four bytes each) */
export const CopyMarkerTexels = 4096;

/**
 * The layers a copy has its own of (`tCopyLayers`), a word each a copy in this order: overpaint, transparency,
 * clipping, emissive, substance, wiggle
 */
export const CopyLayerKinds = 6;
/** A copy's layer word: its value's first byte (the low 30 bits) */
export const CopyLayerValue = 0x80000000;
/** A copy's layer word: its block's first byte (the low 30 bits), a template particle's at its index times the size */
export const CopyLayerBlock = 0x40000000;

/** Texels a pose in `tCopyPoses`: its affine rows */
export const CopyPoseTexels = 3;
/** Poses a row of `tCopyPoses` (its width 1,536 texels: within WebGL2's least `MAX_TEXTURE_SIZE`, 2,048) */
export const CopyPosesPerRow = 512;

/**
 * A template's copies, shared by the objects of its leaf targets (at every coarse level): their poses, and the level
 * each draws its template at by a view (`getCopyLevels`)
 */
export interface ExpansionCopies {
  readonly count: number;
  /** the affine rows of each (`AffineRowsItemSize` floats), in the objects' frame */
  readonly poses: Float32Array;
  /** `poses` as `tCopyPoses` (rows of `CopyPosesPerRow`): one texture for the template's objects (shared) */
  readonly image: TextureImage<Float32Array>;
  /** where each places its template's center (3 floats): what its distance and level are by */
  readonly centers: Float32Array;
  /**
   * its template's particles' spacing (the cube root of its sphere's volume a particle) and coarse levels: a copy draws
   * a level coarser each halving of the spacing's pixels below `ExpansionCoarsePixels` (none without levels)
   */
  readonly spacing: number;
  readonly levelCount: number;
  /** the radius of its template's sphere about `centers` */
  readonly radius: number;
  /**
   * those in the view last asked for (`getCopyBuckets`) by coarse level: `order[offsets[j]]` to `order[offsets[j + 1]]`
   * at level `j`, by their centers' distance from the camera's plane (`distances` theirs, ascending in a level);
   * `levels` each one's (255: out of the view); the view's version (`LodView.version`; -1: none) and the transform
   * moving the copies (`moved` false: none)
   */
  readonly buckets: {
    levels: Uint8Array;
    order: Uint32Array;
    distances: Float32Array;
    offsets: Int32Array;
    version: number;
    moved: boolean;
    readonly transform: Mat4;
  };
  /**
   * `buckets.order` of the views last asked for (`getCopyOrderImage`, by the view and the transform), as `tCopyOrder`'s
   * images: what the template's objects bind (a level's copies a slice of it), one upload a view for all of them
   */
  readonly orders: {
    readonly views: { version: number; moved: boolean; transform: Mat4; image: TextureImage<Float32Array> }[];
    next: number;
  };
}

/** The poses an object's instances are drawn by */
export interface InstanceExpansion {
  readonly count: number;
  /** the copies' (`ExpansionCopies.poses`) */
  readonly poses: Float32Array;
  /** the sphere of the object's instances placed by each (center, radius: 4 floats) */
  readonly spheres: Float32Array;
  readonly copies: ExpansionCopies;
  /** the coarse level of its template the object is (0 the template whole): the copies at it alone drawn */
  readonly level: number;
  /**
   * the poses the visible passes' last cull saw, by level (`SegmentedMdbList.expandedSegment`): a pick's draw slot
   * `v n + k` the pose `order[v]`
   */
  readonly seen: { order: ArrayLike<number>; count: number };
}

/**
 * a template's spacing in CSS pixels below which a copy draws it at a coarse level, a level coarser each halving (the
 * page store's `PagedCoarsePixels`)
 */
export const ExpansionCoarsePixels = 4;

/**
 * The copies of `transforms` (column-major `Mat4`s, 16 floats each), each placing its template's sphere `sphere`, of
 * `levelCount` coarse levels and particles `spacing` apart
 */
export function createExpansionCopies(
  transforms: Float32Array,
  count: number,
  sphere: Sphere3D,
  spacing: number,
  levelCount: number,
): ExpansionCopies {
  // the texture's image the poses (a pose 12 floats, `CopyPosesPerRow` a row)
  const width = CopyPosesPerRow * CopyPoseTexels;
  const height = Math.max(1, Math.ceil(count / CopyPosesPerRow));
  const poses = new Float32Array(width * height * 4);
  const image: TextureImage<Float32Array> = { array: poses, width, height };
  const centers = new Float32Array(count * 3);
  const m = Mat4();
  const c = Vec3();
  for (let i = 0; i < count; ++i) {
    Mat4.fromArray(m, transforms, i * 16);
    for (let r = 0; r < 3; ++r) {
      poses[i * 12 + r * 4] = m[r];
      poses[i * 12 + r * 4 + 1] = m[4 + r];
      poses[i * 12 + r * 4 + 2] = m[8 + r];
      poses[i * 12 + r * 4 + 3] = m[12 + r];
    }
    Vec3.toArray(Vec3.transformMat4(c, sphere.center, m), centers, i * 3);
  }
  const buckets = {
    levels: new Uint8Array(count),
    order: new Uint32Array(count),
    distances: new Float32Array(count),
    offsets: new Int32Array(levelCount + 2),
    version: -1,
    moved: false,
    transform: Mat4(),
  };
  const orders = { views: [], next: 0 };
  return { count, poses, image, centers, spacing, levelCount, radius: sphere.radius, buckets, orders };
}

/**
 * The copies of `copies` in a view by the coarse level each draws its template at (`ExpansionCopies.buckets`): those
 * whose template's sphere reaches into `frustum`, by level (`ExpansionCoarsePixels`: its spacing's pixels at its
 * center's distance from `cameraPlane`, `lod`'s scale; orthographic everywhere the same); `transform` moving them
 * (none: identity). Kept for the view's version (`LodView.version`; not where 0) and the transform: a moved
 * representation's objects (a template's alike) the same buckets, not each its own
 */
export function getCopyBuckets(
  copies: ExpansionCopies,
  cameraPlane: Plane3D,
  frustum: Frustum3D,
  lod: LodView,
  transform?: Mat4,
) {
  const { buckets, levelCount, count, centers, spacing } = copies;
  if (lod.version !== 0 && buckets.version === lod.version && isSameMove(buckets, transform)) return buckets;
  buckets.version = lod.version;
  buckets.moved = !!transform;
  if (transform) Mat4.copy(buckets.transform, transform);
  const ortho = !(lod.depthScale > 0);
  const scale = ortho ? (lod.pixelSize > 0 ? 1 / lod.pixelSize : 0) : 1 / lod.depthScale;
  const radius = copies.radius * (transform ? Mat4.getMaxScaleOnAxis(transform) : 1);
  const { normal, constant } = cameraPlane;
  const { levels, order, distances, offsets } = buckets;
  const c = Vec3();
  offsets.fill(0);
  // by copy first, by order after the sort
  const byCopy = getDistanceScratch(count);
  outer: for (let i = 0; i < count; ++i) {
    Vec3.fromArray(c, centers, i * 3);
    if (transform) Vec3.transformMat4(c, c, transform);
    levels[i] = NotSeen;
    for (let p = 0; p < 6; ++p) {
      const plane = frustum[p as 0];
      if (Vec3.dot(plane.normal, c) + plane.constant < -radius) continue outer;
    }
    const distance = Vec3.dot(normal, c) + constant;
    byCopy[i] = distance;
    let level = 0;
    if (levelCount > 0 && scale > 0) {
      const d = ortho ? 1 : Math.max(distance, 1e-3);
      const px = (spacing * scale) / d;
      if (px < ExpansionCoarsePixels) {
        level = Math.min(levelCount, Math.floor(Math.log2(ExpansionCoarsePixels / Math.max(px, 1e-9))));
      }
    }
    levels[i] = level;
    offsets[level + 1] += 1;
  }
  for (let j = 0; j <= levelCount; ++j) offsets[j + 1] += offsets[j];
  const cursors = offsets.slice(0, levelCount + 1);
  for (let i = 0; i < count; ++i) if (levels[i] !== NotSeen) order[cursors[levels[i]]++] = i;
  // a level by distance: an object's copies at a level of its own a slice (`getCopyRange`)
  const byDistance = (a: number, b: number) => byCopy[a] - byCopy[b];
  for (let j = 0; j <= levelCount; ++j) {
    if (offsets[j + 1] - offsets[j] > 1) order.subarray(offsets[j], offsets[j + 1]).sort(byDistance);
  }
  for (let k = 0, kl = offsets[levelCount + 1]; k < kl; ++k) distances[k] = byCopy[order[k]];
  return buckets;
}

/** whether `a` is moved by `transform` (none: not moved), exactly */
function isSameMove(a: { moved: boolean; transform: Mat4 }, transform: Mat4 | undefined) {
  if (!transform) return !a.moved;
  if (!a.moved) return false;
  for (let i = 0; i < 16; ++i) if (a.transform[i] !== transform[i]) return false;
  return true;
}

/**
 * The views a template's copy order is kept for at once (`getCopyOrderImage`): an object culled by an earlier view (a
 * pick's, another canvas') still draws by its order
 */
const CopyOrderViews = 4;

/** `tCopyOrder`'s images of a template's views (`ExpansionCopies.orders`), not of an object alone */
const sharedCopyOrders = new WeakSet<object>();

/** Whether `image` is a template's copy order (`getCopyOrderImage`), bound by its objects as one texture */
export function isSharedCopyOrder(image: object) {
  return sharedCopyOrders.has(image);
}

/**
 * The image of `copies`' order by the view and the transform their buckets are of (`getCopyBuckets`, a view's version
 * above 0): written
 * where none is of it, into the view's least recent's (its `version` advanced: an object culled by that one draws by
 * it no longer, and is culled again)
 */
export function getCopyOrderImage(copies: ExpansionCopies) {
  const { buckets, orders, levelCount } = copies;
  const { views } = orders;
  for (let i = 0, il = views.length; i < il; ++i) {
    const v = views[i];
    if (v.version === buckets.version && isSameMove(v, buckets.moved ? buckets.transform : undefined)) return v.image;
  }
  const count = buckets.offsets[levelCount + 1];
  let view = views[orders.next];
  if (!view || view.image.height * CopyOrderPerRow < count) {
    const image = createCopyOrderImage(Math.max(1, count));
    sharedCopyOrders.add(image);
    view = views[orders.next] = { version: -1, moved: false, transform: Mat4(), image };
  }
  orders.next = (orders.next + 1) % CopyOrderViews;
  view.version = buckets.version;
  view.moved = buckets.moved;
  Mat4.copy(view.transform, buckets.transform);
  const { image } = view;
  const { array } = image;
  const { order } = buckets;
  for (let i = 0; i < count; ++i) array[i] = order[i];
  image.version = (image.version ?? 0) + 1;
  return image;
}

let distanceScratch = new Float32Array(0);
function getDistanceScratch(count: number) {
  if (distanceScratch.length < count) distanceScratch = new Float32Array(Math.max(count, (count * 3) >> 1));
  return distanceScratch;
}

/**
 * The positions `[from, to)` within `distances`' `[first, last)` (ascending) of the copies within `reach` of the
 * distances `min` to `max`: those whose centers' distance `d` has `d + reach >= min` and `d - reach <= max`
 */
export function getCopyRange(
  distances: Float32Array,
  first: number,
  last: number,
  min: number,
  max: number,
  reach: number,
  out: [number, number],
) {
  out[0] = lowerBound(distances, first, last, min - reach);
  out[1] = Math.max(out[0], upperBound(distances, out[0], last, max + reach));
  return out;
}

/** the first position in `[lo, hi)` whose value is `x` or more (`hi` if none) */
function lowerBound(a: Float32Array, lo: number, hi: number, x: number) {
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (a[mid] < x) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

/** the first position in `[lo, hi)` whose value is more than `x` (`hi` if none) */
function upperBound(a: Float32Array, lo: number, hi: number, x: number) {
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (a[mid] <= x) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

const reaches = new WeakMap<Float32Array, { centers: Float32Array; reach: number }>();
/**
 * How far an object's instances reach from the center its copy places its template's at, the most of any copy
 * (`spheres` about `centers`): an object's sphere's distance within it of its copy's center's; once an object
 */
export function getCopyReach(expansion: InstanceExpansion) {
  const { spheres, copies } = expansion;
  const { centers } = copies;
  const cached = reaches.get(spheres);
  if (cached && cached.centers === centers) return cached.reach;
  let reach = 0;
  for (let i = 0, il = expansion.count; i < il; ++i) {
    const dx = spheres[i * 4] - centers[i * 3];
    const dy = spheres[i * 4 + 1] - centers[i * 3 + 1];
    const dz = spheres[i * 4 + 2] - centers[i * 3 + 2];
    reach = Math.max(reach, Math.sqrt(dx * dx + dy * dy + dz * dz) + spheres[i * 4 + 3]);
  }
  reaches.set(spheres, { centers, reach });
  return reach;
}

/** a copy out of the view (`getCopyBuckets`) */
const NotSeen = 255;

export const InstanceExpansionSchema = {
  dInstanceExpansion: DefineSpec('boolean'),
  /** the object's instances drawn by each pose (the draw's `gl_InstanceID` divided by it) */
  uExpansionCount: UniformSpec('i'),
  /** the first copy of a draw in `tCopyOrder` (its level's, `SegmentedMdbList.expandedSegment`) */
  uCopyBase: UniformSpec('i', undefined, 'draw'),
  /** the template's copies' poses (`ExpansionCopies.image`): one texture for its objects */
  tCopyPoses: TextureSpec('image-float32', 'rgba', 'float', 'nearest', undefined, { shared: true }),
  /** the copies the cull saw, by level (`CopyOrderPerRow`): written by the cull before a pass's draws */
  tCopyOrder: TextureSpec('image-float32', 'rgba', 'float', 'nearest'),
  /**
   * the copies' markers (a representation's, its templates' copies one after another: `uCopyMarkerBase` an object's
   * template's first): two words a copy, then the blocks' bytes (`CopyMarkerTexels` a row, `writeCopyMarkerImage`);
   * one texture for the representation's objects, synced where drawn
   */
  tCopyMarkers: TextureSpec('image-uint8', 'rgba', 'ubyte', 'nearest', undefined, { shared: true }),
  uCopyMarkerBase: UniformSpec('i'),
  /**
   * the copies' own layers (a representation's, by copy as `tCopyMarkers`: `CopyLayerKinds` words a copy, none, a
   * value or a block by template particle, then their bytes): what an instance's copy has replaces its object's
   */
  tCopyLayers: TextureSpec('image-uint8', 'rgba', 'ubyte', 'nearest', undefined, { shared: true }),
  dCopyLayers: DefineSpec('boolean'),
  /** per instance (as stored) its template particle's index: a copy's block's byte (`expansionCopyMarker`) */
  expansionLocals: ValueSpec('unknown'),
  instanceExpansion: ValueSpec('unknown'),
} satisfies RenderableSchema;

export type InstanceExpansionValues = {
  dInstanceExpansion: ValueCell<boolean>;
  uExpansionCount: ValueCell<number>;
  uCopyBase: ValueCell<number>;
  tCopyPoses: ValueCell<TextureImage<Float32Array>>;
  tCopyOrder: ValueCell<TextureImage<Float32Array>>;
  tCopyMarkers: ValueCell<TextureImage<Uint8Array>>;
  uCopyMarkerBase: ValueCell<number>;
  tCopyLayers: ValueCell<TextureImage<Uint8Array>>;
  dCopyLayers: ValueCell<boolean>;
  expansionLocals: ValueCell<Uint32Array>;
  instanceExpansion: ValueCell<InstanceExpansion>;
};

/** Whether `values` are drawn expanded (`InstanceExpansionValues`) */
export function isExpanded(values: object): values is InstanceExpansionValues {
  return !!(values as Partial<InstanceExpansionValues>).dInstanceExpansion?.ref.value;
}

/** An object's expansion by `copies` at its template's coarse level `level`, `sphere` its instances' */
export function createInstanceExpansion(copies: ExpansionCopies, sphere: Sphere3D, level: number): InstanceExpansion {
  const { count, poses } = copies;
  const spheres = new Float32Array(count * 4);
  // rigid poses: the center placed, the radius kept
  const [x, y, z] = sphere.center;
  for (let i = 0; i < count; ++i) {
    const o = i * AffineRowsItemSize;
    spheres[i * 4] = poses[o] * x + poses[o + 1] * y + poses[o + 2] * z + poses[o + 3];
    spheres[i * 4 + 1] = poses[o + 4] * x + poses[o + 5] * y + poses[o + 6] * z + poses[o + 7];
    spheres[i * 4 + 2] = poses[o + 8] * x + poses[o + 9] * y + poses[o + 10] * z + poses[o + 11];
    spheres[i * 4 + 3] = sphere.radius;
  }
  return { count, poses, spheres, copies, level, seen: { order: new Uint32Array(0), count: 0 } };
}

/**
 * The bounding sphere of the object's instances placed by `expansion`'s poses: the corners of their centers' box
 * (`min`, `max`) placed by each, `radius` about them (a list with templates' bounds as `Particle.getBoundary` has them:
 * a copy its template's box's corners, tighter than its sphere)
 */
export function getExpansionBoundingSphere(copies: ExpansionCopies, min: Vec3, max: Vec3, radius: number) {
  const { poses, count } = copies;
  const helper = new BoundaryHelper('98');
  const p = Vec3();
  const corner = (i: number, k: number) => {
    const x = k & 1 ? max[0] : min[0];
    const y = k & 2 ? max[1] : min[1];
    const z = k & 4 ? max[2] : min[2];
    const o = i * AffineRowsItemSize;
    p[0] = poses[o] * x + poses[o + 1] * y + poses[o + 2] * z + poses[o + 3];
    p[1] = poses[o + 4] * x + poses[o + 5] * y + poses[o + 6] * z + poses[o + 7];
    p[2] = poses[o + 8] * x + poses[o + 9] * y + poses[o + 10] * z + poses[o + 11];
    return p;
  };
  for (let i = 0; i < count; ++i) for (let k = 0; k < 8; ++k) helper.includePositionRadius(corner(i, k), radius);
  helper.finishedIncludeStep();
  for (let i = 0; i < count; ++i) for (let k = 0; k < 8; ++k) helper.radiusPositionRadius(corner(i, k), radius);
  return helper.getSphere();
}

/**
 * The values that draw an object of `n` instances by `expansion`'s poses, its copies' markers `markers` from
 * `markerBase` (their layers `layers` by the same copies, read where `hasLayers`), `locals` its instances' template
 * particles
 */
export function createInstanceExpansionValues(
  expansion: InstanceExpansion,
  n: number,
  markers: TextureImage<Uint8Array>,
  markerBase: number,
  locals: Uint32Array,
  layers: TextureImage<Uint8Array>,
  hasLayers: boolean,
): InstanceExpansionValues {
  return {
    tCopyMarkers: ValueCell.create(markers),
    uCopyMarkerBase: ValueCell.create(markerBase),
    tCopyLayers: ValueCell.create(layers),
    dCopyLayers: ValueCell.create(hasLayers),
    expansionLocals: ValueCell.create(locals),
    dInstanceExpansion: ValueCell.create(true),
    uExpansionCount: ValueCell.create(n),
    uCopyBase: ValueCell.create(0),
    tCopyPoses: ValueCell.create(expansion.copies.image),
    tCopyOrder: ValueCell.create(createCopyOrderImage(1)),
    instanceExpansion: ValueCell.create(expansion),
  };
}

function createCopyOrderImage(capacity: number): TextureImage<Float32Array> {
  const width = CopyOrderPerRow / 4;
  const height = Math.max(1, Math.ceil(capacity / CopyOrderPerRow));
  return { array: new Float32Array(width * height * 4), width, height };
}

/** The copies `order` lists (`count` of them) written into `image` (a new one where it holds fewer) */
export function writeCopyOrder(image: TextureImage<Float32Array> | undefined, order: Uint32Array, count: number) {
  const out =
    image && image.height * CopyOrderPerRow >= count ? image : createCopyOrderImage(Math.max(count, (count * 3) >> 1));
  const { array } = out;
  for (let i = 0; i < count; ++i) array[i] = order[i];
  return out;
}

/**
 * The copies' markers' GPU buffer (`words` then `bytes`, `CopyMarkers`'s: `byteLength` of it) as `tCopyMarkers`'s image:
 * `image` written (its version advanced) where it holds them, else a new one
 */
export function writeCopyMarkerImage(
  image: TextureImage<Uint8Array> | undefined,
  words: Uint32Array,
  bytes: Uint8Array,
  byteLength: number,
): TextureImage<Uint8Array> {
  const out = getCopyImage(image, byteLength);
  out.array.set(new Uint8Array(words.buffer, words.byteOffset, words.byteLength), 0);
  out.array.set(bytes.subarray(0, byteLength - words.byteLength), words.byteLength);
  out.version = (out.version ?? 0) + 1;
  return out;
}

/** The copies' layers' buffer (`CopyLayerKinds` words a copy, then the bytes) as `tCopyLayers`'s image, as markers' */
export function writeCopyLayerImage(
  image: TextureImage<Uint8Array> | undefined,
  buffer: Uint8Array,
): TextureImage<Uint8Array> {
  const out = getCopyImage(image, buffer.length);
  out.array.set(buffer, 0);
  out.version = (out.version ?? 0) + 1;
  return out;
}

/** `image` where it holds `byteLength` bytes, else one grown by half again (`CopyMarkerTexels` a row) */
function getCopyImage(image: TextureImage<Uint8Array> | undefined, byteLength: number) {
  const rowBytes = CopyMarkerTexels * 4;
  const height = Math.max(1, Math.ceil(byteLength / rowBytes));
  // in place where it fits (the shared texture the image's: rows uploaded by its version), else grown by half again
  const out: TextureImage<Uint8Array> =
    image && image.height >= height
      ? image
      : {
          array: new Uint8Array(rowBytes * ((height * 3 + 1) >> 1)),
          width: CopyMarkerTexels,
          height: (height * 3 + 1) >> 1,
        };
  return out;
}
