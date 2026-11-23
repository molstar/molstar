/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { InstanceTransforms, type InstancePoses } from '@molstar/core/math/geometry/instance-transforms';
import { Mat3 } from '@molstar/core/math/linear-algebra/3d/mat3';
import type { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Quat } from '@molstar/core/math/linear-algebra/3d/quat';
import type { TransformPart } from './transform-parts.js';
import type { U32Run } from './compute/u32-runs.js';
import { isSharedSerial } from '@molstar/core/util/array';
import { SegmentTexels } from './renderable/schema.js';

/** The segments (a merged object's members, an unmerged object's one) instances belong to, by their first instance */
export interface CulledSegments {
  readonly count: number;
  readonly instanceBases: ArrayLike<number>;
  readonly instanceCounts: ArrayLike<number>;
}

/** The segments' groups and draws (a merged object's, `MergedSegments`) */
export interface CulledSegmentGroups {
  readonly count: number;
  readonly groupCounts: ArrayLike<number>;
  readonly groupInstanceOffsets: ArrayLike<number>;
  readonly groupBases: ArrayLike<number>;
  readonly drawBases: ArrayLike<number>;
}

/** What a writer of `CulledData` is given: a buffer's write (`StorageBuffer.write`) */
export type CulledDataWrite = (data: ArrayBufferView, byteOffset: number) => void;

/** `CulledData` of a set of instances: its size, and its write into a buffer of at least that size */
export interface CulledDataLayout {
  readonly vec4s: number;
  /** the slots' first vec4 */
  readonly slotFirst: number;
  /** the parts of poses, their slots from `slotFirst + base` */
  readonly poseSlots: readonly CulledPoseSlots[];
  /**
   * the words that are runs (`u32-runs.ts`): the instances' slot words of parts in their own order, `aInstance` of
   * the segments where it is serial, the segment indices (a scan of `aInstance` at the first read)
   */
  readonly runs: readonly U32Run[];
  /**
   * `poses` false: the poses' slots left out, written by the GPU (`culled-poses.ts`); `runs` false: the words of
   * `runs` left out, written by the GPU (`createU32RunWriter`)
   */
  write(write: CulledDataWrite, options?: { readonly poses?: boolean; readonly runs?: boolean }): void;
}

/** The slot word's bit of a pose (two vec4s), else the slot holds affine rows (three) */
const PoseBit = 0x80000000;

/** The similarity (`k Rm`, then `m`) a pose's matrix folds into its elements */
export interface CulledSimilarity {
  readonly k: number;
  readonly rotation: Quat;
  readonly matrix: Mat4;
}

/** A part's poses: its ids' elements from `base` (two vec4s an id), its matrix folded in where one */
export interface CulledPoseSlots {
  readonly poses: InstancePoses;
  readonly similarity: CulledSimilarity | undefined;
  /** its ids: all its `order` may hold (an instance count cut to fewer leaves some unread) */
  readonly ids: number;
  readonly base: number;
}

/**
 * The instances' data of a GPU culler (`merged-culler.ts`, `sCulledData`) and of the sphere raster (`transforms`), read
 * by index: vec4s in one buffer (one binding: the vertex stage's storage buffers are few), for `n` instances of
 * `segmentCount` segments, `q = ceil(n / 4)`:
 * - 0: the header, u32s: `n`, `segmentCount`, the slots' first vec4, 0 (`CulledPagedFlag`: a page store's)
 * - from 1: each instance's slot word (u32, four a vec4): its slot's vec4 from the slots' first, `PoseBit` for a pose
 * - from 1 + q: `aInstance` (four a vec4, u32 bits), from 1 + 2 q: each instance's segment index (four a vec4)
 * - from 1 + 3 q: each segment's data, three vec4s (`CulledSegmentVec4s`): its groups (`aSegment.xyz`, the same for its
 *   instances, i32 bits) and its first item (u32 bits, spheres' and points'); its quantized centers' origin (xyz) and scale
 *   (spheres' `dSites`); its first entries of the size table (spheres' `siteType` sizes) and of the second color's
 *   (spheres' `dSiteColor` 'siteType')
 * - the slots: an instance of matrices its affine rows (three vec4s, `affine-rows.ts`); an id of poses (its element of
 *   the particle list, by the part's `indices` or `first`; an instance's slot by its id, `order`) two vec4s:
 *   `t = p - s R c` and `s`, and `R`'s quaternion, its part's matrix folded in (a similarity: `k Rm` and `m`, so
 *   `k s`, `qm q` and `M t`); a matrix that is not one keeps its part's instances on rows. A part's own: the targets
 *   of a list place particles of their own, scattered in it (a list's range for each would repeat it)
 *
 * The reader builds a pose's rows, `[s R | t]` (`InstanceTransforms`): 12 B an instance and 32 B a particle, where rows
 * were 56 B an instance, and a read a vertex more than rows (the slot word; the header is one address for all)
 *
 * A page store's (`paged-cull.ts`, the header's last word `CulledPagedFlag`, `n` 0): no instance arrays, the segments
 * from 1, instance `i` (a pose slot, the culler's regions hold) its composed pose at the slots' first plus `2 i`
 * (`PagedPoseWords`: `t` and the scale f32, the quaternion snorm16, `aInstance` its store index, its segment)
 */
/** the header's last word of a page store's `CulledData` */
export const CulledPagedFlag = 1;
/**
 * a page store's pose's last word: its segment (a member: below 2^16, `PagedSegmentMask`) and, a copy's instance, the
 * copy plus one above it (`PagedCopyShift`: a pick's copy, gathered from the slot it drew, `paged-pick.ts`)
 */
export const PagedCopyShift = 16;
/**
 * a page store's copy's first word of its markers (`sCopyMarkers`): its block's flag, its second word then the block's
 * first byte less its template's first store index; else its marker for all of it
 */
export const CopyMarkerBlock = 0x80000000;
export const PagedSegmentMask = 0xffff;
/** vec4s of a segment's data in `CulledData` */
export const CulledSegmentVec4s = 3;

export const CulledData = {
  /** vec4s of `n` instances' values: four a vec4 */
  quarter: (n: number) => Math.ceil(n / 4),
  /**
   * WGSL readers of instance `i` (`u32`) of `data` (an `array<vec4<f32>>`), named by `prefix`: `<prefix>Rows` (the
   * transform's affine rows, `mat3x4<f32>` a row each), `<prefix>Instance`, `<prefix>Segment` (as `aSegment`: the
   * segment's groups, its index), and of segment `s` (`u32`) `<prefix>SegmentItem` (its first item), `<prefix>SegmentSites`
   * (its quantized centers' origin and scale), `<prefix>SegmentSizeBase` (its first entry of the size table) and
   * `<prefix>SegmentSiteColorBase` (of the second color's)
   */
  wgsl: (data: string, prefix = 'culled') => /* wgsl */ `
fn ${prefix}Quarter() -> u32 { return (bitcast<u32>(${data}[0].x) + 3u) / 4u; }
fn ${prefix}Paged() -> bool { return bitcast<u32>(${data}[0].w) == ${CulledPagedFlag}u; }
/**
 * a page store's pose i: its second vec4's words (the quaternion snorm16, aInstance, the segment in the low 16 bits
 * and a copy's plus one in the high 16: \`PagedCopyShift\`)
 */
fn ${prefix}PagedWords(i: u32) -> vec4<u32> {
    return bitcast<vec4<u32>>(${data}[bitcast<u32>(${data}[0].z) + 2u * i + 1u]);
}
fn ${prefix}Instance(i: u32) -> u32 {
    if (${prefix}Paged()) { return ${prefix}PagedWords(i).z; }
    return bitcast<u32>(${data}[1u + ${prefix}Quarter() + i / 4u][i % 4u]);
}
fn ${prefix}Segment(i: u32) -> vec4<i32> {
    let q = ${prefix}Quarter();
    var s: u32;
    if (${prefix}Paged()) { s = ${prefix}PagedWords(i).w & ${PagedSegmentMask}u; } else { s = u32(${data}[1u + 2u * q + i / 4u][i % 4u]); }
    return vec4<i32>(bitcast<vec3<i32>>(${data}[1u + 3u * q + ${CulledSegmentVec4s}u * s].xyz), i32(s));
}
fn ${prefix}SegmentItem(s: u32) -> u32 {
    return bitcast<u32>(${data}[1u + 3u * ${prefix}Quarter() + ${CulledSegmentVec4s}u * s].w);
}
fn ${prefix}SegmentSites(s: u32) -> vec4<f32> {
    return ${data}[1u + 3u * ${prefix}Quarter() + ${CulledSegmentVec4s}u * s + 1u];
}
fn ${prefix}SegmentSizeBase(s: u32) -> u32 {
    return u32(${data}[1u + 3u * ${prefix}Quarter() + ${CulledSegmentVec4s}u * s + 2u].x);
}
fn ${prefix}SegmentSiteColorBase(s: u32) -> u32 {
    return u32(${data}[1u + 3u * ${prefix}Quarter() + ${CulledSegmentVec4s}u * s + 2u].y);
}
fn ${prefix}Rows(i: u32) -> mat3x4<f32> {
    let paged = ${prefix}Paged();
    // a page store's: no slot word, its pose at 2 i
    var word = 0x80000000u | (2u * i);
    if (!paged) { word = bitcast<u32>(${data}[1u + i / 4u][i % 4u]); }
    let o = bitcast<u32>(${data}[0].z) + (word & 0x7fffffffu);
    if ((word & 0x80000000u) == 0u) { return mat3x4<f32>(${data}[o], ${data}[o + 1u], ${data}[o + 2u]); }
    // the pose: [s R | t], R from the unit quaternion (\`Mat4.fromQuat\`; a page store's snorm16)
    let ts = ${data}[o];
    var r = ${data}[o + 1u];
    if (paged) {
        let w = bitcast<vec4<u32>>(r);
        r = normalize(vec4<f32>(unpack2x16snorm(w.x), unpack2x16snorm(w.y)));
    }
    let s = ts.w;
    let s2 = 2.0 * s;
    let x = r.x;
    let y = r.y;
    let z = r.z;
    let w = r.w;
    return mat3x4<f32>(
        vec4<f32>(s - s2 * (y * y + z * z), s2 * (x * y - w * z), s2 * (x * z + w * y), ts.x),
        vec4<f32>(s2 * (x * y + w * z), s - s2 * (x * x + z * z), s2 * (y * z - w * x), ts.y),
        vec4<f32>(s2 * (x * z - w * y), s2 * (y * z + w * x), s - s2 * (x * x + y * y), ts.z)
    );
}
`,

  /**
   * The data of the first `n` instances of `parts` (`getTransformParts`, each one segment's) with their `aInstance`,
   * `segmentData` the segments' (`CulledSegmentVec4s` vec4s each; zero without), written `block` instances at a time:
   * not composed into one array first (millions of instances)
   */
  layout(
    parts: readonly TransformPart[],
    aInstance: Uint32Array,
    segments: CulledSegments,
    segmentData: Float32Array | undefined,
    n: number,
    block: number,
  ): CulledDataLayout {
    const q = Math.ceil(n / 4);
    const segmentCount = segments.count;
    const segmentBase = 1 + 3 * q;
    const slotFirst = segmentBase + CulledSegmentVec4s * segmentCount;

    // each part's slots: rows of matrices (or of poses whose matrix is not a similarity), else its ids' poses
    const partBases: number[] = [];
    const partPoses: (CulledPoseSlots | undefined)[] = [];
    let slots = 0;
    for (const p of parts) {
      const t = p.transforms;
      partBases.push(slots);
      const similarity = t.kind === 'poses' && t.matrix ? getSimilarity(t.matrix) : undefined;
      if (t.kind === 'matrices' || (t.matrix && !similarity)) {
        partPoses.push(undefined);
        slots += 3 * t.count;
        continue;
      }
      const ids = Math.max(t.count, t.order?.length ?? 0);
      partPoses.push({ poses: t, similarity, ids, base: slots });
      slots += 2 * ids;
    }

    // the words that are runs: a part's slot words but where in an order of its own, a segment's `aInstance` where serial
    // (all of them in one run where they continue), every segment's indices
    let runs: U32Run[] | undefined;
    const slotRuns: boolean[] = [];
    const instanceRuns: boolean[] = [];
    const getRuns = () => {
      if (runs) return runs;
      runs = [];
      for (let k = 0; k < parts.length; ++k) {
        const p = parts[k];
        const t = p.transforms;
        const poses = partPoses[k];
        slotRuns[k] = !poses?.poses.order;
        if (!slotRuns[k] || t.count === 0) continue;
        runs.push(
          poses
            ? { word: 4 + p.first, count: t.count, start: poses.base, step: 2, bits: PoseBit }
            : { word: 4 + p.first, count: t.count, start: partBases[k], step: 3, bits: 0 },
        );
      }
      let instanceRun: { word: number; count: number; start: number; step: number; bits: number } | undefined;
      const serialInstances = isSharedSerial(aInstance);
      for (let s = 0; s < segmentCount; ++s) {
        const first = segments.instanceBases[s];
        const count = Math.min(segments.instanceCounts[s], n - first);
        // a shared serial array's view known serial, not scanned (millions of instances)
        instanceRuns[s] = count <= 0 || serialInstances || isSerial(aInstance, first, count);
        if (count <= 0) continue;
        if (!instanceRuns[s]) {
          instanceRun = undefined;
          continue;
        }
        const word = (1 + q) * 4 + first;
        if (
          instanceRun &&
          instanceRun.word + instanceRun.count === word &&
          (instanceRun.start + instanceRun.count) >>> 0 === aInstance[first]
        ) {
          instanceRun.count += count;
        } else {
          instanceRun = { word, count, start: aInstance[first], step: 1, bits: 0 };
          runs.push(instanceRun);
        }
      }
      for (let s = 0; s < segmentCount; ++s) {
        const first = segments.instanceBases[s];
        const count = Math.min(segments.instanceCounts[s], n - first);
        if (count > 0) runs.push({ word: (1 + 2 * q) * 4 + first, count, start: floatBits(s), step: 0, bits: 0 });
      }
      return runs;
    };

    return {
      vec4s: slotFirst + slots,
      slotFirst,
      poseSlots: partPoses.filter((p): p is CulledPoseSlots => !!p),
      get runs() {
        return getRuns();
      },
      write(write, options) {
        write(new Uint32Array([n, segmentCount, slotFirst, 0]), 0);
        const asRuns = options?.runs === false;
        if (asRuns) getRuns();

        // the instances' slot words
        const words = new Uint32Array(Math.max(1, Math.min(n, block)));
        for (let k = 0; k < parts.length; ++k) {
          const p = parts[k];
          const t = p.transforms;
          const poses = partPoses[k];
          if (asRuns && slotRuns[k]) continue;
          for (let b = 0; b < t.count; b += block) {
            const m = Math.min(block, t.count - b);
            if (poses) {
              const { order } = poses.poses;
              // `|` makes signed 32-bit integers: the typed array keeps their bits
              for (let i = 0; i < m; ++i) words[i] = (poses.base + 2 * (order ? order[b + i] : b + i)) | PoseBit;
            } else {
              const base = partBases[k] + 3 * b;
              for (let i = 0; i < m; ++i) words[i] = base + 3 * i;
            }
            write(words.subarray(0, m), 16 + (p.first + b) * 4);
          }
        }
        if (!asRuns) write(aInstance.subarray(0, n), (1 + q) * 16);
        else {
          // the segments' not serial (instances of no segment: none, as `getRuns` holds them)
          for (let s = 0; s < segmentCount; ++s) {
            if (instanceRuns[s]) continue;
            const first = segments.instanceBases[s];
            const count = Math.min(segments.instanceCounts[s], n - first);
            write(aInstance.subarray(first, first + count), (1 + q) * 16 + first * 4);
          }
        }

        // the instances' segment indices (a float each, exact for 16 M segments)
        const indexBase = (1 + 2 * q) * 16;
        const index = new Float32Array(Math.max(1, Math.min(n, block)));
        for (let s = 0; s < segmentCount && !asRuns; ++s) {
          const first = segments.instanceBases[s];
          const count = Math.min(segments.instanceCounts[s], n - first);
          index.fill(s, 0, Math.min(block, count));
          for (let k = 0; k < count; k += block) {
            const m = Math.min(block, count - k);
            write(index.subarray(0, m), indexBase + (first + k) * 4);
          }
        }
        if (segmentCount > 0) {
          write(segmentData ?? new Float32Array(segmentCount * CulledSegmentVec4s * 4), segmentBase * 16);
        }

        // the slots: rows, the particle lists' poses
        const rows = new Float32Array(Math.max(1, Math.min(n, block)) * 12);
        for (let k = 0; k < parts.length; ++k) {
          if (partPoses[k]) continue;
          const t = parts[k].transforms;
          for (let b = 0; b < t.count; b += block) {
            const m = Math.min(block, t.count - b);
            InstanceTransforms.writeAffineRows(t, b, m, rows);
            write(rows.subarray(0, m * 12), (slotFirst + partBases[k] + 3 * b) * 16);
          }
        }
        if (options?.poses === false) return;
        const poses = new Float32Array(Math.max(1, block) * 8);
        for (const p of partPoses) {
          if (!p) continue;
          for (let b = 0; b < p.ids; b += block) {
            const m = Math.min(block, p.ids - b);
            writePoses(p, b, m, poses);
            write(poses.subarray(0, m * 8), (slotFirst + p.base + 2 * b) * 16);
          }
        }
      },
    };
  },
};

const m3 = Mat3();

const bitsFloat = new Float32Array(1);
const bitsUint = new Uint32Array(bitsFloat.buffer);
/** a float32's bits */
function floatBits(v: number) {
  bitsFloat[0] = v;
  return bitsUint[0];
}

/** whether `count` values of `a` from `first` count up by one */
function isSerial(a: Uint32Array, first: number, count: number) {
  const start = a[first];
  for (let k = 1; k < count; ++k) if (a[first + k] !== (start + k) >>> 0) return false;
  return true;
}

/**
 * `matrix` as a similarity, `k Rm` then a translation (a rotation, a uniform scale `k`, no mirroring, no projection),
 * the form a pose takes in; none where it is not one
 */
function getSimilarity(matrix: Mat4): CulledSimilarity | undefined {
  if (matrix[3] !== 0 || matrix[7] !== 0 || matrix[11] !== 0 || matrix[15] !== 1) return;
  const len = (o: number) => matrix[o] * matrix[o] + matrix[o + 1] * matrix[o + 1] + matrix[o + 2] * matrix[o + 2];
  const dot = (a: number, b: number) =>
    matrix[a] * matrix[b] + matrix[a + 1] * matrix[b + 1] + matrix[a + 2] * matrix[b + 2];
  const k2 = len(0);
  const eps = 1e-6 * k2;
  if (!(k2 > 0) || Math.abs(len(4) - k2) > eps || Math.abs(len(8) - k2) > eps) return;
  if (Math.abs(dot(0, 4)) > eps || Math.abs(dot(0, 8)) > eps || Math.abs(dot(4, 8)) > eps) return;
  Mat3.fromMat4(m3, matrix);
  if (Mat3.determinant(m3) <= 0) return;
  const k = Math.sqrt(k2);
  for (let i = 0; i < 9; ++i) m3[i] /= k;
  const rotation = Quat.normalize(Quat(), Quat.fromMat3(Quat(), m3));
  return { k, rotation, matrix };
}

/**
 * Ids `first` to `first + count` of `slots` as two vec4s each: `t = p - s R c` and `s` (1 not above 0), `R`'s
 * quaternion; the similarity folded in (`M t`, `k s`, `qm q`)
 */
function writePoses(slots: CulledPoseSlots, first: number, count: number, out: Float32Array) {
  const { positions, rotations, scales, center, indices, scale } = slots.poses;
  // the positions' kind once: float32, else the grid's, decoded inline (`GridPositions`)
  const pf = positions instanceof Float32Array ? positions : undefined;
  const g = positions instanceof Float32Array ? undefined : positions;
  const gd = g ? g.data : pf!;
  const gs = g ? g.step : 0;
  const [gox, goy, goz] = g ? g.origin : [0, 0, 0];
  const uniform = scale !== undefined && scale > 0 ? scale : 1;
  const elementFirst = slots.poses.first ?? 0;
  const { similarity } = slots;
  const [cx, cy, cz] = center;
  const m = similarity?.matrix;
  const k = similarity ? similarity.k : 1;
  const [ax, ay, az, aw] = similarity ? similarity.rotation : [0, 0, 0, 1];
  for (let j = 0; j < count; ++j) {
    const id = first + j;
    const e = indices ? indices[id] : elementFirst + id;
    const o = j * 8;
    let sc = scales ? scales[e] : uniform;
    if (!(sc > 0)) sc = 1;
    let x = 0,
      y = 0,
      z = 0,
      w = 1;
    if (rotations) {
      x = rotations[e * 4];
      y = rotations[e * 4 + 1];
      z = rotations[e * 4 + 2];
      w = rotations[e * 4 + 3];
      // unit: rotations of either kind (`InstancePoses`: snorm16 integers), the reader's float32 unit quaternion. Not
      // quantized on the GPU: the reader runs in every raster chunk and vertex, its decoding in their latency
      // (snorm16, 24 B: the yeast cell's sphere raster 5.9 -> 7.9 ms; 16 B, the translation 21 bits an axis and the
      // quaternion smallest-three, one aligned vec4: 7.7 ms)
      const n = Math.sqrt(x * x + y * y + z * z + w * w);
      if (n > 0) {
        x /= n;
        y /= n;
        z /= n;
        w /= n;
      } else {
        w = 1;
      }
    }
    // R c: c + 2 w (q × c) + 2 q × (q × c)
    const tx = 2 * (y * cz - z * cy);
    const ty = 2 * (z * cx - x * cz);
    const tz = 2 * (x * cy - y * cx);
    const rx = cx + w * tx + (y * tz - z * ty);
    const ry = cy + w * ty + (z * tx - x * tz);
    const rz = cz + w * tz + (x * ty - y * tx);
    let px: number, py: number, pz: number;
    if (pf) {
      px = pf[e * 3];
      py = pf[e * 3 + 1];
      pz = pf[e * 3 + 2];
    } else {
      const a = gd[e * 2],
        b = gd[e * 2 + 1];
      px = gox + gs * (a & 0x1fffff);
      py = goy + gs * ((a >>> 21) | ((b & 0x3ff) << 11));
      pz = goz + gs * ((b >>> 10) & 0x1fffff);
    }
    px -= sc * rx;
    py -= sc * ry;
    pz -= sc * rz;
    if (m) {
      const bx = px,
        by = py;
      px = m[0] * bx + m[4] * by + m[8] * pz + m[12];
      py = m[1] * bx + m[5] * by + m[9] * pz + m[13];
      pz = m[2] * bx + m[6] * by + m[10] * pz + m[14];
      // qm q (Hamilton product: R(qm q) = Rm R)
      const qx = aw * x + ax * w + ay * z - az * y;
      const qy = aw * y - ax * z + ay * w + az * x;
      const qz = aw * z + ax * y - ay * x + az * w;
      const qw = aw * w - ax * x - ay * y - az * z;
      x = qx;
      y = qy;
      z = qz;
      w = qw;
    }
    out[o] = px;
    out[o + 1] = py;
    out[o + 2] = pz;
    out[o + 3] = k * sc;
    out[o + 4] = x;
    out[o + 5] = y;
    out[o + 6] = z;
    out[o + 7] = w;
  }
}

/**
 * The segments' data of `CulledData` (`CulledSegmentVec4s` a segment): their groups and first item (u32 bits), their
 * quantized centers' origin and scale (from `tSegment`'s, `SegmentTexels`)
 */
export function getCulledSegmentData(
  segments: CulledSegmentGroups,
  tSegment: { readonly array: Float32Array } | undefined,
) {
  const out = new Float32Array(segments.count * CulledSegmentVec4s * 4);
  const bits = new Uint32Array(out.buffer);
  const signed = new Int32Array(out.buffer);
  const segment = tSegment?.array;
  for (let s = 0; s < segments.count; ++s) {
    const o = s * CulledSegmentVec4s * 4;
    // the groups as i32 bits: exact past 2^24
    signed[o] = segments.groupCounts[s];
    signed[o + 1] = segments.groupInstanceOffsets[s];
    signed[o + 2] = segments.groupBases[s];
    bits[o + 3] = segments.drawBases[s] / 6;
    const t = (s * SegmentTexels + 4) * 4;
    out[o + 4] = segment?.[t] ?? 0;
    out[o + 5] = segment?.[t + 1] ?? 0;
    out[o + 6] = segment?.[t + 2] ?? 0;
    out[o + 7] = segment?.[t + 3] ?? 1;
    out[o + 8] = segment?.[(s * SegmentTexels + 3) * 4 + 2] ?? 0;
    out[o + 9] = segment?.[(s * SegmentTexels + 3) * 4 + 3] ?? 0;
  }
  return out;
}
