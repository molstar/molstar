/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Column } from '@molstar/core/data/db';
import type { PlyFile, PlyTable } from '@molstar/io/reader/ply/schema';
import {
  createShAccessor,
  type GaussianSplat,
  type GaussianSplatData,
  type GaussianSplatReadOptions,
  getShDegree,
  SH_C0,
  setCenterRotation,
  ShBasisCounts,
  sigmoid,
} from './data.js';

/** Splats a chunk of the compressed layout */
const ChunkSize = 256;

function table(plyFile: PlyFile, name: string) {
  return plyFile.getElement(name) as PlyTable | undefined;
}

function hasAll(t: PlyTable | undefined, names: readonly string[]): t is PlyTable {
  return !!t && t.kind === 'table' && names.every((n) => t.propertyNames.includes(n));
}

const StandardProperties = [
  'x',
  'y',
  'z',
  'f_dc_0',
  'f_dc_1',
  'f_dc_2',
  'opacity',
  'scale_0',
  'scale_1',
  'scale_2',
  'rot_0',
  'rot_1',
  'rot_2',
  'rot_3',
] as const;
const PackedProperties = ['packed_position', 'packed_rotation', 'packed_scale', 'packed_color'] as const;
const ChunkProperties = [
  'min_x',
  'min_y',
  'min_z',
  'max_x',
  'max_y',
  'max_z',
  'min_scale_x',
  'min_scale_y',
  'min_scale_z',
  'max_scale_x',
  'max_scale_y',
  'max_scale_z',
] as const;
const ChunkColorProperties = ['min_r', 'min_g', 'min_b', 'max_r', 'max_g', 'max_b'] as const;

/** The compressed layout of PlayCanvas and SuperSplat: a `chunk` element of bounds, `packed_*` vertex properties */
export function isCompressedGaussianSplatPly(plyFile: PlyFile) {
  return hasAll(table(plyFile, 'vertex'), PackedProperties) && hasAll(table(plyFile, 'chunk'), ChunkProperties);
}

/** A 3D Gaussian splatting file: the standard layout (3DGS training output) or the compressed one */
export function isGaussianSplatPly(plyFile: PlyFile) {
  return isCompressedGaussianSplatPly(plyFile) || hasAll(table(plyFile, 'vertex'), StandardProperties);
}

/** The standard layout: logit opacity, log scales, a rotation `rot_0..3` as `w, x, y, z`, the base color by `f_dc` */
function standard(vertex: PlyTable): GaussianSplatData['getSplat'] {
  const [x, y, z, dc0, dc1, dc2, opacity, s0, s1, s2, r0, r1, r2, r3] = StandardProperties.map(
    (n) => vertex.getProperty(n)! as Column<number>,
  );
  return (i: number, out: GaussianSplat) => {
    setCenterRotation(out, x.value(i), y.value(i), z.value(i), r1.value(i), r2.value(i), r3.value(i), r0.value(i));
    out.logScaleX = s0.value(i);
    out.logScaleY = s1.value(i);
    out.logScaleZ = s2.value(i);
    out.opacity = sigmoid(opacity.value(i));
    out.r = 0.5 + SH_C0 * dc0.value(i);
    out.g = 0.5 + SH_C0 * dc1.value(i);
    out.b = 0.5 + SH_C0 * dc2.value(i);
  };
}

function unorm(value: number, bits: number) {
  const max = (1 << bits) - 1;
  return (value & max) / max;
}

/**
 * The compressed layout (PlayCanvas `gsplat-compressed-data.js`): per chunk of 256 splats the bounds of the positions
 * and log scales (and optionally the colors); per splat four `u32`: position and log scale as 11-10-11 bits within
 * the bounds, the rotation as the smallest three (2-10-10-10: the index of the largest component, 0 for `w`; the
 * first field `w` otherwise; the others as `(v - 0.5) * √2`), the color and opacity as bytes
 */
function compressed(vertex: PlyTable, chunk: PlyTable): GaussianSplatData['getSplat'] {
  const [packedPosition, packedRotation, packedScale, packedColor] = PackedProperties.map(
    (n) => vertex.getProperty(n)! as Column<number>,
  );
  if (chunk.rowCount * ChunkSize < vertex.rowCount) throw new Error('fewer chunks than the splats need');
  const hasColorBounds = ChunkColorProperties.every((n) => chunk.propertyNames.includes(n));
  const boundNames = hasColorBounds ? [...ChunkProperties, ...ChunkColorProperties] : ChunkProperties;
  const stride = boundNames.length;
  // the bounds a chunk (256 splats), gathered once: not a copy of the splats'
  const bounds = new Float32Array(chunk.rowCount * stride);
  boundNames.forEach((n, j) => {
    const column = chunk.getProperty(n)!;
    for (let c = 0; c < chunk.rowCount; ++c) bounds[c * stride + j] = column.value(c);
  });
  const lerp = (o: number, t: number) => bounds[o] + (bounds[o + 3] - bounds[o]) * t;
  const sqrt2 = Math.SQRT2;
  return (i: number, out: GaussianSplat) => {
    const b = Math.floor(i / ChunkSize) * stride;
    const pos = packedPosition.value(i);
    const rot = packedRotation.value(i);
    const ra = (unorm(rot >>> 20, 10) - 0.5) * sqrt2;
    const rb = (unorm(rot >>> 10, 10) - 0.5) * sqrt2;
    const rc = (unorm(rot, 10) - 0.5) * sqrt2;
    const rm = Math.sqrt(Math.max(0, 1 - (ra * ra + rb * rb + rc * rc)));
    const largest = rot >>> 30;
    // x, y, z, w by the largest component's index (0 w, 1 x, 2 y, 3 z)
    const qx = largest === 0 ? ra : largest === 1 ? rm : rb;
    const qy = largest === 0 ? rb : largest === 2 ? rm : largest === 1 ? rb : rc;
    const qz = largest === 3 ? rm : rc;
    const qw = largest === 0 ? rm : ra;
    setCenterRotation(
      out,
      lerp(b, unorm(pos >>> 21, 11)),
      lerp(b + 1, unorm(pos >>> 11, 10)),
      lerp(b + 2, unorm(pos, 11)),
      qx,
      qy,
      qz,
      qw,
    );
    const scale = packedScale.value(i);
    // log scales lerped as stored
    out.logScaleX = lerp(b + 6, unorm(scale >>> 21, 11));
    out.logScaleY = lerp(b + 7, unorm(scale >>> 11, 10));
    out.logScaleZ = lerp(b + 8, unorm(scale, 11));
    const color = packedColor.value(i);
    // the opacity a byte, not by the bounds
    out.opacity = unorm(color, 8);
    const r = unorm(color >>> 24, 8);
    const g = unorm(color >>> 16, 8);
    const bl = unorm(color >>> 8, 8);
    out.r = hasColorBounds ? lerp(b + 12, r) : r;
    out.g = hasColorBounds ? lerp(b + 13, g) : g;
    out.b = hasColorBounds ? lerp(b + 14, bl) : bl;
  };
}

/**
 * The view-dependent color of `t` up to `degree` from its `f_rest_*` properties (channel-major: the red coefficients,
 * the green, the blue; compressed: bytes of `v * 8 / 255 - 4`)
 */
function sh(t: PlyTable | undefined, degree: number, isCompressed: boolean) {
  if (!t || degree <= 0) return;
  const names = t.propertyNames.filter((n) => n.startsWith('f_rest_'));
  const d = Math.min(degree, getShDegree(names.length));
  if (d <= 0) return;
  const perChannel = names.length / 3;
  const n = ShBasisCounts[d as 1 | 2 | 3];
  const columns: Column<number>[] = [];
  for (let c = 0; c < 3; ++c) {
    for (let k = 0; k < n; ++k) columns.push(t.getProperty(`f_rest_${c * perChannel + k}`)! as Column<number>);
  }
  const get = isCompressed
    ? (i: number, k: number, c: number) => (columns[c * n + k].value(i) * 8) / 255 - 4
    : (i: number, k: number, c: number) => columns[c * n + k].value(i);
  return createShAccessor(d, get);
}

/** The gaussian splats of `plyFile` (`isGaussianSplatPly`), standard or compressed; stored y-down */
export function gaussianSplatsFromPly(plyFile: PlyFile, options: GaussianSplatReadOptions): GaussianSplatData {
  const vertex = table(plyFile, 'vertex');
  if (!vertex) throw new Error('missing vertex element');
  const isCompressed = isCompressedGaussianSplatPly(plyFile);
  return {
    count: vertex.rowCount,
    upright: false,
    getSplat: isCompressed ? compressed(vertex, table(plyFile, 'chunk')!) : standard(vertex),
    sh: sh(isCompressed ? table(plyFile, 'sh') : vertex, options.shDegree, isCompressed),
  };
}
