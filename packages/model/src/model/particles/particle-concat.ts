/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column } from '@molstar/core/data/db';
import { GridPositions, type Positions } from '@molstar/core/math/geometry/grid-positions';
import { type ParticleAttribute, type ParticleList, ParticleRotations, ParticleValues } from './particle-list.js';

// Particles of lists concatenated (a tile view, a slot list flattened): their segments, each some of a list's particles
// in an order, gathered field by field

/** some particles of a list: `indices` `[start, end)` (a tile's particles of a target: its indices by target) */
export interface Segment {
  readonly tile: number;
  readonly list: ParticleList;
  readonly indices: Int32Array;
  readonly start: number;
  readonly end: number;
}

/** the segment holding view index `i` (`starts` ascending, the count last) */
export function findSegment(starts: Int32Array, i: number) {
  let lo = 0;
  let hi = starts.length - 2;
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    if (starts[mid] <= i) lo = mid;
    else hi = mid - 1;
  }
  return lo;
}

/** the segments' positions: on their grid where all are on one (the words copied: exact), else float32 */
export function concatPositions(segments: readonly Segment[], count: number): Positions {
  const grids = segments.map((s) => s.list.coordinates);
  const g = grids[0];
  const oneGrid =
    GridPositions.is(g) &&
    grids.every(
      (p) =>
        GridPositions.is(p) &&
        p.step === g.step &&
        p.origin[0] === g.origin[0] &&
        p.origin[1] === g.origin[1] &&
        p.origin[2] === g.origin[2],
    );
  if (oneGrid) {
    const data = new Uint32Array(count * 2);
    let o = 0;
    for (const { list, indices, start, end } of segments) {
      const src = (list.coordinates as GridPositions).data;
      for (let i = start; i < end; ++i, ++o) {
        data[o * 2] = src[indices[i] * 2];
        data[o * 2 + 1] = src[indices[i] * 2 + 1];
      }
    }
    return { kind: 'grid', data, origin: g.origin, step: g.step };
  }
  const out = new Float32Array(count * 3);
  let o = 0;
  for (const { list, indices, start, end } of segments) {
    const src = GridPositions.toFloat(list.coordinates);
    for (let i = start; i < end; ++i, ++o) {
      const k = indices[i] * 3;
      out[o * 3] = src[k];
      out[o * 3 + 1] = src[k + 1];
      out[o * 3 + 2] = src[k + 2];
    }
  }
  return out;
}

/**
 * the segments' rotations: as they are stored where all are of one kind (snorm16 or float32), else float32 (a snorm16's
 * components over 32767; none the identity); none where no segment has any
 */
export function concatRotations(segments: readonly Segment[], count: number): ParticleRotations | undefined {
  if (!segments.some((s) => s.list.rotations)) return undefined;
  const snorm = segments.every((s) => s.list.rotations instanceof Int16Array);
  const out = snorm ? new Int16Array(count * 4) : new Float32Array(count * 4);
  let o = 0;
  for (const { list, indices, start, end } of segments) {
    const r = list.rotations;
    const scale = r instanceof Int16Array && !snorm ? 1 / ParticleRotations.Snorm16One : 1;
    for (let i = start; i < end; ++i, ++o) {
      if (!r) {
        out[o * 4 + 3] = snorm ? ParticleRotations.Snorm16One : 1;
        continue;
      }
      const q = indices[i] * 4;
      for (let c = 0; c < 4; ++c) out[o * 4 + c] = r[q + c] * scale;
    }
  }
  return out;
}

/** the segments' values (`get` of each one's list), as runs where fewer than half */
export function concatValues<A extends Int32Array | Float32Array>(
  segments: readonly Segment[],
  count: number,
  get: (list: ParticleList) => ParticleValues<A>,
  ctor: { new (n: number): A },
): ParticleValues<A> {
  const out = new ctor(count);
  let o = 0;
  for (const { list, indices, start, end } of segments) {
    const v = get(list);
    if (!ParticleValues.isRuns(v)) {
      for (let i = start; i < end; ++i, ++o) out[o] = v[indices[i]];
      continue;
    }
    // a segment's indices ascending: its runs walked along them, from its first's (a search)
    const { offsets, values } = v;
    if (end === start) continue;
    let r = ParticleValues.runOf(v, indices[start]);
    for (let i = start; i < end; ++i, ++o) {
      const k = indices[i];
      while (offsets[r + 1] <= k) r += 1;
      out[o] = values[r];
    }
  }
  return ParticleValues.compact(out);
}

/** a list's scale for all its particles (`ParticleList.scale`), 1 where none or not positive */
function getScale(list: ParticleList) {
  return list.scale !== undefined && list.scale > 0 ? list.scale : 1;
}

/** a list's particles' scales (`ParticleList.scales`, else its `scale` a run of all) */
function getScaleValues(list: ParticleList): ParticleValues<Float32Array> {
  return list.scales ?? ParticleValues.ofRuns(Int32Array.of(0, list.count), Float32Array.of(getScale(list)));
}

/**
 * the segments' lists' scales concatenated (a coarse node's particles drawn larger: `ParticleList.scale`, its runs'
 * `scales`), by particle; none where no list has one
 */
export function concatScales(segments: readonly Segment[], count: number) {
  if (!segments.some(({ list }) => list.scales || getScale(list) !== 1)) return undefined;
  return concatValues(segments, count, getScaleValues, Float32Array);
}

/** the attributes every segment's list has, concatenated; their label, range from the first */
export function concatAttributes(
  segments: readonly Segment[],
  count: number,
): Map<string, ParticleAttribute> | undefined {
  const first = segments[0]?.list.attributes;
  if (!first) return undefined;
  const out = new Map<string, ParticleAttribute>();
  for (const [name, a] of first) {
    if (!segments.every((s) => s.list.attributes?.has(name))) continue;
    const values = new Float32Array(count);
    let o = 0;
    let min = a.min;
    let max = a.max;
    for (const { list, indices, start, end } of segments) {
      const attribute = list.attributes!.get(name)!;
      min = Math.min(min, attribute.min);
      max = Math.max(max, attribute.max);
      for (let i = start; i < end; ++i, ++o) values[o] = attribute.column.value(indices[i]);
    }
    out.set(name, { label: a.label, min, max, column: Column.ofFloatArray(values) });
  }
  return out.size ? out : undefined;
}
