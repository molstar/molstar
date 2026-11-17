/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An axis's offsets (_spec/README.md "ViewModel"): a viewport is held in offsets along the axis in display order (0
// at the first position shown), so a range of the system and a cyclic range over the origin read alike. Positions and
// offsets are fractional; these convert between them piece by piece.

import type { Axis } from '@molstar/seqx-core/view';

/** The positions shown, or undefined while the axis waits for its system's length. */
export function axisTotal(axis: Axis): number | undefined {
  const r = axis.ranges;
  if (!r) return undefined;
  let n = 0;
  for (let i = 0; i + 1 < r.length; i += 2) n += r[i + 1] - r[i];
  return n;
}

/** The position at an offset (clamped to the axis); undefined while the axis is pending. */
export function positionAt(axis: Axis, offset: number): number | undefined {
  const r = axis.ranges;
  if (!r || r.length < 2) return undefined;
  let acc = 0;
  for (let i = 0; i + 1 < r.length; i += 2) {
    const w = r[i + 1] - r[i];
    if (offset <= acc + w || i + 2 >= r.length) return r[i] + Math.min(Math.max(offset - acc, 0), w);
    acc += w;
  }
  return undefined;
}

/** The offset of a position, or undefined when the axis does not show it (or is pending). */
export function offsetOf(axis: Axis, position: number): number | undefined {
  const r = axis.ranges;
  if (!r) return undefined;
  let acc = 0;
  for (let i = 0; i + 1 < r.length; i += 2) {
    if (position >= r[i] && position <= r[i + 1]) return acc + position - r[i];
    acc += r[i + 1] - r[i];
  }
  return undefined;
}

/** The offset of a position, the nearest piece's edge when the axis does not show it. */
export function nearestOffsetOf(axis: Axis, position: number): number | undefined {
  const exact = offsetOf(axis, position);
  if (exact !== undefined) return exact;
  const r = axis.ranges;
  if (!r || r.length < 2) return undefined;
  let best = 0;
  let distance = Infinity;
  let acc = 0;
  for (let i = 0; i + 1 < r.length; i += 2) {
    const ds = Math.abs(position - r[i]);
    const de = Math.abs(position - r[i + 1]);
    if (ds < distance) [best, distance] = [acc, ds];
    acc += r[i + 1] - r[i];
    if (de < distance) [best, distance] = [acc, de];
  }
  return best;
}

/** Whether two axes show the same positions of the same system. */
export function axesEqual(a: Axis, b: Axis): boolean {
  if (a.system !== b.system || a.length !== b.length || a.cyclic !== b.cyclic) return false;
  const ra = a.ranges;
  const rb = b.ranges;
  if (!ra || !rb) return ra === rb;
  if (ra.length !== rb.length) return false;
  for (let i = 0; i < ra.length; i++) if (ra[i] !== rb[i]) return false;
  return true;
}

/** A key equal for equal axes: the scroll group of sections on them (EC-113). */
export function axisKey(axis: Axis): string {
  return `${axis.system}|${axis.length ?? '?'}|${axis.ranges?.join(',') ?? '?'}`;
}
