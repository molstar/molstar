/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A follower's window and range from its leader's (_spec/README.md "Resolved view", Follow links): the leader's
// positions carried along the follow path, the hull or the shortest arc on a cyclic follower.

import { HalfOpenRanges } from '@molstar/core/data/int/half-open-ranges';
import { projectEnd, projectRanges } from '@molstar/seqx-core/resolve';
import { type AxisShape, axisTotal, nearestOffsetOf, positionAt } from './axis.js';
import type { ResolvedSection } from './view.js';

/**
 * A leader's continuous positions `[lo, hi]` (the edges of its window, lowest first) carried to the follower's
 * positions, `[lo, hi]` lowest first: the low edge's continuous image (`point`) and the high one's (`projectEnd`, the
 * far edge of the place before it), so a follower moves smoothly and is exact at element edges. Where an edge maps
 * nowhere (an intron, a UTR), the hull of the whole positions' image there (`followRange`). Ends past each other: over a
 * cyclic follower's origin the wrap (`lo > hi`), or its hull with `arc: 'hull'`; on a linear follower the hull (the
 * leader's range wraps the coding sequence's origin). Undefined when nothing maps.
 */
export function followSpan(
  section: ResolvedSection,
  lo: number,
  hi: number,
  arc: 'shortest' | 'hull' = 'shortest',
): [number, number] | undefined {
  const follow = section.follow;
  if (!follow) return undefined;
  let a = follow.point(lo);
  let b = projectEnd(follow.path, hi);
  if (a === undefined || b === undefined) {
    const hull = followRange(section, [Math.floor(lo), Math.ceil(hi)], { arc });
    if (!hull || hull.length < 2) return undefined;
    const first = hull[0];
    const last = hull[hull.length - 1];
    if (a === undefined) a = follow.dir === 1 ? first : last;
    if (b === undefined) b = follow.dir === 1 ? last : first;
  }
  const [x, y] = follow.dir === 1 ? [a, b] : [b, a];
  if (x > y || (x === y && hi > lo)) {
    if (section.axis.cyclic && arc !== 'hull') return [x, y];
    const hull = followRange(section, [Math.floor(lo), Math.ceil(hi)], { arc: 'hull' });
    return hull && hull.length >= 2 ? [hull[0], hull[hull.length - 1]] : undefined;
  }
  return [x, y];
}

/** A follower's window in its own offsets: `start` below `end`, or over a cyclic origin `start > end` (`wraps`). */
export interface FollowWindow {
  readonly start: number;
  readonly end: number;
  readonly wraps: boolean;
}

/**
 * A follower's window from its leader's (`start`, `end`: offsets along the leader's axis, display order): the images of
 * the leader's window edges (`followSpan`) as offsets along the follower's axis, in its display order (a reversed
 * follower's lowest position is its window's end). An offset on a leader's break is read as a window's edge
 * (`positionAt`'s sides). `leader` and `axis` are the axes with their pieces (a pending one placed by the caller).
 */
export function followWindow(
  section: ResolvedSection,
  leader: AxisShape,
  axis: AxisShape,
  start: number,
  end: number,
  arc: 'shortest' | 'hull' = 'shortest',
): FollowWindow | undefined {
  // a cyclic leader's window may pass its total: it wraps the origin (a ring's window, a rotated linear view)
  const total = axisTotal(leader);
  const wraps = leader.cyclic && total !== undefined && end > total && end - start < total;
  const ps = positionAt(leader, start, 'start');
  const pe = positionAt(leader, wraps ? end - total : end, 'end');
  if (ps === undefined || pe === undefined) return undefined;
  // over the origin the positions run from the start's past the length to the end's: low edge past the high
  const lo = wraps ? (leader.dir === -1 ? pe : ps) : Math.min(ps, pe);
  const hi = wraps ? (leader.dir === -1 ? ps : pe) : Math.max(ps, pe);
  const span = followSpan(section, lo, hi, arc);
  if (!span) return undefined;
  const a = nearestOffsetOf(axis, span[0]);
  const b = nearestOffsetOf(axis, span[1]);
  if (a === undefined || b === undefined) return undefined;
  const [s, e] = axis.dir === -1 ? [b, a] : [a, b];
  if (s > e && axis.cyclic) return { start: s, end: e, wraps: true };
  return { start: Math.min(s, e), end: Math.max(s, e), wraps: false };
}

/**
 * A follower's range from its leader's, both as `[start, end, …]` pieces in display order (an axis's shape). A leader
 * piece that wraps its cyclic origin (`start > end`) is split on it first. The image's hull, or on a cyclic follower
 * the shortest arc that holds it (EC-12), as two pieces when it wraps (in the follower's display order: on a reversed
 * axis the piece below the origin first); `arc: 'hull'` takes the hull there too. Undefined when nothing maps.
 */
export function followRange(
  section: ResolvedSection,
  leader: readonly number[],
  options: { readonly arc?: 'shortest' | 'hull' } = {},
): number[] | undefined {
  const follow = section.follow;
  if (!follow) return undefined;
  const pieces: number[] = [];
  const ll = follow.leaderAxis.length;
  for (let i = 0; i + 1 < leader.length; i += 2) {
    const s = Math.floor(leader[i]);
    const e = Math.ceil(leader[i + 1]);
    if (s > e && follow.leaderAxis.cyclic && ll !== undefined) pieces.push(s, ll, 0, e);
    else pieces.push(s, e);
  }
  const image = HalfOpenRanges.normalize(projectRanges(follow.path, pieces));
  const n = image.length;
  if (!n) return undefined;
  const { cyclic, length } = section.axis;
  if (!cyclic || length === undefined || options.arc === 'hull') return [image[0], image[n - 1]];
  // the arc is the circle less its widest gap: around the origin, or between two pieces
  let best = length - image[n - 1] + image[0];
  let arc = [image[0], image[n - 1]];
  for (let i = 2; i < n; i += 2) {
    const gap = image[i] - image[i - 1];
    if (gap > best) {
      best = gap;
      arc = [image[i], length, 0, image[i - 1]];
    }
  }
  return arc.length === 4 && section.axis.dir === -1 ? [arc[2], arc[3], arc[0], arc[1]] : arc;
}
