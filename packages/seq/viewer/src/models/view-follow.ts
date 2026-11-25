/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { inversePath, type Path, projectEnd, projectPoint } from '@molstar/seqx-core/resolve';
import { type Axis, type Follow, followWindow, type ResolvedSection, type ResolvedView } from '@molstar/seqx-core/view';
import { axisTotal, nearestOffsetOf, positionAt } from './axis.js';

// The view's follows (_spec/README.md "ViewModel", Follow): the sections in an order where each leader comes before
// its followers (EC-114); a move asked of a follower carried back to its leader (`follow.followerMoves: lead`), each
// end through the inverse of the forward map; a follower's range, its leader's carried along the follow path and
// clamped like a leader's.

/** How many follows a move asked of a follower goes up at most (a chain of follows; a cycle is cut). */
export const MaxFollowDepth = 16;

/** A record `leaderRange` writes its result into: the caller's own (a view model's), reused so a move makes none. */
export interface LeaderRangeOut {
  start: number;
  end: number;
}

/** A record `followerRange` writes its result into: the caller's own (a view model's), reused so a publish makes none. */
export interface FollowerRangeOut {
  start: number;
  end: number;
  pieces: readonly number[] | undefined;
}

/**
 * The leader's range (offsets, lowest first) a move asked of a follower carries back to (_spec/README.md "ViewModel",
 * Follow): the follower's asked `[start, end)` (offsets) is carried back through the follow, each end to the leader
 * place whose image it is, so the follower shows what was asked and every section moves in step; on a minus strand
 * the follower's ends are the images of the leader's other ends, so its pan moves the leader the other way. `ms`,
 * `me`: which way each end moved. An end that maps nowhere (an intron, a UTR) takes the nearest place that maps the
 * way it moved (`backTo`). A pan (`delta`) moves the follower's unclamped image of the leader, not the window its own
 * `maxPxPerPosition` widened (which would zoom the leader out), and when one end maps nowhere the leader moves by the
 * other's step: the leader keeps its width. A zoom or a range whose end lies past every image outward reaches the
 * leader's axis edge there (a zoom out never stalls); inward, it keeps the leader's end. `lead`: the leader's viewport
 * now. Undefined when the leader cannot move so; the result is `out`, the caller's own record (never shared).
 */
export function leaderRange(
  section: ResolvedSection,
  lead: { readonly start: number; readonly end: number },
  start: number,
  end: number,
  ms: number,
  me: number,
  delta: number | undefined,
  arc: 'shortest' | 'hull',
  out: LeaderRangeOut,
): LeaderRangeOut | undefined {
  const follow = section.follow!;
  const axis = section.axis;
  const total = axisTotal(axis);
  if (total === undefined) return undefined;
  const la = follow.leaderAxis;
  const l0 = positionAt(la, lead.start, 'start');
  const l1 = positionAt(la, lead.end, 'end');
  if (l0 === undefined || l1 === undefined) return undefined;
  // the leader's window in positions, lowest first (a reversed leader's start is its highest)
  const ps = Math.min(l0, l1);
  const pe = Math.max(l0, l1);
  if (delta !== undefined) {
    // the image of the leader's window before the follower's clamp, moved by the pan
    const w = followWindow(section, la, axis, lead.start, lead.end, arc);
    if (w && !w.wraps && w.start < w.end) [start, end] = [w.start + delta, w.end + delta];
  }
  // a cyclic follower's range over its origin runs past the total
  const pos = (o: number, side: 'start' | 'end') => positionAt(axis, axis.cyclic && o > total ? o - total : o, side);
  const f0 = pos(start, 'start');
  const f1 = pos(end, 'end');
  // the follower's ends in positions, lowest first, with the way each moved there (navigation is in display order:
  // on a reversed follower an end moved right moved down)
  const rev = axis.dir === -1;
  const lo = rev ? f1 : f0;
  const hi = rev ? f0 : f1;
  const mLo = rev ? -me : ms;
  const mHi = rev ? -ms : me;
  const inv = inverseOf(follow.path);
  const xLo = lo === undefined ? undefined : backTo(follow, inv, axis, lo, mLo);
  const xHi = hi === undefined ? undefined : backTo(follow, inv, axis, hi, mHi);
  // the leader's lowest position is the preimage of the follower's lowest, or of its highest on a minus strand
  let a = follow.dir === 1 ? xLo : xHi;
  let b = follow.dir === 1 ? xHi : xLo;
  if (delta !== undefined) {
    if (a === undefined && b === undefined) return undefined;
    if (a === undefined) a = ps + (b! - pe);
    if (b === undefined) b = pe + (a - ps);
  } else {
    // an end moved outward past every image: the leader's axis edge on that side (its clamp keeps it in)
    const lowOut = follow.dir === 1 ? mLo < 0 : mHi > 0;
    const highOut = follow.dir === 1 ? mHi > 0 : mLo < 0;
    const lt = axisTotal(la);
    const e0 = lt === undefined ? undefined : positionAt(la, 0, 'start');
    const e1 = lt === undefined ? undefined : positionAt(la, lt, 'end');
    const lowEdge = e0 === undefined || e1 === undefined ? ps : Math.min(e0, e1);
    const highEdge = e0 === undefined || e1 === undefined ? pe : Math.max(e0, e1);
    a ??= lowOut ? lowEdge : ps;
    b ??= highOut ? highEdge : pe;
  }
  const oa = nearestOffsetOf(la, a);
  const ob = nearestOffsetOf(la, b);
  if (oa === undefined || ob === undefined || oa === ob) return undefined;
  // a cyclic follower's range over its end (past the total) is the leader's arc over its own: on from the first end,
  // never the rest of the ring
  const lt = axisTotal(la);
  if (end > total && axis.cyclic && la.cyclic && lt !== undefined && follow.dir === 1 && axis.dir === 1 && ob < oa) {
    out.start = oa;
    out.end = ob + lt;
    return out;
  }
  out.start = Math.min(oa, ob);
  out.end = Math.max(oa, ob);
  return out;
}

/** A follow's path read back, made once per path. */
const inverses = new WeakMap<Path, Path>();
function inverseOf(p: Path): Path {
  let inv = inverses.get(p);
  if (!inv) inverses.set(p, (inv = inversePath(p)));
  return inv;
}

/** Just short of an element's far edge: the limit of its continuous image there, not the next element's start. */
const BeforeEdge = 1 - 1e-9;

/**
 * The leader position whose image (`Follow.point`, what the follower's ends are made of) is the follower's position
 * `y`: the inverse of the forward map itself, so a follower moved through it shows what was asked. The forward map is
 * place by place (a residue's three fine places across a codon step: a codon split by an intron is its nucleotides),
 * each place's image a run of `point(k) + f · slope`, `f` in [0, 1); the places tried are those near `y`'s preimage
 * through `inv` (the follow read back) and near the nearest mapped positions on either side. Where `y` is in no place's
 * image (an intron, a UTR), the nearest image the way the end moved (`m`), so a pan steps over a gap and never sticks.
 * Undefined when nothing lies that way.
 */
function backTo(follow: Follow, inv: Path, axis: Axis, y: number, m: number): number | undefined {
  // a place's width in leader positions: a third of a residue across a codon step
  const u = follow.path.step === 1 ? 1 / 3 : 1;
  const tried = new Set<number>();
  let best: number | undefined;
  let distance = Infinity;
  const near = (x: number | undefined) => {
    if (x === undefined) return undefined;
    const k0 = Math.floor(x / u + 1e-9);
    for (let k = k0 - 1; k <= k0 + 1; k++) {
      if (tried.has(k)) continue;
      tried.add(k);
      const a = follow.point(k * u);
      const half = follow.point((k + 0.5) * u);
      if (a === undefined || half === undefined || half === a) continue;
      const slope = 2 * (half - a);
      const f = (y - a) / slope;
      if (f >= 0 && f < 1) return (k + f) * u;
      if (m === 0) continue;
      // the end of this place's image nearer `y`, when it lies the way the end moved
      const b = a + slope;
      const v = m > 0 ? Math.min(a, b) : Math.max(a, b);
      const d = (v - y) * m;
      if (d > 0 && d < distance) {
        distance = d;
        best = (k + (v === a ? 0 : BeforeEdge)) * u;
      }
    }
    return undefined;
  };
  return (
    near(projectPoint(inv, y)) ??
    near(nearestMapped(inv, axis, y, 1)) ??
    near(nearestMapped(inv, axis, y, -1)) ??
    (m === 0 ? undefined : best)
  );
}

/**
 * The leader place of the follower's nearest mapped position past `y` the way `m` says, within the follower's axis:
 * where the next one starts (`m` > 0), or the far edge of the last one before (`m` < 0).
 */
function nearestMapped(inv: Path, axis: Axis, y: number, m: number): number | undefined {
  const r = axis.ranges;
  if (!r || r.length < 2) return undefined;
  let lo = Infinity;
  let hi = -Infinity;
  for (let i = 0; i + 1 < r.length; i += 2) {
    lo = Math.min(lo, r[i]);
    hi = Math.max(hi, r[i + 1]);
  }
  // the inverse's source is in thirds of a residue when it steps up to nucleotides
  const scale = inv.step === 1 ? 3 : 1;
  if (m > 0) {
    const from = Math.ceil(y);
    if (from >= hi) return undefined;
    const run = inv.linear.forward(from * scale, hi * scale)[0];
    return run ? projectPoint(inv, Math.ceil(run.from / scale)) : undefined;
  }
  const to = Math.floor(y);
  if (to <= lo) return undefined;
  const runs = inv.linear.forward(lo * scale, to * scale);
  const last = runs[runs.length - 1];
  return last ? projectEnd(inv, Math.floor((last.from + last.length - 1) / scale) + 1) : undefined;
}

/**
 * A follower's range: its leader's (`lead`) carried along the follow path at its scale (1, 3 or 1/3). Each end is the
 * continuous image of the leader's (`point` at the start, the far edge of the last element at a whole end), so it
 * moves smoothly within a run; where an end maps nowhere (an intron, a UTR), the hull of the whole positions' image
 * (`followRange`). Over a cyclic follower's origin the range wraps (`pieces`), or with `follow.cyclicArc: 'hull'` is
 * the hull. Clamped by `maxPx` (px a position) like a leader once it has a `width`. Its range before (`prev`) while it
 * cannot be placed. The result is `out`, the caller's own record.
 */
export function followerRange(
  section: ResolvedSection,
  lead: { readonly start: number; readonly end: number },
  prev: { readonly start: number; readonly end: number; readonly pieces?: readonly number[] } | undefined,
  width: number,
  arc: 'shortest' | 'hull',
  maxPx: number,
  out: FollowerRangeOut,
): FollowerRangeOut {
  const follow = section.follow!;
  const axis = section.axis;
  const total = axisTotal(axis);
  let start = prev?.start ?? 0;
  let end = prev?.end ?? total ?? 0;
  let pieces: readonly number[] | undefined = prev?.pieces;
  const w = total === undefined ? undefined : followWindow(section, follow.leaderAxis, axis, lead.start, lead.end, arc);
  if (w && total !== undefined) {
    if (w.wraps) {
      // over the origin: on to the total, then from 0
      start = w.start;
      end = total + w.end;
      pieces = [w.start, total, 0, w.end];
    } else {
      start = w.start;
      end = w.end;
      pieces = undefined;
      if (width > 0) [start, end] = clampRange(start, end, total, width, maxPx);
    }
  }
  out.start = start;
  out.end = end;
  out.pieces = pieces;
  return out;
}

/** At least `width / maxPx` wide (within the total), centered, within the axis. */
function clampRange(start: number, end: number, total: number, width: number, maxPx: number): [number, number] {
  const w = Math.min(Math.max(end - start, Math.min(width / maxPx, total)), total);
  if (w === end - start) return [start, end];
  const c = (start + end) / 2;
  const s = Math.min(Math.max(c - w / 2, 0), total - w);
  return [s, s + w];
}

/** The view's sections with every leader before its followers (a chain of follows, EC-114); a cycle is cut. */
export function followOrder(view: ResolvedView | undefined): string[] {
  const out: string[] = [];
  const done = new Set<string>();
  const byAddress = new Map((view?.sections ?? []).map((s) => [s.address, s] as const));
  const visit = (s: ResolvedSection, seen: Set<string>) => {
    if (done.has(s.address) || seen.has(s.address)) return;
    seen.add(s.address);
    const leader = s.follow && byAddress.get(s.follow.leader);
    if (leader) visit(leader, seen);
    done.add(s.address);
    out.push(s.address);
  };
  for (const s of view?.sections ?? []) visit(s, new Set());
  return out;
}
