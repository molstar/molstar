/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StaticModel } from '@molstar/seqx-core/resolve';
import {
  followRange,
  followWindow,
  offsetOf,
  piecesTotal,
  type ResolvedSection,
  type ResolvedView,
} from '@molstar/seqx-core/view';
import { AxisTransform, type Viewport } from './axis.js';
import type { RenderConfig } from './config.js';
import type { StoreReader } from './layer-type.js';

// A section's axis in the tracks layout (_spec/README.md "Render", Tracks layout): its own pieces, a follower's pieces
// and window from its leader's, a cyclic axis turned past its origin, and which layout and systems a section has.

/**
 * A section's own pieces: its axis's (the array itself: lookups are kept per pieces array), or the whole system once its
 * source gave the length; undefined while unknown.
 */
export function ownPieces(
  section: ResolvedSection,
  model: StaticModel,
  data: StoreReader,
): readonly number[] | undefined {
  if (section.axis.ranges) return section.axis.ranges;
  const system = model.systems.get(section.axis.system);
  if (!system?.source) return undefined;
  const info = data
    .coverage(system.source)
    ?.systems.find((s) => s.name === system.name || s.aliases?.includes(system.name));
  return info?.length === undefined ? undefined : placedPieces(info.length);
}

/** The whole system's one piece for a length, kept per length (lookups are kept per pieces array). */
const placed = new Map<number, readonly number[]>();
function placedPieces(length: number): readonly number[] {
  let p = placed.get(length);
  if (!p) {
    if (placed.size > 64) placed.clear();
    placed.set(length, (p = [0, length]));
  }
  return p;
}

/**
 * A follower's pieces and window from its leader's view (_spec/README.md "Render", Tracks layout): the images of the
 * leader's window edges as offsets along the follower's own pieces (`followWindow`, as the viewer's follow), so a
 * follower moves with its leader and keeps its direction and breaks. A window over a cyclic follower's origin is drawn
 * over the arc's pieces (`followRange`).
 */
export function followAxis(
  section: ResolvedSection,
  own: readonly number[] | undefined,
  leader: AxisTransform,
  window: Viewport,
  config: RenderConfig,
) {
  const axis = { ...section.axis, ...(own ? { ranges: own } : {}) };
  // a cyclic leader's window past its system's end (a linear view across its origin, its pieces twice) wraps
  const turn = leader.cyclic && leader.length !== undefined ? Math.min(leader.total, leader.length) : leader.total;
  if (own && window.end <= turn) {
    const w = followWindow(section, leader, axis, window.start, window.end, config.followArc);
    if (w && !w.wraps) return { pieces: own, viewport: { start: w.start, end: w.end } };
  }
  const pieces = followRange(section, windowPositions(leader, window), { arc: config.followArc });
  return pieces ? { pieces } : undefined;
}

/**
 * The positions a window of a leader's display units shows, `[start, end, …]`: the transform's own window as it is; a
 * cyclic window past the total wraps (two runs, the second from the origin).
 */
function windowPositions(leader: AxisTransform, w: Viewport): number[] {
  if (w.start === leader.u0 && w.end === leader.u1) return leader.visible([]);
  const out: number[] = [];
  const at = (start: number, end: number) =>
    new AxisTransform(leader.system, leader.pieces, 0, 1, 1, { start, end }, { dir: leader.dir }).visible(out);
  const total = leader.total;
  if (leader.cyclic && w.end > total) {
    at(w.start, total);
    at(0, Math.min(w.end - total, w.start));
  } else at(w.start, w.end);
  return out;
}

const Twice = new WeakMap<readonly number[], readonly number[]>();

/**
 * A cyclic axis's pieces under a window past its total (_spec/README.md "Render", Tracks layout: a linear view panned
 * across its origin): the whole system twice, so the window's units go on over the origin as the second turn's (two
 * runs, no break between); else the pieces as they are. Kept per pieces array (a pan makes none).
 */
export function turned(
  section: ResolvedSection,
  own: readonly number[],
  viewport: Viewport | undefined,
): readonly number[] {
  const { cyclic, length } = section.axis;
  if (!cyclic || !viewport || own.length !== 2 || own[0] !== 0 || own[1] !== length || viewport.end <= own[1])
    return own;
  let t = Twice.get(own);
  if (!t) Twice.set(own, (t = [own[0], own[1], own[0], own[1]]));
  return t;
}

/**
 * A following ring's window on its own axis (_spec/README.md "Render", Rings layout: a ring that follows a linear
 * section has the leader's range as its window): the follow's window as it is, else the first piece of the image its
 * offsets on, as long as the image (over the origin: past the total).
 */
export function followedWindow(
  section: ResolvedSection,
  own: readonly number[],
  made: { pieces: readonly number[]; viewport?: Viewport },
): Viewport {
  if (made.viewport) return made.viewport;
  const shape = { ...section.axis, ranges: own };
  const p = made.pieces;
  const len = piecesTotal(p);
  // the image's first piece in display order: its low edge forward, its high edge reversed
  const start = offsetOf(shape, section.axis.dir === -1 ? p[1] : p[0]) ?? 0;
  return { start, end: start + len };
}

/** A section's layout: its own `layout`, else its view's. */
export function sectionLayoutName(section: ResolvedSection, view: ResolvedView): string {
  return section.layout ?? view.layout;
}

/** Whether a section's blocks are on more than one system (_spec/README.md "Render", Blocks: it neither follows nor leads). */
export function multiSystem(section: ResolvedSection): boolean {
  const b = section.blocks;
  if (!b) return false;
  for (let i = 1; i < b.length; i++) if (b[i].axis.system !== b[0].axis.system) return true;
  return false;
}
