/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { NavWindowFlags } from '@molstar/seqx-core/render';

// Where a navigator part's readout goes (_spec/README.md "Composition", Navigators): never over the control in use, a
// window, its grips or the pointer. DOM-free: the host measures the readout's box and hands it in with the avoided boxes.

/** A box in the part's px: `[x0, x1) × [y0, y1)`. */
export interface Box {
  readonly x0: number;
  readonly y0: number;
  readonly x1: number;
  readonly y1: number;
}

/** The gap kept between the readout and what it avoids, and from the part's edge, px. */
export const ReadoutGap = 4;
/** How far a window's box is widened for its grips, px (the least a grip's hit reaches). */
export const GripReach = 8;
/** The pointer's box: its half side, px. */
export const PointerReach = 12;

/** The part of a window as drawn (the frame's `NavigatorWindow`): only the fields the boxes read. */
export interface WindowBounds {
  readonly shown: boolean;
  readonly axes: number;
  readonly x0: number;
  readonly x1: number;
  readonly dx0: number;
  readonly dx1: number;
  readonly wrap: boolean;
  readonly w0: number;
  readonly w1: number;
  readonly y0: number;
  readonly y1: number;
  readonly dy0: number;
  readonly dy1: number;
}

/** The boxes the shown windows cover in a part of `width` × `height`, widened by the grips' reach (a cyclic window's two pieces). */
export function windowBoxes(windows: readonly WindowBounds[], width: number, height: number): Box[] {
  const out: Box[] = [];
  for (const w of windows) {
    if (!w.shown) continue;
    const y0 = w.axes & NavWindowFlags.y ? w.dy0 - GripReach : 0;
    const y1 = w.axes & NavWindowFlags.y ? w.dy1 + GripReach : height;
    const spans: [number, number][] =
      w.axes & NavWindowFlags.x
        ? w.wrap
          ? [
              [w.w0, w.w1],
              [w.x0, w.x1],
            ]
          : [[w.dx0, w.dx1]]
        : [[0, width]];
    for (const [a, b] of spans) {
      const grip = w.axes & NavWindowFlags.x ? GripReach : 0;
      out.push({ x0: a - grip, x1: b + grip, y0, y1 });
    }
  }
  return out;
}

/** The pointer's box, none while it is away. */
export function pointerBox(x: number, y: number): Box | undefined {
  return Number.isNaN(x) || Number.isNaN(y)
    ? undefined
    : { x0: x - PointerReach, x1: x + PointerReach, y0: y - PointerReach, y1: y + PointerReach };
}

const hits = (a: Box, b: Box) => a.x0 < b.x1 && b.x0 < a.x1 && a.y0 < b.y1 && b.y0 < a.y1;

/**
 * The readout's top left in the part's px, for a readout `w` × `h` in a part of `width` × `height` that touches no
 * window box and not the pointer's. In order: the top right corner (where it rests), the top left, beside the windows
 * (right of them, then left of them), above the part, below it; the first clear one (below when none is).
 */
export function placeReadout(
  w: number,
  h: number,
  width: number,
  height: number,
  windows: readonly Box[],
  pointer?: Box,
): readonly [left: number, top: number] {
  let x0 = Number.POSITIVE_INFINITY;
  let x1 = Number.NEGATIVE_INFINITY;
  for (const b of windows) {
    if (b.x0 < x0) x0 = b.x0;
    if (b.x1 > x1) x1 = b.x1;
  }
  const g = ReadoutGap;
  const right = Math.max(g, width - w - g);
  const candidates: (readonly [number, number])[] = [
    [right, g],
    [g, g],
    [x1 + g, g],
    [x0 - g - w, g],
    [right, -h - g],
    [right, height + g],
  ];
  const last = candidates.length - 1;
  for (let i = 0; i <= last; i++) {
    const [left, top] = candidates[i];
    // inside the part: whole within its width (outside it, above and below, only the corner matters)
    if (i < 4 && !(left >= g - 1e-9 && left + w <= width - g + 1e-9)) continue;
    const box: Box = { x0: left, x1: left + w, y0: top, y1: top + h };
    if (pointer && hits(box, pointer)) continue;
    let clear = true;
    for (const b of windows) {
      if (hits(box, b)) {
        clear = false;
        break;
      }
    }
    if (clear) return candidates[i];
  }
  return candidates[last];
}
