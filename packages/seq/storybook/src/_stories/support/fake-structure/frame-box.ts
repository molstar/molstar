/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The camera frame's place (fake-view.ts): the bounding box of the cells it is on, in the coordinates of the element the
// frame sits in (its content box, scrolled): client rectangles against the container's, never the cells' own offsets,
// which are relative to whatever positioned element each is in.

export interface Rect {
  readonly left: number;
  readonly top: number;
  readonly right: number;
  readonly bottom: number;
}

/** The frame, `pad` around the cells' union, as `{ left, top, width, height }` inside the container. */
export function frameBox(
  cells: readonly Rect[],
  container: Rect,
  scroll: { readonly left: number; readonly top: number },
  pad: number,
) {
  const left = Math.min(...cells.map((c) => c.left)) - container.left + scroll.left - pad;
  const top = Math.min(...cells.map((c) => c.top)) - container.top + scroll.top - pad;
  const right = Math.max(...cells.map((c) => c.right)) - container.left + scroll.left + pad;
  const bottom = Math.max(...cells.map((c) => c.bottom)) - container.top + scroll.top + pad;
  return { left, top, width: right - left, height: bottom - top };
}
