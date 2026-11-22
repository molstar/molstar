/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { RangesLay, RowsLay, ViewLayout } from '@molstar/seqx-core/render';
import type { ViewModel } from '@molstar/seqx-viewer/models';

// A rows band's moves the `rows` and `cursor` behaviors share (_spec/README.md "Behaviors", Rows and Cursor): its
// folds, kept with the row at its top, and its pages. Internal: not in the behaviors' entry.

/** New folds on a band, the row at the band's top kept there (a fold above it moves it up, never off). */
export function setBandFolds(view: ViewModel, l: RowsLay | undefined, track: string, names: readonly string[]): void {
  const at = l?.keyAt(l.top);
  const top = l && at ? l.topOf(at.key, at.offset, names) : undefined;
  view.setBand(track, top === undefined ? { folds: names } : { folds: names, bandTop: top });
}

/** The folds a band shows: the session's, else its layer's `collapsed` as laid out. */
export function bandFolds(view: ViewModel, track: string, l: RowsLay): readonly string[] {
  return view.band(track)?.folds ?? l.folds;
}

/** The scrolled rows a page moves (PageDown, PageUp): the rows the band's window shows less one, at least one. */
export function pageRows(l: RowsLay): number {
  return Math.max(1, Math.floor((l.box.height - l.display.pinned * l.rowScale) / l.rowScale) - 1);
}

/**
 * The band's top (display rows) a page on: the scrolled rows in its window (the pinned band's rows are not scrolled)
 * less one, clamped to the band's ends.
 */
export function pagedTop(l: RowsLay, top: number, dir: 1 | -1): number {
  return Math.max(0, Math.min(l.maxTop(l.rowScale, l.display.count), top + dir * pageRows(l)));
}

/**
 * The band's top that shows display row `i` whole from `top` by the least move (a pinned row is always shown): `top`
 * when it is in view already, clamped to the band's ends.
 */
export function revealedTop(l: RowsLay, top: number, i: number): number {
  const P = l.display.pinned;
  if (i < P) return top;
  const area = Math.max(0, l.box.height - P * l.rowScale) / l.rowScale;
  const r = i - P;
  let next = top;
  if (r < top) next = r;
  else if (r + 1 > top + area) next = r + 1 - area;
  return Math.max(0, Math.min(next, l.maxTop(l.rowScale, l.display.count)));
}

/** The ranges layer `key`'s last layout in a layout (a lanes band). */
export function lanesOf(layout: ViewLayout | undefined, key: string): RangesLay | undefined {
  for (const s of layout?.sections ?? [])
    for (const band of s.lays)
      for (const l of band ?? []) if (l.key === key && l.type.name === 'ranges') return l.lay as RangesLay;
  return undefined;
}

/** The lanes a page moves (Page Up, Page Down): the lanes the band shows less one, at least one. */
export function pageLanes(l: RangesLay): number {
  return Math.max(1, Math.floor(l.room) - 1);
}

/** The lanes band's top (lanes) a page on, or at its first (`-Infinity`) or last (`Infinity`) lane: within its lanes. */
export function lanesTop(l: RangesLay, top: number, by: number): number {
  return Math.max(0, Math.min(l.scrollMax, top + by));
}
