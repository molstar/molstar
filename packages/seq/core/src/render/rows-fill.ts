/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedTrack } from '@molstar/seqx-core/view';
import { RowHeight } from './rows-model.js';

// A rows band's height when it fills (_spec/README.md "Render", The rows layer: the band's height).

/** The rows a track's rows layer shows, when its alignment's length is known (an inline one's). */
export function rowsCount(track: ResolvedTrack, model: StaticModel): number | undefined {
  for (const l of track.layers) {
    if (l.type !== 'rows' || !l.data) continue;
    const rows = model.annotations.get(l.data.address)?.rows;
    return rows === undefined ? undefined : model.systems.get(rows)?.length;
  }
  return undefined;
}

/** The display rows a band can show at most: every row, and a pinned row again in its place (`rows.pinInPlace`). */
export function rowsShown(track: ResolvedTrack, model: StaticModel, pinInPlace: boolean): number | undefined {
  const n = rowsCount(track, model);
  if (n === undefined || !pinInPlace) return n;
  const pinned = track.layers.find((l) => l.type === 'rows')?.params.pinned;
  return n + (Array.isArray(pinned) ? pinned.length : 0);
}

/**
 * A filling rows band's height: the room left on the page, at least `least` (its default height), never taller than
 * every display row (`rowsShown`) at `RowHeight` with its title row (`titleRow`) when the rows are counted. Why the
 * cap: on a page as tall as its content (`canvas.height: 'content'`) the room left is the last frame's content, so a
 * band taking it grew by the sections below it on every frame and left an empty well under its rows (the query's
 * alignment, a resize).
 */
export function rowsFillHeight(least: number, room: number, rows: number | undefined, titleRow: number): number {
  const h = Math.max(least, Math.floor(room));
  return rows === undefined ? h : Math.min(h, rows * RowHeight + titleRow);
}
