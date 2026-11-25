/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Key, type Location, rangesLocation } from '@molstar/seqx-core/query';
import type { Hit, HitTarget } from './layer-type.js';
import { TreeHitPx } from './rows-model.js';
import type { RowsLay } from './rows-lay.js';

// The rows layer's hits (_spec/README.md "Render", Lanes and Hits; The rows layer): a row's, a cell's, a summary row's
// or a tree node's, written into the layer's scratch hit (a hover allocates nothing).

/** A hit on display row `i`: a summary row's (its clade), else its data row's, a cell where `position` is given. */
export function rowHit(lay: RowsLay, i: number, position: number | undefined): Hit {
  lay.hitRow = i;
  const e = lay.display.rowAt(i);
  if (e < 0 && lay.clades) {
    const c = -1 - e;
    return writeHit(
      lay,
      position === undefined ? lay.rowsSystem : lay.system,
      position,
      undefined,
      'summary-row',
      cladeHit(lay, c, i),
    );
  }
  if (position === undefined) return writeHit(lay, lay.rowsSystem, e, e, 'row');
  return writeHit(lay, lay.system, position, e, 'cell');
}

/**
 * The inner tree node whose branch is nearest a point within `TreeHitPx` (by the distance to the segment), a folded
 * clade's wedge included; -1 for none.
 */
export function nearestNode(lay: RowsLay, x: number, y: number): number {
  const reach = TreeHitPx;
  // the target stops at the label column: a point over a row's label is the row's (a cherry at the tree's edge would
  // take it, the reach being more than the 4 px between the tree and the labels)
  if (x < lay.treeX0 - reach || x > lay.treeX1 + reach || x >= lay.labelX) return -1;
  const w = lay.wedges;
  for (let k = 0; k + 4 < w.length; k += 5)
    if (x >= w[k] - reach && x <= w[k + 2] && y >= w[k + 1] && y <= w[k + 3]) return w[k + 4];
  let best = -1;
  let bestD = reach;
  const t = lay.tree;
  const owner = lay.treeNode;
  for (let k = 0, s = 0; k + 3 < t.length; k += 4, s++) {
    const node = owner[s];
    if (node < 0) continue;
    const xa = t[k] < t[k + 2] ? t[k] : t[k + 2];
    const xb = t[k] < t[k + 2] ? t[k + 2] : t[k];
    const ya = t[k + 1] < t[k + 3] ? t[k + 1] : t[k + 3];
    const yb = t[k + 1] < t[k + 3] ? t[k + 3] : t[k + 1];
    const dx = x < xa ? xa - x : x > xb ? x - xb : 0;
    const dy = y < ya ? ya - y : y > yb ? y - yb : 0;
    const d = Math.sqrt(dx * dx + dy * dy);
    if (d <= bestD) {
      bestD = d;
      best = node;
    }
  }
  return best;
}

/**
 * The hit written into the layer's scratch (a hover allocates nothing): every field set. A summary row's or a tree
 * node's element is its clade's, its extent the clade's rows (kept by the display).
 */
export function writeHit(
  lay: RowsLay,
  system: string,
  position: number | undefined,
  yPosition: number | undefined,
  target: HitTarget,
  clade?: { element: Key; extent: Location },
): Hit {
  const h = lay.scratch;
  h.layer = lay.key;
  h.section = lay.section;
  h.system = system;
  h.position = position;
  h.ySystem = yPosition === undefined ? undefined : lay.rowsSystem;
  h.yPosition = yPosition;
  h.target = target;
  h.element = clade?.element;
  h.extent = clade?.extent;
  h.state = undefined;
  return h;
}

/** A header's hit naming the column cell it is on (`title: value`, made at layout: a hover allocates nothing). */
export function named(lay: RowsLay, h: Hit, q: number): Hit {
  if (q >= 0 && h === lay.scratch) {
    const c = lay.cols;
    lay.scratch.state = c.name[c.col[q]][c.row[q]];
  }
  return h;
}

/**
 * A clade's element and its rows as a Location on the rows system (its shown rows when folded at display row `i`),
 * made once a display.
 */
export function cladeHit(lay: RowsLay, c: number, i?: number): { element: Key; extent: Location } {
  if (lay.cladeHitsOf !== lay.display.version) {
    lay.cladeHitsOf = lay.display.version;
    lay.cladeHits.clear();
    lay.summaryHits.clear();
  }
  const kept = i === undefined ? lay.cladeHits.get(c) : lay.summaryHits.get(i);
  if (kept) return kept;
  const members = i === undefined ? undefined : lay.display.members.get(i);
  const ranges: number[] = [];
  if (members) for (const r of members) ranges.push(r, r + 1);
  else if (lay.clades)
    for (let k = lay.clades.lo[c]; k < lay.clades.hi[c]; k++) {
      const r = lay.clades.leafRow[k];
      if (r >= 0) ranges.push(r, r + 1);
    }
  const made = {
    element: { owner: lay.owner, key: lay.clades?.name[c] ?? '' },
    extent: rangesLocation(lay.rowsSystem, ranges),
  };
  if (i === undefined) lay.cladeHits.set(c, made);
  else lay.summaryHits.set(i, made);
  return made;
}
