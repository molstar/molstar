/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A liftover's words (_spec/README.md "Hover card", A liftover's places): an item's places through a map read by region
// (the core's `regionProjector`), and a lifted item's strand as carried (turned by an inverted run), read from the
// layer's rows: its far table's strand is the old assembly's.

import type { LayerLayout, RangesLay } from '@molstar/seqx-core/render';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { RegionProjector, ResolvedView } from '@molstar/seqx-core/view';

/** An item's places through a map read by region, when several. */
export function liftPlaces(
  regions: RegionProjector,
  model: StaticModel,
  view: ResolvedView | undefined,
  system: string | undefined,
  start: number,
  end: number,
): string | undefined {
  return view && system !== undefined ? regions.places?.(model, view, system, start, end) : undefined;
}

/** A lifted item's strand as drawn (`+`, `−`, '' when none): its row starting at `start`, else its first. */
export function liftedStrand(l: LayerLayout, key: string | number, start: number): string {
  const rows = (l.lay as RangesLay).rows;
  let at = -1;
  for (let j = 0; j < rows.n; j++)
    if (rows.keys[j] === key) {
      if (at < 0 || rows.starts[j] === start) at = j;
      if (rows.starts[j] === start) break;
    }
  const s = at < 0 ? 0 : rows.strand[at];
  return s < 0 ? '−' : s > 0 ? '+' : '';
}
