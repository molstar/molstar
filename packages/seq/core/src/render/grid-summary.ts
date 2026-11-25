/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Location } from '@molstar/seqx-core/query';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer, ResolvedView } from '@molstar/seqx-core/view';
import { keeps, type Triangle } from './grid.js';

// A summary of a view's grids over a rect (_spec/README.md "Render", The grid layer, A summary over a rect).

/** What `gridSummary` gives: the mean, least and greatest value of the cells counted, and their count (NaN when none). */
export interface GridSummary {
  readonly mean: number;
  readonly min: number;
  readonly max: number;
  readonly count: number;
}

/**
 * A summary of a view's grids over a rect (_spec/README.md "Render", The grid layer, A summary over a rect): the cells
 * of every grid a grid layer of the view draws (each grid once) inside a rect of `rect` on the grid's two systems, each
 * cell once however many rects hold it, within the layer's `triangle`, missing values left out. A host's readout reads it (a selected rect's mean PAE); a library function, no state.
 */
export function gridSummary(
  view: { readonly model: StaticModel; readonly view: ResolvedView },
  rect: Location,
): GridSummary {
  const seen = new Set<string>();
  let sum = 0;
  let count = 0;
  let min = Number.NaN;
  let max = Number.NaN;
  const visit = (layer: ResolvedLayer) => {
    const owner = layer.data?.address;
    if (layer.type !== 'grid' || !owner || seen.has(owner)) return;
    seen.add(owner);
    const a = view.model.annotations.get(owner);
    const grid = a?.shape === 'grid' ? a.grid : undefined;
    if (!a || !grid) return;
    const t = layer.params.triangle;
    const triangle: Triangle = t === 'upper' || t === 'lower' ? t : 'both';
    // a cell once however many rects hold it (rects added by a drag overlap): a mark a cell of this grid
    let counted: Uint8Array | undefined;
    for (const r of rect.rects) {
      if (r.a.system !== a.target || r.b.system !== a.ySystem) continue;
      counted ??= new Uint8Array(grid.nx * grid.ny);
      const xs = r.a.ranges;
      const ys = r.b.ranges;
      for (let p = 0; p + 1 < xs.length; p += 2)
        for (let i = Math.max(0, xs[p]); i < Math.min(grid.nx, xs[p + 1]); i++)
          for (let q = 0; q + 1 < ys.length; q += 2)
            for (let j = Math.max(0, ys[q]); j < Math.min(grid.ny, ys[q + 1]); j++) {
              const k = i * grid.ny + j;
              if (counted[k] || !keeps(triangle, i, j)) continue;
              counted[k] = 1;
              const v = grid.values[k];
              if (Number.isNaN(v)) continue;
              sum += v;
              count++;
              if (!(v >= min)) min = v;
              if (!(v <= max)) max = v;
            }
    }
  };
  for (const s of view.view.sections) {
    if (s.layout !== 'grid') continue;
    for (const t of s.tracks) for (const l of t.layers) visit(l);
    for (const b of s.blocks ?? []) for (const ls of b.layers) for (const l of ls) visit(l);
  }
  return { mean: count ? sum / count : Number.NaN, min, max, count };
}
