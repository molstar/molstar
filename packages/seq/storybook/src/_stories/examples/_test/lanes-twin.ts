/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { LayerType } from '@molstar/seqx-core/render';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { vi } from 'vitest';

// The canvas and the SVG of one settled view lay their lanes out alike (core _spec/README.md "Render", Lanes and Hits:
// sticky lanes give way to an item packed first once the view rests): the SVG export's own layout is read through
// the layer types it is given (each type's `layout` wrapped to keep the lay it returns, by the layer's key), and each
// item's lane compared with the canvas's frame's. A lay with lanes is one with `ids` and `lanes` (the ranges types').

interface LanesLay {
  readonly items: number;
  readonly ids: readonly string[];
  readonly lanes: Int32Array;
}

const lanesOf = (lay: unknown): LanesLay | undefined => {
  const l = lay as Partial<LanesLay> | undefined;
  return l && Array.isArray(l.ids) && l.lanes instanceof Int32Array && typeof l.items === 'number'
    ? (l as LanesLay)
    : undefined;
};

/** Each item's lane, by its id (only items in a lane: −1 is out of the band or past the lanes). */
function laneMap(lay: LanesLay): Map<string, number> {
  const out = new Map<string, number>();
  for (let i = 0; i < lay.items; i++) if (lay.lanes[i] >= 0) out.set(lay.ids[i], lay.lanes[i]);
  return out;
}

/** What differs between the canvas's lanes and the SVG export's, and how many items were compared. */
export interface LanesTwin {
  readonly wrong: string[];
  readonly compared: number;
  /** The canvas's layers with an item in a lane, and those of them with no item in the SVG's lanes (nothing compared). */
  readonly laid: number;
  readonly uncompared: string[];
}

/** The settled viewer's lanes against its SVG export's (the export run on the fake clock). */
export async function lanesTwin(v: SeqViewer, id: string): Promise<LanesTwin> {
  // the SVG's lays, the last each layer's layout returned (the settled frame's)
  const svgLays = new Map<string, unknown>();
  const types: Record<string, LayerType<any, any>> = {};
  for (const [name, type] of Object.entries(v.context.types))
    types[name] = {
      ...type,
      layout: (prep, ctx, data, overlay, prev) => {
        const lay = type.layout(prep, ctx, data, overlay, prev);
        svgLays.set(ctx.key, lay);
        return lay;
      },
    };
  // the clock run until the export is done (a story's staggered answers), a turn of real I/O between (a gzip's inflate)
  let done = false;
  const svg = v.renderSvg({ idPrefix: 'lanes-', types }).finally(() => {
    done = true;
  });
  for (let k = 0; !done && k < 200; k++) {
    await vi.advanceTimersByTimeAsync(250);
    await new Promise((resolve) => setImmediate(resolve));
  }
  await svg;
  const wrong: string[] = [];
  let compared = 0;
  let laid = 0;
  const uncompared: string[] = [];
  for (const s of v.frame?.layout.sections ?? [])
    for (const t of s.tracks)
      for (const l of t.layers) {
        const canvas = lanesOf(l.lay);
        const svg = lanesOf(svgLays.get(l.key));
        if (!canvas || !svg) continue;
        const a = laneMap(canvas);
        const before = compared;
        const b = laneMap(svg);
        for (const [item, lane] of a) {
          const other = b.get(item);
          if (other === undefined) continue;
          compared++;
          if (other !== lane) wrong.push(`${id} ${l.key} ${item}: lane ${lane} on the canvas, ${other} in the SVG`);
        }
        if (a.size > 0) {
          laid++;
          if (compared === before) uncompared.push(l.key);
        }
      }
  // a row each: the first few and a count say enough
  return {
    wrong: wrong.length > 8 ? [...wrong.slice(0, 8), `and ${wrong.length - 8} more`] : wrong,
    compared,
    laid,
    uncompared,
  };
}
