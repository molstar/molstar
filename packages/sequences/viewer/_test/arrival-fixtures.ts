/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Arrivals under a test's control, for per-region redraws (arrivals.test.ts on the recording backend,
// canvas/band-copy.test.ts in Chromium): a lazy alignment's rows layer, its stream wrapped so the fine level's tiles
// over some tile columns are held back (pending) and then released one column at a time, each release logged with its
// columns as the stream manager logs a tile's arrival, so the next frame redraws those columns alone. The rest is the
// real stream's. Hold columns before any composer draws them: a drawn tile's image is kept by the renderer.

import { TILE } from "@molstar/seq-core/format";
import type { SeqViewerContext } from "../src/index.ts";
import { blockLevels, lanesInView, RegionState, type StreamView } from "../src/streams.ts";

export function holdBack(ctx: SeqViewerContext, section: string) {
  const m = ctx.graph.assembled.value?.view;
  const av = ctx.layout.axis(section);
  const l = ctx.layout.get(section)?.current;
  if (!m || !av || !l) throw new Error("no view");
  const k = m.prepared.findIndex((p) => p.renderer === "rows");
  const ref = m.set.strings.get(m.resolved.layers.ref[k]);
  const streams = ctx.streams as unknown as { get(ref: string): StreamView | null };
  const get = streams.get.bind(streams);
  const real = get(ref);
  if (!real) throw new Error("no stream");
  // The fine level the rows draw: the stream manager's (`blockLevels` over the lanes in view).
  const lanes = new Float64Array(8);
  lanesInView(l, av, lanes);
  const lv = new Int32Array(4);
  blockLevels(av, 0, real.budget, lanes, lv);
  const held = new Set<number>();
  let bumps = 0;
  const span = [0, 0];
  const hides = (ly: number, lx: number, tx: number) => ly === lv[0] && lx === lv[1] && held.has(tx);
  const proxy = Object.create(real, { version: { get: () => real.version + bumps } }) as StreamView;
  Object.assign(proxy, {
    state: (ly: number, lx: number, ty: number, tx: number) =>
      hides(ly, lx, tx) ? RegionState.pending : real.state(ly, lx, ty, tx),
    tile: (ly: number, lx: number, ty: number, tx: number) => (hides(ly, lx, tx) ? null : real.tile(ly, lx, ty, tx)),
    // One change since the version drawn: the held column's; else unknown (everything).
    changed: (since: number, out: Int32Array) => {
      if (since !== proxy.version - 1) return -1;
      out[0] = span[0];
      out[1] = span[1];
      return 1;
    },
  });
  streams.get = (r: string) => (r === ref ? proxy : get(r));
  return {
    /** The fine level (ly, lx) and the coarse one drawn in place of a held tile. */
    levels: lv,
    /** The fine tile column under the view's center. */
    center: () => av.indexOf(0, ((av.visLo[0] + av.visHi[0]) >> 1) as never) >> (8 + lv[1]),
    hold: (tx: number) => held.add(tx),
    /** Column tx arrives: one change, its columns. */
    release: (tx: number) => {
      held.delete(tx);
      span[0] = (tx * TILE) << lv[1];
      span[1] = ((tx + 1) * TILE) << lv[1];
      bumps++;
    },
  };
}
