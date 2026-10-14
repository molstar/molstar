/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST3b's fixtures (streams-progress.test.ts and its heap guards): the streamed genome composed at DPR 2, 1,000 × 300,
// with the tracks asked for, its providers gated (regions answered on `release`), failing on `fail`, and a frame clock
// of its own (`clock[0]`, ms: the overlay's shimmer and fades read the scheduler's `now`, set here allocation-free).

import type { Loaders, StreamRequest } from "@molstar/seq-core/format";
import { GRCH38, hg38Streams } from "../../core/_test/format/stream-provider.js";
import { SkinToken } from "../src/draw/backend.js";
import { Dirty } from "../src/draw/composer.js";
import { RegionState, SeqViewerConfig } from "../src/index.js";
import type { Streams } from "../src/streams.js";
import type { RecordingBackend } from "../src/testing/index.js";
import { composed, recording } from "./composer-fixtures.js";
import { GENOME, streamedDoc } from "./stream-fixtures.js";
export interface ProgressOptions {
  items?: boolean;
  sequence?: boolean;
  values?: boolean;
  /** The composer is a surface's: skeletons shimmer, arrivals fade. */
  live?: boolean;
}

/** The genome at [start, stop) of chr2 (gated: nothing answered until `release`); `failing` answers with errors. */
export async function progressGenome(start: number, stop: number, o: ProgressOptions = {}) {
  /** Every data request fails while `on`; or those `when` picks (ST4b). */
  const failing = { on: false, when: null as ((r: StreamRequest) => boolean) | null };
  const s = hg38Streams({ chromosomes: GRCH38.slice(0, 3), letters: true, gated: true });
  // Every provider's data requests fail at once while `failing` is on.
  const wrap = (p: (typeof s)["genome"]["provider"]): typeof p => ({
    fetch: (r, sig) =>
      r.system !== null && (failing.on || failing.when?.(r))
        ? Promise.reject(new Error("503 from the fixture"))
        : p.fetch(r, sig),
  });
  const loaders: Loaders = {
    load: () => Promise.reject(new Error("none")),
    stream: (src) => {
      const p = s.stream(src);
      return p ? wrap(p) : null;
    },
  };
  const t = await composed(streamedDoc({ items: o.items, sequence: o.sequence, values: o.values }), {
    width: 1000,
    height: 300,
    dpr: 2,
    loaders,
  });
  t.composer.live = o.live ?? false;
  // The scheduler's clock, read through an array (a double set on a closure would be boxed).
  const clock = new Float64Array(1);
  (t.scheduler as { now: () => number }).now = () => clock[0];
  const m = t.ctx.streams as Streams;
  const vp = t.vp;
  if (!vp) throw new Error("no viewport");
  vp.zoomTo(1, start, stop);
  t.layout(GENOME);
  const be = recording();
  const flush = () => new Promise((ok) => setTimeout(ok, 0));
  const providers = [s.genome, s.phylop, s.genes, s.variants];
  return {
    ...t,
    s,
    m,
    vp,
    clock,
    backend: be,
    failing,
    /** A pass (what the view needs asked), then every held request answered, and their patches applied. */
    async release() {
      m.pass();
      for (const p of providers) p.release();
      await flush();
      await flush();
    },
    asked: () => providers.reduce((n, p) => n + p.requests.length, 0),
    /** Pans of ±9 px with every region answered, so a pan's pass finds the margin held (the guards' warm state). */
    async settle() {
      for (let k = 0; k < 8; k++)
        for (const dx of [9, -9]) {
          vp.panBy(dx, 0, "pointer");
          t.layout(GENOME);
          await this.release();
        }
    },
    /** The frame drawn whole (tracks and overlay) into the recording. */
    frame(flags: number = Dirty.all) {
      be.reset();
      t.composer.frame(be, flags);
      return be;
    },
  };
}

export type ProgressCase = "pending" | "landed" | "failed";
/**
 * The heap guards' overlay frames (streams-progress*.heap.ts): a live composer over the genome with genes, variants,
 * letters and a signal at 100 kb in state `c` (landed: after a second of pending frames, as a page sees it; the clock
 * then held, so the fade never ends); `next` counts the frames the surface would ask for (`animating`).
 */
export async function overlayFrame(c: ProgressCase) {
  const g = await progressGenome(100_000_000, 100_100_000, { items: true, sequence: true, values: true, live: true });
  g.ctx.config.set(SeqViewerConfig.Stream.MaxInFlight, 64);
  g.clock[0] = 1000;
  if (c === "failed") g.failing.on = true;
  g.m.pass();
  await new Promise((ok) => setTimeout(ok, 0));
  if (c === "landed") {
    for (let k = 0; k < 64; k++) {
      g.clock[0] += 16;
      g.frame(Dirty.overlay);
    }
    await g.release();
  }
  g.frame(Dirty.tracks);
  const next = new Int32Array(1);
  const frame = () => {
    // The shimmer's clock moves a frame; a landed region's stays.
    if (c === "pending") g.clock[0] += 16;
    g.frame(Dirty.overlay);
    if (g.composer.animating) next[0]++;
  };
  return { g, frame, next };
}

/** The rects of `color` (a skin token), and their alphas, in a recording. */
export function rectsOf(r: RecordingBackend, color: number): { n: number; alphas: number[] } {
  const ri = r.rectItems;
  const alphas: number[] = [];
  for (let i = 0; i < ri.n; i++) if (ri.color[i] === color) alphas.push(ri.alpha[i]);
  return { n: alphas.length, alphas };
}
export { Dirty, RegionState, SkinToken };
