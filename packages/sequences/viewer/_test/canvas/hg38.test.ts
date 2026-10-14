/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST1's check (ST1's Perf): an hg38-sized genome streamed through the port, GRCh38's 24 chromosomes (3.1 G positions)
// from a letters stream's manifest and a values stream over them whose bins are sampled on demand
// (core/_test/format/stream-provider.ts: nothing held by the provider), on a real surface in Chromium at 1,320 × 300 CSS
// px, DPR 1 and 2. At the whole view, 1 Mb, 10 kb and 100 positions, once the view's regions have arrived
// (`streams.settled`), a full frame stays under G1's 16 ms; a tour of 10 kb windows over every chromosome keeps the
// cache's bytes under a small `Stream.CacheBytes`. ST3a: with a sequence track, its letters streamed (cells thinner
// than a device px at 10 kb, letters at 100 positions; none asked at 1 Mb or the whole view). `HG38=1` prints the
// numbers to __debug__/hg38.txt. ST3a-2: genes and variants streamed too (density from `count` bins zoomed out, items
// region by region below their first level).

import { afterAll, expect, it } from "vitest";
import { commands } from "vitest/browser";
import { GRCH38, hg38Streams, streamedDoc } from "../../../core/_test/format/stream-provider.js";
import { mountSection } from "../../src/canvas/index.js";
import {
  createSeqViewer,
  DefaultSeqViewerSpec,
  type ManualScheduler,
  manualScheduler,
  SeqViewerConfig,
  type SeqViewerContext,
} from "../../src/index.js";
import type { Streams } from "../../src/streams.js";

const SECTION = "browser.main";
const median = (xs: number[]) => xs.slice().sort((a, b) => a - b)[xs.length >> 1];
const out: string[] = [];
afterAll(() => commands.writeFile("__debug__/hg38.txt", out.join("\n")));

/** The clock run until every region in view has arrived; the frames it took. */
async function settle(ctx: SeqViewerContext, manual: ManualScheduler): Promise<number> {
  let done = false;
  void ctx.streams.settled().then(() => {
    done = true;
  });
  let k = 0;
  for (; k < 400 && !done; k++) {
    manual.advance(16);
    await new Promise((ok) => setTimeout(ok, 0));
  }
  if (!done) throw new Error("the streams did not settle");
  return k;
}

for (const dpr of [1, 2])
  it(`a streamed hg38: whole view, 1 Mb and 10 kb under 16 ms a frame, the cache bounded (DPR ${dpr})`, async () => {
    const manual = manualScheduler(0);
    const s = hg38Streams({ letters: true });
    const ctx = createSeqViewer(DefaultSeqViewerSpec(), {
      scheduler: manual,
      loaders: { load: () => Promise.reject(new Error("the bench loads nothing")), stream: s.stream },
    });
    const t0 = performance.now();
    const r = await ctx.setDocument(streamedDoc({ sequence: true, items: true }));
    if (!r.ok) throw new Error(JSON.stringify(r.issues));
    const host = document.createElement("div");
    host.style.cssText = "position:absolute;left:0;top:0;width:1320px;height:300px";
    document.body.appendChild(host);
    const surface = mountSection(ctx, SECTION, host, { dpr });
    await ctx.graph.update();
    const vp = ctx.viewports.get(SECTION);
    if (!vp) throw new Error("no viewport");
    vp.fit();
    await settle(ctx, manual);
    surface.renderNow();
    const load = performance.now() - t0;
    const canvas = host.querySelector("canvas") as HTMLCanvasElement;
    const full = () => {
      surface.composer.redraw();
      surface.invalidate(7);
      const t = performance.now();
      surface.renderNow();
      canvas.getContext("2d")?.getImageData(0, 0, 1, 1);
      return performance.now() - t;
    };
    const manager = ctx.streams as Streams;
    const rows = [`load ${load.toFixed(0)} ms`];
    for (const [label, span] of [
      ["whole", 0],
      ["1 Mb", 1_000_000],
      ["10 kb", 10_000],
      ["100 b", 100],
    ] as const) {
      // 100 positions clear of the synthetic genome's N runs (each 100 kb's first 1 kb).
      const at = span === 100 ? 100_040_100 : 100_000_000;
      if (span > 0) vp.zoomTo(1, at, at + span);
      else vp.fit();
      const frames = await settle(ctx, manual);
      for (let i = 0; i < 5; i++) full();
      const times = Array.from({ length: 20 }, full);
      const p0 = performance.now();
      for (let i = 0; i < 100; i++) manager.pass();
      const pass = (performance.now() - p0) / 100;
      rows.push(
        `${label}: full ${median(times).toFixed(1)} ms (max ${Math.max(...times).toFixed(1)}), pass ${pass.toFixed(3)} ms, settled in ${frames} frames`,
      );
      expect(median(times), `${label} at DPR ${dpr}`).toBeLessThan(16);
    }
    // A tour: 10 kb windows over every chromosome, each settled, under an 8 MiB cache.
    const cap = 8 * 2 ** 20;
    ctx.config.set(SeqViewerConfig.Stream.CacheBytes, cap);
    let most = 0;
    const tour = performance.now();
    for (let k = 0; k < 96; k++) {
      const b = k % GRCH38.length;
      const at = Math.floor(((k * 7919) % 1000) * (GRCH38[b][1] / 1000));
      vp.zoomTo(b, at, at + 10_000);
      await settle(ctx, manual);
      most = Math.max(most, manager.stats.current.bytes);
    }
    const st = manager.stats.current;
    const items = s.genes.requests.length + s.variants.requests.length;
    rows.push(
      `tour of 96 windows: ${((performance.now() - tour) / 96).toFixed(1)} ms each, ${st.requests} requests ` +
        `(${s.genome.requests.length} of letters, ${items} of genes and variants), ${st.cancelled} cancelled, ` +
        `at most ${(most / 2 ** 20).toFixed(2)} MiB held (bound 8 MiB)`,
    );
    expect(most).toBeLessThanOrEqual(cap);
    expect(s.phylop.requests.length + s.genome.requests.length + items).toBe(st.requests);
    out.push(`DPR ${dpr}: ${rows.join("; ")}`);
    surface.dispose();
    ctx.dispose();
    host.remove();
  }, 120_000);
