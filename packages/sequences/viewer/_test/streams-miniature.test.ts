/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A miniature over streamed 1-D layers (GM-04; seqspec:sources#the-viewers-side): the navigator's map draws every
// block's whole shown range, so the manager asks for it too, at the level whose bin is about a map pixel (at most
// `MAX_REGIONS` regions a block, else coarser), after the view's regions and before the margin. Once settled, the map
// image holds no placeholder band; without a miniature nothing more is asked.

import type { Loaders } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import { hg38Streams, streamedDoc } from "../../core/_test/format/stream-provider.ts";
import { SkinToken } from "../src/draw/backend.ts";
import { Dirty, SectionComposer } from "../src/draw/composer.ts";
import {
  createSeqViewer,
  DefaultPresets,
  DefaultRegistry,
  DefaultSeqViewerSpec,
  manualScheduler,
} from "../src/index.ts";
import type { Streams } from "../src/streams.ts";
import { fixtureMeasurer, RecordingBackend } from "../src/testing/index.ts";

const flush = () => new Promise((ok) => setTimeout(ok, 0));
const PENDING = SkinToken.pending;

/** The pending bands in each image the frame built (the map's). */
function pendingIn(be: RecordingBackend): number[] {
  const out: number[] = [];
  for (const c of be.images.content)
    if (c instanceof RecordingBackend) {
      let n = 0;
      for (let i = 0; i < c.rectItems.n; i++) if (c.rectItems.color[i] === PENDING) n++;
      out.push(n);
    }
  return out;
}

async function run(nav: string) {
  const s = hg38Streams({
    letters: true,
    chromosomes: [
      ["chr1", 248_956_422],
      ["chr2", 242_193_529],
    ],
  });
  const loaders: Loaders = { load: (r) => Promise.reject(new Error(`no ${r.uri}`)), stream: s.stream };
  const presets = DefaultPresets.map((p) =>
    p.geometry === "linear" && p.name !== "strip" ? { ...p, navigator: nav } : p,
  );
  const scheduler = manualScheduler();
  const ctx = createSeqViewer(
    { ...DefaultSeqViewerSpec(), registry: [{ ...DefaultRegistry, presets }] },
    { scheduler, loaders, text: fixtureMeasurer() },
  );
  const size = () => {
    for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(1000, 300, 1);
  };
  const done = ctx.setDocument(streamedDoc({ items: true }));
  size();
  await done;
  await ctx.graph.update();
  const sec = ctx.viewports.sections.value[0];
  ctx.viewports.get(sec)?.zoomTo(1, 100_000_000, 100_200_000);
  const comp = new SectionComposer(ctx, sec);
  const be = new RecordingBackend();
  for (let k = 0; k < 30; k++) {
    size();
    scheduler.advance(16);
    await flush();
    await ctx.graph.update();
    be.reset();
    comp.frame(be, Dirty.all);
  }
  await (ctx.streams as Streams).settled();
  const out = { pending: pendingIn(be), stats: { ...ctx.streams.stats.current } };
  ctx.dispose();
  return out;
}

describe("a miniature over streamed 1-D layers (GM-04)", () => {
  it("the map's regions are asked and land: no placeholder band in the map once settled", {
    timeout: 60_000,
  }, async () => {
    const strip = await run("strip");
    expect(strip.pending.length).toBeGreaterThan(0);
    expect(strip.pending.every((n) => n === 0)).toBe(true);
    expect(strip.stats.inFlight).toBe(0);
    // Without a miniature, fewer regions: the map's are asked only for it.
    const rail = await run("rulerRail");
    expect(rail.pending).toEqual([]);
    expect(strip.stats.requests).toBeGreaterThan(rail.stats.requests);
  });
});
