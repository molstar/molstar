/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST4b: the JBrowse 2 demo's navigation (David, 2026-10-12: after the whole chromosome, TP53 drew empty genes and
// ClinVar's "zoom in"): TP53, then chr17 whole, then TP53 again over the demo's sources (recorded slices,
// jbrowse-sources.ts): each view asks for its regions and draws them, and the second TP53 draws as the first.

import { afterAll, expect, it } from "vitest";
import { Dirty } from "../src/draw/composer.ts";
import { SeqViewerConfig } from "../src/index.ts";
import type { Streams } from "../src/streams.ts";
import { composed, recording } from "./composer-fixtures.ts";
import { demoSources, RECORD } from "./jbrowse-sources.ts";
import { JBROWSE_SECTION, jbrowseDoc } from "./jbrowse-story-fixtures.ts";

const src = demoSources("nav");
afterAll(src.save);
const TP53: [number, number] = [7_661_779, 7_687_538];

it("TP53, chr17 whole, TP53: each view asks for and draws its regions; the second TP53 draws as the first", async () => {
  const t = await composed(jbrowseDoc("chr17", ...TP53), { width: 800, height: 400, dpr: 1, loaders: src.loaders });
  const m = t.ctx.streams as Streams;
  // The recorded slices hold what the views ask for, not what a pan would prefetch.
  t.ctx.config.set(SeqViewerConfig.Stream.Prefetch, false);
  const settle = async () => {
    // Until every region in view is in (`settled`), frames running the passes as a page's would.
    t.layout(JBROWSE_SECTION);
    let done = false;
    void m.settled().then(() => {
      done = true;
    });
    for (let k = 0; k < 2000 && !done; k++) {
      t.scheduler.advance(16);
      await new Promise((ok) => setTimeout(ok, RECORD ? 20 : 0));
    }
    expect(done).toBe(true);
  };
  const draw = () => {
    const be = t.frame(recording(), Dirty.tracks);
    return {
      rects: be.rectItems.n,
      glyphs: be.glyphItems.n,
      texts: be.textItems.text.slice(0, be.textItems.n).sort(),
    };
  };
  const go = async (a: number, e: number) => {
    t.vp?.zoomTo(0, a, e);
    await settle();
    return draw();
  };
  const first = await go(...TP53);
  expect(first.texts).toContain("TP53");
  expect(first.rects).toBeGreaterThan(100);
  // The navigator's miniature asked for the whole chromosome at a map pixel's level already (GM-04): the whole view
  // needs nothing coarser.
  const asked = src.asked.length;
  expect(src.asked.some((q) => q.start <= 0 && q.stop >= 40_000_000)).toBe(true);
  const whole = await go(0, 83_257_441);
  expect(src.asked.length).toBeGreaterThanOrEqual(asked);
  expect(whole.texts).toContain("zoom in to see variants");
  expect(whole.texts).toContain("zoom in to see letters");
  const again = await go(...TP53);
  expect(again).toEqual(first);
  expect(t.ctx.issues.current.value.map((i) => i.message)).toEqual([]);
  expect(m.stats.current.errors).toBe(0);
  t.ctx.dispose();
}, 600_000);

const genome = demoSources("nav-genome", ["chr17", "chr21"]);
afterAll(genome.save);

it("over a genome of blocks (chr17, chr21): TP53, the whole genome, chr17 whole, TP53 again", async () => {
  const t = await composed(jbrowseDoc("chr17", ...TP53, { genome: true }), {
    width: 800,
    height: 400,
    dpr: 1,
    loaders: genome.loaders,
  });
  const m = t.ctx.streams as Streams;
  // The recorded slices hold what the views ask for, not what a pan would prefetch.
  t.ctx.config.set(SeqViewerConfig.Stream.Prefetch, false);
  const settle = async () => {
    // Until every region in view is in (`settled`), frames running the passes as a page's would.
    t.layout(JBROWSE_SECTION);
    let done = false;
    void m.settled().then(() => {
      done = true;
    });
    for (let k = 0; k < 2000 && !done; k++) {
      t.scheduler.advance(16);
      await new Promise((ok) => setTimeout(ok, RECORD ? 20 : 0));
    }
    expect(done).toBe(true);
  };
  const draw = () => {
    const be = t.frame(recording(), Dirty.tracks);
    return { rects: be.rectItems.n, glyphs: be.glyphItems.n, texts: be.textItems.text.slice(0, be.textItems.n).sort() };
  };
  const at = async (go: () => void) => {
    go();
    await settle();
    return draw();
  };
  const chr17 = () => {
    const ax = t.ctx.layout.axis(JBROWSE_SECTION)?.axis;
    const set = t.ctx.graph.assembled.value?.view?.set;
    for (let b = 0; ax && set && b < ax.count; b++)
      if (set.strings.get(set.units.key[ax.unit[b]]).endsWith("/chr17")) return b;
    return -1;
  };
  const b = chr17();
  expect(b).toBeGreaterThanOrEqual(0);
  const first = await at(() => t.vp?.zoomTo(b, ...TP53));
  expect(first.texts).toContain("TP53");
  await at(() => t.vp?.fit());
  await at(() => t.vp?.zoomTo(b, 0, 83_257_441));
  const again = await at(() => t.vp?.zoomTo(b, ...TP53));
  expect(again).toEqual(first);
  expect(t.ctx.issues.current.value.map((i) => i.message)).toEqual([]);
  t.ctx.dispose();
}, 600_000);
