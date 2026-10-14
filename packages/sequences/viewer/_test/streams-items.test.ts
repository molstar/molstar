/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST3a-2: a streamed items layer drawn from the region cache (seqspec:sources#the-viewers-side, `StreamView.items`;
// seqspec:sources#levels-and-the-zoom-gate), over the streamed genome with genes and variants
// (core/_test/format/stream-provider.ts: synthetic genes of 2 to 6 exons a
// row each, one `id`; variants about one a kb). Regions of 65,536 positions of data below the first `count` level,
// bins above; a gene crossing regions keeps one lane and is drawn once; exons joined; density from bins; a pending
// region's placeholder band, gone when it arrives; an arrival logs its columns (an item past its region's included);
// the letters gate's hatched band and hint, and `Stream.LettersPerPx`.

import { describe, expect, it } from "vitest";
import { GENE_LEVELS, GRCH38, geneAt, hg38Streams, variantAt } from "../../core/_test/format/stream-provider.ts";
import { SkinToken } from "../src/draw/backend.ts";
import { RegionState, SeqViewerConfig } from "../src/index.ts";
import type { Streams } from "../src/streams.ts";
import { composed, recording, textsOf } from "./composer-fixtures.ts";
import { GENOME, MaxInFlight, streamedDoc } from "./stream-fixtures.ts";

const flush = () => new Promise((ok) => setTimeout(ok, 0));
const THREE = GRCH38.slice(0, 3);
const CHR2 = 1;
const SPAN = 65_536;

/** The streamed genome with genes and variants (and letters), composed at 1,000 × 300, its providers gated. */
async function genome(o: { letters?: boolean; variantsAs?: "variants" | "points" } = {}) {
  const s = hg38Streams({ chromosomes: THREE, letters: true, gated: true });
  const loaders = { load: () => Promise.reject(new Error("none")), stream: s.stream };
  const doc = streamedDoc({ items: true, sequence: o.letters, variantsAs: o.variantsAs });
  const t = await composed(doc, { width: 1000, height: 300, loaders });
  t.ctx.config.set(MaxInFlight, 256);
  const m = t.ctx.streams as Streams;
  const ref = (type: string) => {
    const v = t.ctx.graph.assembled.value?.view;
    for (let k = 0; v && k < v.resolved.layers.count; k++)
      if (v.set.strings.get(v.resolved.layers.type[k]) === type) return v.set.strings.get(v.resolved.layers.ref[k]);
    throw new Error(`no ${type} layer`);
  };
  const show = (start: number, stop: number) => {
    if (stop > start) t.vp?.zoomTo(CHR2, start, stop);
    else t.vp?.fit();
    t.layout(GENOME);
    m.pass();
  };
  const drain = async () => {
    for (let k = 0; k < 6; k++) {
      for (const p of [s.genome, s.phylop, s.genes, s.variants]) p.release();
      await flush();
      m.pass();
    }
  };
  return { t, s, m, ref, show, drain, frame: () => t.frame(recording()) };
}

describe("requests and the read", () => {
  it("bins above the first level, data regions of 65,536 below; each region in order, an item's lane kept", async () => {
    const g = await genome();
    g.show(100_000_000, 105_000_000);
    expect(g.s.genes.requests.every((r) => r.level >= 0)).toBe(true);
    g.show(100_000_000, 100_300_000);
    const data = g.s.genes.requests.filter((r) => r.level === -1);
    expect(data.length).toBeGreaterThan(4);
    expect(data.every((r) => r.start % SPAN === 0 && r.stop - r.start === SPAN)).toBe(true);
    await g.drain();
    const st = g.m.get(g.ref("ranges"));
    const out = new Int32Array(6);
    out[0] = 300;
    // Every region over [100.0, 100.3) Mb, left to right; the first owns what starts before it.
    const lanes = new Map<string, [number, number]>();
    let k = 0;
    let crossing = 0;
    let joins = 0;
    for (; ; k++) {
      const r = st?.items?.(CHR2, 100_000_000, 100_300_000, k, out);
      if (!r) break;
      expect(out[4]).toBe(RegionState.loaded);
      expect(out[3]).toBe(k === 0 ? -1 : out[1]);
      const f = r.features;
      if (!f) throw new Error("no features");
      for (let q = 0; q < f.count; q++) {
        // A gene seen in an earlier region (it crosses their edge) has the lane it had there.
        const key = `${r.first[q]}`;
        const was = lanes.get(key);
        if (was !== undefined && was[0] < k) {
          crossing++;
          expect(f.lane[q]).toBe(was[1]);
        }
        if (was === undefined) lanes.set(key, [k, f.lane[q]]);
        // A gap between exons is joined, from an exon of the same gene ending before it.
        const j = f.joinFrom[q];
        if (j >= 0) {
          joins++;
          expect(r.first[j]).toBe(r.first[q]);
          expect(f.stop[j]).toBeLessThan(f.start[q]);
        }
      }
    }
    expect(k).toBe(Math.floor(100_299_999 / SPAN) - Math.floor(100_000_000 / SPAN) + 1);
    expect(crossing).toBeGreaterThan(0);
    expect(joins).toBeGreaterThan(0);
  });

  it("an arrival logs its columns, what its items reach included", async () => {
    const g = await genome();
    g.show(100_000_000, 100_300_000);
    await g.drain();
    const st = g.m.get(g.ref("ranges"));
    const before = st?.version ?? 0;
    g.show(150_000_000, 150_100_000);
    await g.drain();
    const out = new Int32Array(64);
    const n = st?.changed?.(before, out) ?? -1;
    expect(n).toBeGreaterThan(0);
    // Each at least its region; some reach past it (a gene crossing its edge).
    let past = 0;
    for (let q = 0; q < n; q++) {
      expect(out[2 * q + 1] - out[2 * q]).toBeGreaterThanOrEqual(SPAN);
      if (out[2 * q] % SPAN !== 0 || out[2 * q + 1] % SPAN !== 0) past++;
    }
    expect(past).toBeGreaterThan(0);
  });
});

describe("drawn", () => {
  const PENDING = SkinToken.pending;
  const pendingRects = (r: ReturnType<typeof recording>) => {
    let n = 0;
    for (let k = 0; k < r.rectItems.n; k++) if (r.rectItems.color[k] === PENDING) n++;
    return n;
  };

  it("a gene is labeled once, however many regions hold it; placeholders until its region arrives", async () => {
    const g = await genome();
    g.show(100_050_000, 100_150_000);
    const waiting = g.frame();
    expect(pendingRects(waiting)).toBeGreaterThan(0);
    await g.drain();
    const done = g.frame();
    expect(pendingRects(done)).toBe(0);
    const labels = textsOf(done).filter((x) => /^G\d+$/.test(x));
    expect(labels.length).toBeGreaterThan(0);
    expect(new Set(labels).size).toBe(labels.length);
    // Every gene over the view is labeled where its label fits: the one crossing 100,073,728 (a region edge) too.
    const edge = Math.ceil(100_050_000 / SPAN) * SPAN;
    const crossing = Array.from({ length: 4000 }, (_, q) => geneAt(CHR2, q + 2700)).find(
      (x) => x.start < edge && x.stop > edge && x.start > 100_050_000,
    );
    if (crossing) expect(labels).toContain(crossing.label);
  });

  it("variants as glyphs, one a variant in view, a second at its position stacked above the first", async () => {
    const g = await genome();
    g.show(100_000_000, 100_008_000);
    await g.drain();
    const f = g.frame();
    let n = 0;
    for (let x = 100_000_000; x < 100_008_000; x++)
      n += +(variantAt(CHR2, x) !== "") + +(variantAt(CHR2, x, true) !== "");
    expect(n).toBeGreaterThan(2);
    expect(f.glyphItems.n).toBe(n);
  });

  it("points over streamed items: a glyph per item in view, by region, as variants draw them unstacked", async () => {
    const g = await genome({ variantsAs: "points" });
    g.show(100_000_000, 100_008_000);
    await g.drain();
    expect(g.s.variants.requests.length).toBeGreaterThan(0);
    const f = g.frame();
    let n = 0;
    for (let x = 100_000_000; x < 100_008_000; x++)
      n += +(variantAt(CHR2, x) !== "") + +(variantAt(CHR2, x, true) !== "");
    expect(n).toBeGreaterThan(2);
    expect(f.glyphItems.n).toBe(n);
  });

  it("density from bins when zoomed out, no labels", async () => {
    const g = await genome();
    g.show(0, 0);
    await g.drain();
    const f = g.frame();
    expect(textsOf(f).filter((x) => /^G\d+$/.test(x))).toEqual([]);
    expect(f.rectItems.n).toBeGreaterThan(200);
    // A level whose bins are at most a device px's positions: never the data, nothing coarser first.
    expect(g.s.genes.requests.every((r) => r.level >= 1 && r.level < GENE_LEVELS.length)).toBe(true);
  });

  it("past the letters gate a hatched band and its hint; inside it none; `Stream.LettersPerPx` moves it", async () => {
    const g = await genome({ letters: true });
    g.show(0, 0);
    expect(textsOf(g.frame())).toContain("zoom in to see letters");
    g.show(100_000_000, 100_010_000);
    expect(textsOf(g.frame())).not.toContain("zoom in to see letters");
    expect(g.s.genome.requests.length).toBeGreaterThan(0);
    const asked = g.s.genome.requests.length;
    g.t.ctx.config.set(SeqViewerConfig.Stream.LettersPerPx, 2);
    g.show(100_100_000, 100_110_000);
    expect(textsOf(g.frame())).toContain("zoom in to see letters");
    expect(g.s.genome.requests.length).toBe(asked);
  });
});
