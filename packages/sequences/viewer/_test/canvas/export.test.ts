/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The genome story rendered headless in the page to SVG (case 14 imports Node's fs through its fixtures, so it is in
// _test/svg/export.test.ts), and the SVG rasterized to a PNG by the page itself (SV1b), with the browser's own
// `CompressionStream` for images. Bounded time and size; the PNG is the SVG's size and not blank: dark ink on white.
// An export over a mounted, zoomed and scrolled section leaves its canvases as they were, pixel for pixel.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import { genomeDoc, synthGenome } from "../../../core/_test/format/genome-fixture.js";
import { FONT_FAMILIES, FONT_WEIGHTS, probeFonts } from "../../scripts/font-probe.js";
import { canvasMeasurer } from "../../src/canvas/index.js";
import { browserExportAction, canvasEmBox, compressionDeflate } from "../../src/dom/index.js";
import { createSeqViewer, DefaultRegistry, manualScheduler } from "../../src/index.js";
import { hasMetrics, type MetricsMeasurer, renderSvg } from "../../src/svg/index.js";
import { canvasesOf, close, diff, opened, pixelsOf, settle } from "./open.js";

/** The SVG drawn on white at its own size and written as PNG; the pixels too. */
async function rasterize(svg: string) {
  const size = /<svg[^>]* width="([\d.]+)" height="([\d.]+)"/.exec(svg);
  const [w, h] = [Math.ceil(Number(size?.[1])), Math.ceil(Number(size?.[2]))];
  const url = URL.createObjectURL(new Blob([svg], { type: "image/svg+xml" }));
  const img = new Image();
  img.src = url;
  await img.decode();
  const c = new OffscreenCanvas(w, h);
  const g = c.getContext("2d", { willReadFrequently: true }) as OffscreenCanvasRenderingContext2D;
  g.fillStyle = "#fff";
  g.fillRect(0, 0, w, h);
  g.drawImage(img, 0, 0, w, h);
  URL.revokeObjectURL(url);
  const png = new Uint8Array(await (await c.convertToBlob({ type: "image/png" })).arrayBuffer());
  return { w, h, png, px: g.getImageData(0, 0, w, h).data };
}
const darkPixels = (px: Uint8ClampedArray) => {
  let n = 0;
  for (let i = 0; i < px.length; i += 4) if (px[i] + px[i + 1] + px[i + 2] < 300) n++;
  return n;
};

describe("a headless render as PNG", () => {
  it("the genome story: 3,000,000 nucleotides across 24 chromosomes", async () => {
    const doc = genomeDoc(synthGenome());
    const t = performance.now();
    const svg = await renderSvg(doc, { width: 1200, background: true, ports: { deflate: compressionDeflate } });
    const ms = performance.now() - t;
    const r = await rasterize(svg);
    expect(ms).toBeLessThan(30_000); // about 1.5 s
    expect(svg.length).toBeLessThan(2_000_000);
    expect(r.w).toBe(1200);
    expect(r.png.length).toBeLessThan(2_000_000);
    expect(darkPixels(r.px)).toBeGreaterThan(2000);
  }, 120_000);
});

describe("the browser's export", () => {
  it("canvasEmBox is the em box a fresh probe reads, to a hundredth of an em, for every family and weight", () => {
    const box = canvasEmBox();
    for (const p of probeFonts({ families: FONT_FAMILIES, weights: FONT_WEIGHTS })) {
      const b = box(`${p.weight} 13px ${FONT_FAMILIES[p.family]}`);
      expect(b?.ascent, `${p.family} ${p.weight}`).toBeCloseTo(p.ascent / 1000, 1);
      expect(b?.descent, `${p.family} ${p.weight}`).toBeCloseTo(p.descent / 1000, 1);
      expect(Math.abs((b?.ascent ?? 0) - p.ascent / 1000)).toBeLessThanOrEqual(0.011);
      expect(Math.abs((b?.descent ?? 0) - p.descent / 1000)).toBeLessThanOrEqual(0.011);
      // The ink heights at 100 px are the probe's at 1,000 px (to its hundredth); at 13 px, hinted, a px of them.
      const big = box(`${p.weight} 100px ${FONT_FAMILIES[p.family]}`);
      expect(Math.abs((big?.cap ?? 0) - p.cap / 1000), `${p.family} ${p.weight}`).toBeLessThanOrEqual(0.011);
      expect(Math.abs((big?.xHeight ?? 0) - p.xHeight / 1000), `${p.family} ${p.weight}`).toBeLessThanOrEqual(0.011);
      expect(Math.abs((b?.cap ?? 0) * 13 - (p.cap / 1000) * 13)).toBeLessThanOrEqual(1);
    }
    // The em box is the family's and weight's at any size; the ink heights are each size's own.
    const [a, b] = [box("12px sans-serif"), box("24px sans-serif")];
    expect([a?.ascent, a?.descent]).toEqual([b?.ascent, b?.descent]);
  });

  it("view.exportSvg on a live context: the file offered has dy offsets and no dominant-baseline, and its images are deflated", async () => {
    const b = createSeqViewSpec({ title: "gene" });
    const s = b.data.sequence({ ref: "s", alphabet: "dna", value: "GAATTCGCCACCATGGCTTCCAAGGAACTGACCCCACGTGACGGCC" });
    const v = b.view({ ref: "v", label: "V" }).block({ sequence: s });
    v.track({ label: "" }).ruler();
    v.track({ label: "DNA" }).letters();
    const offers: { svg: string; name: string }[] = [];
    const ctx = createSeqViewer(
      {
        registry: [
          DefaultRegistry,
          { actions: [browserExportAction({ offer: (svg, name) => void offers.push({ svg, name }) })] },
        ],
        behaviors: [],
      },
      // The context's measurer is the page's own (a mounted surface sets it); the export takes its em boxes.
      { scheduler: manualScheduler(), text: canvasMeasurer() },
    );
    const r = await ctx.setDocument(b.getDocument());
    expect(r.ok).toBe(true);
    expect(ctx.actions.run("view.exportSvg", { width: 600, scale: 2 })).toBe(true);
    while (offers.length === 0) await new Promise((res) => setTimeout(res, 10));
    expect(offers[0].name).toBe("v.svg");
    expect(offers[0].svg).toMatch(/<svg [^>]* width="1200"/);
    expect(offers[0].svg).toMatch(/<text [^>]*dy="/);
    expect(offers[0].svg).not.toContain("dominant-baseline");
    ctx.dispose();
  });
});

describe("the canvas measurer's em boxes", () => {
  /** A small document of letters, whose texts are drawn at `middle`. */
  const gene = () => {
    const b = createSeqViewSpec({ title: "gene" });
    const s = b.data.sequence({ ref: "s", alphabet: "dna", value: "GAATTCGCCACCATGGCTTCC" });
    const v = b.view({ ref: "v", label: "V" }).block({ sequence: s });
    v.track({ label: "" }).ruler();
    v.track({ label: "DNA" }).letters();
    return b.getDocument();
  };

  it("a plain renderSvg over a mounted section (the canvas measurer, no metrics given) writes dy, not dominant-baseline", async () => {
    const t = await opened(gene(), { width: 500, height: 120 });
    expect(hasMetrics(t.ctx.text.port)).toBe(true);
    const svg = await renderSvg(t.ctx, { width: 500 });
    expect(svg).toMatch(/<text [^>]*dy="/);
    expect(svg).not.toContain("dominant-baseline");
    close(t);
  });

  it("are read again after a font load: the context's text version moves (its listener), the probe runs again", async () => {
    const t = await opened(gene(), { width: 500, height: 120 });
    const m = t.ctx.text.port as MetricsMeasurer;
    const proto = OffscreenCanvasRenderingContext2D.prototype;
    const read = proto.getImageData;
    let reads = 0;
    proto.getImageData = function (this: OffscreenCanvasRenderingContext2D, ...a: [number, number, number, number]) {
      reads++;
      return read.apply(this, a);
    };
    try {
      const box = m.metrics("400 13px sans-serif");
      expect(box).not.toBeNull();
      const once = reads;
      expect(once).toBeGreaterThan(0);
      expect(m.metrics("400 13px sans-serif")).toBe(box);
      expect(reads).toBe(once);
      // Another size: its ink heights read (two probes), the em box not.
      const other = m.metrics("400 20px sans-serif");
      expect([other?.ascent, other?.descent]).toEqual([box?.ascent, box?.descent]);
      expect(reads).toBe(once + 2);
      document.fonts.dispatchEvent(new Event("loadingdone"));
      expect(m.metrics("400 13px sans-serif")).toEqual(box);
      expect(reads).toBe(2 * once + 2);
    } finally {
      proto.getImageData = read;
      close(t);
    }
  });
});

describe("an export over a mounted section", () => {
  it("leaves its canvases as they were: no resize, the same pixels, zoomed in and scrolled to the far end", async () => {
    const b = createSeqViewSpec({ title: "gene" });
    const s = b.data.sequence({
      ref: "s",
      alphabet: "dna",
      value: "GAATTCGCCACCATGGCTTCCAAGGAACTGACCCCACGTGACGGCC".repeat(8),
    });
    const v = b.view({ ref: "v", label: "V" }).block({ sequence: s });
    v.track({ label: "" }).ruler();
    v.track({ label: "DNA" }).letters();
    const t = await opened(b.getDocument(), { width: 500, height: 120 });
    const vp = t.ctx.viewports.get(t.section);
    vp?.zoomAt(6, 200);
    vp?.scrollTo(1e9, 0);
    settle(t);
    const [main] = canvasesOf(t);
    const before = { size: [main.width, main.height], px: pixelsOf(main), resizes: t.surface.stats.resizes };
    // Wider: the scroll would be clamped to a nearer end; narrower: the zoom raised to the narrower fit.
    for (const width of [1400, 200]) await renderSvg(t.ctx, { width, scale: 2 });
    settle(t);
    expect([main.width, main.height]).toEqual(before.size);
    expect(t.surface.stats.resizes).toBe(before.resizes);
    expect(diff(pixelsOf(main), before.px, main.width, 1).count).toBe(0);
    close(t);
  });
});
