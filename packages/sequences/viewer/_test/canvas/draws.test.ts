/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// RN1's draws as pixels in Chromium (seqspec:viewer#the-render-backend): each RN1 layer drawn on the Canvas2D backend
// equals its recording replayed into another, at every level of detail and DPR 1 and 2; `color` layers ink the covered
// cells alone (bands) or every cell (values); letters' cells carry their theme's alpha; emphasis underlines only its
// own
// positions on each wrapped line and its bold letters cover the regular ones.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import { Canvas2DBackend, Canvas2DShared } from "../../src/canvas/index.js";
import { type DrawResources, SkinToken } from "../../src/draw/backend.js";
import {
  createSeqViewer,
  DefaultSeqViewerSpec,
  DefaultSkin,
  type LayerRenderer,
  manualScheduler,
  type SeqViewerContext,
  type Skin,
} from "../../src/index.js";
import { fixtureMeasurer, LayerFrame, RecordingBackend } from "../../src/testing/index.js";
import { logoDoc } from "../logo-fixtures.js";
import { allTypes, colorDoc } from "../prepare-fixtures.js";

/** A context showing `doc` (inline data only), its viewports `width` × `height` at `dpr`. */
async function open(doc: ReturnType<typeof allTypes>, width: number, height: number, dpr: number, render = {}) {
  const ctx = createSeqViewer(
    { ...DefaultSeqViewerSpec(), render },
    { scheduler: manualScheduler(), text: fixtureMeasurer() },
  );
  const done = ctx.setDocument(doc);
  for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(width, height, dpr);
  const r = await done;
  if (!r.ok) throw new Error(JSON.stringify(r.issues));
  for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(width, height, dpr);
  await ctx.graph.update();
  return ctx;
}
/** The context's skin tokens, fonts and glyphs, as a surface gives a backend. */
function resources(ctx: SeqViewerContext): DrawResources {
  const tokens = new Map<number, string>();
  for (const [name, id] of Object.entries(SkinToken))
    tokens.set(id, (DefaultSkin as Skin)[name as keyof Skin] as string);
  return {
    font: (id) => ctx.text.fontString(id),
    fontPx: (id) => Number.parseFloat(/([0-9.]+)px/.exec(ctx.text.fontString(id))?.[1] ?? "10"),
    token: (t) => tokens.get(t) ?? "#ff00ff",
    glyph: () => null,
  };
}
/** A layer of a track (by label) and its frame. */
function frameOf(ctx: SeqViewerContext, label: string, k0 = 0) {
  const m = ctx.graph.assembled.value?.view;
  if (!m) throw new Error("no view");
  const r = m.resolved;
  for (let k = 0; k < r.layers.count; k++) {
    const t = r.layers.track[k];
    if (m.set.strings.get(r.tracks.label[t]) !== label) continue;
    const section = m.set.strings.get(r.sections.ref[r.tracks.section[t]]);
    const entry = ctx.registry.get("renderers", m.prepared[k + k0].renderer) as LayerRenderer;
    return { f: new LayerFrame(ctx, entry, section, k + k0), section };
  }
  throw new Error(`no track "${label}"`);
}
/** An offscreen canvas and a Canvas2D backend drawing straight to it. */
function canvas(ctx: SeqViewerContext, w: number, h: number) {
  const c = new OffscreenCanvas(w, h);
  const g = c.getContext("2d", { willReadFrequently: true }) as OffscreenCanvasRenderingContext2D;
  const be = new Canvas2DBackend([g], new Canvas2DShared(resources(ctx)), 0);
  return { be, px: () => g.getImageData(0, 0, w, h).data };
}

/** Each label's layer at each cell size drawn on Canvas2D and replayed from its recording: the same pixels, some ink. */
async function replays(ctx: SeqViewerContext, labels: string[], cells: number[], w: number, h: number) {
  for (const cell of cells) {
    for (const label of labels) {
      const { f, section } = frameOf(ctx, label);
      ctx.viewports.get(section)?.setView(cell, 0, "api");
      await ctx.graph.update();
      const direct = canvas(ctx, w, h);
      expect(f.draw(direct.be), `${label} at ${cell}`).toBeGreaterThan(0);
      const rec = new RecordingBackend();
      f.draw(rec);
      const replayed = canvas(ctx, w, h);
      rec.replay(replayed.be);
      const [a, b] = [direct.px(), replayed.px()];
      let diff = 0;
      let ink = 0;
      for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) diff++;
      for (let i = 3; i < a.length; i += 4) if (a[i] > 40) ink++;
      expect(diff, `${label || "ruler"} at ${cell}`).toBe(0);
      expect(ink, `${label || "ruler"} at ${cell}`).toBeGreaterThan(0);
    }
  }
}

describe("a frame drawn on Canvas2D equals its recording replayed", () => {
  it.each([1, 2])("every RN1 layer at each level of detail, at DPR %i", async (dpr) => {
    const ctx = await open(allTypes({ length: 400 }), 900, 420, dpr);
    await replays(ctx, ["", "Letters", "Emphasis", "Values", "Color"], [14, 4, 1.5, 0.3], 900 * dpr, 420 * dpr);
    ctx.dispose();
  });

  it.each([1, 2])("every RN2 layer at each level of detail, at DPR %i", async (dpr) => {
    const ctx = await open(allTypes({ length: 400 }), 900, 420, dpr);
    await replays(
      ctx,
      ["Contacts", "Features", "Sites", "Variants", "Secondary"],
      [14, 4, 1.5, 0.3],
      900 * dpr,
      420 * dpr,
    );
    ctx.dispose();
  });

  it.each([1, 2])("the logo, its stacks as stretched characters and its thin columns, at DPR %i", async (dpr) => {
    const wide = await open(logoDoc(12, 3000), 900, 300, dpr);
    await replays(wide, ["Logo"], [0.24, 16], 900 * dpr, 300 * dpr);
    wide.dispose();
    const ctx = await open(logoDoc(), 900, 300, dpr);
    await replays(ctx, ["Logo"], [24], 900 * dpr, 300 * dpr);
    ctx.dispose();
  });
});

describe("color layers as pixels", () => {
  it.each([1, 2])("bands ink the covered cells alone, values every cell, at DPR %i", async (dpr) => {
    const ctx = await open(colorDoc(), 900, 200, dpr);
    const inked = (label: string, i: number) => {
      const { f, section } = frameOf(ctx, label);
      ctx.viewports.get(section)?.setView(10, 0, "api");
      const c = canvas(ctx, 900 * dpr, 200 * dpr);
      f.draw(c.be);
      const x = (f.v.xd(0, i) + f.v.xe(0, i)) >> 1;
      const y = Math.round((f.v.top + f.v.height / 2) * dpr);
      return c.px()[(y * 900 * dpr + x) * 4 + 3] > 0;
    };
    await ctx.graph.update();
    expect([15, 35, 25, 5].map((i) => inked("Bands", i))).toEqual([true, true, false, false]);
    expect([5, 25, 59].map((i) => inked("Shade", i))).toEqual([true, true, true]);
    ctx.dispose();
  });
});

describe("emphasis as pixels", () => {
  it.each([1, 2])(
    "wrapped: the accent underlines only the emphasized positions of each line, at DPR %i",
    async (dpr) => {
      const b = createSeqViewSpec();
      const a = b.data.sequence({ ref: "a", value: "MKVLATGHEK".repeat(8) });
      const d = b.data.annotation({ ref: "d", geometry: "ranges", target: a, value: { start: [12], stop: [70] } });
      const track = b.view({ ref: "v", label: "V", layout: "text" }).track({ label: "Letters" });
      track.letters();
      track.emphasis({ data: d });
      const ctx = await open(b.getDocument(), 300, 600, dpr, { wrappedGroupLabels: "none" });
      const [w, h] = [300 * dpr, 600 * dpr];
      const { f } = frameOf(ctx, "Letters", 1);
      const c = canvas(ctx, w, h);
      f.draw(c.be);
      const px = c.px();
      // The skin's accent, as the canvas parsed it.
      const probe = new OffscreenCanvas(1, 1).getContext("2d") as OffscreenCanvasRenderingContext2D;
      probe.fillStyle = DefaultSkin.accent;
      probe.fillRect(0, 0, 1, 1);
      const accent = [...probe.getImageData(0, 0, 1, 1).data.slice(0, 3)];
      const l = ctx.layout.get(f.section)?.current;
      if (!l) throw new Error("no layout");
      let lines = 0;
      for (let e = 0; e < l.entryCount; e++) {
        if (l.entryTrack[e] < 0) continue;
        const line = l.entryRow[e];
        const [k0, k1] = [Math.max(12, l.lines.start[line]), Math.min(70, l.lines.stop[line])];
        const top = Math.round((l.scrolledTop + l.entryTop[e]) * dpr);
        const bottom = Math.round((l.scrolledTop + l.entryTop[e] + l.entryHeight[e]) * dpr);
        // The accent's columns on this line's band.
        const xs = new Set<number>();
        for (let y = top; y < bottom; y++)
          for (let x = 0; x < w; x++) {
            const i = (y * w + x) * 4;
            if (px[i + 3] === 255 && px[i] === accent[0] && px[i + 1] === accent[1] && px[i + 2] === accent[2])
              xs.add(x);
          }
        if (k0 >= k1) {
          expect(xs.size, `line ${line}`).toBe(0);
          continue;
        }
        lines++;
        f.v.from(ctx.layout.axis(f.section) as never);
        f.v.lineStart = l.lines.start[line];
        f.v.lineStop = l.lines.stop[line];
        f.v.onlyBlock = 0;
        const [x0, x1] = [f.v.xd(0, k0), f.v.xe(0, k1 - 1)];
        expect(Math.min(...xs), `line ${line}`).toBeGreaterThanOrEqual(x0);
        expect(Math.max(...xs), `line ${line}`).toBeLessThan(x1);
        expect(xs.size, `line ${line}`).toBeGreaterThan(0);
      }
      expect(lines).toBeGreaterThan(1);
      ctx.dispose();
    },
  );
});
