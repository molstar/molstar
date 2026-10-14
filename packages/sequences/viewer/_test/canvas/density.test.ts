/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// RN2's density as pixels, on Canvas2D in Chromium: a device column holding items of two colors shows both, stacked in
// palette order (the first at the bottom) with heights by share; under `majority`, the plurality color alone at its
// share. The band is drawn by the shared helper through the test `points` entry.

import { describe, expect, it } from "vitest";
import { Canvas2DBackend, Canvas2DShared } from "../../src/canvas/index.ts";
import { type DrawResources, SkinToken } from "../../src/draw/backend.ts";
import {
  createSeqViewer,
  DefaultRegistry,
  DefaultSeqViewerSpec,
  DefaultSkin,
  type LayerRenderer,
  manualScheduler,
  type SeqViewerContext,
  type Skin,
} from "../../src/index.ts";
import { fixtureMeasurer, LayerFrame } from "../../src/testing/index.ts";
import { densityDoc, densityItems, withDensity } from "../rn2-fixtures.ts";

async function open(doc: ReturnType<typeof densityDoc>, dpr: number) {
  const spec = DefaultSeqViewerSpec();
  const ctx = createSeqViewer(
    { ...spec, registry: [{ ...DefaultRegistry, renderers: withDensity("points") }] },
    { scheduler: manualScheduler(), text: fixtureMeasurer() },
  );
  const done = ctx.setDocument(doc);
  for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(900, 300, dpr);
  const r = await done;
  if (!r.ok) throw new Error(JSON.stringify(r.issues));
  for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(900, 300, dpr);
  await ctx.graph.update();
  return ctx;
}
function resources(): DrawResources {
  const tokens = new Map<number, string>();
  for (const [name, id] of Object.entries(SkinToken))
    tokens.set(id, (DefaultSkin as Skin)[name as keyof Skin] as string);
  return {
    font: () => "10px monospace",
    fontPx: () => 10,
    token: (t) => tokens.get(t) ?? "#ff00ff",
    glyph: () => null,
  };
}
/** The track's frame, its palette (legend order) and its items' categories by data row. */
function frameOf(ctx: SeqViewerContext) {
  const m = ctx.graph.assembled.value?.view;
  if (!m) throw new Error("no view");
  const k = 0;
  const section = m.set.strings.get(m.resolved.sections.ref[m.resolved.tracks.section[m.resolved.layers.track[k]]]);
  const entry = ctx.registry.get("renderers", m.prepared[k].renderer) as LayerRenderer;
  const colors: number[] = [];
  for (const t of m.prepared[k].themes)
    t?.legend?.({ item: (_l, c) => colors.push(c), ramp: () => {}, size: () => {} });
  return { f: new LayerFrame(ctx, entry, section, k), colors, state: m.prepared[k].state as never as Features };
}
interface Features {
  count: number;
  block: Int32Array;
  start: Int32Array;
  row: Int32Array;
}

describe("density composed per device column, as pixels", () => {
  for (const density of ["stacked", "majority"] as const)
    it.each([1, 2])(`${density}, at DPR %i`, async (dpr) => {
      const ctx = await open(densityDoc({ density }), dpr);
      const { f, colors, state } = frameOf(ctx);
      const [w, h] = [900 * dpr, 300 * dpr];
      const c = new OffscreenCanvas(w, h);
      const g = c.getContext("2d", { willReadFrequently: true }) as OffscreenCanvasRenderingContext2D;
      f.draw(new Canvas2DBackend([g], new Canvas2DShared(resources()), 0));
      const px = g.getImageData(0, 0, w, h).data;
      const v = f.v;
      const y0 = Math.round(v.top * dpr);
      const y1 = Math.round((v.top + v.height) * dpr);
      // Per device column: the items by category, brute force (each at its center's column).
      const category = densityItems({ density }).category;
      const counts = new Map<number, [number, number]>();
      for (let i = 0; i < state.count; i++) {
        const q = v.shownOf(0, state.start[i]);
        const x = q >= 0 ? v.xd(0, q) : -1;
        if (x < Math.round(v.left * dpr) || x >= Math.round(v.right * dpr)) continue;
        const e = counts.get(x) ?? [0, 0];
        e[category[state.row[i]] === "k1" ? 1 : 0]++;
        counts.set(x, e);
      }
      const rgb = (y: number, x: number) => {
        const i = (y * w + x) * 4;
        return px[i + 3] === 0 ? -1 : (px[i] << 16) | (px[i + 1] << 8) | px[i + 2];
      };
      let both = 0;
      for (const [x, [n0, n1]] of counts) {
        // The column's runs of color from the bottom up.
        const runs: [number, number][] = [];
        for (let y = y1 - 1; y >= y0; y--) {
          const col = rgb(y, x);
          if (col < 0) break;
          if (runs.length > 0 && runs[runs.length - 1][0] === col) runs[runs.length - 1][1]++;
          else runs.push([col, 1]);
        }
        const bar = runs.reduce((s, r) => s + r[1], 0);
        expect(bar, `column ${x}`).toBeGreaterThan(0);
        if (density === "majority") {
          expect(runs.map((r) => r[0])).toEqual([n1 > n0 ? colors[1] : colors[0]]);
          continue;
        }
        if (n0 === 0 || n1 === 0) continue;
        both++;
        expect(runs.map((r) => r[0])).toEqual([colors[0], colors[1]]);
        expect(Math.abs(runs[0][1] - (bar * n0) / (n0 + n1))).toBeLessThanOrEqual(1);
      }
      if (density === "stacked") expect(both).toBeGreaterThan(50);
      ctx.dispose();
    });
});
