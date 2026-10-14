/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hover emphasis in Chromium (HV1; seqspec:interaction#hover): a range hovered through a mounted surface at DPR 2,
// light and dark. The overlay canvas with the emphasis against the same frame without it: every pixel that changes lies
// inside the hovered item's own bars (David, 2026-10-12: never over its edge), some do, and the outline's pixels are
// the hover ink of the bar's color, exactly.

import { describe, expect, it } from "vitest";
import { SectionSurface } from "../../src/canvas/index.ts";
import { Pass } from "../../src/draw/backend.ts";
import { hoverInk } from "../../src/draw/index.ts";
import { createSeqViewer, DarkSkin, DefaultSeqViewerSpec, DefaultSkin, Hit, manualScheduler } from "../../src/index.ts";
import { fixtureMeasurer } from "../../src/testing/index.ts";
import { hoverDoc } from "../hover-fixtures.ts";

const DPR = 2;

describe.each(["light", "dark"] as const)("a hovered range's outline in Chromium, %s", (skin) => {
  it("changes only pixels inside the hovered item's bars, in the hover ink of their color", async () => {
    const ctx = createSeqViewer(DefaultSeqViewerSpec(), { scheduler: manualScheduler(), text: fixtureMeasurer() });
    const r = await ctx.setDocument(hoverDoc(100));
    if (!r.ok) throw new Error(JSON.stringify(r.issues));
    const sk = skin === "dark" ? DarkSkin : DefaultSkin;
    ctx.setSkin(sk);
    await ctx.graph.update();
    const host = document.createElement("div");
    host.style.cssText = "position:absolute;left:0;top:0;width:1000px;height:600px";
    document.body.appendChild(host);
    const section = ctx.viewports.sections.value[0];
    const surface = new SectionSurface(ctx, section, host, { dpr: DPR });
    try {
      ctx.behaviors.get("hoverLabel")?.remove();
      surface.renderNow();
      const ax = ctx.layout.axis(section);
      const l = ctx.layout.get(section)?.current;
      if (!ax || !l) throw new Error("no layout");
      // The ranges track's entry and its first lane (laneRows: 2 px padding, lanes of (22 − 4) / 2, 1 px inset).
      const rv = ctx.graph.assembled.value?.view?.resolved;
      const set = ctx.graph.assembled.value?.set;
      let e = 0;
      while (set?.strings.get(rv?.tracks.ref[l.entryTrack[e]] ?? -1) !== "t.ranges") e++;
      const top = l.scrolledTop + l.entryTop[e];
      const y0 = Math.round((top + 2 + 1) * DPR);
      const y1 = Math.round((top + 2 + 9 - 1) * DPR);
      const box = (a: number, b: number) => [
        Math.round(ax.x(0, a as never) * DPR),
        y0,
        Math.round(ax.x(0, b as never) * DPR),
        y1,
      ];
      const bars = [box(10, 20), box(30, 45)];
      // g1's second fragment, alone at 42.
      const hit = Hit.create();
      expect(ctx.layout.hitTest(section, ax.x(0, 42 as never) + 2, top + 6, hit)).toBe(true);
      ctx.hover.setFromHit(hit, "pointer");
      const canvases = host.querySelectorAll("canvas");
      const read = (pass: number) => {
        const c = canvases[pass] as HTMLCanvasElement;
        return (c.getContext("2d") as CanvasRenderingContext2D).getImageData(0, 0, c.width, c.height);
      };
      surface.renderNow();
      const on = read(Pass.overlay);
      ctx.setRenderProps({ hoverEmphasis: false });
      await ctx.graph.update();
      surface.renderNow();
      const off = read(Pass.overlay);
      expect([on.width, on.height]).toEqual([off.width, off.height]);
      let changed = 0;
      let outside = 0;
      for (let y = 0; y < on.height; y++)
        for (let x = 0; x < on.width; x++) {
          const i = 4 * (y * on.width + x);
          if (
            on.data[i] === off.data[i] &&
            on.data[i + 1] === off.data[i + 1] &&
            on.data[i + 2] === off.data[i + 2] &&
            on.data[i + 3] === off.data[i + 3]
          )
            continue;
          changed++;
          if (!bars.some((b) => x >= b[0] && x < b[2] && y >= b[1] && y < b[3])) outside++;
        }
      expect(changed).toBeGreaterThan(50);
      expect(outside).toBe(0);
      // The outline's left edge in the second bar: the hover ink of the bar's color (read on the main canvas).
      const main = read(Pass.main);
      const at = (img: ImageData, x: number, y: number) => {
        const i = 4 * (y * img.width + x);
        return [img.data[i], img.data[i + 1], img.data[i + 2], img.data[i + 3]];
      };
      const mid = (y0 + y1) >> 1;
      const fill = at(main, bars[1][0] + 8, mid);
      const ink = hoverInk(sk, (fill[0] << 16) | (fill[1] << 8) | fill[2]);
      expect(at(on, bars[1][0], mid)).toEqual([(ink >> 16) & 255, (ink >> 8) & 255, ink & 255, 255]);
    } finally {
      surface.dispose();
      ctx.dispose();
      host.remove();
    }
  });
});
