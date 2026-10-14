/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// GenomeScale against its targets (CS2's verify step): the synthetic genome (3M nucleotides in 24 chromosome blocks,
// 300 proteins, every track) on a surface through the timed renderers of genome-probe.ts, at 1,320 × 640 CSS px, in
// Chromium. A whole view under 16 ms (a frame, from the full redraw to the pixels), a hover step under 2 ms (the
// pointer move, hit test, label and overlay frame), the numbers printed by `GENOME=1` for the progress file; the heap
// of the loaded genome stays under 400 MB in Node.

import { afterAll, expect, it } from "vitest";
import { commands } from "vitest/browser";
import { mountSection } from "../../src/canvas/index.ts";
import { Inputs, manualScheduler, PointerType } from "../../src/index.ts";
import { genomeContext, genomeDocument, PROBE, SECTION } from "./genome-probe.ts";

const median = (xs: number[]) => xs.slice().sort((a, b) => a - b)[xs.length >> 1];
const out: string[] = [];
// The numbers go to __debug__/genome.txt (gitignored) for the progress file.
afterAll(() => commands.writeFile("__debug__/genome.txt", out.join("\n")));

for (const dpr of [1, 2])
  it(`the genome's whole view draws under 16 ms and a hover step under 2 ms (DPR ${dpr})`, async () => {
    const manual = manualScheduler(0);
    const { ctx, ready } = genomeContext(manual, genomeDocument());
    await ready;
    const host = document.createElement("div");
    host.style.cssText = "position:absolute;left:0;top:0;width:1320px;height:640px";
    document.body.appendChild(host);
    const surface = mountSection(ctx, SECTION, host, { dpr });
    await ctx.graph.update();
    manual.advance(16);
    surface.renderNow();
    const vp = ctx.viewports.get(SECTION);
    const full = () => {
      surface.composer.redraw();
      surface.invalidate(7);
      const t = performance.now();
      surface.renderNow();
      // The deferred raster is flushed by a one-pixel read, as a visible page's compositor would.
      (host.querySelector("canvas") as HTMLCanvasElement).getContext("2d")?.getImageData(0, 0, 1, 1);
      return performance.now() - t;
    };
    const rows: string[] = [];
    for (const [label, cell] of [
      ["whole", 0],
      ["1 Mb", 0.001],
      ["10 kb", 0.1],
    ] as const) {
      if (cell > 0) vp?.setView(cell, 0, "api");
      else vp?.fit();
      manual.advance(16);
      for (let i = 0; i < 5; i++) full();
      const times = Array.from({ length: 20 }, full);
      PROBE.reset();
      full();
      rows.push(`${label}: full ${median(times).toFixed(1)} ms, max ${Math.max(...times).toFixed(1)}`);
      expect(median(times), `${label} at DPR ${dpr}`).toBeLessThan(16);
    }
    // Hover: a pointer moving over the tracks, each move its overlay frame.
    vp?.fit();
    manual.advance(16);
    surface.renderNow();
    const l = ctx.layout.get(SECTION)?.current;
    if (!l) throw new Error("no layout");
    const p = Inputs.pointer();
    p.section = ctx.layout.sectionIndex(SECTION);
    p.type = PointerType.move;
    const hovers: number[] = [];
    for (let i = 0; i < 60; i++) {
      p.x = 200 + (i % 40) * 17;
      p.y = l.scrolledTop + 4 + (i % 9) * 25;
      const t = performance.now();
      ctx.input.pointer(p);
      manual.advance(16);
      surface.renderNow();
      hovers.push(performance.now() - t);
    }
    rows.push(`hover: ${median(hovers).toFixed(2)} ms, max ${Math.max(...hovers).toFixed(2)}`);
    expect(median(hovers), `hover at DPR ${dpr}`).toBeLessThan(2);
    out.push(`DPR ${dpr}: ${rows.join("; ")}`);
    surface.dispose();
    ctx.dispose();
    host.remove();
  }, 120_000);
