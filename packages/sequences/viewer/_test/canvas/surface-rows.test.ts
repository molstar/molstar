/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The rows draw now (RN3), and a vertical scroll copies the scrolled band (CS1a): no pixel of a scrolled row reaches
// the pinned band above it, after random vertical scrolls and row zooms, on the surface (Chromium;
// seqspec:viewer#the-render-backend). And a lazy alignment's row labels, arriving after its tiles, are drawn in the
// header (an arrival once redrew only the data box).

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { tilesStream } from "@molstar/seq-core/stream/tiles";
import { tileServer } from "@molstar/seq-core/stream/tiles/server";
import { afterEach, expect, it } from "vitest";
import { Dirty } from "../../src/draw/index.js";
import { syntheticFamily } from "../families.js";
import { canvasesOf, close, diff, type Opened, opened, pixelsOf } from "./open.js";

let t: Opened | null = null;
afterEach(() => {
  if (t) close(t);
  t = null;
});

/** A seeded generator, so a failure reproduces. */
function random(seed: number) {
  let s = seed >>> 0;
  return () => {
    s = (s * 1664525 + 1013904223) >>> 0;
    return s / 2 ** 32;
  };
}

for (const dpr of [1, 2])
  it(`no pixel of a scrolled row reaches the pinned band, after random vertical scrolls and row zooms (DPR ${dpr})`, async () => {
    t = await opened(syntheticFamily(120, 8), { dpr, width: 800, height: 260 });
    const o = t;
    const vp = o.ctx.viewports.get(o.section);
    const l = o.ctx.layout.get(o.section)?.current;
    if (!vp || !l) throw new Error("no viewport");
    const [main] = canvasesOf(o);
    const g = main.getContext("2d") as CanvasRenderingContext2D;
    const y0 = Math.round(l.navHeight * dpr);
    const band = () => g.getImageData(0, y0, main.width, Math.round(l.scrolledTop * dpr) - y0).data;
    const redraw = () => {
      o.surface.composer.redraw();
      o.surface.invalidate(Dirty.all);
      o.surface.renderNow();
    };
    redraw();
    const pinned = band();
    expect(pinned.some((v, i) => i % 4 === 3 && v > 0)).toBe(true);
    const rnd = random(9);
    for (let k = 0; k < 60; k++) {
      if (k % 20 === 19) vp.setRowScale(0.5 + rnd(), 0);
      vp.panBy(0, (rnd() - 0.45) * 90, "pointer");
      o.advance(16);
      o.surface.renderNow();
      const copied = band();
      // The band as copies left it equals the band redrawn whole, and the first one: no scrolled row leaked into it.
      redraw();
      expect(diff(copied, band(), main.width).count, `scroll ${k} against a redraw`).toBe(0);
      expect(diff(copied, pinned, main.width).count, `scroll ${k} against the first`).toBe(0);
    }
    expect(o.surface.composer.stats.copies).toBeGreaterThan(30);
  });

it("a lazy alignment's row labels that arrive after its tiles are drawn in the header, as a whole redraw draws them", async () => {
  const uri = "https://example.org/held";
  const server = tileServer({ rowCount: 400, columnCount: 300 });
  // The labels' chunks wait until released; the tiles come at once.
  let release = () => {};
  const held = new Promise<void>((r) => {
    release = r;
  });
  const b = createSeqViewSpec({ title: "held" });
  const aln = b.data.alignment({
    ref: "held",
    row_count: 400,
    column_count: 300,
    source: { uri, protocol: "seqstar-tiles/1" },
  } as never);
  b.view({ ref: "v", label: "V", style: "msa", height: "box" })
    .block({ alignment: aln })
    .track({ label: "Rows" })
    .rows();
  t = await opened(b.getDocument(), {
    width: 900,
    height: 300,
    loaders: {
      load: (r) => Promise.reject(new Error(`no fixture for ${r.uri ?? r.hostFile}`)),
      stream: tilesStream((u) =>
        u === uri
          ? async (path, signal) => {
              if (path.startsWith("labels/")) await held;
              return server.transport(path, signal);
            }
          : null,
      ),
    },
  });
  const o = t;
  const [main] = canvasesOf(o);
  const frames = async () => {
    for (let k = 0; k < 6; k++) {
      await new Promise((r) => setTimeout(r, 0));
      o.advance(16);
      o.surface.renderNow();
    }
  };
  await frames();
  const unlabelled = pixelsOf(main);
  release();
  await frames();
  const arrived = pixelsOf(main);
  o.surface.composer.redraw();
  o.surface.invalidate(Dirty.all);
  o.surface.renderNow();
  const whole = pixelsOf(main);
  expect(diff(unlabelled, whole, main.width, 1).count, "the labels show in a whole redraw").toBeGreaterThan(100);
  expect(diff(arrived, whole, main.width, 1).count, "the frames after the labels arrived").toBe(0);
});
