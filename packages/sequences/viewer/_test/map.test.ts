/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// RN3's `seqv:map` and the navigator's miniature (seqspec:viewer#extension-points): the tracks in miniature as an image
// kept under a key (a pan or a zoom builds none, a width change builds another), the view's window washed and outlined
// over it; the navigator's map image through the same host, none on the ruler.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import { Miniature } from "../src/draw/miniature.ts";
import { drawNavigator, NavView } from "../src/draw/nav.ts";
import { DefaultNavigators, KernelRenderers, type LayerRenderer } from "../src/index.ts";
import { LayerFrame, Op, RecordingBackend } from "../src/testing/index.ts";
import { miniatureDoc } from "./conformance-fixtures.ts";
import { layerOf } from "./draw-fixtures.ts";
import { opened } from "./layout-fixtures.ts";

const open = (d = miniatureDoc()) => opened(d, { renderers: [...KernelRenderers] as LayerRenderer[] });
const layers = (r: RecordingBackend) =>
  r.images.content.filter((c): c is RecordingBackend => c instanceof RecordingBackend);

describe("the seqv:map miniature", () => {
  it("is built once across pans and zooms, again on a width change; its window is outlined", async () => {
    const t = await open();
    const { k, section } = layerOf(t.ctx, "Map");
    const f = new LayerFrame(t.ctx, t.ctx.registry.get("renderers", "seqv:map") as LayerRenderer, section, k);
    const vp = t.ctx.viewports.get(section);
    f.draw(f.direct);
    expect(layers(f.direct)).toHaveLength(1);
    // Hit the miniature's tracks: the ranges' and the values' rects are in the layer, in lanes.
    const layer = layers(f.direct)[0];
    expect(layer.rectItems.n).toBeGreaterThan(3);
    expect(Array.from(f.direct.rectItems.stroke.subarray(0, f.direct.rectItems.n)).some((s) => s > 0)).toBe(true);
    expect(f.direct.count(Op.image)).toBe(1);
    for (const cell of [8, 4, 12]) {
      vp?.setView(cell, 0, "api");
      vp?.panBy(40, 0, "api");
      await t.ctx.graph.update();
      f.draw(f.direct);
    }
    expect(layers(f.direct)).toHaveLength(1);
    vp?.setSize(700, 400, 1);
    await t.ctx.graph.update();
    f.draw(f.direct);
    expect(f.direct.images.content.filter((c) => c !== null)).toHaveLength(1);
    expect(f.direct.images.w.length).toBe(2);
  });

  it("draws the tracks it names, else every track whose renderer is in miniatures", async () => {
    const all = await open();
    const only = await open(miniatureDoc(["t.hot"]));
    const sizes = [all, only].map((t) => {
      const { k, section } = layerOf(t.ctx, "Map");
      const f = new LayerFrame(t.ctx, t.ctx.registry.get("renderers", "seqv:map") as LayerRenderer, section, k);
      f.draw(f.direct);
      return layers(f.direct)[0].rectItems.n;
    });
    expect(sizes[0]).toBeGreaterThan(sizes[1]);
    expect(sizes[1]).toBeGreaterThan(0);
  });
});

describe("the navigator's map image", () => {
  it("is the miniature in the band, kept across calls; none on the ruler", async () => {
    const t = await open();
    const { section } = layerOf(t.ctx, "Map");
    const av = t.ctx.layout.axis(section);
    const rec = new RecordingBackend();
    const n = new NavView().from(av as never);
    n.backend = rec;
    n.text = t.ctx.text;
    n.parts = DefaultNavigators[2].parts;
    n.height = n.gripBottom = 20;
    const mini = new Miniature(t.ctx, section);
    n.map = mini.image(n);
    expect(n.map).toBeGreaterThanOrEqual(0);
    expect(mini.image(n)).toBe(n.map);
    expect(rec.images.content.filter((c) => c !== null)).toHaveLength(1);
    drawNavigator(n);
    expect([rec.count(Op.image), n.overMap]).toEqual([1, true]);
    n.parts = DefaultNavigators[0].parts;
    expect(mini.image(n)).toBe(-1);
    mini.dispose();
  });

  it("with no track to draw is none, kept under its key: no image built frame after frame", async () => {
    const b = createSeqViewSpec();
    b.view({ ref: "v", label: "V" })
      .block({ sequence: b.data.sequence({ ref: "a", value: "MKTAYIAKQR".repeat(20) }) })
      .track({ label: "L" })
      .letters();
    const t = await opened(b.getDocument(), { renderers: [...KernelRenderers] as LayerRenderer[] });
    const av = t.ctx.layout.axis("v.main");
    const rec = new RecordingBackend();
    const n = new NavView().from(av as never);
    n.backend = rec;
    n.parts = DefaultNavigators[2].parts;
    n.height = 20;
    const mini = new Miniature(t.ctx, "v.main");
    for (let i = 0; i < 3; i++) expect(mini.image(n)).toBe(-1);
    expect(rec.images.w.length).toBe(0);
  });
});
