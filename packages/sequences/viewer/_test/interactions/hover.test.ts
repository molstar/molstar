/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The fingerprint's hover emphasis (HV1; seqspec:interaction#hover) on the recording backend: the hovered stripe
// outlined inside its own pixels in the hover ink of its kind's color, in both skins, and nothing else; a cell's ground
// (no kind under the point) outlined inside the ground; nothing over a cell without contacts.

import { describe, expect, it } from "vitest";
import { SkinToken } from "../../src/draw/backend.js";
import { Dirty, hoverInk, Pass, SectionComposer } from "../../src/draw/index.js";
import { DarkSkin, DefaultSkin, Hit } from "../../src/index.js";
import { FingerprintRenderer } from "../../src/interactions/viewer/index.js";
import { LayerFrame, RecordingBackend } from "../../src/testing/index.js";
import { BoxKind, passBoxes, within } from "../hover-fixtures.js";
import { layerOf, open, unitScene } from "./fixtures.js";

describe.each(["light", "dark"] as const)("the fingerprint's hover emphasis, %s", (skin) => {
  it("the hovered stripe alone, inside its own rect, in its kind's hover ink; a ground inside the ground", async () => {
    const scene = unitScene();
    const t = await open(scene, { dpr: 2, size: [1000, 600], render: { maxCell: 400 } });
    t.ctx.behaviors.get("hoverLabel")?.remove();
    const sk = skin === "dark" ? DarkSkin : DefaultSkin;
    if (skin === "dark") t.ctx.setSkin(DarkSkin);
    t.ctx.viewports.get(scene.section)?.setView(24, 0, "api");
    await t.ctx.graph.update();
    const { k } = layerOf(t.ctx);
    const main = new RecordingBackend();
    new LayerFrame(t.ctx, FingerprintRenderer as never, scene.section, k).draw(main);
    const b = main.rectItems;
    const rects = Array.from({ length: b.n }, (_, i) => [b.x[i], b.y[i], b.x[i] + b.w[i], b.y[i] + b.h[i], b.color[i]]);
    const composer = new SectionComposer(t.ctx, scene.section);
    const hit = Hit.create();
    const outlines = (r: number[]) => {
      t.ctx.layout.hitTest(scene.section, (r[0] + r[2]) / 4, (r[1] + r[3]) / 4, hit);
      t.ctx.hover.setFromHit(hit, "pointer");
      const rec = new RecordingBackend();
      composer.frame(rec, Dirty.overlay);
      return passBoxes(rec, Pass.overlay).filter((x) => x[5] === BoxKind.outline || x[4] === hoverInk(sk, r[4]));
    };
    const stripes = rects.filter((r) => r[4] !== SkinToken.grid && r[3] - r[1] >= 6);
    expect(stripes.length).toBeGreaterThan(3);
    for (const r of stripes.slice(0, 6)) {
      const o = outlines(r);
      expect(o.length, `stripe ${r.join(",")}`).toBe(1);
      expect(within(o[0], r)).toBe(true);
      expect(o[0][4]).toBe(hoverInk(sk, r[4]));
    }
    // A ground: the cell where no stripe is under the point.
    const ground = rects.find(
      (g) => g[4] === SkinToken.grid && !stripes.some((s) => within([g[0], g[1], g[0] + 1, g[1] + 1], s)),
    );
    if (ground) {
      const o = outlines([ground[0], ground[1], ground[0] + 2, ground[1] + 2, ground[4]]);
      expect(o.length).toBe(1);
      expect(within(o[0], ground)).toBe(true);
    }
  });
});
