/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ideogram navigator as pixels (ST5; seqspec:viewer#extension-points), in Chromium: the navigator pass's canvas
// holds the strip the draw describes, in the skin's colors (made distinct to read them exactly): a band's ink at its
// stain's alpha over the ground, a negative band the ground, the centromere's triangles in red meeting at the middle
// (the corners beyond them the ground), a chromosome's outline, and a hairline across a block without bands. At 1 and 2
// device px per CSS px.

import { afterEach, describe, expect, it } from "vitest";
import type { SectionAxisView } from "../../src/axis.ts";
import { type CytoBands, ideogramNavigator } from "../../src/index.ts";
import { blocksDoc } from "../axis-fixtures.ts";
import { canvasesOf, close, type Opened, opened, pixelsOf, settle } from "./open.ts";

let t: Opened | null = null;
afterEach(() => {
  if (t) close(t);
  t = null;
});

const SKIN =
  "--seq-surface:rgb(0 255 0);--seq-surface-alt:rgb(255 255 0);--seq-ink:rgb(0 0 0);--seq-ink-quiet:rgb(50 50 50);" +
  "--seq-ink-faint:rgb(200 200 200);--seq-grip:rgb(255 0 0);--seq-grip-hover:rgb(0 0 255)";
const BANDS: CytoBands = {
  c0: [
    [0, 100, "p2", "gneg"],
    [100, 180, "p1", "gpos50"],
    [180, 220, "p11", "acen"],
    [220, 260, "q11", "acen"],
    [260, 400, "q1", "gpos100"],
  ],
  c1: [[0, 400, "p1", "gneg"]],
};

async function mount(dpr: number) {
  t = await opened(blocksDoc(3, 400), {
    width: 840,
    height: 200,
    dpr,
    style: SKIN,
    registry: [
      {
        navigators: [ideogramNavigator(BANDS)],
        presets: [{ name: "karyo", geometry: "linear", render: {}, navigator: "ideogram" }],
      },
    ],
  });
  t.ctx.ui.setStyle(t.section, "karyo");
  await t.ctx.graph.update();
  settle(t);
  return t;
}

describe("the ideogram navigator's pixels", () => {
  for (const dpr of [1, 2])
    it(`bands in ink at their stain's alpha, red triangles, an outline, a hairline (${dpr} device px per CSS px)`, async () => {
      const o = await mount(dpr);
      const nav = canvasesOf(o)[2];
      const px = pixelsOf(nav);
      const v = o.ctx.layout.axis(o.section) as SectionAxisView;
      const at = (b: number, p: number, y: number) => {
        const x =
          Math.round((v.navX0s[b] - v.left + (p * (v.navX1s[b] - v.navX0s[b])) / 400) * dpr) + Math.round(v.left * dpr);
        const i = (Math.round(y * dpr) * nav.width + x) * 4;
        return [px[i], px[i + 1], px[i + 2], px[i + 3]];
      };
      const near = (a: number[], b: number[], tol = 2) => a.every((c, k) => Math.abs(c - b[k]) <= tol);
      // The strip is 14 px in the 24 px row, from y 5; its middle is y 12.
      expect(near(at(0, 50, 12), [0, 255, 0, 255]), "negative band: the ground").toBe(true);
      expect(near(at(0, 140, 12), [0, 143, 0, 255]), "gpos50: ink at 112/255 over the ground").toBe(true);
      expect(near(at(0, 330, 12), [0, 19, 0, 255]), "gpos100: ink at 236/255").toBe(true);
      expect(near(at(0, 182, 12), [181, 72, 63, 255]), "the short arm's triangle near its base").toBe(true);
      expect(near(at(0, 218, 6), [0, 255, 0, 255]), "beyond the triangle's slope: the ground").toBe(true);
      expect(near(at(0, 238, 12), [181, 72, 63, 255]), "the long arm's triangle near its base").toBe(true);
      // The hairline is the strip's middle device row(s).
      expect(near(at(2, 200, dpr === 1 ? 11 : 11.5), [200, 200, 200, 255]), "a block without bands: a hairline").toBe(
        true,
      );
    });
});
