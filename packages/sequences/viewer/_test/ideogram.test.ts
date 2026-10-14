/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ideogram navigator (ST5; seqspec:viewer#extension-points): each block's axis as a chromosome's bands, drawn from
// named data, found by the block's name; its window and gestures are any navigator's. Bands as ink at an alpha by
// stain, the centromere's `acen` pair as two triangles, a block without bands a hairline; hg38's committed cytobands
// (_test/fixtures/hg38-cytobands.json) in whole, and over the streamed genome's chromosomes.

import { readFileSync } from "node:fs";
import { describe, expect, it } from "vitest";
import { GRCH38, hg38Streams, streamedDoc } from "../../core/_test/format/stream-provider.js";
import type { SectionAxisView } from "../src/axis.js";
import { SkinToken } from "../src/draw/backend.js";
import { Dirty } from "../src/draw/composer.js";
import { type CytoBands, ideogramNavigator, type PresetEntry } from "../src/index.js";
import { WindowKind } from "../src/windows.js";
import { blocksDoc } from "./axis-fixtures.js";
import { composed, recording } from "./composer-fixtures.js";
import { driven } from "./input-fixtures.js";
import { GENOME } from "./stream-fixtures.js";

const HG38 = JSON.parse(readFileSync(new URL("./fixtures/hg38-cytobands.json", import.meta.url), "utf8")) as {
  chromosomes: CytoBands;
};
const KARYO: PresetEntry = { name: "karyo", geometry: "linear", render: {}, navigator: "ideogram" };
const entry = (bands: CytoBands) => ({ navigators: [ideogramNavigator(bands)], presets: [KARYO] });

/** Three sequences of 400, named c0..c2 (`c0/unit`): bands for the first two. */
const BANDS: CytoBands = {
  c0: [
    [0, 100, "p2", "gneg"],
    [100, 180, "p1", "gpos50"],
    [180, 220, "p11", "acen"],
    [220, 260, "q11", "acen"],
    [260, 400, "q1", "gpos100"],
  ],
  // Named without "chr": a block "chr1" or "1" finds it too.
  c1: [
    [0, 200, "p1", "gvar"],
    [200, 400, "q1", "stalk"],
  ],
};

async function karyo(doc = blocksDoc(3, 400), bands = BANDS, width = 840) {
  const t = await composed(doc, { width, height: 200, dpr: 1, registry: [entry(bands)] });
  t.ctx.ui.setStyle(t.section, "karyo");
  await t.ctx.graph.update();
  return t;
}
const rects = (b: ReturnType<typeof recording>, color: number) => {
  const r = b.rectItems;
  const out: number[][] = [];
  for (let i = 0; i < r.n; i++) if (r.color[i] === color) out.push([r.x[i], r.y[i], r.w[i], r.h[i], r.alpha[i]]);
  return out;
};

describe("the ideogram navigator draws", () => {
  it("bands by stain, the centromere as two triangles, a block without bands as a hairline", async () => {
    const t = await karyo();
    const be = t.frame(recording(), Dirty.nav);
    const av = t.ctx.layout.axis(t.section) as SectionAxisView;
    expect(t.ctx.layout.navigatorOf(0)?.name).toBe("ideogram");
    const x = (b: number, p: number) => Math.round(av.navX0s[b] + (p * (av.navX1s[b] - av.navX0s[b])) / 400);
    const ink = rects(be, SkinToken.ink);
    // Row [0, 24), the strip 14 px centered: y 5. Block 0's gpos50 and gpos100, block 1's gvar and stalk.
    expect(ink).toEqual([
      [x(0, 100), 5, x(0, 180) - x(0, 100), 14, 112],
      [x(0, 260), 5, x(0, 400) - x(0, 260), 14, 236],
      [x(1, 0), 5, x(1, 200) - x(1, 0), 14, 80],
      [x(1, 200), 5, x(1, 400) - x(1, 200), 14, 40],
    ]);
    // Two triangles on the centromere's bands: the short arm's points right, the long arm's left.
    const p = be.pathItems;
    expect(p.n).toBe(1);
    expect(p.color[0]).toBe(0xb5483f);
    expect(Array.from(p.verbs.subarray(0, p.nv))).toEqual([0, 1, 1, 6, 0, 1, 1, 6]);
    const c = Array.from(p.coords.subarray(0, p.nc));
    expect(c).toEqual([x(0, 180), 5, x(0, 220), 12, x(0, 180), 19, x(0, 260), 5, x(0, 220), 12, x(0, 260), 19]);
    // An outline per chromosome, a hairline across the third block, which has none.
    expect(rects(be, SkinToken.inkQuiet).length).toBe(2);
    expect(rects(be, SkinToken.inkFaint)).toEqual([
      [Math.round(av.navX0s[2]), 11, Math.round(av.navX1s[2]) - Math.round(av.navX0s[2]), 1, 255],
    ]);
    // The window's wash over it, as over any map (the window is the whole fit here: no wash outside it).
    expect(t.composer.stats.frames).toBeGreaterThan(0);
    t.ctx.dispose();
  });

  it("hg38's committed cytobands over the streamed genome's chromosomes: every band, both centromere halves", async () => {
    const three = GRCH38.slice(0, 3);
    const s = hg38Streams({ chromosomes: three, letters: true });
    const loaders = { load: () => Promise.reject(new Error("none")), stream: s.stream };
    const t = await composed(streamedDoc({ values: false }), {
      width: 1000,
      height: 300,
      dpr: 1,
      registry: [entry(HG38.chromosomes)],
      loaders,
    });
    t.ctx.ui.setStyle(GENOME, "karyo");
    await t.ctx.graph.update();
    const be = t.frame(recording(), Dirty.nav);
    const av = t.ctx.layout.axis(GENOME) as SectionAxisView;
    expect(av.axis?.names).toEqual(three.map(([k]) => `hg38/${k}`));
    const drawn = three.flatMap(([k]) => HG38.chromosomes[k]);
    const plain = drawn.filter(([, , , stain]) => stain !== "gneg" && stain !== "acen").length;
    expect(rects(be, SkinToken.ink).length).toBe(plain);
    expect(be.pathItems.nv).toBe(8 * 3);
    expect(rects(be, SkinToken.inkQuiet).length).toBe(3);
    // Chromosome 1's centromere: the short arm's triangle starts where its first `acen` band does.
    const [acen] = HG38.chromosomes.chr1.filter((r) => r[3] === "acen");
    const at = Math.round(av.navX0s[0] + (acen[0] * (av.navX1s[0] - av.navX0s[0])) / three[0][1]);
    expect(be.pathItems.coords[0]).toBe(at);
    t.ctx.dispose();
  });
});

describe("the committed hg38 cytobands", () => {
  it("are contiguous bands of the primary chromosomes, ending at their lengths, with a centromere's two halves", () => {
    const names = Object.keys(HG38.chromosomes);
    expect(names).toEqual([...Array.from({ length: 22 }, (_, k) => `chr${k + 1}`), "chrX", "chrY"]);
    for (const [k, length] of GRCH38) {
      const rows = HG38.chromosomes[k];
      expect(rows[0][0], k).toBe(0);
      expect(rows.at(-1)?.[1], k).toBe(length);
      for (let i = 1; i < rows.length; i++) expect(rows[i][0], `${k} ${rows[i][2]}`).toBe(rows[i - 1][1]);
      expect(rows.filter((r) => r[3] === "acen").length, k).toBe(2);
    }
  });
});

describe("the ideogram's gestures are the navigator behavior's", () => {
  it("a press jumps the window to the pointer, an end dragged zooms, the body pans", async () => {
    const t = await driven(blocksDoc(3, 400), { size: {}, registry: [entry(BANDS)] });
    const sec = t.ctx.viewports.sections.value[0];
    t.ctx.ui.setStyle(sec, "karyo");
    await t.ctx.graph.update();
    t.ctx.viewports.get(sec)?.zoomTo(0, 100, 160);
    const nav = () => t.frame(sec, WindowKind.navigator);
    const av = () => t.ctx.layout.axis(sec) as SectionAxisView;
    const y = 10;
    const x = av().navX0s[1] + (av().navX1s[1] - av().navX0s[1]) / 4;
    t.pointer("down", sec, x, y);
    t.pointer("up", sec, x, y);
    expect(Math.abs((nav()[8] + nav()[9]) / 2 - x)).toBeLessThan(1);
    // Its right end dragged: the left end held, the right one under the pointer.
    const [x0, x1] = [nav()[8], nav()[9]];
    t.drag(sec, x1 - 2, y, x1 + 38, y);
    expect(Math.abs(nav()[8] - x0)).toBeLessThan(0.5);
    expect(nav()[9] - x1).toBeGreaterThan(30);
    // The body dragged back by 30 px.
    const m = (nav()[8] + nav()[9]) / 2;
    const before = nav()[8];
    t.drag(sec, m, y, m - 30, y);
    expect(Math.abs(nav()[8] - (before - 30))).toBeLessThan(1);
  });
});
