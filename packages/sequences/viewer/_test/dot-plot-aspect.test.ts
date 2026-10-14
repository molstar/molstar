/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A dot plot's `aspect_ratio` (seqspec:format#layer-types): its track's height is the data area's width ÷ the ratio,
// asked again when the section's width changes (one layout a resize), within the room a box shows its rows in; `null`
// keeps the fixed height.

import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import { KernelRenderers, type LayerRenderer } from "../src/index.ts";
import { dotPlotDoc } from "./dot-plot-fixtures.ts";
import { driven } from "./input-fixtures.ts";

const SECTION = "v.s";
/** The fixture with every dot plot layer's `aspect_ratio` set (undefined: the default). */
function withRatio(ratio: number | null | undefined): SeqViewSpecDoc {
  const d = structuredClone(dotPlotDoc());
  type N = { kind: string; params?: Record<string, unknown>; children?: N[] };
  const walk = (n: N) => {
    if (n.kind === "layer" && n.params?.type === "dot_plot" && ratio !== undefined)
      n.params = { ...n.params, aspect_ratio: ratio };
    n.children?.forEach(walk);
  };
  walk(d.root as N);
  return d;
}
async function plot(ratio: number | null | undefined, width: number, height = 3000) {
  const t = await driven(withRatio(ratio), {
    renderers: KernelRenderers as LayerRenderer[],
    size: { [SECTION]: [width, height] },
  });
  const sizeOf = () => {
    const l = t.layout(SECTION);
    const e = l.entryCount - 1;
    return { plot: l.entryHeight[e], data: l.dataRight - l.dataLeft };
  };
  return { t, sizeOf };
}

describe("a dot plot's aspect ratio", () => {
  it("square by default: the data area's width, at two widths, asked again on a resize", async () => {
    const { t, sizeOf } = await plot(undefined, 800);
    const a = sizeOf();
    expect(a.plot).toBe(Math.round(a.data));
    const builds = t.layouts.stats.builds;
    t.vp(SECTION).setSize(1200, 3000, 1);
    t.advance(16);
    const b = sizeOf();
    expect(b.plot).toBe(Math.round(b.data));
    expect(b.plot).toBeGreaterThan(a.plot);
    // One layout for the resize.
    expect(t.layouts.stats.builds - builds).toBe(1);
  });

  it("width ÷ ratio at another ratio; a box's room caps it; null keeps the fixed height", async () => {
    const two = (await plot(2, 1200)).sizeOf();
    expect(two.plot).toBe(Math.round(two.data / 2));
    // In a box 600 px tall the plot is at most the rows' room less its margin (or the stretch's floor, 1200 ÷ 2.5).
    const box = (await plot(1, 1200, 600)).sizeOf();
    expect(box.plot).toBeLessThanOrEqual(Math.max(540, Math.ceil(box.data / 2.5)));
    expect(box.plot).toBeGreaterThanOrEqual(120);
    const fixed = (await plot(null, 1200)).sizeOf();
    expect(fixed.plot).toBe(320);
    expect((await plot(null, 800)).sizeOf().plot).toBe(320);
  });

  it("never flatter than the ratio stretched by dotPlotMaxStretch: at 1000, 2000 and 4000 px of data, a box's room caps the square until the stretch's floor wins", async () => {
    const { t, sizeOf } = await plot(1, 1200, 700);
    const off = 1200 - sizeOf().data;
    const heights: number[] = [];
    for (const data of [1000, 2000, 4000]) {
      const builds = t.layouts.stats.builds;
      t.vp(SECTION).setSize(data + off, 700, 1);
      t.advance(16);
      // One layout for the resize.
      expect(t.layouts.stats.builds - builds).toBe(1);
      const l = t.layout(SECTION);
      const room = Math.floor(700 - l.scrolledTop);
      // The room less its margin: a tenth of it, at least two rows.
      const cap = room - Math.max(Math.round(room / 10), 2 * t.ctx.layout.props(SECTION).rowHeight);
      const s = sizeOf();
      expect(Math.round(s.data)).toBe(data);
      // The square capped by the room less its margin, raised to 120; then at least width ÷ (1 × 2.5), which wins.
      expect(s.plot).toBe(Math.max(Math.max(120, Math.min(cap, data)), Math.ceil(data / 2.5)));
      if (data === 1000) expect(s.plot).toBe(cap);
      heights.push(s.plot);
    }
    // The room caps the first; the stretch's floor the others, past the visible rows.
    expect(heights[1]).toBe(800);
    expect(heights[2]).toBe(1600);
    expect(heights[0]).toBeLessThan(700);
  });

  it("settles: a bounce (a page scrollbar coming and going) keeps the taller height after one round; a resize lands exactly", async () => {
    const { t, sizeOf } = await plot(1, 1200);
    const at = (w: number) => {
      t.vp(SECTION).setSize(w, 3000, 1);
      t.advance(16);
      return sizeOf();
    };
    const a = sizeOf().plot;
    // The scrollbar comes (narrower), goes (the plot grows back), comes again: held at the taller height from then on.
    const b = at(1185).plot;
    expect(b).toBeLessThan(a);
    expect(at(1200).plot).toBe(a);
    expect(at(1185).plot).toBe(a);
    expect(at(1185).plot).toBe(a);
    // A slow resize narrower follows its target exactly (no 31 px left over).
    for (const w of [1180, 1176, 1172, 1168]) {
      const s = at(w);
      expect(s.plot).toBe(Math.round(s.data));
    }
  });
});
