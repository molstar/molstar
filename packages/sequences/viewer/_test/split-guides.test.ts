/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Split rows' guides (seqspec:viewer#render-props): every row of an expanded split track but the last has a hairline
// at its bottom edge in the grid tone, drawn whether or not `rowBorders` is on, under the data, and only where the row
// is at least 4 device px tall; an unsplit track has none, and a real border (a row over 8 px) takes its place.

import { describe, expect, it } from "vitest";
import { case7, case13 } from "../../core/_test/format/cases.ts";
import { Dirty, SkinToken } from "../src/draw/index.ts";
import type { LayerRenderer } from "../src/index.ts";
import { VariantsRenderer } from "../src/renderers/variants.ts";
import { composed, opNames } from "./composer-fixtures.ts";
import { rects } from "./draw-fixtures.ts";

const ROWS = 22;

async function guides(o: { dpr?: number; rowHeight?: number; rowBorders?: boolean; doc?: ReturnType<typeof case13> }) {
  const dpr = o.dpr ?? 2;
  const t = await composed(o.doc ?? case13(), {
    width: 1000,
    height: 900,
    dpr,
    render: { rowBorders: o.rowBorders ?? false },
    replace: o.rowHeight ? { variants: { ...VariantsRenderer, splitHeight: o.rowHeight } as LayerRenderer } : undefined,
  });
  const rec = t.frame(undefined, Dirty.all);
  const all = rects(rec).filter(([, , , , c]) => c === SkinToken.grid);
  // The ruler's baseline is a full-strength grid line too (in the first 60 px); a guide is fainter.
  const grid = all.filter(([, y]) => y > 60 * dpr);
  return { t, grid, rec };
}

describe("a split track's row guides", () => {
  it("case 13: a hairline at the bottom of each of the 22 rows but the last, in the grid tone, half-strength", async () => {
    const { grid, t } = await guides({});
    const dpr = 2;
    // The rows are 8 px tall: the guides are 2 device px (a CSS px's hairline at DPR 2) wide across the data area.
    expect(grid).toHaveLength(ROWS - 1);
    const l = t.layout(t.section);
    for (const [x, , w, h, , alpha] of grid) {
      expect([x, w, h]).toEqual([Math.round(l.dataLeft * dpr), Math.round((l.dataRight - l.dataLeft) * dpr), 2]);
      expect(alpha).toBeLessThan(255);
    }
    // One per row edge: 16 device px apart.
    for (let i = 1; i < grid.length; i++) expect(grid[i][1] - grid[i - 1][1]).toBe(16);
  });

  it("with `rowBorders` on the 8 px rows still get the guides, and the last row the track's border (full strength)", async () => {
    const { grid } = await guides({ rowBorders: true });
    expect(grid.filter((r) => r[5] < 255)).toHaveLength(ROWS - 1);
    expect(grid.filter((r) => r[5] === 255).length).toBeGreaterThanOrEqual(1);
  });

  it("rows 4 device px tall have guides, 3 none (at DPR 1)", async () => {
    expect((await guides({ dpr: 1, rowHeight: 4 })).grid).toHaveLength(ROWS - 1);
    expect((await guides({ dpr: 1, rowHeight: 3 })).grid).toHaveLength(0);
    // The threshold is in device px: 2 CSS px at DPR 2.
    expect((await guides({ dpr: 2, rowHeight: 2 })).grid).toHaveLength(ROWS - 1);
  });

  it("rows over 8 px with `rowBorders` have borders, not guides; without it, guides", async () => {
    const borders = (await guides({ rowHeight: 12, rowBorders: true })).grid;
    expect(borders.every((r) => r[5] === 255)).toBe(true);
    expect(borders.length).toBeGreaterThanOrEqual(ROWS);
    const plain = (await guides({ rowHeight: 12 })).grid;
    expect(plain).toHaveLength(ROWS - 1);
    expect(plain.every((r) => r[5] < 255)).toBe(true);
  });

  it("an unsplit document has no guides", async () => {
    expect((await guides({ doc: case7() })).grid).toHaveLength(0);
  });

  it("the guides are under the data: the first guide is drawn before the variants' glyphs", async () => {
    const { rec } = await guides({});
    const ops = opNames(rec);
    const rectsAt = ops.indexOf("rects");
    const glyphsAt = ops.indexOf("glyphs");
    expect(rectsAt).toBeGreaterThanOrEqual(0);
    expect(glyphsAt).toBeGreaterThan(rectsAt);
  });
});
