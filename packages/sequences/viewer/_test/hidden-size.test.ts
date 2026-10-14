/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A size that leaves the data area no width (the view stays as it is) is still a size the layout is told of:
// `Seq*/Views and layouts/Section style › Data floor at a narrow width`'s first card (`data_floor: 0`, four data
// columns) was laid out at the height its box had before its rows (40 px) and kept it, showing one row, until some
// input bumped the viewport again.

import { describe, expect, it } from "vitest";
import { fixtureLoaders } from "../../core/_test/format/load-fixtures.ts";
import { createSeqViewer, DefaultRegistry, DefaultSeqViewerSpec, manualScheduler } from "../src/index.ts";
import { fixtureMeasurer } from "../src/testing/index.ts";
import { familyWithTree } from "./families.ts";

const STYLE = { preset: "msa", data_floor: 0 };
function styled() {
  const doc = familyWithTree(60, 12, { columns: true });
  const walk = (n: typeof doc.root): typeof doc.root => ({
    ...n,
    params: n.ref === "fam.main" ? { ...n.params, style: STYLE } : n.params,
    children: n.children?.map(walk),
  });
  return { ...doc, root: walk(doc.root) };
}

describe("a viewport that is hidden for want of a data width", () => {
  it("lays out for the height that came after the first layout: every row that fits is visible, no interaction", async () => {
    const scheduler = manualScheduler();
    const ctx = createSeqViewer(
      { ...DefaultSeqViewerSpec(), registry: [DefaultRegistry] },
      { loaders: fixtureLoaders(), scheduler, text: fixtureMeasurer() },
    );
    await ctx.setDocument(styled());
    const section = "fam.main";
    const vp = ctx.viewports.get(section);
    // The box is 40 px tall before its rows are known; then it grows to hold them. 300 px leaves the letters nothing.
    vp?.setSize(300, 40, 1);
    scheduler.advance(50);
    const first = ctx.layout.get(section)?.current;
    expect(first?.columnCount).toBe(4);
    expect((first?.dataRight ?? 1) - (first?.dataLeft ?? 0)).toBeLessThanOrEqual(0);
    vp?.setSize(300, 242, 1);
    scheduler.advance(50);
    const l = ctx.layout.get(section)?.current;
    const rows = ctx.layout.visibleRows(section)?.current;
    expect(l?.rows.count).toBe(60);
    // The rows fill the 242 px under the pinned ruler (18 px each, the last partly): not the one the 40 px box held.
    const fit = Math.floor((242 - (l?.scrolledTop ?? 0)) / 18) + 1;
    expect(rows?.count).toBe(fit);
    expect(fit).toBeGreaterThan(10);
  });
});
