/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A section's inline style (seqspec:format#sections; seqspec:viewer#render-props): its keys over the preset's render
// props, the host's below both and its clamp over everything; a bad key a warning in the issues, the section drawn with
// the rest.

import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import { case3 } from "../../core/_test/format/cases.ts";
import { DefaultPresets, type RenderProps, SeqViewerConfig } from "../src/index.ts";
import { hoverDoc } from "./hover-fixtures.ts";
import { opened } from "./layout-fixtures.ts";

const S = "family.main";
/** Case 3 with its section's `style` set. */
function styled(style: unknown): SeqViewSpecDoc {
  const doc = case3();
  const walk = (n: SeqViewSpecDoc["root"]) => {
    if (n.ref === S) n.params = { ...n.params, style };
    n.children?.forEach(walk);
  };
  walk(doc.root);
  return doc;
}
/** [pinned, columns, the last pinned entry's bottom, the heads' top and height, the rows' top] (layout.test.ts). */
const band = async (style: unknown, render?: Partial<RenderProps>) => {
  const t = await opened(styled(style), { size: { [S]: [1000, 300] }, render });
  const l = t.layout(S);
  const last = l.pinnedCount - 1;
  return [
    l.pinnedCount,
    l.columnCount,
    l.entryTop[last] + l.entryHeight[last],
    l.headsTop,
    l.headsHeight,
    l.scrolledTop,
  ];
};

describe("an inline style's `column_heads`", () => {
  it('`"own_row"` reserves the heads a row of their own, as the render prop does; `"shared"` keeps the line', async () => {
    const own = [5, 1, 110, 110, 18, 128];
    expect(await band({ column_heads: "own_row" })).toEqual(own);
    expect(await band({ preset: "msa", column_heads: "own_row" })).toEqual(own);
    expect(await band(null, { dataColumnHeads: "ownRow" })).toEqual(own);
    expect(await band({ column_heads: "shared" }, { dataColumnHeads: "ownRow" })).toEqual([5, 1, 110, 92, 18, 110]);
  });
});

describe("precedence", () => {
  const props = async (style: unknown, render?: Partial<RenderProps>) => {
    const t = await opened(styled(style), { render });
    return Object.assign(() => t.ctx.layout.props(S), { t });
  };

  it("the document's key over its preset's value, the preset's over the host's", async () => {
    // The fixture `msa` sets 18 px rows; the host asks for 40.
    expect((await props("msa", { rowHeight: 40 }))().rowHeight).toBe(18);
    expect((await props({ preset: "msa", row_height: 30 }, { rowHeight: 40 }))().rowHeight).toBe(30);
    // No preset named: the default for a columns axis (msa) is the base.
    const p = (await props({ data_floor: 0, header_width: 90, row_borders: true }))();
    expect([p.rowHeight, p.dataFloor, p.headerWidth, p.rowBorders]).toEqual([18, 0, 90, true]);
  });

  it("a reader's preset switch keeps the document's keys over the new base", async () => {
    const p = await props({ preset: "msa", column_heads: "own_row", row_height: 30 });
    p.t.ctx.ui.setStyle(S, "bare");
    await p.t.ctx.graph.update();
    expect([p().dataColumnHeads, p().rowHeight, p().verticalNav]).toEqual(["ownRow", 30, "scrollbar"]);
  });

  it("the overview asked for where the preset gives it no width gets msa's 72 px", async () => {
    const p = await props({ preset: "bare", vertical_nav: "overview" });
    expect([p().verticalNav, p().verticalOverview]).toEqual(["overview", 72]);
    // The fixture msa's own width is kept.
    expect((await props({ preset: "msa", vertical_nav: "overview" }))().verticalOverview).toBe(24);
  });

  it("the host's clamp bounds the numbers over presets and styles, read when the layout next updates", async () => {
    const p = await props("msa");
    p.t.ctx.config.set(SeqViewerConfig.Render.Clamp, { min: { dataFloor: 96 }, max: { rowHeight: 20 } });
    await p.t.ctx.setDocument(styled({ preset: "msa", data_floor: 0, row_height: 30 }));
    expect([p().dataFloor, p().rowHeight]).toEqual([96, 20]);
    expect(p.t.vp(S).axis?.floor).toBe(96);
  });
});

describe("a bad key or value", () => {
  it("is a warning in the issues naming what is allowed; the section is drawn with the rest", async () => {
    const t = await opened(styled({ preset: "msa", colum_heads: "own_row", data_floor: -4, row_height: 30 }));
    const listed = t.ctx.issues.current.value
      .filter((i) => i.code === "svs.section.style")
      .map((i) => [i.code, i.severity, i.path.split("/").pop()]);
    expect(listed).toEqual([
      ["svs.section.style", "warning", "colum_heads"],
      ["svs.section.style", "warning", "data_floor"],
    ]);
    expect(t.ctx.issues.current.value.find((i) => i.path.endsWith("colum_heads"))?.message).toContain(
      'unknown style key "colum_heads"; the keys: preset, column_heads, data_floor',
    );
    const p = t.ctx.layout.props(S);
    expect([p.rowHeight, p.dataFloor, p.dataColumnHeads]).toEqual([30, 64, "shared"]);
    // The next document clears them.
    await t.ctx.setDocument(styled("msa"));
    expect(t.ctx.issues.current.value.filter((i) => i.code === "svs.section.style")).toEqual([]);
  });
});

describe("row borders", () => {
  it("the built-in tracks preset draws none; a document's `row_borders` or a host's render props turn them on", async () => {
    const borders = async (style?: Record<string, unknown>, render?: Partial<RenderProps>) => {
      const doc = hoverDoc(40, style ? { style } : {});
      const t = await opened(doc, { presets: false, registry: [{ presets: DefaultPresets }], render });
      const S0 = t.ctx.viewports.sections.value[0];
      return [t.ctx.layout.presetOf(0)?.name, t.ctx.layout.props(S0).rowBorders];
    };
    expect(await borders()).toEqual(["tracks", false]);
    expect(await borders({ row_borders: true })).toEqual(["tracks", true]);
    expect(await borders(undefined, { rowBorders: true })).toEqual(["tracks", true]);
  });
});
