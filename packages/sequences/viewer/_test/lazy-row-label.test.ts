/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A streamed (lazy) row's hover label reads what the region cache holds for it: "loading" before its tile arrived,
// its label and residue after (case 9's rows are streamed in tiles; the label provider reads no rows in memory).

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { proceduralCode } from "@molstar/seq-core/stream/tiles/server";
import { describe, expect, it } from "vitest";
import { HoverKind, KernelRenderers, type LayerRenderer } from "../src/index.ts";
import type { Streams } from "../src/streams.ts";
import { driven } from "./input-fixtures.ts";
import { lazyLoaders } from "./layout-fixtures.ts";

const [uri, ROWS, COLUMNS] = ["https://example.org/lazy-label", 3000, 400];

describe("a streamed row's hover label", () => {
  it("says it is loading before the tile arrived, then reads the row's label and residue from it", async () => {
    const b = createSeqViewSpec({ title: "lazy" });
    const aln = b.data.alignment({
      ref: "lazy",
      row_count: ROWS,
      column_count: COLUMNS,
      source: { uri, protocol: "seqstar-tiles/1" },
    });
    b.view({ ref: "v", label: "V", style: "msa", height: "box" })
      .block({ alignment: aln })
      .track({ label: "Rows" })
      .rows();
    const t = await driven(b.getDocument(), {
      renderers: [...KernelRenderers] as LayerRenderer[],
      loaders: lazyLoaders(uri, ROWS, COLUMNS),
      size: { "v.main": [1000, 400] },
    });
    t.vp("v.main").setView(20, 0, "api");
    await t.ctx.graph.update();
    const l = t.layout("v.main");
    // A cell with a residue (the procedural rows have a gap in about one cell in six).
    const h = t.ctx.hover.live.current;
    let [row, column] = [-1, -1];
    for (let dx = 30; proceduralCode(row, column) === 45 || row < 0; dx += 20) {
      t.pointer("move", "v.main", l.dataLeft + dx, l.scrolledTop + 5);
      expect(h.kind).toBe(HoverKind.aligned);
      [row, column] = [h.aligned.row as number, h.aligned.column as number];
    }
    const label = () => t.ctx.labels.current.value;
    // Nothing arrived: its tile is in flight, which is not a gap.
    expect(label()?.lines).toContainEqual({ key: "residue", value: "loading" });
    for (let i = 0; i < 4; i++) {
      (t.ctx.streams as Streams).pass();
      await new Promise((ok) => setTimeout(ok, 0));
    }
    // Arrived (the label was made again for it, the pointer never moved): the row's own label and its residue.
    const code = proceduralCode(row, column);
    expect(label()?.title).toMatch(new RegExp(`^SEQ${row}/1-\\d+$`));
    expect(label()?.lines).toContainEqual({ key: "residue", value: String.fromCharCode(code) });
    expect(label()?.lines).toContainEqual({ key: "column", value: `${column + 1}` });
    t.ctx.dispose();
  });
});
