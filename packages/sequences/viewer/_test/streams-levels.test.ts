/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A streamed alignment whose manifest lists some of the grid's levels (seqspec:sources#levels-and-the-zoom-gate; VT-03):
// the level the view picks is snapped to one the source serves (the coarsest no wider, else the nearest coarser), so
// it is asked and lands, and the draws read the level the manager asked (`listedLevels`, one rule for both).

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { DataPatch, Loaders, StreamSource } from "@molstar/seq-core/format";
import { tilesStream } from "@molstar/seq-core/stream/tiles";
import { tileServer } from "@molstar/seq-core/stream/tiles/server";
import { describe, expect, it } from "vitest";
import { RegionState } from "../src/index.js";
import type { Streams } from "../src/streams.js";
import { listedLevels } from "../src/streams.js";
import { layer, opened, preparedRenderers } from "./layout-fixtures.js";

const URI = "https://example.org/lazy";
const flush = () => new Promise((ok) => setTimeout(ok, 0));

/** A tiles provider whose manifest keeps the levels `keep` passes; the data levels it was asked. */
async function sparse(keep: (ly: number, lx: number) => boolean) {
  const server = tileServer({ rowCount: 40_000, columnCount: 3000 });
  const inner = tilesStream((u) => (u === URI ? server.transport : null));
  const asked: number[] = [];
  const loaders: Loaders = {
    load: (r) => Promise.reject(new Error(`no fixture ${r.uri}`)),
    stream: (src: StreamSource) => {
      const p = inner(src);
      if (!p) return null;
      return {
        fetch: (r, signal) =>
          p.fetch(r, signal).then((x: DataPatch) => {
            if (typeof r.level === "number" && r.level >= 0) asked.push(r.level);
            if (x.kind !== "manifest") return x;
            const levels = (x.manifest.levels as [number, number][]).filter(([a, b]) =>
              keep(Math.log2(a), Math.log2(b)),
            );
            return { ...x, manifest: { ...x.manifest, levels } };
          }),
      };
    },
  };
  const b = createSeqViewSpec({ title: "lazy" });
  const aln = b.data.alignment({
    ref: "lazy",
    row_count: 40_000,
    column_count: 3000,
    source: { uri: URI, protocol: "seqstar-tiles/1" },
  });
  const v = b.view({ ref: "v", label: "A", style: "bare", height: "box" });
  v.block({ alignment: aln });
  v.track({ label: "Rows" }).rows();
  const t = await opened(b.getDocument(), { renderers: preparedRenderers(), loaders, size: { "v.main": [1000, 600] } });
  const s = t.ctx.streams as Streams;
  const vp = t.vp("v.main");
  vp.zoomTo(0, 0, 3000);
  vp.setRowScale(0.01, 0);
  t.layout("v.main");
  s.pass();
  await flush();
  for (let i = 0; i < 8; i++) {
    t.scheduler.advance(16);
    s.pass();
    await flush();
  }
  const { k, m } = layer(t.ctx, "Rows");
  const ref = m.set.strings.get(m.resolved.layers.ref[k]);
  return { t, s, ref, asked };
}

describe("an alignment's levels the manifest does not list (VT-03)", () => {
  it("zoomed out over fine levels only, the coarsest served is asked and lands", async () => {
    const { s, ref, asked } = await sparse((ly, lx) => ly <= 2 && lx <= 2);
    expect(asked.length).toBeGreaterThan(0);
    expect(s.state(ref, 0, 1500)).toBe(RegionState.loaded);
    const view = s.get(ref);
    const lv = Int32Array.of(9, 7, 12, 11);
    listedLevels(view?.served, lv);
    expect(Array.from(lv)).toEqual([2, 2, 2, 2]);
  });

  it("a level between served ones takes the coarsest no wider, else the nearest coarser", async () => {
    const { s, ref } = await sparse((ly, lx) => ly % 4 === 0 && lx % 4 === 0 && ly >= 4);
    const served = s.get(ref)?.served;
    const at = (ly: number, lx: number) => {
      const lv = Int32Array.of(ly, lx, ly, lx);
      listedLevels(served, lv);
      return [lv[0], lv[1]];
    };
    expect(at(6, 5)).toEqual([4, 4]);
    expect(at(8, 8)).toEqual([8, 8]);
    expect(at(2, 9)).toEqual([4, 8]);
    expect(at(13, 2)).toEqual([12, 0]);
    expect(s.state(ref, 0, 1500)).toBe(RegionState.loaded);
  });

  it("a local stream, or a manifest listing every level, keeps the levels chosen", async () => {
    const { s, ref } = await sparse(() => true);
    const lv = Int32Array.of(9, 7, 12, 11);
    listedLevels(s.get(ref)?.served, lv);
    expect(Array.from(lv)).toEqual([9, 7, 12, 11]);
    listedLevels(null, lv);
    expect(Array.from(lv)).toEqual([9, 7, 12, 11]);
  });
});
