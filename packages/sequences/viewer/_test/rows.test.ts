/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// RN3's rows draws on the recording backend (seqspec:viewer#extension-points; seqspec:viewer#the-render-backend,
// seqspec:sources#levels-and-the-zoom-gate): rows at each level of detail and row height (dense rows in one `rows`
// call, single rows with gaps and letters), the tree header and a hovered clade, the data columns drawn through their
// renderers (an unknown type as `text`, a field provider's values, a long name cut with an ellipsis), lazy rows from
// tiles as they arrive (coarse under fine), the vertical overview's map, the `seqv:map` miniature as a keyed image. The
// old test's pixel parts are CS2's (copies equal a redraw) and CS1b's (the selection of cells on one row: the overlay's
// row rects).

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { TILE } from "@molstar/seq-core/format";
import { NotServed, tilesStream } from "@molstar/seq-core/stream/tiles";
import { describe, expect, it } from "vitest";
import { ChannelMarks } from "../src/channels.ts";
import { SkinToken, TextAlign } from "../src/draw/backend.ts";
import { Dirty, VerticalNav } from "../src/draw/index.ts";
import { MapView } from "../src/draw/miniature.ts";
import { DrawView } from "../src/draw/view.ts";
import { KernelRenderers, type LayerRenderer } from "../src/index.ts";
import { letterLines } from "../src/renderer-helpers.ts";
import { RowsRenderer } from "../src/renderers/rows.ts";
import type { RowsState } from "../src/rows.ts";
import type { Streams } from "../src/streams.ts";
import { AffinityField, LayerFrame, Op, RecordingBackend, StarsFormat } from "../src/testing/index.ts";
import { TREE_INSET } from "../src/tree.ts";
import { set } from "./axis-fixtures.ts";
import { drawn, layerOf, rects, texts } from "./draw-fixtures.ts";
import { familyWithTree, wideFamily } from "./families.ts";
import { layer, lazyLoaders, opened } from "./layout-fixtures.ts";

const count = (r: RecordingBackend, op: number) => {
  let n = 0;
  for (let k = 0; k < r.opCount; k++) if (r.ops[4 * k] === op) n++;
  return n;
};
/** The region cache's passes until its local tiles for the view are made. */
async function settle(ctx: { streams: unknown }) {
  for (let i = 0; i < 4; i++) {
    (ctx.streams as Streams).pass();
    await new Promise((ok) => setTimeout(ok, 0));
  }
}
/** The Members track's rows layer drawn through the kit's composer: its recording and entries drawn. */
function frame(ctx: Awaited<ReturnType<typeof drawn>>["ctx"], label = "Members") {
  const { k, section } = layerOf(ctx, label);
  const entry = ctx.registry.get("renderers", "rows") as LayerRenderer;
  const f = new LayerFrame(ctx, entry, section, k);
  const n = f.draw(f.direct);
  return { rec: f.direct, n, k, section };
}
/** A view on the rows layer as the composer sets it for one row (`row`) at (`top`, `height`). */
function view(ctx: Awaited<ReturnType<typeof drawn>>["ctx"], rec: RecordingBackend, label = "Members") {
  const { k, section } = layerOf(ctx, label);
  const p = ctx.graph.assembled.value?.view?.prepared[k];
  const v = new DrawView().from(ctx.layout.axis(section) as never);
  v.backend = rec;
  v.layer = k;
  v.setThemes(p?.themes ?? []);
  v.text = ctx.text;
  v.props = ctx.layout.props(section);
  v.skin = ctx.skin.value;
  return { v, s: p?.state as RowsState, section, k };
}

describe("rows, the tree, columns draw", () => {
  it.each([1, 2])("at each level of detail and row height, DPR %i; dense rows in one call", async (dpr) => {
    const t = await drawn(wideFamily(160, 80), { dpr, width: 1000, height: 420 });
    const { section } = layerOf(t.ctx, "Members");
    for (const cell of [16, 4, 1.2, 0.3])
      for (const scale of [1, 0.15, 0.01]) {
        t.ctx.viewports.get(section)?.setView(cell, 0, "api");
        t.ctx.viewports.get(section)?.setRowScale(scale, 0);
        await t.ctx.graph.update();
        await settle(t.ctx);
        const at = `cell ${cell}, row scale ${scale}`;
        const { rec, n } = frame(t.ctx);
        const l = t.ctx.layout.get(section)?.current;
        const dense = (l?.rows.height ?? 0) < t.ctx.layout.props(section).denseRowHeight;
        const images = count(rec, Op.image);
        expect(images, at).toBeGreaterThan(0);
        // Dense: one `rows` call, a few tile parts (odd sample counts of 512 device px at most: 18 at most here); else
        // a row's own images, at least one a row.
        if (dense) expect([n, images <= 24], at).toEqual([1, true]);
        else expect(images, at).toBeGreaterThanOrEqual(n);
        // Gaps (the ground over the cells' right edges) between wide cells of a readable row only; letters over them.
        const cells = t.ctx.viewports.get(section)?.live.current.cell ?? 0;
        const wide = cells > 6;
        const gaps = rects(rec).filter((r) => r[2] === 1 && r[4] === SkinToken.surface).length;
        expect(gaps > 0, at).toBe(wide && !dense);
        const letters = texts(rec).length;
        expect(letters > 0, at).toBe(cells >= 8 && !dense && (l?.rows.height ?? 0) >= 10);
      }
  });

  it("letters are inked to read on their cell, a gap quiet; single rows carry their row's residues", async () => {
    const t = await drawn(familyWithTree(40, 12), { width: 1000, height: 420 });
    const { section } = layerOf(t.ctx, "Members");
    t.ctx.viewports.get(section)?.setView(16, 0, "api");
    await t.ctx.graph.update();
    const { rec } = frame(t.ctx);
    const ts = texts(rec);
    expect(ts.length).toBeGreaterThan(50);
    for (const x of ts) expect(x.text.length).toBe(1);
    for (const x of ts) expect([0x181a1e, 0xffffff, SkinToken.inkQuiet]).toContain(x.color);
    expect(ts.filter((x) => x.text === "-").every((x) => x.color === SkinToken.inkQuiet)).toBe(true);
    expect(ts.filter((x) => x.text !== "-").every((x) => x.color >= 0)).toBe(true);
  });

  it.each([1, 2])("the tree header and a hovered clade, DPR %i", async (dpr) => {
    const t = await drawn(familyWithTree(40, 10), { dpr, width: 1000, height: 600 });
    const rec = new RecordingBackend();
    const { v, s } = view(t.ctx, rec);
    const l = t.ctx.layout.get("fam.main")?.current;
    if (!l) throw new Error("no layout");
    v.top = l.scrolledTop + l.entryTop[l.rows.entry];
    v.height = l.rows.height;
    v.left = 0;
    const last = l.rows.count - 1;
    RowsRenderer.header?.draw?.(v, s, {}, 0, last);
    const whole = rec.pathItems;
    expect([whole.n, whole.color[0], whole.mode[0], whole.width[0]]).toEqual([1, SkinToken.inkQuiet, 1, 2 * dpr]);
    const verbs = whole.nv;
    // Every branch within the tree's box.
    for (let i = 0; i < whole.nc; i += 2) expect(whole.coords[i]).toBeLessThanOrEqual(TREE_INSET * dpr);
    rec.reset();
    RowsRenderer.header?.hover?.(v, s, {}, 2, 0, last);
    expect([rec.pathItems.color[0], rec.pathItems.width[0]]).toEqual([SkinToken.accent, 4 * dpr]);
    expect(rec.pathItems.nv).toBeGreaterThan(0);
    expect(rec.pathItems.nv).toBeLessThan(verbs);
    // Rows thinner than a device pixel draw no leaves' branches.
    rec.reset();
    v.height = 0.3 / dpr;
    RowsRenderer.header?.draw?.(v, s, {}, 0, last);
    expect(rec.pathItems.nv).toBeLessThan(verbs);
  });

  it("a band of rows: a connector cut at its edges runs a line width past them (no seam where bands meet)", async () => {
    const t = await drawn(familyWithTree(40, 10), { dpr: 2, width: 1000, height: 600 });
    const rec = new RecordingBackend();
    const { v, s } = view(t.ctx, rec);
    const l = t.ctx.layout.get("fam.main")?.current;
    if (!l) throw new Error("no layout");
    v.top = l.scrolledTop + l.entryTop[l.rows.entry];
    v.height = l.rows.height;
    v.left = 0;
    const [first, last] = [10, 19];
    RowsRenderer.header?.draw?.(v, s, {}, first, last);
    const p = rec.pathItems;
    const [top, bottom] = [Math.round((v.top + first * v.height) * 2), Math.round((v.top + (last + 1) * v.height) * 2)];
    let [above, below] = [false, false];
    // Move and line pairs: a vertical one is a connector.
    for (let i = 0; i + 3 < p.nc; i += 4) {
      if (p.coords[i] !== p.coords[i + 2]) continue;
      if (Math.min(p.coords[i + 1], p.coords[i + 3]) < top) above = true;
      if (Math.max(p.coords[i + 1], p.coords[i + 3]) > bottom) below = true;
    }
    expect([above, below]).toEqual([true, true]);
  });

  it("collapsed clades draw as wedges: a filled triangle each", async () => {
    const t = await drawn(familyWithTree(40, 10), { width: 1000, height: 600 });
    const rec = new RecordingBackend();
    const { v, s, section } = view(t.ctx, rec);
    const wedges = () => {
      const l = t.ctx.layout.get(section)?.current;
      if (!l) throw new Error("no layout");
      rec.reset();
      v.top = l.scrolledTop + l.entryTop[l.rows.entry];
      v.height = l.rows.height;
      RowsRenderer.header?.draw?.(v, s, {}, 0, l.rows.count - 1);
      return Array.from({ length: rec.pathItems.n }, (_, i) => rec.pathItems.mode[i]).filter((m) => m === 0).length;
    };
    expect(wedges()).toBe(0);
    t.ctx.ui.toggleClade("family", 2);
    await t.ctx.graph.update();
    expect(wedges()).toBe(1);
  });
});

/** The family's data columns drawn for alignment row 3 into a fresh recording, one cell each. */
async function cells(width: number, dpr = 1) {
  const t = await drawn(familyWithTree(40, 12, { columns: true }), { dpr });
  const rec = new RecordingBackend();
  const { v, s } = view(t.ctx, rec);
  v.top = 50;
  v.height = 14;
  v.columnX = 100;
  v.columnWidth = width;
  return {
    t,
    rec,
    v,
    s,
    draw: (c: number, row = 3) => {
      rec.reset();
      const dc = s.dataColumns[c];
      t.ctx.registry.get("columns", dc.renderer)?.draw?.(v, dc.state as never, row);
    },
  };
}

describe("data columns", () => {
  it("number right-aligned, fraction a bar and a percentage, text left, category a swatch and a label", async () => {
    const { rec, draw, s } = await cells(80);
    expect(s.dataColumns.map((c) => c.renderer)).toEqual(["number", "fraction", "text", "category"]);
    expect(s.dataColumns.map((c) => c.label)).toEqual(["Len", "Id.", "Species", "Clade"]);
    draw(0);
    expect(texts(rec).map((x) => [x.text, x.align, x.x, x.color])).toEqual([
      ["28", TextAlign.right, 174, SkinToken.ink],
    ]);
    draw(1);
    // The bar's track and its share of it, then the percentage, quiet.
    expect(rects(rec).map((r) => [r[0], r[2], r[4]])).toEqual([
      [104, 36, SkinToken.inkFaint],
      [104, Math.round(36 * 0.11), SkinToken.accent],
    ]);
    expect(texts(rec).map((x) => [x.text, x.align, x.color])).toEqual([["11%", TextAlign.right, SkinToken.inkQuiet]]);
    draw(2);
    expect(texts(rec).map((x) => [x.text, x.align, x.x])).toEqual([["Species 3", TextAlign.left, 106]]);
    draw(3);
    expect(rects(rec).length).toBe(1);
    expect(rects(rec)[0][4]).toBeGreaterThanOrEqual(0);
    expect(texts(rec).map((x) => x.text)).toEqual(["clade 3"]);
    // A row without a value draws nothing; a row past the annotation too.
    draw(0, -1);
    expect(rec.opCount).toBe(0);
  });

  it("a long name is cut with an ellipsis at the column's width", async () => {
    const full = await cells(110);
    full.draw(2);
    expect(texts(full.rec)[0].text).toBe("Species 3");
    const narrow = await cells(36);
    narrow.draw(2);
    const text = texts(narrow.rec)[0].text;
    expect([text.endsWith("…"), text.length < "Species 3".length]).toEqual([true, true]);
    const tiny = await cells(12);
    tiny.draw(2);
    expect(tiny.rec.opCount).toBe(0);
  });

  it("an unknown type is drawn as `text`, reported once with svs.render.fallback", async () => {
    const t = await drawn(familyWithTree(40, 12, { columns: true }));
    const issues: string[] = [];
    t.ctx.events.issue.subscribe((i) => issues.push(i.code));
    const m = t.ctx.graph.assembled.value?.view;
    const { k, section } = layerOf(t.ctx, "Members");
    if (!m) throw new Error("no view");
    const cols = m.resolved.layers.params[k].data_columns as Record<string, unknown>[];
    await t.ctx.applyPatch(
      set(m.set.strings.get(m.resolved.layers.ref[k]), {
        data_columns: [{ ...cols[2], type: "org.example:sparkline" }],
      }),
    );
    t.ctx.layout.get(section);
    const { v, s } = view(t.ctx, new RecordingBackend());
    void v;
    expect(s.dataColumns.map((c) => c.renderer)).toEqual(["text"]);
    expect(issues.filter((c) => c === "svs.render.fallback")).toHaveLength(1);
  });

  it("a field provider's values fill a column wherever a field is (org.example:affinity)", async () => {
    const b = createSeqViewSpec({ title: "Affinity" });
    const aln = b.data.alignment({
      ref: "aln",
      source: { value: "# STOCKHOLM 1.0\nA   MKTAYI\nB   MKTAYV\n//\n", format: "stockholm" },
    });
    const meta = b.data.annotation({
      ref: "meta",
      geometry: "rows",
      target: aln,
      value: { row: [0, 1], kd_nm: [1, 3] },
      custom: { "org.example": { unit: "pKd" } },
    } as never);
    b.view({ ref: "v", label: "V" })
      .block({ alignment: aln })
      .track({ label: "Members" })
      .rows({ data_columns: [{ data: meta, field: "org.example:affinity", label: "pKd", type: "number" }] } as never);
    const t = await opened(b.getDocument(), {
      renderers: [...KernelRenderers] as LayerRenderer[],
      registry: [{ svs: StarsFormat }],
    });
    const rec = new RecordingBackend();
    const { v, s } = view(t.ctx, rec);
    v.top = 0;
    v.height = 14;
    v.columnWidth = 64;
    const dc = s.dataColumns[0];
    const draw = t.ctx.registry.get("columns", "number")?.draw;
    draw?.(v, dc.state as never, 0);
    draw?.(v, dc.state as never, 1);
    expect(texts(rec).map((x) => x.text)).toEqual(["9", "8.52"]);
    void AffinityField;
  });
});

describe("the rows' stream key", () => {
  it("two lane orders of the same rows are two keys: a reorder drops the old order's tiles (GV-04)", async () => {
    // 40 rows, three kept (0, 31, 32), ordered by `p` (31, 32, 0) or by `q` (32, 0, 31): a base-31 hash's collision.
    const N = 40;
    const rows = Array.from({ length: N }, (_, i) =>
      Array.from({ length: 12 }, (_, j) => "ACDEFGHIKLMNPQRSTVWY"[(i * 7 + j * (i + 1)) % 20]).join(""),
    );
    const b = createSeqViewSpec({ title: "R" });
    const aln = b.data.alignment({ ref: "aln", value: rows, labels: rows.map((_, i) => `r${i}`) });
    const rank = (order: number[]) => rows.map((_, i) => (order.includes(i) ? order.indexOf(i) : 9));
    const keep = rows.map((_, i) => i === 0 || i === 31 || i === 32);
    b.data.annotation({
      ref: "kt",
      geometry: "rows",
      target: aln,
      value: { row: rows.map((_, i) => i), p: rank([31, 32, 0]), q: rank([32, 0, 31]), keep },
    } as never);
    const filter = [{ data: "kt", field: "keep", values: [true] }];
    b.view({ ref: "v", label: "V", style: "msa" })
      .block({ alignment: aln })
      .track({ label: "Rows" })
      .rows({ ref: "members", order: { data: "kt", field: "p" }, filter } as never);
    const t = await opened(b.getDocument(), { renderers: [...KernelRenderers] as LayerRenderer[] });
    const state = () => layer(t.ctx, "Rows").state as { rowOf: Int32Array | null; source: { key: string } };
    const before = state();
    expect(Array.from(before.rowOf ?? [])).toEqual([31, 32, 0]);
    await t.ctx.applyPatch(set("members", { order: { data: "kt", field: "q" }, filter }));
    await t.ctx.graph.update();
    const after = state();
    expect(Array.from(after.rowOf ?? [])).toEqual([32, 0, 31]);
    expect(after.source.key).not.toBe(before.source.key);
    t.ctx.dispose();
  });
});

describe("gutters and lazy rows", () => {
  it("a row's first and last residue numbers sit in the gutters", async () => {
    const t = await drawn(familyWithTree(40, 12), { width: 1000, height: 420 });
    const rec = new RecordingBackend();
    const { v, s } = view(t.ctx, rec);
    v.top = 100;
    v.height = 14;
    v.dataRow = 3;
    v.gutterLeft = 56;
    RowsRenderer.gutter?.(v, s, {});
    const [left, right] = texts(rec);
    expect([left.align, right.align]).toEqual([TextAlign.right, TextAlign.left]);
    expect([left.x, right.x]).toEqual([Math.round(v.left - 4), Math.round(v.right + 4)]);
    // 30 columns, a row's own numbering: the first is near 1, the last near 30.
    expect(Number(left.text)).toBeLessThan(4);
    expect(Number(right.text)).toBeGreaterThan(26);
    // Centered on the digits' height in the row, as the row's letters: an alphabetic baseline under the middle (107).
    const lines = new Float64Array(2);
    letterLines(v, v.text.font("ui", 10), 100, 114, lines);
    expect([left.y, right.y]).toEqual([lines[0], lines[0]]);
    expect(left.y).toBeGreaterThan(110);
    rec.reset();
    v.height = 8;
    RowsRenderer.gutter?.(v, s, {});
    expect(rec.opCount).toBe(0);
  });

  it("lazy rows draw from tiles as they arrive; the overview's map pending until its tiles do, in any box", async () => {
    const uri = "https://example.org/lazy";
    const b = createSeqViewSpec({ title: "lazy" });
    const aln = b.data.alignment({
      ref: "lazy",
      row_count: 3000,
      column_count: 400,
      source: { uri, protocol: "seqstar-tiles/1" },
    });
    b.view({ ref: "v", label: "V", style: "msa", height: "box" })
      .block({ alignment: aln })
      .track({ label: "Rows" })
      .rows();
    const t = await opened(b.getDocument(), {
      renderers: [...KernelRenderers] as LayerRenderer[],
      loaders: lazyLoaders(uri, 3000, 400),
      size: { "v.main": [1000, 400] },
    });
    t.vp("v.main").setView(1.5, 0, "api");
    await t.ctx.graph.update();
    // The vertical overview's map (`MapView`) in the column's own box and in a navigator's band: the level the region
    // cache loads for the section's overview column either way (`overviewLevels`), so it fills once its tiles arrive.
    const { v: v0, s, k } = view(t.ctx, new RecordingBackend(), "Rows");
    const m = t.ctx.graph.assembled.value?.view;
    const av = t.ctx.layout.axis("v.main");
    if (!m || !av?.axis) throw new Error("no view");
    expect(av.vnavWidth).toBeGreaterThan(0);
    const map = new MapView();
    const cache = RowsRenderer.cache?.create(s, v0.skin) ?? null;
    const overview = (w: number, h: number) => {
      const rec = new RecordingBackend();
      const v = map.of(av, new Float64Array(1), 0, w / 400, w);
      Object.assign(v, { backend: rec, layer: k, text: v0.text, props: v0.props, skin: v0.skin, cache, top: 0 });
      v.height = h;
      v.setThemes(v0.themes);
      v.stream = t.ctx.streams.get(m.set.strings.get(m.resolved.layers.ref[k]));
      expect(RowsRenderer.overview?.(v, s, {})).toBe(true);
      return [count(rec, Op.image) > 0, rects(rec).some((r) => r[4] === SkinToken.pending)];
    };
    const boxes = [
      [av.vnavWidth, av.regionBottom - av.regionTop],
      [700, 28],
    ];
    const before = frame(t.ctx, "Rows");
    // Nothing arrived: no image, nothing drawn but what the draw has of its own; the maps pending.
    expect(count(before.rec, Op.image)).toBe(0);
    for (const [w, h] of boxes) expect(overview(w, h), `${w} × ${h}`).toEqual([false, true]);
    await settle(t.ctx);
    const after = frame(t.ctx, "Rows");
    expect(count(after.rec, Op.image)).toBeGreaterThan(0);
    for (const [w, h] of boxes) expect(overview(w, h), `${w} × ${h}`).toEqual([true, false]);
    // No overview column: the map at `MAP_WIDTH`'s levels (ST3b; the manager loads them where a miniature names the
    // track), pending here since nothing asked them.
    t.ctx.layout.axis("v.main");
    const v = map.of(av, new Float64Array(1), 0, 1, 400);
    v.vnavWidth = 0;
    v.cache = cache;
    expect(RowsRenderer.overview?.(v, s, {})).toBe(true);
  });

  it("the vertical overview's image keys on its rows layer's regions: another layer's arrival does not rebuild it", async () => {
    const uri = "https://example.org/lazy";
    const b = createSeqViewSpec({ title: "lazy" });
    const aln = b.data.alignment({
      ref: "lazy",
      row_count: 3000,
      column_count: 400,
      source: { uri, protocol: "seqstar-tiles/1" },
    });
    b.view({ ref: "v", label: "V", style: "msa", height: "box" })
      .block({ alignment: aln })
      .track({ label: "Rows" })
      .rows();
    const t = await opened(b.getDocument(), {
      renderers: [...KernelRenderers] as LayerRenderer[],
      loaders: lazyLoaders(uri, 3000, 400),
      size: { "v.main": [1000, 400] },
    });
    await t.ctx.graph.update();
    const { graph, selection, highlight, systems, layout } = t.ctx;
    const marks = new ChannelMarks({ graph, hover: t.ctx.hover, selection, highlight, systems, layout });
    const vnav = new VerticalNav(t.ctx, "v.main", marks);
    const rec = new RecordingBackend();
    const builds = () => {
      rec.reset();
      vnav.frame(rec, Dirty.nav);
      return (vnav as unknown as { key: { builds: number } }).key.builds;
    };
    expect(builds()).toBe(1);
    expect(builds()).toBe(1);
    // Another layer's region arriving bumps the manager's version, not this layer's.
    const version = t.ctx.streams.version as unknown as { value: number; set(v: number): void };
    version.set(version.value + 1);
    expect(builds()).toBe(1);
    // Its own regions arriving rebuild it.
    await settle(t.ctx);
    expect(builds()).toBeGreaterThan(1);
    vnav.dispose();
    marks.dispose();
    t.ctx.dispose();
  });

  it("10,000,000 lazy rows: the layout, the region cache's passes and frames hold nothing per row", async () => {
    // A transport that serves one tile and labels at once: the procedural tile server (the format's tests) reduces every
    // cell of a sample and keeps every row's string, so at a million rows its overview tile alone runs out of memory.
    const uri = "https://example.org/huge";
    const tile = new Uint8Array(2 * TILE * TILE).fill(65, 0, TILE * TILE).fill(200, TILE * TILE);
    const calls: string[] = [];
    const transport = async (path: string) => {
      calls.push(path);
      const [kind, a, b] = path.split("/");
      if (kind === "tiles") return tile;
      if (kind !== "labels") throw new NotServed(`404 ${path}`);
      const labels = Array.from({ length: Number(b) - Number(a) }, (_, i) => `S${Number(a) + i}`);
      return new TextEncoder().encode(JSON.stringify(labels));
    };
    const N = 10_000_000;
    const before = process.memoryUsage();
    const b = createSeqViewSpec({ title: "huge" });
    const aln = b.data.alignment({
      ref: "huge",
      row_count: N,
      column_count: 5000,
      source: { uri, protocol: "seqstar-tiles/1" },
    });
    b.view({ ref: "v", label: "V", style: "msa", height: "box" })
      .block({ alignment: aln })
      .track({ label: "Rows" })
      .rows();
    const t = await opened(b.getDocument(), {
      renderers: [...KernelRenderers] as LayerRenderer[],
      loaders: {
        load: (r) => Promise.reject(new Error(`no fixture for ${r.uri}`)),
        stream: tilesStream((u) => (u === uri ? (transport as never) : null)),
      },
      size: { "v.main": [1000, 400] },
    });
    for (const scale of [0.0001, 0.01, 1]) {
      t.vp("v.main").setView(1.5, 0, "api");
      t.vp("v.main").setRowScale(scale, 0);
      await t.ctx.graph.update();
      await settle(t.ctx);
      const { rec, n } = frame(t.ctx, "Rows");
      expect([n, count(rec, Op.image) > 0], `row scale ${scale}`).toEqual([1, true]);
    }
    const after = process.memoryUsage();
    // A 4-byte array per row is 40 MB (the alignment's parent rows); anything per row larger shows here.
    const grown = after.heapUsed + after.arrayBuffers - before.heapUsed - before.arrayBuffers;
    expect(grown / N).toBeLessThan(16);
    expect(calls.length).toBeLessThan(400);
  });
});
