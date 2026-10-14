/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// RN3's hot-path frames for the heap guards (rn3.heap.ts), in a module of their own: the guard imports it afresh per
// case (`vi.resetModules`), so each frame's code meets one copy of the viewer's classes. Frames built in the test
// module are shared by the cases, and a load there polymorphic over two copies' double fields boxes (RN1).

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { TILE } from "@molstar/seq-core/format";
import { AxisGeometry, SectionAxisView } from "../src/axis.js";
import * as api from "../src/draw/api.js";
import { Dirty } from "../src/draw/index.js";
import * as nav from "../src/draw/nav.js";
import { navigatorFrom } from "../src/index.js";
import { SectionMarks } from "../src/marks.js";
import { RowsRenderer } from "../src/renderers/rows.js";
import * as tiles from "../src/renderers/rows-tiles.js";
import { lanesInView, type Streams } from "../src/streams.js";
import { fixtureMeasurer, RecordingBackend } from "../src/testing/index.js";
import { TextManager } from "../src/text.js";
import { axisOf, resolved } from "./axis-fixtures.js";
import { COPY_BACK, composed, recording } from "./composer-fixtures.js";
import { drawn, layerOf } from "./draw-fixtures.js";
import { familyWithTree, syntheticFamily, wideFamily } from "./families.js";

/** A pan's frame: the axis view synced (its window mapped), the navigator's axis and window over 24 blocks. */
export async function navigatorFrame(stacked: boolean, labels: boolean) {
  const b = createSeqViewSpec();
  const view = b.view({ ref: "v", label: "V", arrangement: stacked ? "stacked" : "side" });
  for (let i = 0; i < 24; i++)
    view.block({ sequence: b.data.sequence({ ref: `c${i}`, value: "M".repeat(200 + i * 40) }) });
  view.track({ ref: "v.l" }).letters();
  const r = await resolved(b.getDocument(), "v");
  const ax = axisOf(r);
  const av = new SectionAxisView();
  const L = { geometry: AxisGeometry.linear, dataLeft: 40, dataRight: 1040, scrolledTop: 20, height: 400 };
  const gaps = ax.stacked ? 0 : (ax.count - 1) * ax.gap;
  const minCell = (1000 - gaps) / ax.cells;
  const s = { cell: minCell * 20, minCell, scrollX: 0, scrollY: 0, dpr: 2, lod: 0, contentHeight: 0 };
  const st = { ...s, visLo: new Int32Array(ax.count), visHi: new Int32Array(ax.count) };
  av.sync(ax, L, st as never);
  const marks = new SectionMarks();
  for (const set of [marks.selection, marks.highlight, marks.accent, marks.hover]) {
    set.reset(ax.count);
    for (let bl = 0; bl < ax.count; bl += 3) {
      set.begin(bl);
      for (let k = 0; k < 10; k++) set.add(k * 17, k * 17 + 5);
      set.close(0, ax.shown[bl]);
    }
  }
  marks.updatePrimary(av, ax.count);
  const parts = navigatorFrom({
    lanes: stacked,
    height: 36,
    thickness: 28,
    axis: labels ? ["labels"] : ["hairline"],
  });
  const n = new nav.NavView().from(av);
  const rec = new RecordingBackend();
  n.backend = rec;
  n.text = new TextManager(fixtureMeasurer(), () => 20_000);
  n.parts = parts.parts;
  n.names = ax.names;
  n.marks = marks;
  n.height = 36;
  n.gripBottom = 36;
  n.dpr = 2;
  let i = 0;
  return () => {
    // A pan: the axis view synced (its window mapped), then the axis and the window drawn.
    st.scrollX = (i++ % 500) * 7.5;
    av.sync(ax, L, st as never);
    n.x0 = av.navX0;
    n.x1 = av.navX1;
    n.hot = (i & 1) === 0;
    n.lit = i % 3;
    // Now and then idle (nowhere to go): the quiet draw, with no grips.
    n.idle = (i & 7) === 0;
    rec.reset();
    nav.drawNavigator(n);
    nav.drawNavigatorWindow(n);
  };
}

/** A frame of the family's rows: dense lanes, single rows, the overview's map, coarse under fine when `remote`. */
export async function rowsFrame(remote: boolean) {
  const t = await drawn(syntheticFamily(600, 40), { width: 1000, height: 600 });
  const { k, section } = layerOf(t.ctx, "Members");
  const p = t.ctx.graph.assembled.value?.view?.prepared[k];
  const s0 = p?.state as api.RowsState;
  const s = remote ? ({ ...s0, codes: null, slot: null } as typeof s0) : s0;
  const v = new api.DrawView().from(t.ctx.layout.axis(section) as never);
  const rec = new RecordingBackend();
  v.backend = rec;
  v.layer = k;
  v.setThemes(p?.themes ?? []);
  const cls = new Uint8Array(TILE * TILE).fill(65);
  const tile = { cls, share: new Uint8Array(TILE * TILE).fill(200) };
  // Remote: the fine level's tile 0 arrived, tile 1 not (coarse under it); the overview's level pending past tile 0.
  v.stream = {
    source: s.source,
    version: 0,
    budget: 0,
    tile: (ly, _lx, ty) => (ly === 0 && ty > 0 ? null : tile),
    state: (ly, _lx, ty) => (ly === 0 && ty > 0 ? 0 : 1),
    label: () => null,
  };
  const c = new tiles.RowTiles();
  const lv = new Int32Array([0, 0, 2, 0]);
  return () => {
    rec.reset();
    tiles.beginTiles(v, c);
    // Dense lanes in one call, then single rows with their inset, then the overview's map.
    tiles.rowGeo[0] = 30;
    tiles.rowGeo[1] = 0.8;
    tiles.rowGeo[2] = 0;
    tiles.paintLevels(v, s, c, 0, lv, 0, 600, v.visLo[0], v.visHi[0]);
    tiles.rowGeo[1] = 14;
    tiles.rowGeo[2] = 1;
    for (let r = 250; r < 262; r++) {
      tiles.rowGeo[0] = 30 + (r - 250) * 14;
      tiles.paintLevels(v, s, c, 0, lv, r, 1, v.visLo[0], v.visHi[0]);
    }
    tiles.rowGeo[0] = 0;
    tiles.rowGeo[1] = 0.5;
    tiles.rowGeo[2] = 0;
    tiles.paintTiles(v, s, c, 0, 1, 0, 0, 600, v.visLo[0], v.visHi[0], api.SkinToken.pending);
  };
}

/**
 * A frame of the rows entry's own draws (seqspec:viewer#extension-points): single rows with their gaps and letters, the
 * gutters' numbers, the tree header and a hovered clade, the data columns' cells, or (`dense`) the dense band in one
 * call over reduced tiles; `remote`: the same single rows from the region cache's tiles, none held.
 */
export async function rowsDrawFrame(kind: "single" | "dense" | "remote", parts = 31) {
  const doc = kind === "dense" ? wideFamily(400, 80) : familyWithTree(400, 40, { columns: true });
  const t = await drawn(doc, { width: 1000, height: 600, dpr: 2 });
  const { k, section } = layerOf(t.ctx, "Members");
  const vp = t.ctx.viewports.get(section);
  if (kind === "dense") {
    vp?.setView(0.5, 0, "api");
    vp?.setRowScale(0.04, 0);
  } else vp?.setView(16, 0, "api");
  await t.ctx.graph.update();
  for (let i = 0; i < 4; i++) {
    (t.ctx.streams as Streams).pass();
    await new Promise((ok) => setTimeout(ok, 0));
  }
  const m = t.ctx.graph.assembled.value?.view;
  const p = m?.prepared[k];
  const s0 = p?.state as api.RowsState;
  const s = kind === "remote" ? ({ ...s0, codes: null, slot: null } as typeof s0) : s0;
  const l = t.ctx.layout.get(section)?.current;
  const v = new api.DrawView().from(t.ctx.layout.axis(section) as never);
  const rec = new RecordingBackend();
  v.backend = rec;
  v.layer = k;
  v.setThemes(p?.themes ?? []);
  v.text = t.ctx.text;
  v.props = t.ctx.layout.props(section);
  v.skin = t.ctx.skin.value;
  v.cache = RowsRenderer.cache?.create(s, v.skin) ?? null;
  const tile = { cls: new Uint8Array(TILE * TILE).fill(65), share: new Uint8Array(TILE * TILE).fill(200) };
  v.stream =
    kind === "remote"
      ? {
          source: s.source,
          version: 0,
          budget: 0,
          tile: () => tile,
          state: () => 1,
          label: () => null,
        }
      : t.ctx.streams.get(m?.set.strings.get(m.resolved.layers.ref[k]) ?? "");
  if (l) lanesInView(l, v, v.lanes);
  const h = l?.rows.height ?? 18;
  const top = (l?.scrolledTop ?? 0) + (l?.entryTop[l.rows.entry] ?? 0);
  const P = {};
  const cols = s.dataColumns.map((c) => t.ctx.registry.get("columns", c.renderer));
  v.height = h;
  v.gutterLeft = 56;
  v.columnWidth = 80;
  return () => {
    rec.reset();
    if (kind === "dense") {
      v.top = top;
      v.row = v.dataRow = -1;
      RowsRenderer.rows?.(v, s, P, 0, 300);
      return;
    }
    for (let r = 0; r < 12; r++) {
      v.row = v.dataRow = r;
      v.top = top + r * h;
      if (parts & 1) RowsRenderer.draw?.(v, s, P);
      if (parts & 2) RowsRenderer.gutter?.(v, s, P);
      for (let c = 0; parts & 4 && c < cols.length; c++) {
        v.columnX = 100 + 90 * c;
        cols[c]?.draw?.(v, s.dataColumns[c].state, r);
      }
    }
    v.top = top;
    if (parts & 8) RowsRenderer.header?.draw?.(v, s, P, 0, 40);
    if (parts & 16) RowsRenderer.header?.hover?.(v, s, P, 2, 0, 40);
  };
}

/**
 * A pan's frames through the composer (every pass): `seqv:map` and a `band` navigator's map over a ranges track, their
 * images kept under their keys; `bare`: no track the miniature draws (letters only), so neither has an image.
 */
export async function mapFrame(bare = false) {
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: "a", value: "MKTAYIAKQR".repeat(60) });
  const hot = b.data.annotation({
    ref: "hot",
    geometry: "ranges",
    target: a,
    value: { start: [10, 120, 300], stop: [60, 180, 320], category: ["a", "b", "a"] },
  } as never);
  const v = b.view({ ref: "v", label: "V" }).block({ sequence: a });
  v.track({ label: "L" }).letters();
  if (!bare) {
    v.track({ label: "Hot" }).ranges({ data: hot } as never);
    v.track({ label: "Map" }).layer({ type: "seqv:map", fallback: "omit" } as never);
  }
  const preset = { name: "lab:band", geometry: "linear", navigator: "band" } as const;
  const t = await composed(b.getDocument(), { width: 1000, height: 300, dpr: 2, registry: [{ presets: [preset] }] });
  t.ctx.ui.setStyle(t.section, preset.name);
  await t.ctx.graph.update();
  const be = recording(COPY_BACK);
  t.vp?.setView(8, 200, "api");
  t.frame(be, Dirty.all);
  let i = 0;
  return {
    composer: t.composer,
    frame: () => {
      // A pan (small ints: a double passed to a call is boxed), every pass drawn.
      t.vp?.panBy((i++ & 1) === 0 ? 9 : -9, 0, "pointer");
      t.frame(be, Dirty.all);
    },
  };
}
