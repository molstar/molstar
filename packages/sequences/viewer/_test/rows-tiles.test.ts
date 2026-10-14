/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The rows layer's tile images (RN3 hot paths; seqspec:viewer#extension-points, seqspec:viewer#the-render-backend,
// seqspec:sources#levels-and-the-zoom-gate), on the recording backend: cells colored by the layer's theme with its
// alpha (case 3's conservation shading), reduced tiles by their share, gaps quiet; images made once and kept (builds
// counted), dropped on a theme change, the least recently drawn evicted; coarse in place of a fine tile not arrived
// (only there), the overview's pending shade; tiles at their whole scale.

import { TILE, type Tile } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import {
  DrawView,
  NO_COLOR,
  type RowsState,
  SkinToken,
  type StreamView,
  type Theme,
  type ThemeLocation,
} from "../src/draw/api.js";
import { paintLevels, paintTiles, RowTiles, rowGeo, tileImage } from "../src/renderers/rows-tiles.js";
import { Op, RecordingBackend } from "../src/testing/index.js";
import { drawn, layerOf } from "./draw-fixtures.js";
import { syntheticFamily } from "./families.js";

/** A theme by letter, its alpha by column: what `rules` with a column alpha gives (case 3). */
const shading: Theme<"color"> = {
  kind: "color",
  granularity: 3,
  reduction: "classes",
  value: (l: ThemeLocation) => (l.letter === 45 ? NO_COLOR : l.letter * 0x010203),
  alpha: (l: ThemeLocation) => (l.column * 8) & 255,
  contextHash: 1,
};

/** The family's rows layer drawn by a `DrawView` from its section's axis view, the theme given. */
async function family(theme: Theme<"color"> = shading) {
  const t = await drawn(syntheticFamily(120, 12));
  const { k, section } = layerOf(t.ctx, "Members");
  const m = t.ctx.graph.assembled.value?.view;
  const s = m?.prepared[k].state as RowsState;
  const v = new DrawView().from(t.ctx.layout.axis(section) as never);
  const rec = new RecordingBackend();
  v.backend = rec;
  v.layer = k;
  v.setThemes([theme]);
  return { t, s, v, rec, c: new RowTiles() };
}
/** An image's pixel (x, y) as [0xRRGGBB, alpha]. */
function pixelOf(rec: RecordingBackend, img: number, x: number, y: number) {
  const p = (rec.images.content[img] as Uint32Array)[y * TILE + x];
  return [((p & 255) << 16) | (p & 0xff00) | ((p >>> 16) & 255), p >>> 24];
}
/** A region cache's view of fixed tiles: `tiles` by "ly/lx/ty/tx", the rest pending. */
function stream(source: RowsState, tiles: Record<string, Tile>): StreamView {
  return {
    source: source.source,
    version: 0,
    budget: 0,
    tile: (ly, lx, ty, tx) => tiles[`${ly}/${lx}/${ty}/${tx}`] ?? null,
    state: (ly, lx, ty, tx) => (tiles[`${ly}/${lx}/${ty}/${tx}`] ? 1 : 0),
    label: () => null,
  };
}
const tile = (cls: number, share: number): Tile => ({
  cls: new Uint8Array(TILE * TILE).fill(cls),
  share: new Uint8Array(TILE * TILE).fill(share),
});

describe("the tile image build", () => {
  it("held rows at the finest level: each cell its theme's color at its theme's alpha", async () => {
    const { s, v, rec, c } = await family();
    const img = tileImage(v, s, c, 0, 0, 0, 0, 0, false);
    expect([img >= 0, c.builds, rec.images.w[img], rec.images.h[img]]).toEqual([true, 1, TILE, TILE]);
    const codes = s.codes as Uint8Array;
    for (const [k, col] of [
      [0, 0],
      [5, 7],
      [119, 29],
    ]) {
      const code = codes[(s.slot as Int32Array)[k] * s.columns + col];
      expect(pixelOf(rec, img, col, k)).toEqual(code === 45 ? [0, 0] : [code * 0x010203, (col * 8) & 255]);
    }
    // Past the alignment's columns and rows: nothing.
    expect([pixelOf(rec, img, 30, 0), pixelOf(rec, img, 0, 120)]).toEqual([
      [0, 0],
      [0, 0],
    ]);
    // Kept: asked again, no build.
    expect([tileImage(v, s, c, 0, 0, 0, 0, 0, false), c.builds]).toEqual([img, 1]);
  });

  it("a reduced tile: opaque by its share times the theme's alpha; gap samples a quiet gray at half; -1 when absent", async () => {
    const { s, v, rec, c } = await family();
    v.stream = stream(s, { "1/1/0/0": tile(65, 200), "1/2/0/0": tile(45, 200) });
    const img = tileImage(v, s, c, 0, 1, 1, 0, 0, false);
    // Column of sample 3 at lx 1: 6, alpha 48; 200 × 48 / 255.
    expect(pixelOf(rec, img, 3, 0)).toEqual([65 * 0x010203, Math.round((200 * 48) / 255)]);
    expect(pixelOf(rec, tileImage(v, s, c, 0, 1, 2, 0, 0, false), 0, 0)).toEqual([0x9aa0a8, 100]);
    expect(tileImage(v, s, c, 0, 2, 2, 0, 0, false)).toBe(-1);
    expect(c.builds).toBe(2);
  });

  it("a new theme or another backend's images drops them; past the slots, the least recently drawn goes", async () => {
    const { s, v, rec, c } = await family();
    const { beginTiles } = await import("../src/renderers/rows-tiles.js");
    beginTiles(v, c);
    const img = tileImage(v, s, c, 0, 0, 0, 0, 0, false);
    beginTiles(v, c);
    tileImage(v, s, c, 0, 0, 0, 0, 0, false);
    expect(c.builds).toBe(1);
    v.setThemes([{ ...shading, contextHash: 2 }]);
    beginTiles(v, c);
    tileImage(v, s, c, 0, 0, 0, 0, 0, false);
    expect([c.builds, rec.images.content[img]]).toEqual([2, null]);
    // 320 slots: a wide source's tiles drawn over frames; the oldest frame's image is the one disposed.
    const wide = { ...s, codes: null, slot: null, columns: TILE * 400 } as RowsState;
    const t = tile(65, 255);
    v.stream = { ...stream(wide, {}), tile: () => t, state: () => 1 };
    const ids: number[] = [];
    for (let tx = 0; tx < 330; tx++) {
      beginTiles(v, c);
      ids.push(tileImage(v, wide, c, 0, 1, 0, 0, tx, false));
    }
    expect(c.used).toBe(320);
    expect(rec.images.content[ids[0]]).toBeNull();
    expect(rec.images.content[ids[329]]).not.toBeNull();
  });
});

/** The image ops of a recording: image id, its rect (sx sy sw sh dx dy dw dh) and what of it shows (the clip's). */
function imageParts(rec: RecordingBackend) {
  const out: { img: number; r: number[]; shown: number[] }[] = [];
  const clips: number[][] = [];
  for (let k = 0; k < rec.opCount; k++) {
    const a = rec.ops[4 * k + 3];
    const code = rec.ops[4 * k];
    if (code === Op.pushClip) clips.push(Array.from(rec.args.subarray(a, a + 4)));
    if (code === Op.popClip) clips.pop();
    if (code !== Op.image) continue;
    const r = Array.from(rec.args.subarray(a, a + 8));
    const [x0, y0, x1, y1] = [r[4], r[5], r[4] + r[6], r[5] + r[7]];
    const c = clips[clips.length - 1] ?? [x0, y0, x1 - x0, y1 - y0];
    const shown = [Math.max(x0, c[0]), Math.max(y0, c[1]), Math.min(x1, c[0] + c[2]), Math.min(y1, c[1] + c[3])];
    out.push({ img: rec.ops[4 * k + 1], r, shown });
  }
  return out;
}

describe("painting the rows", () => {
  it("a tile in parts at its whole scale, odd sample counts placed through xd from its first column; clipped to the view", async () => {
    const { s, v, rec, c } = await family();
    const { beginTiles } = await import("../src/renderers/rows-tiles.js");
    rowGeo[0] = 50;
    rowGeo[1] = 12;
    rowGeo[2] = 1;
    const lv = new Int32Array([0, 0, 0, 0]);
    for (let f = 0; f < 2; f++) {
      rec.reset();
      beginTiles(v, c);
      paintLevels(v, s, c, 0, lv, 3, 1, v.visLo[0], v.visHi[0]);
    }
    expect(c.builds).toBe(1);
    const [lo, hi] = [v.visLo[0], v.visHi[0]];
    const at = (u: number) => v.xd(0, lo - v.indexOf(0, lo as never) + u);
    const parts = imageParts(rec);
    expect(parts.length).toBeGreaterThan(0);
    for (const [k, { r }] of parts.entries()) {
      // Row 3's sample row over the row inset by 1 px; odd sample counts, consecutive, never smoothed.
      expect([r[1], r[3], r[5], r[7], r[2] % 2]).toEqual([3, 1, 51, 10, 1]);
      if (k > 0) expect(r[0]).toBe(parts[k - 1].r[0] + parts[k - 1].r[2]);
      expect([r[4], r[4] + r[6]]).toEqual([at(r[0]), at(r[0] + r[2])]);
      expect(r[6]).toBeLessThanOrEqual(512 + at(1) - at(0));
    }
    // What shows: the columns in view, no more.
    expect([parts[0].shown[0], parts[parts.length - 1].shown[2]]).toEqual([v.xd(0, lo), v.xd(0, hi)]);
  });

  it("coarse in place of a fine tile not arrived; the overview's pending shade where none has", async () => {
    const { s, v, rec, c } = await family();
    const remote = { ...s, codes: null, slot: null } as RowsState;
    v.stream = stream(remote, { "2/2/0/0": tile(65, 255) });
    rowGeo[0] = 0;
    rowGeo[1] = 1;
    rowGeo[2] = 0;
    paintLevels(v, remote, c, 0, new Int32Array([0, 0, 2, 2]), 0, 120, v.visLo[0], v.visHi[0]);
    // The fine level (0, 0) has nothing: the coarse tile drawn alone (in parts).
    const parts = imageParts(rec);
    expect([parts.length > 0, new Set(parts.map((p) => p.img)).size]).toEqual([true, 1]);
    rec.reset();
    v.stream = stream(remote, {});
    expect(paintTiles(v, remote, c, 0, 2, 2, 0, 120, v.visLo[0], v.visHi[0], SkinToken.pending)).toBe(1);
    expect([rec.count(Op.image), rec.rectItems.n, rec.rectItems.color[0]]).toEqual([0, 1, SkinToken.pending]);
  });

  it("the coarse level only in place of a fine tile not arrived, clipped to its part", async () => {
    const { s, v, rec, c } = await family();
    const remote = { ...s, codes: null, slot: null } as RowsState;
    // Lanes 0 to 399 at 1 px: fine tile row 0 arrived, row 1 not; the coarse level's one tile holds both.
    v.stream = stream(remote, { "0/0/0/0": tile(65, 255), "2/0/0/0": tile(66, 255) });
    rowGeo[0] = 10;
    rowGeo[1] = 1;
    rowGeo[2] = 0;
    expect(paintTiles(v, remote, c, 0, 0, 0, 0, 400, v.visLo[0], v.visHi[0], NO_COLOR, (2 << 4) | 0)).toBe(1);
    const parts = imageParts(rec);
    const ids = [...new Set(parts.map((p) => p.img))];
    // Rows shown per image: the fine tile its 256 lanes, the coarse one the missing tile's 144 only.
    const rows = (img: number) => {
      const ys = parts.filter((p) => p.img === img && p.shown[3] > p.shown[1]).map((p) => p.shown);
      return [Math.min(...ys.map((y) => y[1])), Math.max(...ys.map((y) => y[3]))];
    };
    expect(ids.map(rows)).toEqual([
      [10, 266],
      [266, 410],
    ]);
    const [lo, hi] = [v.visLo[0], v.visHi[0]];
    for (const p of parts) expect([p.shown[0] >= v.xd(0, lo), p.shown[2] <= v.xd(0, hi)]).toEqual([true, true]);
  });
});
