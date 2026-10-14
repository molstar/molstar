/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The fingerprint's draws on the recording backend (seqspec:viewer#extension-points): the rects against an oracle
// written in the old doubles' terms (`Math.round` rows and stripes), at DPR 1, 2 and 3 and several zooms; a stripe's
// center hits its row as `state` and its kind as `part`; a row's gap hits the row above, outside the rows none; 20
// poses × 8 kinds keep every stripe; a short track goes compact with no kind missing from drawing or hit; grounds are
// drawn before stripes; below a device pixel the cells merge by device column and cover what the cells would; the
// blocks without the property are the `matrix` renderer's tiles.

import { describe, expect, it } from "vitest";
import { case6 } from "../../../core/_test/format/cases.js";
import { SkinToken } from "../../src/draw/backend.js";
import { Hit, KernelRenderers } from "../../src/index.js";
import { Kind, KindStyle, Kinds } from "../../src/interactions/index.js";
import {
  FingerprintRenderer,
  type FingerprintState,
  InteractionsViewerRegistry,
} from "../../src/interactions/viewer/index.js";
import { LayerFrame, Op, type RecordedImages, RecordingBackend } from "../../src/testing/index.js";
import { opened } from "../layout-fixtures.js";
import {
  absentScene,
  layerOf,
  manyScene,
  mixedScene,
  open,
  type Scene,
  ShortFingerprint,
  unitScene,
  wideScene,
} from "./fixtures.js";

const COLOR = Int32Array.from(Kinds, (k) => Number.parseInt(KindStyle[k][1].slice(1), 16));
const kindOfColor = (c: number) => COLOR.indexOf(c);
const GROUND = SkinToken.grid;
const CASE6: Scene = { doc: case6(), section: "poses-view.main", track: "poses-view.fingerprint" };

type Opened = Awaited<ReturnType<typeof open>>;
interface Rect {
  x: number;
  y: number;
  w: number;
  h: number;
  color: number;
}

/** The scene drawn at `cell` (null: the fit) and `dpr`: the entry's box, its rects, and the recording. */
async function draw(
  scene: Scene,
  o: { dpr?: number; cell?: number | null; registry?: Parameters<typeof open>[1] } = {},
) {
  const dpr = o.dpr ?? 1;
  const t: Opened =
    scene === CASE6
      ? Object.assign(
          await opened(CASE6.doc, {
            renderers: [...KernelRenderers],
            view: "poses-view",
            registry: [InteractionsViewerRegistry],
            render: { maxCell: 400 },
          }),
          {},
        )
      : await open(scene, { dpr, size: [1000, 600], render: { maxCell: 400 }, ...o.registry });
  if (scene === CASE6) {
    t.ctx.viewports.get(CASE6.section)?.setSize(1000, 600, dpr);
    await t.ctx.graph.update();
  }
  const vp = t.ctx.viewports.get(scene.section);
  if (o.cell !== undefined && o.cell !== null) vp?.setView(o.cell, 0, "api");
  await t.ctx.graph.update();
  const { k, state } = layerOf(t.ctx);
  const rec = new RecordingBackend();
  const f = new LayerFrame(t.ctx, FingerprintRenderer as never, scene.section, k);
  f.draw(rec);
  const l = t.layout(scene.section);
  const e = Array.from(l.entryTrack.subarray(0, l.entryCount)).indexOf(t.trackIndex(scene.track));
  const top = e < l.pinnedCount ? l.entryTop[e] : l.scrolledTop + l.entryTop[e];
  const height = l.entryHeight[e];
  const b = rec.rectItems;
  const rects: Rect[] = Array.from({ length: b.n }, (_, i) => ({
    x: b.x[i],
    y: b.y[i],
    w: b.w[i],
    h: b.h[i],
    color: b.color[i],
  }));
  return {
    t,
    k,
    s: state as FingerprintState,
    rec,
    rects,
    top,
    height,
    dpr,
    ax: t.ctx.layout.axis(scene.section),
    scene,
  };
}
type Drawn = Awaited<ReturnType<typeof draw>>;

/** The old rows' and stripes' edges: `Math.round` of the doubles, in device pixels. */
function oracle(top: number, height: number, dpr: number, S: number, K: number) {
  const a = Math.round((top + 1) * dpr);
  const H = Math.round((top + height - 1) * dpr) - a;
  const G = H >= 6 * S ? 2 : 0;
  const rows = Array.from({ length: S }, (_, q) => {
    const y0 = a + Math.round((q * (H + G)) / S);
    const y1 = a + Math.round(((q + 1) * (H + G)) / S) - G;
    return {
      y0,
      y1,
      stripes: Array.from({ length: K + 1 }, (_, j) => y0 + Math.round((j * (y1 - y0)) / Math.max(1, K))),
    };
  });
  return { a, H, G, rows, compact: Math.floor((H + G) / S) - G < 2 * K };
}

/** Each contacted position's device x extent [x0, max(x1, x0 + 1)) in block 0, the ones in view (a thin view: all). */
function extents(d: Drawn, b = 0) {
  const base = d.ax?.axis?.base[b] ?? 0;
  const thin = (d.ax?.cell ?? 1) * (d.ax?.span[b] ?? 1) * d.dpr < 1;
  const shown = Array.from(d.s.contacted[b], (j) => ({ j, k: d.ax?.shownOf(b, base + j) ?? -1 }));
  return shown
    .filter(({ k }) => k >= 0 && (thin || (k >= (d.ax?.visLo[b] ?? 0) && k < (d.ax?.visHi[b] ?? 0))))
    .map(({ j, k }) => {
      const x0 = Math.round((d.ax?.x(b, k as never) ?? 0) * d.dpr);
      const x1 = Math.round((d.ax?.x(b, (k + 1) as never) ?? 0) * d.dpr);
      return { j, col: d.s.colOf[b][j], x0, x1: Math.max(x1, x0 + 1) };
    });
}
const maskAt = (s: FingerprintState, b: number, q: number, col: number) => s.cells[b][q * s.width[b] + col];

/** The hit at a device point of the drawn entry. */
function hitAt(d: Drawn, x: number, y: number, hit = Hit.create()) {
  d.t.ctx.layout.hitTest(d.scene.section, x / d.dpr, y / d.dpr, hit);
  return hit;
}

describe("case 6's rects against the old rows and stripes, and hits on every stripe", () => {
  it.each([
    [1, null],
    [2, null],
    [3, 60],
    [2, 150],
  ])("DPR %i at cell %s: a ground and a stripe per kind held, on the rounded rows and stripes", async (dpr, cell) => {
    const d = await draw(CASE6, { dpr, cell });
    const S = d.s.matrix.states;
    const K = d.s.kinds.length;
    const o = oracle(d.top, d.height, dpr, S, K);
    expect(o.compact).toBe(false);
    const want: Rect[] = [];
    for (const { col, x0, x1 } of extents(d))
      for (let q = 0; q < S; q++) {
        const m = maskAt(d.s, 0, q, col);
        if (m === 0) continue;
        const inset = x1 - x0 >= 4 ? 1 : 0;
        const w = Math.max(1, x1 - x0 - 2 * inset);
        const r = o.rows[q];
        want.push({ x: x0 + inset, y: r.y0, w, h: r.y1 - r.y0, color: GROUND });
        for (let j = 0; j < K; j++)
          if (m & (1 << d.s.kinds[j]))
            want.push({
              x: x0 + inset,
              y: r.stripes[j],
              w,
              h: r.stripes[j + 1] - r.stripes[j],
              color: COLOR[d.s.kinds[j]],
            });
      }
    const key = (r: Rect) => `${r.color},${r.x},${r.y},${r.w},${r.h}`;
    expect(d.rects.map(key).sort()).toEqual(want.map(key).sort());
    expect(want.length).toBeGreaterThan(40);
  });

  it("grounds are drawn before every stripe, whatever a backend's color grouping", async () => {
    const d = await draw(CASE6);
    const colors = d.rects.map((r) => r.color);
    const lastGround = colors.lastIndexOf(GROUND);
    const firstStripe = colors.findIndex((c) => c !== GROUND);
    expect(lastGround).toBeGreaterThanOrEqual(0);
    expect(lastGround).toBeLessThan(firstStripe);
    // Two rects ops: the grounds', then the stripes'.
    const ops = Array.from({ length: d.rec.opCount }, (_, i) => d.rec.ops[4 * i]).filter((c) => c === Op.rects);
    expect(ops).toHaveLength(2);
  });

  it.each([1, 2, 3])("DPR %i: every stripe's center hits its row as `state` and its kind as `part`", async (dpr) => {
    const d = await draw(CASE6, { dpr });
    const o = oracle(d.top, d.height, dpr, d.s.matrix.states, d.s.kinds.length);
    const hit = Hit.create();
    let stripes = 0;
    for (const r of d.rects) {
      if (r.color === GROUND) continue;
      hitAt(d, r.x + r.w / 2, r.y + r.h / 2, hit);
      const q = o.rows.findIndex((row) => r.y >= row.y0 && r.y < row.y1);
      expect([hit.state, hit.part]).toEqual([q, kindOfColor(r.color)]);
      stripes++;
    }
    expect(stripes).toBeGreaterThan(20);
    // A ground (the cell's, where a stripe of a kind it lacks would be) hits the row and no kind.
    const g = d.rects.find((r) => r.color === GROUND);
    if (!g) throw new Error("no ground");
    for (let y = g.y; y < g.y + g.h; y++) {
      hitAt(d, g.x + g.w / 2, y + 0.5, hit);
      expect(hit.state).toBe(o.rows.findIndex((row) => g.y >= row.y0 && g.y < row.y1));
    }
  });

  it("a row's gap hits the row above; outside the rows nothing", async () => {
    const d = await draw(CASE6, { dpr: 2 });
    const o = oracle(d.top, d.height, 2, 5, 4);
    expect(o.G).toBe(2);
    const x = (extents(d)[3].x0 + extents(d)[3].x1) / 2;
    const hit = Hit.create();
    for (let q = 0; q < 4; q++) {
      for (let y = o.rows[q].y1; y < o.rows[q + 1].y0; y++) {
        hitAt(d, x, y + 0.5, hit);
        expect([y, hit.state]).toEqual([y, q]);
        expect(hit.part).toBe(-1);
      }
    }
    for (const y of [o.a - 1.5, o.a + o.H + 0.5, 0.5, d.top * 2 + 0.5]) {
      hitAt(d, x, y, hit);
      expect(hit.state).toBe(-1);
    }
    // The last row ends where the track's bottom border begins.
    expect(o.rows[4].y1).toBe(o.a + o.H);
  });
});

describe("the hit at every device pixel row equals what the draw painted there (DPR 1, 1.5 and 2)", () => {
  // A brute force over the recorded rects: per contacted cell in view and device row, the painted stripe (the last rect
  // covering the pixel) is the hit's `part`, the row (its gap included) its `state`, outside the rows none. Points at the
  // top, middle and bottom of each device row, so the hit's one-step corrections meet every rounding edge.
  const SCENES = [
    ["case 6", () => CASE6],
    ["20 × 8", () => manyScene()],
    ["compact", () => unitScene({ hint: 1 })],
  ] as const;
  for (const [name, scene] of SCENES)
    it.each([1, 1.5, 2])(`${name} at DPR %f`, async (dpr) => {
      const d = await draw(scene(), { dpr });
      const S = d.s.matrix.states;
      const o = oracle(d.top, d.height, dpr, S, d.s.kinds.length);
      const stateAt = (y: number) => {
        if (y < o.a || y >= o.a + o.H) return -1;
        let q = S - 1;
        while (q > 0 && o.rows[q].y0 > y) q--;
        return q;
      };
      const hit = Hit.create();
      let checked = 0;
      const parts = new Set<number>();
      // Disagreements collected and compared once: an expect per point was most of the test.
      const wrong: [number, number, number, number, number][] = [];
      for (const { x0, x1 } of extents(d)) {
        const x = (x0 + x1) / 2;
        for (let y = Math.floor(d.top * dpr) - 2; y < Math.ceil((d.top + d.height) * dpr) + 2; y++) {
          let part = -1;
          for (const r of d.rects)
            if (r.color !== GROUND && r.x <= x && x < r.x + r.w && r.y <= y && y < r.y + r.h)
              part = kindOfColor(r.color);
          for (const f of [0.02, 0.5, 0.98]) {
            hitAt(d, x, y + f, hit);
            if (hit.state !== stateAt(y) || hit.part !== part)
              wrong.push([y + f, hit.state, hit.part, stateAt(y), part]);
            checked++;
          }
          parts.add(part);
        }
      }
      // [y, the hit's state and part, the painted state and part] where they differ.
      expect(wrong).toEqual([]);
      expect(checked).toBeGreaterThan(300);
      expect(parts.size).toBeGreaterThan(2);
    });
});

describe("20 poses with 8 kinds keep every stripe", () => {
  it.each([1, 2])("DPR %i: each kind drawn and hit at least once per row that holds it", async (dpr) => {
    const d = await draw(manyScene(), { dpr });
    const S = 20;
    const K = 8;
    const o = oracle(d.top, d.height, dpr, S, K);
    expect(o.compact).toBe(false);
    expect(o.rows[0].y1 - o.rows[0].y0).toBeGreaterThanOrEqual(16 * dpr);
    const hit = Hit.create();
    for (let q = 0; q < S; q++) {
      const row = o.rows[q];
      const drawn = new Set(
        d.rects.filter((r) => r.color !== GROUND && r.y >= row.y0 && r.y < row.y1).map((r) => kindOfColor(r.color)),
      );
      const held = new Set<number>();
      for (const { col } of extents(d)) for (const k of d.s.kinds) if (maskAt(d.s, 0, q, col) & (1 << k)) held.add(k);
      expect(held.size).toBe(8);
      expect(drawn).toEqual(held);
      // Every stripe at least 2 device px: the rects' height.
      for (const r of d.rects)
        if (r.color !== GROUND && r.y >= row.y0 && r.y < row.y1) expect(r.h).toBeGreaterThanOrEqual(2);
    }
    // Hit: every cell's every stripe gives its kind.
    let n = 0;
    for (const { col, x0, x1 } of extents(d))
      for (let q = 0; q < S; q++)
        for (let j = 0; j < K; j++) {
          const m = maskAt(d.s, 0, q, col);
          if (!(m & (1 << d.s.kinds[j]))) continue;
          const y = (o.rows[q].stripes[j] + o.rows[q].stripes[j + 1]) / 2;
          hitAt(d, (x0 + x1) / 2, y, hit);
          expect([hit.state, hit.part]).toEqual([q, d.s.kinds[j]]);
          n++;
        }
    expect(n).toBe(S * K);
  });
});

describe("a short track draws compact: the first kind in the stripes' order, none missing from drawing or hit", () => {
  it.each([1, 2])("DPR %i: one rect over the row in the lead kind's color, no ground", async (dpr) => {
    // 4 poses in 22 px (a hint of one row): under 2 device px a kind.
    const d = await draw(unitScene({ hint: 1 }), { dpr });
    const S = 4;
    const K = d.s.kinds.length;
    const o = oracle(d.top, d.height, dpr, S, K);
    expect(o.compact).toBe(true);
    expect(d.rects.some((r) => r.color === GROUND)).toBe(false);
    const lead = (m: number) => d.s.kinds.find((k) => (m & (1 << k)) !== 0) ?? -1;
    const want: Rect[] = [];
    for (const { col, x0, x1 } of extents(d))
      for (let q = 0; q < S; q++) {
        const m = maskAt(d.s, 0, q, col);
        if (m === 0) continue;
        const inset = 0;
        want.push({
          x: x0 + inset,
          y: o.rows[q].y0,
          w: Math.max(1, x1 - x0),
          h: o.rows[q].y1 - o.rows[q].y0,
          color: COLOR[lead(m)],
        });
      }
    const key = (r: Rect) => `${r.color},${r.x},${r.y},${r.w},${r.h}`;
    expect(d.rects.map(key).sort()).toEqual(want.map(key).sort());
    // The cell with two kinds (pose 2 at residue 20: π stacking and hydrophobic) is π stacking's; the hit says so too.
    const hit = Hit.create();
    for (const { col, x0, x1 } of extents(d))
      for (let q = 0; q < S; q++) {
        const m = maskAt(d.s, 0, q, col);
        const y = (o.rows[q].y0 + o.rows[q].y1) / 2;
        hitAt(d, (x0 + x1) / 2, y, hit);
        expect([hit.state, hit.part]).toEqual([q, m === 0 ? -1 : lead(m)]);
      }
    expect(want.some((r) => r.color === COLOR[Kind.pi_stacking])).toBe(true);
  });

  // `minTrackHeight` is DPR-free (S + 2 CSS px): at any DPR from 1 a row is at least a device px, fractional ones too.
  it.each([1, 1.25, 1.5, 2])(
    "DPR %f: every kind of a many-kinds layer survives in a track of a CSS px a row",
    async (dpr) => {
      const d = await draw(manyScene(20, 8, 1), { dpr });
      expect(d.height).toBe(22);
      expect(new Set(d.rects.map((r) => kindOfColor(r.color)))).toEqual(new Set(Array.from(d.s.kinds)));
      const o = oracle(d.top, d.height, dpr, 20, 8);
      for (const row of o.rows) expect(row.y1 - row.y0).toBeGreaterThanOrEqual(1);
      const hit = Hit.create();
      const seen = new Set<number>();
      for (const { col, x0, x1 } of extents(d))
        for (let q = 0; q < 20; q++) {
          const m = maskAt(d.s, 0, q, col);
          hitAt(d, (x0 + x1) / 2, (o.rows[q].y0 + o.rows[q].y1) / 2, hit);
          if (m !== 0) seen.add(hit.part);
          expect(hit.state).toBe(q);
        }
      expect(seen).toEqual(new Set(Array.from(d.s.kinds)));
    },
  );
});

describe("thin cells merge by device column and cover what the cells would", () => {
  it.each([null, 0.03])(
    "at cell %s of 20,000 positions: the stripes' coverage equals the brute force's",
    async (cell) => {
      const d = await draw(wideScene(), { cell: cell ?? null });
      const o = oracle(d.top, d.height, 1, d.s.matrix.states, d.s.kinds.length);
      const px = d.ax?.cell ?? 1;
      expect(px).toBeLessThan(1);
      // Per row and kind, the device columns covered.
      const covered = (rects: Rect[], q: number, kind: number) => {
        const out = new Set<number>();
        const [y0, y1] = [o.rows[q].stripes[d.s.kinds.indexOf(kind)], o.rows[q].stripes[d.s.kinds.indexOf(kind) + 1]];
        for (const r of rects)
          if (r.color === COLOR[kind] && r.y === y0 && r.y + r.h === y1)
            for (let x = r.x; x < r.x + r.w; x++) out.add(x);
        return out;
      };
      const want = new Map<string, Set<number>>();
      for (const { col, x0, x1 } of extents(d))
        for (let q = 0; q < d.s.matrix.states; q++)
          for (const kind of d.s.kinds)
            if (maskAt(d.s, 0, q, col) & (1 << kind)) {
              const set = want.get(`${q}:${kind}`) ?? new Set<number>();
              for (let x = x0; x < x1; x++) set.add(x);
              want.set(`${q}:${kind}`, set);
            }
      let compared = 0;
      for (const [key, cols] of want) {
        const [q, kind] = key.split(":").map(Number);
        expect([...covered(d.rects, q, kind)].sort((a, b) => a - b)).toEqual([...cols].sort((a, b) => a - b));
        compared++;
      }
      expect(compared).toBeGreaterThan(10);
      // Bounded by device columns: no more rects than a ground and a stripe per kind per row per column.
      expect(d.rects.length).toBeLessThanOrEqual(1000 * d.s.matrix.states * (1 + d.s.kinds.length));
      // Fewer rects than cells: the merge happened.
      const cells = extents(d).length * d.s.matrix.states;
      expect(d.rects.filter((r) => r.color === GROUND).length).toBeLessThanOrEqual(cells);
    },
  );
});

describe("the blocks without the property are the `matrix` renderer's tiles", () => {
  /** The pixels of the images a frame drew. */
  const images = (rec: RecordingBackend) => {
    const out: Uint32Array[] = [];
    for (let k = 0; k < rec.opCount; k++)
      if (rec.ops[4 * k] === Op.image) {
        const c = (rec.images as RecordedImages).content[rec.ops[4 * k + 1]];
        if (c instanceof Uint32Array) out.push(c);
      }
    return out;
  };

  it("the property absent: tiles equal a kernel `matrix` layer's over the same data, in the counts' Blues", async () => {
    const mine = await draw(absentScene());
    const doc = structuredClone(absentScene().doc);
    const layer = JSON.stringify(doc).includes("interactions:fingerprint");
    expect(layer).toBe(true);
    const text = JSON.stringify(doc).replace('"type":"interactions:fingerprint"', '"type":"matrix"');
    const t = await opened(JSON.parse(text), { renderers: [...KernelRenderers], view: "v" });
    const k = t.ctx.graph.assembled.value?.view?.prepared.findIndex((p) => p.renderer === "matrix") ?? -1;
    const rec = new RecordingBackend();
    const entry = t.ctx.registry.get("renderers", "matrix");
    if (!entry) throw new Error("no matrix");
    new LayerFrame(t.ctx, entry, "v.main", k).draw(rec);
    const a = images(mine.rec);
    const b = images(rec);
    expect(a.length).toBeGreaterThan(0);
    expect(a.map((x) => Array.from(x))).toEqual(b.map((x) => Array.from(x)));
    // Nothing else of its own: no ground, no stripe.
    expect(mine.rects).toEqual([]);
  });

  it("mixed data: rows of each fingerprint, bounded by its own rows, the matrix block as an image beside them", async () => {
    const d = await draw(mixedScene());
    const o = oracle(d.top, d.height, 1, 5, d.s.kinds.length);
    // Block 1 is cut to three rows: no rect of it below row 3's top.
    const ax = d.ax;
    const x1 = Math.round((ax?.x(1, 0 as never) ?? 0) * d.dpr);
    const x2 = Math.round((ax?.x(2, 0 as never) ?? 0) * d.dpr);
    const b1 = d.rects.filter((r) => r.x >= x1 && r.x < x2);
    expect(b1.length).toBeGreaterThan(0);
    expect(Math.max(...b1.map((r) => r.y + r.h))).toBeLessThanOrEqual(o.rows[2].y1);
    const b0 = d.rects.filter((r) => r.x < x1);
    expect(Math.max(...b0.map((r) => r.y + r.h))).toBeGreaterThan(o.rows[3].y0);
    expect(images(d.rec).length).toBeGreaterThan(0);
    // The tiles are inset a pixel inside the track, as the rows are.
    for (let k = 0; k < d.rec.opCount; k++)
      if (d.rec.ops[4 * k] === Op.image)
        expect(d.rec.args[d.rec.ops[4 * k + 3] + 5]).toBe(Math.round((d.top + 1) * d.dpr));
  });

  it("a throw inside the tiles leaves the view as it was: its cache, color, top and height", async () => {
    const d = await draw(mixedScene());
    const f = new LayerFrame(d.t.ctx, FingerprintRenderer as never, d.scene.section, d.k);
    const rec = new RecordingBackend();
    f.draw(rec);
    const v = f.v;
    const before = [v.cache, v.color, v.top, v.height];
    expect(before[0]).not.toBeNull();
    // The tiles' image call throws: the composer catches it and omits the layer, the view goes on to the next one.
    v.backend = Object.assign(new RecordingBackend(), {
      image: () => {
        throw new Error("image");
      },
    });
    expect(() => FingerprintRenderer.draw?.(v, d.s, null)).toThrow("image");
    expect([v.cache, v.color, v.top, v.height]).toEqual(before);
    expect(() => FingerprintRenderer.tiny?.(v, d.s, null)).toThrow("image");
    expect([v.cache, v.color, v.top, v.height]).toEqual(before);
  });

  it("the format half alone draws `matrix`, with one `svs.render.fallback` for the layer", async () => {
    const t = await open(mixedScene(), {
      registry: [{ svs: InteractionsViewerRegistry.svs }, ShortFingerprint] as never,
    });
    const fb = t.ctx.issues.current.value.filter((i) => i.code === "svs.render.fallback");
    expect(fb).toHaveLength(1);
    expect(t.ctx.graph.assembled.value?.view?.prepared.some((p) => p.renderer === "matrix")).toBe(true);
  });
});
