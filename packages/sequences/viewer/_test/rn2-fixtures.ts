/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Fixtures for RN2's hot paths (K4): documents of two-color points and ranges (density composed), many pairs with one
// long pair (the interval tree), a 40,000-wide matrix (its tiles); test entries that draw a feature layer's density
// through the shared helper (`drawDensity`, `densityMode`), as `ranges`, `points` and `variants` will (RN2's Sonnet
// half), and the readers of a recorded band.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import {
  DensityKind,
  type DrawView,
  densityColumns,
  densityMode,
  drawDensity,
  type Features,
  itemsInView,
  NO_COLOR,
  RectBatch,
} from "../src/draw/api.ts";
import { densityKindOf } from "../src/draw/density.ts";
import { KernelRenderers, type LayerRenderer } from "../src/index.ts";
import type { RecordingBackend } from "../src/testing/index.ts";

/** A seeded generator (the same documents every run). */
export function rng(seed = 7) {
  let s = seed;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 2 ** 32;
  };
}

type DensityOptions = {
  n?: number;
  items?: number;
  blocks?: number;
  type?: "points" | "ranges";
  len?: number;
  density?: string;
};
/** The items of `densityDoc`: per row its chain, start, stop and category (`k0` first, about 60% of them). */
export function densityItems(o: DensityOptions = {}) {
  const n = o.n ?? 20_000;
  const rnd = rng();
  const t = { target: [] as string[], start: [] as number[], stop: [] as number[], category: [] as string[] };
  for (let r = 0; r < (o.items ?? 3000); r++) {
    const len = o.type === "ranges" ? 1 + Math.floor(rnd() * (o.len ?? 200)) : 1;
    const at = Math.floor(rnd() * (n - len));
    t.target.push(`s${Math.floor(rnd() * (o.blocks ?? 1))}`);
    t.start.push(at);
    t.stop.push(at + len);
    t.category.push(r === 0 || rnd() < 0.6 ? "k0" : "k1");
  }
  return t;
}
/**
 * `n` positions on one chain (or `blocks` chains) and `items` items of two categories (`densityItems`): points one
 * position wide, or ranges of up to `len`. A track "D" draws them as `type` with its `density`, colored by category.
 */
export function densityDoc(o: DensityOptions = {}) {
  const n = o.n ?? 20_000;
  const b = createSeqViewSpec({ title: "Density" });
  const view = b.view({ ref: "v", label: "V", height: "box" });
  for (let c = 0; c < (o.blocks ?? 1); c++)
    view.block({
      sequence: b.data.sequence({ ref: `s${c}`, label: `S${c}`, value: "ACDEFGHIKLMNPQRSTVWY".repeat(n / 20) }),
    });
  const data = b.data.annotation({ ref: "d", geometry: "ranges", value: densityItems(o) } as never);
  const color = { theme: "field_palette", params: { field: "category" } };
  const layer = { data, color, ...(o.density ? { density: o.density } : {}) };
  (view.track({ label: "D" }) as never as Record<string, (x: unknown) => void>)[o.type ?? "points"](layer);
  return b.getDocument();
}

/** What the test entries drew last: 1 density, 0 items, and the items in view counted then. */
export const lastDraw = { mode: -1, count: 0, columns: 0 };
const rects = new RectBatch();
/** Set by a test: per block, the filter the density reads (`v.positions`), to draw a block through its items. */
export const forced: { positions: readonly unknown[] | null } = { positions: null };
/**
 * A feature entry drawn through the shared density (the shape of RN2's `ranges`, `points`, `variants`): one mode per
 * frame (`ranges`: over the data area's device columns; points: over 64), the band over the whole entry; in items, one
 * rect per item in view (enough for these tests).
 */
export function densityEntry(name: "points" | "ranges"): LayerRenderer {
  const base = KernelRenderers.find((r) => r.name === name) as LayerRenderer<Features>;
  const point = name !== "ranges";
  return {
    ...base,
    draw(v: DrawView, f: Features, p: unknown) {
      const threshold = point ? 64 : densityColumns(v);
      const count = itemsInView(v, f, threshold + 1);
      const on = densityMode(v, count, threshold);
      lastDraw.mode = on ? 1 : 0;
      lastDraw.count = count;
      lastDraw.columns = threshold;
      rects.reset();
      const y0 = Math.round(v.top * v.dpr) | 0;
      const y1 = Math.round((v.top + v.height) * v.dpr) | 0;
      if (forced.positions) v.positions = forced.positions as never;
      if (on) drawDensity(v, f, f, point, densityKindOf(p), NO_COLOR, rects, y0, y1);
      else
        for (let k = 0; k < f.count; k++) {
          rects.begin(0x123456);
          rects.push(v.xd(f.block[k], 0), y0 + 1, 1, 1);
        }
      if (rects.n > 0) v.backend.rects(rects);
    },
  } as LayerRenderer;
}
export const withDensity = (name: "points" | "ranges") =>
  [...KernelRenderers.filter((r) => r.name !== name), densityEntry(name)] as LayerRenderer[];
export { DensityKind };

/** A recording's rects per device column: [y, h, color] bottom up, for columns [x0, x1). */
export function bandColumns(r: RecordingBackend, x0: number, x1: number): [number, number, number][][] {
  const out: [number, number, number][][] = Array.from({ length: x1 - x0 }, () => []);
  const b = r.rectItems;
  for (let i = 0; i < b.n; i++)
    if (b.x[i] >= x0 && b.x[i] < x1 && b.w[i] === 1) out[b.x[i] - x0].push([b.y[i], b.h[i], b.color[i]]);
  for (const c of out) c.sort((a, b) => b[0] - a[0]);
  return out;
}

/**
 * `count` short pairs within one chain of `n`, the first of them a long one over nearly the whole chain (`long: false`:
 * short too); `colorBy`: colored by that field through a palette ("category": three colors and the first's; a field
 * they lack: no color).
 */
export function pairsDoc(o: { n?: number; count?: number; arcs?: string; colorBy?: string; long?: boolean } = {}) {
  const n = o.n ?? 20_000;
  const rnd = rng(11);
  const b = createSeqViewSpec({ title: "Pairs" });
  const a = b.data.sequence({ ref: "a", label: "A", value: "ACDEFGHIKLMNPQRSTVWY".repeat(n / 20) });
  const start: number[] = [0];
  const bStart: number[] = [o.long === false ? 30 : n - 2];
  const category: string[] = ["long"];
  for (let k = 1; k < (o.count ?? 5000); k++) {
    const at = Math.floor(rnd() * (n - 40));
    start.push(at);
    bStart.push(at + 3 + Math.floor(rnd() * 30));
    category.push(`k${k % 3}`);
  }
  const pairs = b.data.annotation({
    ref: "pp",
    geometry: "pairs",
    target: a,
    value: { start, b_target: start.map(() => "a"), b_start: bStart, category },
  } as never);
  const v = b.view({ ref: "v", label: "V", height: "box" }).block({ sequence: a });
  v.track({ label: "Letters" }).letters();
  const color = o.colorBy ? { color: { theme: "field_palette", params: { field: o.colorBy } } } : {};
  v.track({ label: "Contacts" }).pairs({ data: pairs, arcs: o.arcs ?? "always", ...color } as never);
  return b.getDocument();
}

/** A matrix of `n` positions over 3 states (a fingerprint over titin's length), one chain. */
export function wideMatrix(n = 40_000) {
  const b = createSeqViewSpec({ title: "Wide matrix" });
  const a = b.data.sequence({ ref: "a", label: "A", value: "ACDEFGHIKLMNPQRSTVWY".repeat(n / 20) });
  const index = Array.from({ length: n }, (_, i) => i);
  const value = { index, value: index.map((i) => [(i % 7) / 7, (i % 11) / 11, i % 2]) };
  const m = b.data.annotation({ ref: "m", geometry: "matrix", target: a, states: ["x", "y", "z"], value } as never);
  b.view({ ref: "v", label: "V", height: "box" }).block({ sequence: a }).track({ label: "M" }).matrix({ data: m });
  return b.getDocument();
}
