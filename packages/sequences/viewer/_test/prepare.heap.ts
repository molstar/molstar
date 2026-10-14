/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for LY2's per-frame and per-move reads (seqspec:viewer#performance), on `heapScaling` (a step is one
// read): segment windows and the item count that picks items or density, the hit's `featureAt` over overlapping
// fragments, `ofPartAt` stepping and jumping from its cursor, the level choice, a column's reads of class, value and
// count pyramids, the
// cell key and the mode's hysteresis. Each is what a draw or a hit runs per column, item or move. The positive control:
// the windows' reads allocating an object each (`leaking`), measured alike, are caught.

import { createSeqViewSpec, translate } from "@molstar/seq-core/builder";
import { beforeAll, describe, expect, it } from "vitest";
import {
  type Features,
  featureAt,
  type LettersState,
  type OfLetters,
  ofPartAt,
  type PerPosition,
  positionKey,
} from "../src/prepare.ts";
import {
  binClass,
  binCounts,
  binValues,
  type CountPyramid,
  levelFor,
  segWindow,
  type ValuePyramid,
} from "../src/pyramids.ts";
import { type CellClasses, cellClasses, densityOn, ItemDensity } from "../src/reductions.ts";
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from "./heap.ts";
import { layer, opened, preparedRenderers } from "./layout-fixtures.ts";
import { itemsIn, scaleDoc } from "./prepare-fixtures.ts";

const N = 50_000;
const P = { renderers: preparedRenderers() };
let sink = 0;

let genes: Features;
let variants: Features;
let of: OfLetters;
let values: PerPosition;
let bands: PerPosition;
let classes: CellClasses;
let counts: CountPyramid;
const LENGTH = 200_000;

beforeAll(async () => {
  const t = await opened(scaleDoc({ length: LENGTH, ranges: 3000, variants: 4000 }), P);
  genes = layer(t.ctx, "Genes").state as Features;
  variants = layer(t.ctx, "Variants").state as Features;
  values = layer(t.ctx, "Values").state as PerPosition;
  bands = layer(t.ctx, "Bands").state as PerPosition;
  const b = createSeqViewSpec();
  const dna = Array.from({ length: 60_000 }, (_, i) => "ACGT"[(i * 7 + (i >> 3)) % 4]).join("");
  const gene = b.data.sequence({ ref: "gene", alphabet: "dna", value: dna });
  const cds: [number, number][] = [];
  for (let s = 30; s + 600 < 60_000; s += 900) cds.push([s, s + 600]);
  const p = translate(b.data, { ref: "p", from: gene, cds });
  b.view({ ref: "v", label: "V" }).block({ sequence: gene }).track({ label: "P" }).letters({ of: p });
  const u = await opened(b.getDocument(), P);
  of = (layer(u.ctx, "P").state as LettersState).of as OfLetters;
  const c = cellClasses(
    {
      color: (bb, i) => (bands.row[bb]?.[i] ?? -1) & 7,
      alpha: () => 255,
      key: (bb, i) => positionKey(bands, bb, i, true),
    },
    0,
    0,
    LENGTH,
  );
  if (!c) throw new Error("no classes");
  classes = c;
  counts = new ItemDensity(genes, (k) => genes.annotation[k] & 3, false).pyramid(0, 0, LENGTH) as CountPyramid;
}, 120_000);

const window2 = new Int32Array(2);
/** Segment windows and the count in view (items or density) as the view pans and zooms. */
const windows = (i: number) => {
  const from = (i * 7919) % LENGTH;
  const to = from + 50 + (i & 4095) * 16;
  const g = i % (genes.seg.length - 1);
  segWindow(genes.start, genes.reach, genes.seg[g], genes.seg[g + 1], from, to, window2);
  sink = (sink + window2[1] - window2[0] + itemsIn(genes, i & 1, from, to, -1, 2000)) & 0xffff;
  sink = (sink + itemsIn(variants, i & 1, from, to, -1, 64)) & 0xffff;
  sink = (sink + levelFor(((i & 1023) + 1) * 0.37, 16)) & 0xffff;
  sink = (sink + (densityOn(i & 255, 100, (i & 1) === 1) ? 1 : 0)) & 0xffff;
};

describe("LY2's reads allocate nothing", () => {
  it("segment windows and the count in view (items or density) as the view pans and zooms", () => {
    expect(itemsIn(genes, 0, 0, LENGTH, -1, 1e9)).toBeGreaterThan(100);
    const r = heapScaling(N, windows);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it("the hit's `featureAt` per move, every lane, and `positionKey`", () => {
    const r = heapScaling(N, (i) => {
      const at = (i * 104729) % LENGTH;
      sink = (sink + featureAt(genes, i & 1, i % genes.lanes, at)) & 0xffff;
      sink = (sink + featureAt(variants, i & 1, 0, at)) & 0xffff;
      sink = (sink + positionKey(bands, i & 1, at, (i & 2) === 0)) & 0xffff;
    });
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it("`ofPartAt` stepping along a frame and jumping between frames", () => {
    expect(of.count).toBeGreaterThan(10_000);
    const r = heapScaling(N, (i) => {
      const at = i & 1 ? i % 60_000 : (i * 7919) % 60_000;
      sink = (sink + ofPartAt(of, 0, at)) & 0xffff;
    });
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it("a device column's reads: its plurality class, its values, its counts", () => {
    const out = new Float64Array(8);
    const v = values.pyramid[0] as ValuePyramid;
    const r = heapScaling(N, (i) => {
      const level = 2 + (i % 12);
      const j0 = (i * 31) % Math.max(1, LENGTH >> level);
      sink = (sink + binClass(classes.pyramid, level, j0, j0 + 1 + (i & 3), out) + out[0]) & 0xffff;
      binValues(v, level, j0, j0 + 2, out);
      sink = (sink + out[3]) & 0xffff;
      const cl = Math.min(6 + (i % 10), 6 + counts.counts.length - 1);
      binCounts(counts, cl, j0 >> 4, (j0 >> 4) + 2, (i & 1) === 1, out);
      sink = (sink + out[0]) & 0xffff;
    });
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it("the windows' reads allocating an object each are caught (positive control)", () => {
    expect(controlCaught(() => heapScaling(N, leaking(windows)))).toBe(true);
  });
});
