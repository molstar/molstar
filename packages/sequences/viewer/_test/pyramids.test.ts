/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The index for scale (G1; seqspec:sources#levels-and-the-zoom-gate): segment windows, the pyramids and a column's
// reads against brute force; and the reductions by theme kind (LY2, K1a's headless half): a `classes` theme's cells by
// their plurality class per device column, a `mean` theme's by the mean value, at three zooms; past 255 classes no
// folding; a dense rowless layer colored per position; items counted once at their center as points.

import { createSeqViewSpec, translate } from '@molstar/seq-core/builder';
import { describe, expect, it } from 'vitest';
import { type Features, type PerPosition, positionKey, rowValues } from '../src/prepare.js';
import {
  binClass,
  binCounts,
  binValues,
  classPyramid,
  classPyramidFrom,
  countPyramid,
  levelFor,
  PYRAMID_BASE,
  predecessor,
  prefixReach,
  segWindow,
  valuePyramid,
} from '../src/pyramids.js';

import { type CellColoring, cellClasses, densityOn, ItemDensity } from '../src/reductions.js';
import { layer, opened, preparedRenderers } from './layout-fixtures.js';

const P = { renderers: preparedRenderers() };

function random(seed: number) {
  let s = seed >>> 0;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 2 ** 32;
  };
}

/** Items sorted by start, overlapping or not, with boundaries (start === stop). */
function items(rnd: () => number, n: number, maxLen: number) {
  const list = Array.from({ length: n }, () => {
    const s = Math.floor(rnd() * 1000);
    return [s, rnd() < 0.1 ? s : s + Math.floor(rnd() * maxLen) + 1] as const;
  }).sort((a, b) => a[0] - b[0]);
  const start = Int32Array.from(list, (x) => x[0]);
  const stop = Int32Array.from(list, (x) => x[1]);
  const reach = new Int32Array(n);
  prefixReach(start, stop, 0, n, reach);
  return { start, stop, reach };
}
/** Brute force: items overlapping [from, to) (a boundary: at a position in it). */
const overlaps = (s: number, e: number, from: number, to: number) =>
  s === e ? s >= from && s < to : s < to && e > from;

/** Runs of classes (as cytobands), with noise and gaps (0). */
function bands(n: number, seed: number, classes = 5) {
  const rnd = random(seed);
  const out = new Uint8Array(n);
  for (let i = 0, c = 1; i < n; i++) {
    if (rnd() < 0.01) c = 1 + Math.floor(rnd() * classes);
    out[i] = rnd() < 0.15 ? 0 : rnd() < 0.2 ? 1 + Math.floor(rnd() * classes) : c;
  }
  return out;
}

describe('segments', () => {
  it('predecessor: the last index at or below x, lo - 1 when none', () => {
    const a = Int32Array.from([1, 3, 3, 7]);
    expect([0, 1, 3, 4, 7, 9].map((x) => predecessor(a, 0, 4, x))).toEqual([-1, 0, 2, 2, 3, 3]);
    expect(predecessor(a, 2, 2, 5)).toBe(1);
  });

  it.each([
    ['overlapping', 400, 200],
    ['sparse, short', 40, 3],
    ['empty', 0, 10],
  ])('a window holds every item overlapping it (%s), from two searches', (_, n, maxLen) => {
    const rnd = random(n + maxLen);
    const t = items(rnd, n, maxLen);
    const out = new Int32Array(2);
    for (let q = 0; q < 300; q++) {
      const from = Math.floor(rnd() * 1100) - 50;
      const to = from + Math.floor(rnd() * 120);
      segWindow(t.start, t.reach, 0, n, from, to, out);
      for (let k = 0; k < n; k++) {
        const inside = k >= out[0] && k < out[1];
        if (overlaps(t.start[k], t.stop[k], from, to)) expect(inside, `${k} in [${from}, ${to})`).toBe(true);
      }
      // Nothing before the window reaches it, nothing after starts in it.
      if (out[0] > 0) expect(t.reach[out[0] - 1]).toBeLessThanOrEqual(from);
      if (out[1] < n) expect(t.start[out[1]]).toBeGreaterThanOrEqual(to);
    }
  });

  it("edges: an item ending at the window's start is out, one starting at its end is out, a boundary at it is in", () => {
    const start = Int32Array.from([0, 10, 20, 30]);
    const stop = Int32Array.from([10, 20, 20, 40]);
    const reach = new Int32Array(4);
    prefixReach(start, stop, 0, 4, reach);
    const out = new Int32Array(2);
    segWindow(start, reach, 0, 4, 20, 30, out);
    expect([...out]).toEqual([2, 3]);
    segWindow(start, reach, 0, 4, 40, 50, out);
    expect(out[1] - out[0]).toBe(0);
  });
});

describe('pyramids against brute force', () => {
  const n = 5000;
  const classes = bands(n, 11);
  const rnd = random(11);

  it('classes: the base is the exact plurality; a strict majority found at every level; shares a lower bound', () => {
    const p = classPyramid(classes);
    expect(p.cls.length).toBe(Math.ceil(Math.log2(n / 4)) + 1);
    let strict = 0;
    for (let L = 0; L < p.cls.length; L++) {
      const size = 2 ** (L + PYRAMID_BASE);
      for (let j = 0; j < p.cls[L].length; j++) {
        const counts = new Map<number, number>();
        for (let i = j * size; i < Math.min(n, (j + 1) * size); i++)
          if (classes[i]) counts.set(classes[i], (counts.get(classes[i]) ?? 0) + 1);
        const top = Math.max(0, ...counts.values());
        const c = p.cls[L][j];
        if (L === 0) expect(counts.get(c) ?? 0).toBe(top);
        if (top > size / 2) {
          strict++;
          expect(c).toBe([...counts].find(([, v]) => v === top)?.[0]);
        }
        if (L === 0 && c > 0) expect(p.share[L][j]).toBe(Math.round((top / size) * 255));
        else if (c > 0) expect(p.share[L][j]).toBeLessThanOrEqual(Math.round(((counts.get(c) ?? 0) / size) * 255));
        if (counts.size === 0) expect(c).toBe(0);
      }
    }
    expect(strict).toBeGreaterThan(100);
  });

  it('classes: no fixed preference amplified (four equal classes, as nucleotides, all win at the top)', () => {
    const r = random(5);
    const dna = Uint8Array.from({ length: 1 << 16 }, () => 1 + Math.floor(r() * 4));
    const p = classPyramid(dna);
    const L = p.cls.length - 7;
    const wins = new Set(p.cls[L]);
    expect(wins.size).toBe(4);
  });

  it("the base level filled 64 KiB at a time, equal to the whole array's; a fill can stop the build", () => {
    const big = bands(300_001, 3, 9);
    const whole = classPyramid(big);
    const starts: number[] = [];
    const chunked = classPyramidFrom(big.length, (s, out) => {
      starts.push(s);
      out.set(big.subarray(s, s + out.length));
      return true;
    });
    expect(starts).toEqual([0, 65536, 131072, 196608, 262144]);
    expect(chunked).toEqual(whole);
    expect(classPyramidFrom(big.length, (s) => s < 65536)).toBeNull();
  });

  it('values: min, max, sum and count per bin; NaN left out; without min and max, sums and counts alike', () => {
    const v = Float64Array.from({ length: n }, (_, i) => (i % 97 === 0 ? Number.NaN : Math.sin(i / 50) * 10 + rnd()));
    const p = valuePyramid(v);
    const q = valuePyramid(v, false);
    expect([q.min.length, q.max.length]).toEqual([0, 0]);
    expect(q.sum).toEqual(p.sum);
    expect(q.count).toEqual(p.count);
    for (let L = 0; L < p.min.length; L++) {
      const size = 2 ** (L + PYRAMID_BASE);
      for (let j = 0; j < p.min[L].length; j++) {
        const xs = Array.from(v.subarray(j * size, Math.min(n, (j + 1) * size))).filter((x) => !Number.isNaN(x));
        expect(p.count[L][j]).toBe(xs.length);
        expect(p.min[L][j]).toBe(Math.min(...xs));
        expect(p.max[L][j]).toBe(Math.max(...xs));
        expect(p.sum[L][j]).toBeCloseTo(
          xs.reduce((a, b) => a + b, 0),
          1,
        );
      }
    }
  });

  it("counts: per bin and slot the items overlapping it, or (points) each once in its center's bin; an offset", () => {
    const t = items(rnd, 300, 150);
    const slot = Uint8Array.from(t.start, (_, k) => k % 3);
    for (const offset of [0, 100])
      for (const point of [false, true]) {
        const p = countPyramid(1200, 3, 3, { ...t, slot, from: 0, to: 300, offset, point });
        for (let L = 0; L < p.counts.length; L++) {
          const size = 2 ** (L + 3);
          let total = 0;
          for (let j = 0; j < p.counts[L].length / 3; j++)
            for (let s = 0; s < 3; s++) {
              let want = 0;
              for (let k = 0; k < 300; k++) {
                const a = t.start[k] - offset;
                const e = Math.max(t.stop[k], t.start[k] + 1) - offset;
                const c = (a + e - 1) >> 1;
                const inBin = point ? c >= j * size && c < (j + 1) * size : overlaps(a, e, j * size, (j + 1) * size);
                if (slot[k] === s && inBin) want++;
              }
              expect(p.counts[L][j * 3 + s]).toBe(want);
              total += want;
            }
          if (point && offset === 0) expect(total).toBe(300);
        }
      }
  });

  it('the level for a column: the widest bin not wider; none below the base', () => {
    expect([levelFor(3, 10), levelFor(4, 10), levelFor(1000, 10), levelFor(1e9, 10), levelFor(100, 10, 6)]).toEqual([
      -1, 2, 9, 11, 6,
    ]);
  });

  it("a column's reads: values' min, max, sum and count and counts' sums and maxima over its bins", () => {
    const v = Float64Array.from({ length: n }, (_, i) => (i % 13 === 0 ? Number.NaN : (i * 7919) % 101));
    const p = valuePyramid(v);
    const t = items(rnd, 200, 60);
    const slot = Uint8Array.from(t.start, (_, k) => k % 2);
    const c = countPyramid(1000, 3, 2, { ...t, slot, from: 0, to: 200 });
    const out = new Float64Array(4);
    for (const level of [2, 5, 8])
      for (let j0 = 0; j0 < n >> level; j0 += 3) {
        const j1 = j0 + 3;
        binValues(p, level, j0, j1, out);
        const xs = Array.from(v.subarray(j0 << level, Math.min(n, j1 << level))).filter((x) => !Number.isNaN(x));
        expect([out[0], out[1], out[3]]).toEqual([Math.min(...xs), Math.max(...xs), xs.length]);
        expect(out[2]).toBeCloseTo(
          xs.reduce((a, b) => a + b, 0),
          6,
        );
      }
    for (let j0 = 0; j0 < 1000 >> 4; j0 += 2) {
      binCounts(c, 4, j0, j0 + 2, false, out);
      const sum = [0, 1].map((s) => c.counts[1][j0 * 2 + s] + (c.counts[1][(j0 + 1) * 2 + s] ?? 0));
      expect([out[0], out[1]]).toEqual(sum);
      binCounts(c, 4, j0, j0 + 2, true, out);
      const max = [0, 1].map((s) => Math.max(c.counts[1][j0 * 2 + s], c.counts[1][(j0 + 1) * 2 + s] ?? 0));
      expect([out[0], out[1]]).toEqual(max);
    }
  });
});

/** A coloring over a block's class array: color per class, alpha 255, the class as its key. */
function coloring(cls: Uint8Array, colors: (c: number) => number, alpha = (_c: number) => 255): CellColoring {
  return {
    color: (_b, i) => (cls[i] === 0 ? -1 : colors(cls[i])),
    alpha: (_b, i) => alpha(cls[i]),
    key: (_b, i) => cls[i],
  };
}

describe('reductions by theme kind (K1a, the headless half)', () => {
  const n = 40_000;
  const cls = bands(n, 21, 6);
  const PALETTE = [0, 0xe41a1c, 0x377eb8, 0x4daf4a, 0x984ea3, 0xff7f00, 0xffff33];

  it('classes: per device column the plurality class of its bins, as brute force, at three zooms', () => {
    const c = cellClasses(
      coloring(cls, (k) => PALETTE[k]),
      0,
      0,
      n,
    );
    if (!c) throw new Error('no classes');
    expect(c.count).toBe(6);
    const out = new Float64Array(1);
    for (const per of [6, 40, 300]) {
      const L = levelFor(per, c.pyramid.cls.length);
      let columns = 0;
      let equal = 0;
      for (let col = 0; col * per < n; col++) {
        // The bins whose start lies in the column; brute force over their positions.
        const j0 = Math.ceil((col * per) / 2 ** L);
        const j1 = Math.ceil(((col + 1) * per) / 2 ** L);
        if (j1 <= j0) continue;
        const counts = new Map<number, number>();
        for (let i = j0 << L; i < Math.min(n, j1 << L); i++)
          if (cls[i]) counts.set(PALETTE[cls[i]], (counts.get(PALETTE[cls[i]]) ?? 0) + 1);
        const top = Math.max(0, ...counts.values());
        const got = c.color[binClass(c.pyramid, L, j0, j1, out)];
        columns++;
        if (counts.get(got) === top) equal++;
        // A strict majority of the cells is always found.
        if (top > 0.5 * ((j1 - j0) << L)) expect(counts.get(got), `${per}: column ${col}`).toBe(top);
      }
      // Where no class holds most of a column the tournament may differ (measured: 99.5, 99.1 and 93.3% equal).
      expect(equal / columns, `${per} positions a column`).toBeGreaterThan(0.9);
    }
  });

  it('mean: per device column palette(mean) from the sums and counts, exactly as brute force, at three zooms', () => {
    const rnd = random(4);
    const v = Float64Array.from({ length: n }, (_, i) => (i % 31 === 0 ? Number.NaN : i / n + 0.05 * rnd()));
    const p = valuePyramid(v, false);
    const out = new Float64Array(4);
    for (const per of [6, 40, 300]) {
      const L = levelFor(per, p.sum.length);
      for (let col = 0; col * per < n; col++) {
        const j0 = Math.ceil((col * per) / 2 ** L);
        const j1 = Math.ceil(((col + 1) * per) / 2 ** L);
        if (j1 <= j0) continue;
        binValues(p, L, j0, j1, out);
        const xs = Array.from(v.subarray(j0 << L, Math.min(n, j1 << L))).filter((x) => !Number.isNaN(x));
        expect(out[2] / out[3]).toBeCloseTo(xs.reduce((a, b) => a + b, 0) / xs.length, 9);
      }
    }
  });

  it('past 255 (color, alpha) classes a block has none (the raw path), never one class for the rest', () => {
    const many = Uint16Array.from({ length: 3000 }, (_, i) => 1 + (i % 300));
    const at = (k: number) => 0x010101 * k;
    const c300: CellColoring = { color: (_b, i) => at(many[i]), alpha: () => 255, key: () => -1 };
    expect(cellClasses(c300, 0, 0, 3000)).toBeNull();
    const c255: CellColoring = { color: (_b, i) => at(1 + (i % 255)), alpha: () => 255, key: () => -1 };
    const ok = cellClasses(c255, 0, 0, 3000);
    expect([ok?.count, ok?.color[1], ok?.color[255]]).toEqual([255, at(1), at(255)]);
    // A block's classes are its own: the same block after a block that overflowed is unchanged.
    const again = cellClasses(c255, 1, 0, 3000);
    expect(again?.pyramid).toEqual(ok?.pyramid);
  });

  it("a theme's alpha is its own byte: one color at two alphas is two classes", () => {
    const two = Uint8Array.from({ length: 64 }, (_, i) => 1 + (i >> 5));
    const c = cellClasses(
      coloring(
        two,
        () => 0x336699,
        (k) => (k === 1 ? 255 : 64),
      ),
      0,
      0,
      64,
    );
    expect([c?.count, c?.color[1], c?.color[2], c?.alpha[1], c?.alpha[2]]).toEqual([2, 0x336699, 0x336699, 255, 64]);
  });

  it("a dense rowless color layer: every position its own key, so its cells take every value's color", async () => {
    const b = createSeqViewSpec({ title: 'GC' });
    const length = 20_000;
    const dna = Array.from({ length }, (_, i) => 'ACGT'[(i * 7) % 4]).join('');
    const a = b.data.sequence({ ref: 'a', alphabet: 'dna', value: dna });
    const value = Float64Array.from({ length }, (_, i) => i / length);
    const gc = b.data.annotation({ ref: 'gc', geometry: 'values', target: a, value: { value } });
    b.view({ ref: 'v', label: 'V' }).block({ sequence: a }).track({ label: 'GC' }).color({ data: gc });
    const t = await opened(b.getDocument(), P);
    const s = layer(t.ctx, 'GC').state as PerPosition;
    expect([s.row[0], s.value[0]?.length, s.pyramid[0]?.sum.length]).toEqual([null, length, 14]);
    const values = s.value[0] as Float64Array;
    // A theme by value in ten colors (a `classes` theme over values): keys never stand for a value.
    const decile = (i: number) => Math.min(9, Math.floor(values[i] * 10));
    const c = cellClasses(
      {
        color: (_b, i) => 0x111111 * (1 + decile(i)),
        alpha: () => 255,
        key: (bb, i) => positionKey(s, bb, i, true),
      },
      0,
      0,
      length,
    );
    expect(c?.count).toBe(10);
    // And a `mean` theme: the column means rise along the block.
    const out = new Float64Array(4);
    const pyr = s.pyramid[0];
    if (!pyr) throw new Error('no pyramid');
    const means = [0, 1, 2, 3].map((q) => {
      binValues(pyr, 10, q * 4, q * 4 + 4, out);
      return out[2] / out[3];
    });
    expect(means).toEqual([...means].sort((x, y) => x - y));
    expect(new Set(means).size).toBe(4);
  });

  it("rowValues: a `mean` theme over a row field reads each covering row's value", async () => {
    const b = createSeqViewSpec({ title: 'Scores' });
    const a = b.data.sequence({ ref: 'a', value: 'MKTAYIAKQR'.repeat(10) });
    const value = { start: [10, 40], stop: [20, 60], score: [0.25, 0.75] };
    const d = b.data.annotation({ ref: 'd', geometry: 'ranges', target: a, value });
    b.view({ ref: 'v', label: 'V' }).block({ sequence: a }).track({ label: 'Scores' }).color({ data: d });
    const t = await opened(b.getDocument(), P);
    const s = layer(t.ctx, 'Scores').state as PerPosition;
    const set = t.ctx.graph.assembled.value?.set;
    const v = rowValues(s, 0, (ann, row) => Number(set?.annotations[ann]?.fields[0]?.data[row] ?? Number.NaN));
    expect([v[15], v[45], v[30]]).toEqual([0.25, 0.75, Number.NaN]);
  });
});

describe('item density', () => {
  it('a codon-long point among 500 counts once, at its center, at every level; ranges count as coverage', async () => {
    const b = createSeqViewSpec({ title: 'Points' });
    const length = 6000;
    const dna = Array.from({ length }, (_, i) => 'ACGT'[(i * 5 + (i >> 2)) % 4]).join('');
    const gene = b.data.sequence({ ref: 'gene', alphabet: 'dna', value: dna });
    const protein = translate(b.data, { ref: 'p', from: gene, cds: [[62, 3062]] });
    const at = Array.from({ length: 500 }, (_, k) => 7 + 11 * k);
    const nt = b.data.annotation({
      ref: 'nt',
      geometry: 'ranges',
      target: gene,
      value: { start: at, stop: at.map((x) => x + 1) },
    });
    // Residue 0 over [62, 65): across the bin edge at 64.
    const aa = b.data.annotation({ ref: 'aa', geometry: 'ranges', target: protein, value: { start: [0], stop: [1] } });
    b.view({ ref: 'v', label: 'V' })
      .block({ sequence: gene })
      .track({ label: 'Sites' })
      .points({ data: [nt, aa] });
    const t = await opened(b.getDocument(), P);
    const f = layer(t.ctx, 'Sites').state as Features;
    expect(f.count).toBe(501);
    const k = Array.from(f.stop).findIndex((e, i) => e - f.start[i] === 3);
    expect([f.start[k], f.stop[k]]).toEqual([62, 65]);
    const point = new ItemDensity(f, () => 0x1f77b4, true).pyramid(0, 0, length);
    const cover = new ItemDensity(f, () => 0x1f77b4, false).pyramid(0, 0, length);
    if (!point || !cover) throw new Error('no pyramids');
    for (const c of point.counts) expect(c.reduce((x, y) => x + y, 0)).toBe(501);
    // The codon's center, 63, in bin 0 of 64 positions; as coverage it is in bins 0 and 1.
    expect(cover.counts[0].reduce((x, y) => x + y, 0)).toBe(502);
  });

  it("slots: a layer's first 8 colors, later ones sharing the last; none on a rows track", () => {
    const n = 20;
    const f = {
      count: n,
      start: Int32Array.from({ length: n }, (_, k) => k * 10),
      stop: Int32Array.from({ length: n }, (_, k) => k * 10 + 5),
      lanes: 1,
      seg: Int32Array.of(0, n),
      hasRows: false,
    } as unknown as Features;
    const d = new ItemDensity(f, (k) => 0x10 * (k % 10), false);
    expect([d.slots, Array.from(d.color), d.slot[9], d.slot[17]]).toEqual([
      8,
      [0, 0x10, 0x20, 0x30, 0x40, 0x50, 0x60, 0x70],
      7,
      7,
    ]);
    expect(d.pyramid(0, 0, 200)?.counts[0].length).toBe(4 * 8);
    expect(new ItemDensity({ ...f, hasRows: true }, () => 0, false).pyramid(0, 0, 200)).toBeNull();
  });

  it('one mode per layer: density past the threshold, items again only below 0.8 of it', () => {
    let on = false;
    const modes: boolean[] = [];
    for (const count of [90, 101, 95, 85, 81, 79, 95, 101]) {
      on = densityOn(count, 100, on);
      modes.push(on);
    }
    expect(modes).toEqual([false, true, true, true, true, false, false, true]);
  });
});
