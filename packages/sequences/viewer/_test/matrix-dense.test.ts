/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// KS2: dense matrix states (seqspec:format#keyed-sets, seqspec:format#layer-types): a genotype matrix's states drawn in
// their keyed order and kept by the filter, a row each while bands are a device pixel or more, and below that a tile
// row per device pixel (⌊height·dpr⌋ rows) holding the plurality of its bands' colors (ties: the color met first),
// opaque by the share of its bands with a color (at DPR 1 and 1.5); the hit names the state whose color the draw put
// there; a patch of the layer's order and filter reorders the bands.

import { describe, expect, it } from 'vitest';
import { Hit, HoverKind, type MatrixState } from '../src/draw/api.js';
import { Dirty } from '../src/draw/index.js';
import { matrixFirstBand } from '../src/draw/matrix-tiles.js';
import { renderSvg } from '../src/svg/index.js';
import { set } from './axis-fixtures.js';
import { composed, recording } from './composer-fixtures.js';
import { drawn, entryOf, layerOf, record } from './draw-fixtures.js';
import { genotype, genotypeDoc } from './genotype-fixtures.js';

/**
 * A reduced row's pixel and drawing band (its offset in the row; -1: none) from its bands' pixels (-1: no color),
 * counted independently of the draw: colors by first appearance, the most votes (the first met on a tie) as the hue,
 * opaque by the colored share.
 */
function reduce(pixels: readonly number[]): { pixel: number; band: number } {
  const votes = new Map<number, number>();
  for (const c of pixels) if (c >= 0) votes.set(c, (votes.get(c) ?? 0) + 1);
  let hue = -1;
  for (const [c, n] of votes) if (hue < 0 || n > (votes.get(hue) as number)) hue = c;
  if (hue < 0) return { pixel: 0, band: -1 };
  const colored = pixels.filter((c) => c >= 0).length;
  const alpha = Math.floor(((hue >>> 24) * colored) / pixels.length + 0.5);
  return { pixel: ((hue & 0xffffff) | (alpha << 24)) >>> 0, band: pixels.indexOf(hue) };
}

type Opened = Awaited<ReturnType<typeof drawn>>;
const Op = { image: 13 };

/** The layer's prepared state, its color theme, and its one recorded tile with its image op's rects. */
function tileOf(t: Opened) {
  const { k } = layerOf(t.ctx, 'G');
  const p = t.ctx.graph.assembled.value?.view?.prepared[k];
  const s = p?.state as MatrixState;
  const theme = p?.themes[0];
  const { rec, frame } = record(t.ctx, 'G');
  const ops: number[][] = [];
  for (let j = 0; j < rec.opCount; j++)
    if (rec.ops[4 * j] === Op.image) {
      const a = rec.ops[4 * j + 3];
      ops.push([rec.ops[4 * j + 1], ...Array.from({ length: 8 }, (_, m) => rec.args[a + m])]);
    }
  expect(ops.length).toBe(1);
  const [img] = ops[0];
  const px = rec.images.content[img] as Uint32Array;
  if (!theme) throw new Error('no theme');
  /** The pixel band q's state has at position i, through the layer's theme. */
  const pixel = (i: number, q: number, at = genotype) => {
    const loc = { ...frame.v.loc, value: at(i, s.stateOf[q]), state: s.stateOf[q], index: i };
    const c = theme.value(loc as never);
    if (c < 0) return -1;
    return ((255 << 24) | ((c & 255) << 16) | (c & 0xff00) | ((c >> 16) & 255)) >>> 0;
  };
  return { s, px, w: rec.images.w[img], h: rec.images.h[img], op: ops[0].slice(1), frame, pixel };
}

describe("a dense matrix's per-device-row reduction", () => {
  it.each([
    [1, false],
    [1.5, false],
    [1.5, true],
  ])(
    "a tile row per device pixel, the plurality of its bands' colors by its share, DPR %f, a categorical palette %s",
    async (dpr, holes) => {
      const t = await drawn(genotypeDoc({ samples: 3000, positions: 40, holes }), { width: 900, height: 400, dpr });
      const { s, px, w, h, op, pixel } = tileOf(t);
      const e = entryOf(t.ctx, 'G');
      const H = Math.floor(e.height * dpr);
      expect(s.states).toBe(3000);
      expect(H).toBeLessThan(3000);
      // The tile is the entry's device rows, drawn 1:1 from its top.
      expect([h, op[1], op[3], op[5], op[7]]).toEqual([H, 0, H, Math.round(e.top * dpr), H]);
      expect(w).toBe(40);
      let unlikeFirst = 0;
      for (let y = 0; y < H; y++)
        for (let i = 0; i < 40; i += 3) {
          const [q0, q1] = [matrixFirstBand(y, H, 3000), matrixFirstBand(y + 1, H, 3000)];
          const row = Array.from({ length: q1 - q0 }, (_, k) => pixel(i, q0 + k));
          const want = reduce(row);
          expect(px[y * w + i], `row ${y}, position ${i}`).toBe(want.pixel);
          if (want.band > 0) unlikeFirst++;
        }
      // A reduction, not a sample of the first band: some rows' plurality is another band's color.
      expect(unlikeFirst).toBeGreaterThan(0);
    },
  );

  it("a tie goes to the color met first; the alpha is the colored share, whatever the hue's own count", async () => {
    // 640 states over 160 device rows: four bands a row, by q mod 4.
    const run = async (pattern: number[]) => {
      const value = (_i: number, q: number) => pattern[q % 4];
      const t = await drawn(genotypeDoc({ samples: 640, positions: 4, holes: true, value }), { width: 400 });
      const { px, w, h, pixel } = tileOf(t);
      expect(h).toBe(160);
      return { px: px[w + 1], first: pixel(1, 4, value), second: pixel(1, 5, value) };
    };
    // Amber (1), red (2), red, amber: two each, amber met first.
    const tie = await run([1, 2, 2, 1]);
    expect(tie.px).toBe(tie.first);
    // Red, amber, then two uncolored: red met first, half the bands colored.
    const half = await run([2, 1, 0, 0]);
    expect(half.px).toBe(((half.first & 0xffffff) | (128 << 24)) >>> 0);
    expect(half.second).not.toBe(half.first);
  });

  it('the bands partition the states: each state in exactly one device row, the rows in band order', () => {
    for (const [H, B] of [
      [160, 3000],
      [240, 2501],
      [7, 8],
    ]) {
      let next = 0;
      for (let y = 0; y < H; y++) {
        const q0 = matrixFirstBand(y, H, B);
        const q1 = matrixFirstBand(y + 1, H, B);
        expect(q0).toBe(next);
        expect(q1).toBeGreaterThan(q0);
        for (let q = q0; q < q1; q++) expect(((q * H) / B) | 0).toBe(y);
        next = q1;
      }
      expect(next).toBe(B);
    }
  });

  it("density maximum: each reduced row the strongest band's own color, the hit its state; majority by contrast", async () => {
    // Four bands a row: amber (1), none, red (2), none.
    const value = (_i: number, q: number) => [1, 0, 2, 0][q % 4];
    const at = async (density: string) => {
      const doc = genotypeDoc({ samples: 640, positions: 4, holes: true, value, layer: { density } });
      const t = await drawn(doc, { width: 400 });
      const { px, w, pixel, s, frame } = tileOf(t);
      const { section } = layerOf(t.ctx, 'G');
      const e = entryOf(t.ctx, 'G');
      const hit = Hit.create();
      t.ctx.layout.hitTest(section, (frame.v.xd(0, 1) + frame.v.xe(0, 1)) / 2, e.top + 1.5, hit);
      return { px: px[w + 1], amber: pixel(1, 4, value), red: pixel(1, 6, value), state: hit.state, s };
    };
    const max = await at('maximum');
    expect(max.px).toBe(max.red);
    expect(max.state).toBe(max.s.stateOf[6]);
    const maj = await at('majority');
    expect(maj.px).toBe(((maj.amber & 0xffffff) | (128 << 24)) >>> 0);
    expect(maj.state).toBe(maj.s.stateOf[4]);
  });

  it("density maximum over the genotypes: the first band with the row's largest value", async () => {
    const t = await drawn(genotypeDoc({ samples: 3000, positions: 40, layer: { density: 'maximum' } }), {
      width: 900,
      height: 400,
      dpr: 1.5,
    });
    const { s, px, w, pixel } = tileOf(t);
    const H = Math.floor(entryOf(t.ctx, 'G').height * 1.5);
    for (let y = 0; y < H; y += 7)
      for (let i = 0; i < 40; i += 5) {
        const [q0, q1] = [matrixFirstBand(y, H, 3000), matrixFirstBand(y + 1, H, 3000)];
        let best = q0;
        for (let q = q0; q < q1; q++) if (genotype(i, s.stateOf[q]) > genotype(i, s.stateOf[best])) best = q;
        const want = pixel(i, best);
        expect(px[y * w + i], `row ${y}, position ${i}`).toBe(want < 0 ? 0 : want);
      }
  });

  it("the hit names the row's state whose color is drawn (its first band without one); the frame at a fraction", async () => {
    const dpr = 1.5;
    const t = await drawn(genotypeDoc({ samples: 3000, positions: 40, holes: true }), { width: 900, height: 400, dpr });
    const { s, frame, pixel, px, w } = tileOf(t);
    const { section } = layerOf(t.ctx, 'G');
    const e = entryOf(t.ctx, 'G');
    const y0 = Math.round(e.top * dpr);
    const H = Math.floor(e.height * dpr);
    const hit = Hit.create();
    const x = (frame.v.xd(0, 5) + frame.v.xe(0, 5)) / 2 / dpr;
    let other = 0;
    for (let y = 0; y < H; y++) {
      t.ctx.layout.hitTest(section, x, (y0 + y + 0.5) / dpr, hit);
      const [q0, q1] = [matrixFirstBand(y, H, 3000), matrixFirstBand(y + 1, H, 3000)];
      const want = reduce(Array.from({ length: q1 - q0 }, (_, k) => pixel(5, q0 + k)));
      const band = q0 + Math.max(0, want.band);
      expect([hit.kind, hit.state], `row ${y}`).toEqual([HoverKind.seq, s.stateOf[band]]);
      // The named state's color is the pixel's hue.
      if (want.band >= 0) expect(pixel(5, band) & 0xffffff).toBe(px[y * w + 5] & 0xffffff);
      if (want.band > 0) other++;
    }
    expect(other).toBeGreaterThan(0);
  });
});

describe("a matrix's states in their keyed order", () => {
  const byDepth = { data: 'samples', field: 'depth', direction: 'desc' };

  it("bands follow the order, a filter leaves states out; cells, hits and band labels name the band's state", async () => {
    const layer = { order: byDepth, filter: [{ data: 'samples', field: 'population', values: ['AFR', 'EAS'] }] };
    const t = await drawn(genotypeDoc({ samples: 9, positions: 30, layer }), { width: 900, height: 400 });
    const { s, px, w, h, pixel } = tileOf(t);
    // Depth is (37q mod 101); AFR and EAS are q mod 3 of 0 and 2.
    const kept = [0, 2, 3, 5, 6, 8].sort((a, b) => ((b * 37) % 101) - ((a * 37) % 101));
    expect(Array.from(s.stateOf)).toEqual(kept);
    expect(s.states).toBe(6);
    expect(Array.from(s.bandOf)).toEqual(Array.from({ length: 9 }, (_, q) => kept.indexOf(q)));
    expect(h).toBe(6);
    for (let q = 0; q < 6; q++) for (let i = 0; i < 30; i++) expect(px[q * w + i]).toBe(pixel(i, q));
    const { k, section } = layerOf(t.ctx, 'G');
    const bands = t.ctx.registry.get('renderers', 'matrix')?.bands;
    expect(
      Array.from({ length: 6 }, (_, q) => bands?.label(t.ctx.graph.assembled.value?.view?.prepared[k].state, q)),
    ).toEqual(kept.map((q) => `n${q}`));
    const e = entryOf(t.ctx, 'G');
    const hit = Hit.create();
    for (let q = 0; q < 6; q++) {
      t.ctx.layout.hitTest(section, 300.5, e.top + ((q + 0.5) * e.height) / 6, hit);
      expect(hit.state).toBe(kept[q]);
    }
  });

  it("a patch of the layer's order and filter reorders the bands, and back to data order", async () => {
    const t = await drawn(genotypeDoc({ samples: 9, positions: 30 }), { width: 900, height: 400 });
    const order = () => Array.from(tileOf(t).s.stateOf);
    expect(order()).toEqual([0, 1, 2, 3, 4, 5, 6, 7, 8]);
    await t.ctx.applyPatch(set('geno', { order: [{ data: 'samples', field: 'population' }, byDepth] }));
    const pop = (q: number) => ['AFR', 'EUR', 'EAS'][q % 3];
    const depth = (q: number) => (q * 37) % 101;
    const want = [0, 1, 2, 3, 4, 5, 6, 7, 8].sort((a, b) => pop(a).localeCompare(pop(b)) || depth(b) - depth(a));
    expect(order()).toEqual(want);
    await t.ctx.applyPatch(set('geno', { filter: [{ data: 'samples', field: 'depth', min: 50 }] }));
    expect(order()).toEqual(want.filter((q) => depth(q) >= 50));
    await t.ctx.applyPatch(set('geno', { order: null, filter: [] }));
    expect(order()).toEqual([0, 1, 2, 3, 4, 5, 6, 7, 8]);
  });
});

describe('VL-04: the hit reduces its row again, whichever cache drew last', () => {
  it("after the navigator's miniature and an SVG export drew the same state, every device row hits as before", async () => {
    const dpr = 1.5;
    const t = await composed(genotypeDoc({ samples: 3000, positions: 40, holes: true }), {
      width: 900,
      height: 400,
      dpr,
    });
    const { section } = layerOf(t.ctx, 'G');
    t.frame(recording(), Dirty.tracks);
    const e = entryOf(t.ctx, 'G');
    const y0 = Math.round(e.top * dpr);
    const H = Math.floor(e.height * dpr);
    const av = t.ctx.layout.axis(section);
    if (!av) throw new Error('no axis');
    const x = (av.x(0, 5 as never) + av.x(0, 6 as never)) / 2;
    const hit = Hit.create();
    const states = () =>
      Array.from({ length: H }, (_, y) => {
        t.ctx.layout.hitTest(section, x, (y0 + y + 0.5) / dpr, hit);
        return hit.state;
      });
    const before = states();
    // Some rows name a state other than their first band's (a hit that fell back to it would pass blindly).
    expect(before.filter((st, y) => st !== matrixFirstBand(y, H, 3000)).length).toBeGreaterThan(0);
    t.frame(recording(), Dirty.nav);
    expect(states()).toEqual(before);
    await renderSvg(t.ctx, { width: 600 });
    expect(states()).toEqual(before);
  });
});
