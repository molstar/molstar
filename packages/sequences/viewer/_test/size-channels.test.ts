/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The size channels in the draws (seqspec:format#layer-types, seqspec:format#themes): `pairs` width and end radius
// (arcs batched by color then width, the widest first), `ranges` thickness (a bar's share of its lane) and a boundary's
// radius, `points` radius, `values` thickness (a bar's share of its cell); a size theme on a channel its renderer does
// not draw is reported, once.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { describe, expect, it } from 'vitest';
import { drawn, record, rects, texts } from './draw-fixtures.js';

const SEQ = 'ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY';
const b0 = () => createSeqViewSpec({ title: 'Sizes' });
const scale = (field: string, range: [number, number], more: Record<string, unknown> = {}) => ({
  theme: 'field_scale',
  params: { field, range, ...more },
});

describe('pairs', () => {
  it("arcs are stroked at their width theme's, grouped by color then widest first; the ends take the radius theme's", async () => {
    const b = b0();
    const a = b.data.sequence({ ref: 'a', value: SEQ });
    const count = [1, 8, 4, 8, 2, 4];
    const pairs = b.data.annotation({
      ref: 'pp',
      geometry: 'pairs',
      target: a,
      value: {
        start: [0, 2, 4, 6, 8, 10],
        b_target: count.map(() => 'a'),
        b_start: [20, 22, 24, 26, 28, 30],
        category: ['x', 'y', 'x', 'y', 'x', 'y'],
        count,
      },
    } as never);
    const v = b.view({ ref: 'v', label: 'V', height: 'box' }).block({ sequence: a });
    v.track({ label: 'Pairs' }).pairs({
      data: pairs,
      display: 'both',
      width: scale('count', [1, 8]),
      radius: scale('count', [2, 6]),
    } as never);
    const t = await drawn(b.getDocument(), { width: 800, height: 200 });
    const { rec } = record(t.ctx, 'Pairs');
    const p = rec.pathItems;
    const arcs = Array.from({ length: p.n }, (_, i) => [p.color[i], p.width[i]]);
    // Half device px: a count of 1 is 1 px (2), 8 is 8 px (16); two colors, each widest first.
    expect(arcs.map(([, w]) => w).sort((x, y) => x - y)).toEqual([2, 4, 8, 8, 16, 16].map((w) => w));
    const byColor = new Map<number, number[]>();
    for (const [c, w] of arcs) byColor.set(c, [...(byColor.get(c) ?? []), w]);
    expect(byColor.size).toBe(2);
    for (const ws of byColor.values()) expect(ws).toEqual([...ws].sort((x, y) => y - x));
    // Ends: radius 2 px at the smallest count, 6 at the largest.
    const r = Array.from(rec.glyphItems.radius.subarray(0, rec.glyphItems.n));
    expect([Math.min(...r), Math.max(...r)]).toEqual([2, 6]);
    t.ctx.dispose();
  });

  it("without a width theme the arcs keep the viewer's width", async () => {
    const b = b0();
    const a = b.data.sequence({ ref: 'a', value: SEQ });
    const pairs = b.data.annotation({
      ref: 'pp',
      geometry: 'pairs',
      target: a,
      value: { start: [0, 2], b_target: ['a', 'a'], b_start: [20, 22], category: ['x', 'y'] },
    } as never);
    b.view({ ref: 'v', label: 'V', height: 'box' })
      .block({ sequence: a })
      .track({ label: 'Pairs' })
      .pairs({ data: pairs });
    const t = await drawn(b.getDocument(), { width: 800, height: 200 });
    const p = record(t.ctx, 'Pairs').rec.pathItems;
    expect(Array.from({ length: p.n }, (_, i) => p.width[i])).toEqual([2, 2]);
    t.ctx.dispose();
  });
});

describe('ranges and points', () => {
  function doc() {
    const b = b0();
    const a = b.data.sequence({ ref: 'a', value: SEQ });
    const items = b.data.annotation({
      ref: 'f',
      geometry: 'ranges',
      target: a,
      value: { start: [0, 10, 20], stop: [8, 18, 28], conf: [0.25, 0.5, 1], category: ['x', 'x', 'x'] },
    } as never);
    const v = b.view({ ref: 'v', label: 'V', height: 'box' }).block({ sequence: a });
    v.track({ label: 'Bars' }).ranges({
      data: items,
      thickness: scale('conf', [0.25, 1], { domain: [0.25, 1] }),
    } as never);
    v.track({ label: 'Dots' }).points({ data: items, radius: scale('conf', [2, 6]) } as never);
    return b.getDocument();
  }
  it("a bar is its thickness theme's share of its lane, centered; the full lane at a share of 1", async () => {
    const t = await drawn(doc(), { width: 800, height: 300 });
    const r = rects(record(t.ctx, 'Bars').rec)
      .filter((x) => x[2] > 10)
      .sort((p, q) => p[0] - q[0]);
    expect(r.length).toBe(3);
    const full = r[2][3];
    expect(r[0][3]).toBe(Math.max(1, Math.round(full * 0.25)));
    expect(r[1][3]).toBe(Math.max(1, Math.round(full * (((0.5 - 0.25) / 0.75) * 0.75 + 0.25))));
    // Centered on the lane: the thinner bars sit inside the full one.
    expect(r[0][1]).toBeGreaterThan(r[2][1]);
    expect(r[0][1] + r[0][3]).toBeLessThan(r[2][1] + r[2][3] + 1);
    t.ctx.dispose();
  });
  it("a glyph's radius is its radius theme's, in CSS px times the device ratio", async () => {
    const t = await drawn(doc(), { width: 800, height: 300, dpr: 2 });
    const g = record(t.ctx, 'Dots').rec.glyphItems;
    const radii = Array.from(g.radius.subarray(0, g.n)).sort((x, y) => x - y);
    expect([radii[0], radii.at(-1)]).toEqual([4, 12]);
    t.ctx.dispose();
  });
});

describe('values', () => {
  function doc(thickness?: unknown) {
    const b = b0();
    const a = b.data.sequence({ ref: 'a', value: SEQ });
    const d = b.data.annotation({
      ref: 'd',
      geometry: 'values',
      target: a,
      value: { index: [0, 1, 2, 3], value: [1, 2, 3, 4] },
    });
    b.view({ ref: 'v', label: 'V', height: 'box' })
      .block({ sequence: a })
      .track({ label: 'Bars' })
      .values({ data: d, ...(thickness ? { thickness } : {}) } as never);
    return b.getDocument();
  }
  it("a bar is its thickness theme's share of its cell, centered", async () => {
    const wide = await drawn(doc(), { width: 400, height: 120 });
    const [w0] = rects(record(wide.ctx, 'Bars').rec).sort((p, q) => p[0] - q[0]);
    const half = await drawn(doc({ theme: 'uniform', params: { size: 0.5 } }), { width: 400, height: 120 });
    const [w1] = rects(record(half.ctx, 'Bars').rec).sort((p, q) => p[0] - q[0]);
    expect(w1[2]).toBe(Math.round(w0[2] / 2));
    expect(w1[0]).toBe(w0[0] + ((w0[2] - w1[2]) >> 1));
    wide.ctx.dispose();
    half.ctx.dispose();
  });
});

describe('a channel its renderer does not draw', () => {
  it("is reported once and drawn at the viewer's size: a matrix's height (its bands are image tiles)", async () => {
    const b = b0();
    const a = b.data.sequence({ ref: 'a', value: SEQ });
    const index = Array.from({ length: 20 }, (_, i) => i);
    const value = { index, value: index.map((i) => [i % 3, 1]) };
    const m = b.data.annotation({ ref: 'm', geometry: 'matrix', target: a, states: ['x', 'y'], value } as never);
    b.view({ ref: 'v', label: 'V', height: 'box' })
      .block({ sequence: a })
      .track({ label: 'M' })
      .matrix({ data: m, height: { theme: 'uniform', params: { size: 0.5 } } } as never);
    const t = await drawn(b.getDocument(), { width: 400, height: 160 });
    const said = t.ctx.issues.current.value.filter((i) => i.code === 'svs.render.theme');
    expect(said.map((i) => i.message)).toEqual(['"height" is drawn at the viewer\'s size here']);
    t.ctx.dispose();
  });
});

describe('EN2b: pairs height and labels, variants and secondary thickness', () => {
  function pairsDoc(more: Record<string, unknown>) {
    const b = b0();
    const a = b.data.sequence({ ref: 'a', value: SEQ });
    const pairs = b.data.annotation({
      ref: 'pp',
      geometry: 'pairs',
      target: a,
      value: { start: [0, 10], b_target: ['a', 'a'], b_start: [30, 12], count: [1, 4], label: [null, 'x'] },
    } as never);
    b.view({ ref: 'v', label: 'V', height: 'box' })
      .block({ sequence: a })
      .track({ label: 'Pairs', height: 6 })
      .pairs({ data: pairs, display: 'arcs', side: 'above', ...more } as never);
    return b.getDocument();
  }
  /** Each arc's apex height over its base (device px), from its quad's control point. */
  function apexes(t: Awaited<ReturnType<typeof drawn>>) {
    const p = record(t.ctx, 'Pairs').rec.pathItems;
    return Array.from({ length: p.n }, (_, q) => {
      const c = p.cstart[q];
      return Math.round((p.coords[c + 1] - p.coords[c + 3]) / 2);
    });
  }
  it("an arc's apex is its height theme's share of the band, whatever its span", async () => {
    const own = apexes(await drawn(pairsDoc({}), { width: 800, height: 300 }));
    // By span: the wide arc is the band's height, the narrow one a small share.
    expect(own[0]).toBeGreaterThan(4 * own[1]);
    const t = await drawn(pairsDoc({ height: scale('count', [0.25, 1]) }), { width: 800, height: 300 });
    const sized = apexes(t);
    // Count 1 (the wide arc) a quarter of the band, count 4 (the narrow one) all of it.
    expect(sized[1]).toBe(own[0]);
    expect(Math.abs(sized[0] - own[0] / 4)).toBeLessThanOrEqual(1);
    expect(t.ctx.issues.current.value.filter((i) => i.code === 'svs.render.theme')).toEqual([]);
    t.ctx.dispose();
  });
  it("labels: a row's label, else the width theme's field's value, under the apex where it fits", async () => {
    const t = await drawn(pairsDoc({ labels: true, width: scale('count', [1, 4]) }), { width: 800, height: 300 });
    const said = texts(record(t.ctx, 'Pairs').rec).map((x) => x.text);
    // The wide arc (no label) shows its count; the narrow one (label "x") is too low for text.
    expect(said).toEqual(['1']);
    t.ctx.dispose();
  });
  it("the sashimi preset: the + strand's junctions rise from the line, the - strand's hang below it, labeled by count", async () => {
    const b = b0();
    const g = b.data.sequence({ ref: 'g', alphabet: 'dna', value: 'ACGT'.repeat(30) });
    const j = b.data.annotation({
      ref: 'j',
      geometry: 'pairs',
      target: g,
      value: { start: [5, 40], b_target: ['g', 'g'], b_start: [60, 100], count: [12, 3], strand: ['+', '-'] },
    } as never);
    const s = b.view({ ref: 'v', label: 'V', height: 'box' }).block({ sequence: g });
    b.preset.sashimi(s, j);
    const t = await drawn(b.getDocument(), { width: 800, height: 300 });
    const drawnBy = (track: string) => {
      const rec = record(t.ctx, track).rec;
      const p = rec.pathItems;
      return { said: texts(rec), arcs: p.n, base: p.coords[1], apex: p.coords[3] };
    };
    const [plus, minus] = [drawnBy('Junctions (+)'), drawnBy('Junctions (-)')];
    expect([plus.said.map((x) => x.text), minus.said.map((x) => x.text)]).toEqual([['12'], ['3']]);
    expect([plus.arcs, minus.arcs]).toEqual([1, 1]);
    // Above: the arc rises from its base (the control point higher up); below: it hangs.
    expect(plus.apex).toBeLessThan(plus.base);
    expect(minus.apex).toBeGreaterThan(minus.base);
    t.ctx.dispose();
  });
  function featuresDoc(type: 'variants' | 'secondary', more: Record<string, unknown>) {
    const b = b0();
    const a = b.data.sequence({ ref: 'a', value: SEQ });
    const d = b.data.annotation({
      ref: 'd',
      geometry: 'ranges',
      target: a,
      value: { start: [2, 20], stop: [3, 32], alt: ['A', 'C'], kind: ['helix', 'helix'], w: [0.5, 1] },
    } as never);
    b.view({ ref: 'v', label: 'V', height: 'box' })
      .block({ sequence: a })
      .track({ label: 'T', height: 2 })
      [type]({ data: d, ...more } as never);
    return b.getDocument();
  }
  it("a variant's marker is its thickness theme's share of its radius", async () => {
    const plain = record((await drawn(featuresDoc('variants', {}), { width: 800, height: 200 })).ctx, 'T').rec;
    const t = await drawn(featuresDoc('variants', { thickness: scale('w', [0.5, 1], { domain: [0.5, 1] }) }), {
      width: 800,
      height: 200,
    });
    const g = record(t.ctx, 'T').rec.glyphItems;
    const r0 = plain.glyphItems.radius[0];
    expect(Array.from(g.radius.subarray(0, 2))).toEqual([r0 / 2, r0]);
    expect(t.ctx.issues.current.value.filter((i) => i.code === 'svs.render.theme')).toEqual([]);
    t.ctx.dispose();
  });
  it("a secondary structure's shape is its thickness theme's share of its height", async () => {
    const tall = (doc: ReturnType<typeof featuresDoc>) =>
      drawn(doc, { width: 800, height: 200 }).then((t) => Math.max(...rects(record(t.ctx, 'T').rec).map((x) => x[3])));
    // The helix's half height h at a share of 0.5, rounded.
    const full = await tall(featuresDoc('secondary', {}));
    const sized = await tall(featuresDoc('secondary', { thickness: { theme: 'uniform', params: { size: 0.5 } } }));
    expect(sized).toBe(2 * (((full / 2) * 128 + 128) >> 8));
  });
});
