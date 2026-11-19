/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The variants layer (_spec/README.md "Render", Layer types): a stem and a head a site, the heads of one site stacked,
// drawn against the data (a head for every variant in view, every head hit as its variant), a count strip at the bars
// and density steps, the same on a reversed, collapsed axis and through a codon map, the fallback to ranges.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { hitView, type VariantsLay } from '../index.js';
import { BigRuntime, GENE_EXONS, geneLetters, P53, VariantRows } from './fixtures.js';
import { settle } from './helpers.js';
import { frameWithout, layerLay } from './small-helpers.js';

const SIZE = { width: 800, height: 40 } as const;

type Sites = {
  start: number[];
  end: number[];
  id: string[];
  category: string[];
};

/** `n` variants a site over the given positions (ids `v0`, `v1`, …). */
function sites(positions: number[]): Sites {
  return {
    start: positions,
    end: positions.map((p) => p + 1),
    id: positions.map((_, i) => `v${i}`),
    category: positions.map((_, i) => ['missense', 'nonsense', 'synonymous'][i % 3]),
  };
}

function spec(rows: Sites, length: number, section: Record<string, unknown> = {}, pieces?: Sites): Node {
  const b = createSeqSpec({ title: 'variants' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(Math.ceil(length / 4)).slice(0, length) });
  const vs = b.data.ranges('vs', { target: s, rows });
  const pc = pieces ? b.data.ranges('pieces', { target: s, rows: pieces }) : undefined;
  b.views.view('main', (v) => {
    v.section('s', { of: s, ...(pc ? { pieces: pc } : {}), ...section })
      .track('vs')
      .layer('vs', { type: 'variants', data: vs });
  });
  return b.build();
}

/** The keys of the heads in a layout, with their centers. */
function heads(lay: VariantsLay) {
  const out: { key: string | number; x: number; y: number }[] = [];
  for (let j = 0; j < lay.m; j++) out.push({ key: lay.rows.keys[lay.mrow[j]], x: lay.mx[j], y: lay.my[j] });
  return out;
}

describe('variants: drawn against the data', () => {
  it('p53 at the fit view: a head for each variant, each head hit as its variant, three stacked at one site', async () => {
    const b = createSeqSpec({ title: 'p53 variants' });
    const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
    const vs = b.data.ranges('vs', { target: p53, rows: VariantRows });
    b.views.view('main', (v) => {
      v.section('protein', { of: p53 }).track('vs').layer('vs', { type: 'variants', data: vs });
    });
    const s = await settle(b.build());
    const f = await s.frame({ size: SIZE });
    const lay = layerLay<VariantsLay>(f, 0);
    const h = heads(lay);
    expect(h.map((x) => x.key).sort()).toEqual([...VariantRows.id].sort());
    for (const head of h) expect(hitView(f.layout, head.x, head.y)?.element?.key).toBe(head.key);
    // the heads of site 175 share a stem and stand one above the other, inside the track
    const site = h.filter((x) => ['v4', 'v5', 'v6'].includes(String(x.key)));
    expect(new Set(site.map((x) => x.x)).size).toBe(1);
    expect(new Set(site.map((x) => x.y)).size).toBe(3);
    const box = f.layout.sections[0].tracks[0].box;
    for (const x of h) {
      expect(x.y - 4).toBeGreaterThanOrEqual(box.y);
      expect(x.y + 4).toBeLessThanOrEqual(box.y + box.height);
    }
    // every head is a polygon and every stem a rect: 12 of each, and the track's own label
    expect(f.lists.passes.data.paths.n).toBe(12);
    expect(f.lists.passes.data.rects.n).toBe(12);
  });

  it('forty variants at one site are squeezed into the track, each still hit as itself', async () => {
    const s = await settle(spec(sites(new Array(40).fill(50)), 100));
    const f = await s.frame({ size: { width: 800, height: 80 } });
    const lay = layerLay<VariantsLay>(f, 0);
    const box = f.layout.sections[0].tracks[0].box;
    const h = heads(lay);
    expect(h.length).toBe(40);
    for (const x of h) {
      expect(x.y - 4).toBeGreaterThanOrEqual(box.y - 1e-6);
      expect(x.y + 4).toBeLessThanOrEqual(box.y + box.height);
      expect(hitView(f.layout, x.x, x.y)?.element?.key).toBe(x.key);
    }
    expect(new Set(h.map((x) => x.y)).size).toBe(40);
  });

  it('across a tile edge: a variant over the edge is one head, however the view meets it', async () => {
    const b = createSeqSpec({ title: 'tiles' });
    const src = b.data.source('synthetic', { format: 'synthetic', uri: 'synthetic:chr' });
    const chr = b.data.sequence('chr', { alphabet: 'dna', source: src });
    const items = b.data.ranges('items', { target: chr, source: src });
    b.views.view('main', (v) => {
      v.section('chr', { of: chr }).track('vs').layer('vs', { type: 'variants', data: items });
    });
    const s = await settle(b.build(), [BigRuntime]);
    for (const [start, end] of [
      [990, 1060],
      [1000, 1100],
      [1010, 1030],
    ]) {
      const f = await s.frame({ size: SIZE, viewport: { start, end } });
      const keys = heads(layerLay<VariantsLay>(f, 0)).map((x) => x.key);
      expect(new Set(keys).size).toBe(keys.length);
      const want: string[] = [];
      for (let p = 0; p < end; p += 50) if (p + 70 > start) want.push(`r${p}`);
      expect(keys.sort()).toEqual(want.sort());
    }
  });

  it('at the bars and density steps a strip counts the variants: no heads, one image', async () => {
    const b = createSeqSpec({ title: 'strip' });
    const src = b.data.source('synthetic', { format: 'synthetic', uri: 'synthetic:chr' });
    const chr = b.data.sequence('chr', { alphabet: 'dna', source: src });
    const items = b.data.ranges('items', { target: chr, source: src });
    b.views.view('main', (v) => {
      v.section('chr', { of: chr }).track('vs').layer('vs', { type: 'variants', data: items });
    });
    const s = await settle(b.build(), [BigRuntime]);
    for (const end of [150_000, 3_000_000]) {
      const f = await s.frame({ size: { width: 1440, height: 40 }, viewport: { start: 0, end } });
      const lay = layerLay<VariantsLay>(f, 0);
      expect(lay.density).toBeDefined();
      expect(lay.m).toBe(0);
      expect(f.lists.passes.data.images.n).toBe(1);
      expect(lay.density!.data.some((v) => v > 0)).toBe(true);
    }
  });
});

describe('variants: a reversed, collapsed axis and a codon map', () => {
  it('every variant in a shown piece has its head at its cell’s center, counting down on a reversed axis', async () => {
    const pieces = { start: [0, 40, 70], end: [20, 60, 90], id: ['a', 'b', 'c'], category: ['x', 'x', 'x'] };
    const positions = [5, 15, 19, 25, 45, 59, 65, 75, 89];
    for (const direction of ['genomic', 'reversed']) {
      const s = await settle(spec(sites(positions), 100, { direction }, pieces));
      const f = await s.frame({ size: SIZE });
      const lay = layerLay<VariantsLay>(f, 0);
      const h = heads(lay);
      // 25, 65 are in the pieces' gaps: the 7 others are shown
      expect(h.map((x) => x.key).sort()).toEqual(
        positions
          .map((p, i) => [p, `v${i}`] as const)
          .filter(([p]) => [5, 15, 19, 45, 59, 75, 89].includes(p))
          .map(([, k]) => k)
          .sort(),
      );
      for (const x of h) {
        const p = positions[Number(String(x.key).slice(1))];
        expect(lay.axis.cell(p)).toBe(true);
        expect(x.x).toBeCloseTo(lay.axis.snap((lay.axis.sx0 + lay.axis.sx1) / 2), 6);
        expect(lay.axis.positionAt(x.x)).toBe(p);
      }
      // left to right, positions rise forward and fall reversed
      const ps = [...h].sort((a, b) => a.x - b.x).map((x) => positions[Number(String(x.key).slice(1))]);
      const sorted = [...ps].sort((a, b) => (direction === 'genomic' ? a - b : b - a));
      expect(ps).toEqual(sorted);
    }
  });

  it('through the codon map, a protein’s variant is a head over its codon’s cell, on the gene', async () => {
    const b = createSeqSpec({ title: 'variants on the gene' });
    const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
    const gene = b.data.sequence('tp53', { alphabet: 'dna', letters: geneLetters() });
    b.data.map('p53-on-gene', { from: p53, to: gene, method: 'codon', exons: GENE_EXONS, strand: '-' });
    const vs = b.data.ranges('vs', { target: p53, rows: VariantRows });
    b.views.view('main', (v) => {
      const protein = v.section('protein', { of: p53 });
      protein.track('vs').layer('vs', { type: 'variants', data: vs });
      v.section('gene', { of: gene, follow: protein }).track('vs').layer('vs', { type: 'variants', data: vs });
    });
    const s = await settle(b.build());
    const f = await s.frame({ size: { width: 1400, height: 120 }, viewport: { start: 165, end: 185 } });
    const onProtein = heads(layerLay<VariantsLay>(f, 0, 0, 0));
    const onGene = heads(layerLay<VariantsLay>(f, 0, 0, 1));
    expect(onProtein.map((x) => x.key).sort()).toEqual(['v4', 'v5', 'v6']);
    expect(onGene.map((x) => x.key).sort()).toEqual(['v4', 'v5', 'v6']);
    // the gene reads along the protein: the heads stand in the same order, each over the codon of its residue
    const order = (hs: typeof onGene) => [...hs].sort((a, b) => a.x - b.x || a.y - b.y).map((x) => x.key);
    expect(order(onGene)).toEqual(order(onProtein));
    expect(f.layout.issues).toEqual([]);
  });
});

describe('variants: the fallback', () => {
  it('a target without the type draws ranges', async () => {
    const sp = spec(sites([5, 50]), 100);
    const f = await frameWithout(sp, ['variants'], { size: SIZE });
    expect(f.layout.issues).toEqual([]);
    expect(f.layout.sections[0].tracks[0].layers[0].type.name).toBe('ranges');
  });
});
