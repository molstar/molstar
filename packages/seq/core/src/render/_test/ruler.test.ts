/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ruler layer (_spec/README.md "Render", Layer types): labels at steps of 1, 2, 5 × 10ⁿ that never overlap and are
// each a position of the numbering, insertion labels as minor ticks with room only, a reversed axis counting down, the
// numbering a spec param, the unknown one reported and the ordinal one used; `units: bases` in bp, kb or Mb by the step,
// the canvas's labels the export's.

import { SeqRuntime, seqToValue } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import type { Numbering } from '@molstar/seqx-core/query';
import { describe, expect, it } from 'vitest';
import {
  AxisTransform,
  bundledTextMetrics,
  DefaultRenderConfig,
  LightTokens,
  type LayoutContext,
  labelStep,
  type RulerLay,
  RulerType,
  renderConfig,
  renderSvg,
} from '../index.js';
import { GENE_EXONS, GENE_LENGTH, geneLetters, P53 } from './fixtures.js';
import { settle } from './helpers.js';
import { layerLay, rng } from './small-helpers.js';

const GAP = DefaultRenderConfig.ruler.labelGap;

/** The labels drawn, with their x, in the order of x. */
function shown(lay: RulerLay): { label: string; x: number; w: number; major: boolean }[] {
  const out: { label: string; x: number; w: number; major: boolean }[] = [];
  for (let j = 0; j < lay.n; j++)
    if (lay.shown[j]) out.push({ label: lay.labels[j], x: lay.x[j], w: lay.widths[j], major: lay.major[j] === 1 });
  return out.sort((a, b) => a.x - b.x);
}

function p53Ruler(numbering?: string) {
  const b = createSeqSpec({ title: 'ruler' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  b.views.view('main', (v) => {
    v.section('protein', { of: p53 })
      .track('r')
      .layer('r', { type: 'ruler', ...(numbering ? { numbering } : {}) });
  });
  return b.build();
}

describe('the ruler: steps', () => {
  it('a step is 1, 2 or 5 × 10ⁿ, the least whose ticks are a label gap apart', () => {
    const r = rng(5);
    const steps = new Set<number>();
    for (let i = 0; i < 2000; i++) {
      const scale = 10 ** (r() * 6 - 4);
      const s = labelStep(scale, GAP);
      const e = Math.floor(Math.log10(s));
      expect([1, 2, 5]).toContain(Math.round(s / 10 ** e));
      expect(s * scale).toBeGreaterThanOrEqual(GAP * (1 - 1e-9));
      // the step below it (the previous of 1, 2, 5, 10) is too close
      const below = s === 10 ** e ? 5 * 10 ** (e - 1) : s === 2 * 10 ** e ? 10 ** e : 2 * 10 ** e;
      expect(below * scale).toBeLessThan(GAP);
      steps.add(s);
    }
    expect(steps.size).toBeGreaterThan(10);
  });
});

describe('the ruler: labels', () => {
  it('never overlap, each is the label of the position under its tick, over seeded views and widths', async () => {
    const s = await settle(p53Ruler());
    const r = rng(11);
    let labels = 0;
    for (let i = 0; i < 60; i++) {
      const w = 5 + r() * 390;
      const start = r() * (393 - w);
      const size = { width: 300 + Math.floor(r() * 1000), height: 60 };
      const f = await s.frame({ size, viewport: { start, end: start + w }, dpr: [1, 2][i & 1] });
      const lay = layerLay<RulerLay>(f, 0);
      const out = shown(lay);
      labels += out.length;
      for (let k = 0; k < out.length; k++) {
        const p = lay.axis.positionAt(out[k].x);
        // a label is its position's number, and its tick is in the data box
        expect(out[k].label).toBe(String(p + 1));
        expect(out[k].x).toBeGreaterThanOrEqual(lay.box.x);
        expect(out[k].x).toBeLessThanOrEqual(lay.box.x + lay.box.width);
        if (k) expect(out[k].x - out[k].w / 2 - (out[k - 1].x + out[k - 1].w / 2)).toBeGreaterThanOrEqual(4);
      }
    }
    expect(labels).toBeGreaterThan(100);
  });

  it('at the fit view the labels are multiples of a step, drawn once each', async () => {
    const s = await settle(p53Ruler());
    const f = await s.frame({ size: { width: 700, height: 60 } });
    const out = shown(layerLay<RulerLay>(f, 0));
    expect(out.map((l) => l.label)).toEqual(['50', '100', '150', '200', '250', '300', '350']);
    expect(f.layout.issues).toEqual([]);
  });

  it('`numbering: index` labels from 0', async () => {
    const s = await settle(p53Ruler('index'));
    const f = await s.frame({ size: { width: 700, height: 60 }, viewport: { start: 0, end: 30 } });
    expect(shown(layerLay<RulerLay>(f, 0)).map((l) => l.label)).toEqual(['0', '5', '10', '15', '20', '25']);
  });

  it('a numbering the system lacks is reported once and the ordinal one is used', async () => {
    const s = await settle(p53Ruler('kabat'));
    const f = await s.frame({ size: { width: 700, height: 60 } });
    expect(f.layout.issues.map((i) => i.code)).toEqual(['seq.render.numbering']);
    expect(shown(layerLay<RulerLay>(f, 0)).map((l) => l.label)[0]).toBe('50');
  });

  it('draws a tick a label and the labels in the labels pass; a hit is the position under the pointer', async () => {
    const s = await settle(p53Ruler());
    const f = await s.frame({ size: { width: 700, height: 60 } });
    const lay = layerLay<RulerLay>(f, 0);
    expect(f.lists.passes.data.rects.n).toBe(lay.n);
    // the labels, and the track's own label
    expect(f.lists.passes.labels.texts.n).toBe(shown(lay).length + 1);
    const x = lay.x[0];
    const h = RulerType.hit(lay, Float64Array.of(x, lay.box.y + 2));
    expect(h?.position).toBe(lay.axis.positionAt(x));
  });
});

/** A chromosome-long DNA system (its length a `seq` source's coverage, no letters) under a ruler. */
function chromosomeRuler(params: Record<string, unknown>) {
  const b = createSeqSpec({ title: 'ruler' });
  const value = seqToValue({
    coverage: { systems: [{ name: 'chr', length: 83_257_441 }], shapes: [], levels: [1] },
    pairs: [],
  });
  const source = b.data.source('chr-source', { format: 'seq', value: value as never });
  const chr = b.data.sequence('chr', { alphabet: 'dna', source, name: 'chr' });
  b.views.view('main', (v) => {
    v.section('chr', { of: chr })
      .track('r')
      .layer('r', { type: 'ruler', ...params });
  });
  return b.build();
}

describe('the ruler: units', () => {
  const size = { width: 700, height: 60 };
  const labelsAt = async (params: Record<string, unknown>, viewport: { start: number; end: number }) => {
    const spec = chromosomeRuler(params);
    const s = await settle(spec, [SeqRuntime]);
    const canvas = shown(layerLay<RulerLay>(await s.frame({ size, viewport }), 0)).map((l) => l.label);
    // the export's SVG, a second consumer of the same spec: its text the canvas's labels
    const r = await renderSvg(spec, { size, viewport, runtimes: [SeqRuntime], legend: 'none' });
    const texts = [...r.svg.matchAll(/<text[^>]*>([^<]*)<\/text>/g)].map((m) => m[1]);
    for (const l of canvas) expect(texts).toContain(l);
    return canvas;
  };

  it('`bases` labels an ordinal ruler in bp, kb or Mb by its step; the canvas’s labels are the SVG’s', async () => {
    const fit = await labelsAt({ units: 'bases' }, { start: 0, end: 83_257_441 });
    expect(fit.slice(0, 2)).toEqual(['10 Mb', '20 Mb']);
    const mb = await labelsAt({ units: 'bases' }, { start: 7_600_000, end: 7_700_000 });
    expect(mb.length).toBeGreaterThan(2);
    for (const l of mb) expect(l).toMatch(/^7\.\d\d Mb$/);
    const kb = await labelsAt({ units: 'bases' }, { start: 7_668_000, end: 7_678_000 });
    for (const l of kb) expect(l).toMatch(/^7,6\d\d kb$/);
    const bp = await labelsAt({ units: 'bases' }, { start: 7_668_400, end: 7_668_500 });
    expect(bp.length).toBeGreaterThan(2);
    for (const l of bp) expect(l).toMatch(/^7,668,4\d\d$/);
  });

  it('`none` (and absent) labels the number; another numbering than ordinal ignores `bases`', async () => {
    const at = { start: 7_600_000, end: 7_700_000 };
    const plain = await labelsAt({}, at);
    expect(plain[0]).toMatch(/^76\d0000$/);
    expect(await labelsAt({ units: 'none' }, at)).toEqual(plain);
    const index = await labelsAt({ units: 'bases', numbering: 'index' }, at);
    for (const l of index) expect(l).toMatch(/^\d+$/);
  });
});

describe('the ruler: a reversed, collapsed axis', () => {
  function gene() {
    const b = createSeqSpec({ title: 'reversed ruler' });
    const g = b.data.sequence('tp53', { alphabet: 'dna', letters: geneLetters() });
    const exons = b.data.ranges('exons', { target: g, rows: GENE_EXONS });
    b.views.view('main', (v) => {
      v.section('gene', { of: g, direction: 'reversed', pieces: exons }).track('r').layer('r', { type: 'ruler' });
    });
    return b.build();
  }

  it('counts down on the screen, every label a position of a piece in view', async () => {
    const s = await settle(gene());
    const r = rng(3);
    let checked = 0;
    for (let i = 0; i < 20; i++) {
      const total = GENE_EXONS.end.reduce((a, e, k) => a + e - GENE_EXONS.start[k], 0);
      const w = 40 + r() * (total - 40);
      const start = r() * (total - w);
      const f = await s.frame({ size: { width: 800, height: 60 }, viewport: { start, end: start + w } });
      const lay = layerLay<RulerLay>(f, 0);
      const out = shown(lay);
      const nums = out.map((l) => Number(l.label));
      for (let k = 1; k < nums.length; k++) expect(nums[k]).toBeLessThan(nums[k - 1]);
      for (const l of out) {
        const p = lay.axis.positionAt(l.x);
        expect(p).toBeGreaterThanOrEqual(0);
        expect(l.label).toBe(String(p + 1));
        // the position is inside an exon
        expect(GENE_EXONS.start.some((st, k) => p >= st && p < GENE_EXONS.end[k])).toBe(true);
      }
      checked += out.length;
    }
    expect(checked).toBeGreaterThan(30);
    expect(GENE_LENGTH).toBeGreaterThan(0);
  });
});

/** A numbering with insertion labels: 1 to 120, 100A and 100B after 100, 111.1 and 111.2 after 111. */
function kabat(): Numbering {
  const labels: string[] = [];
  for (let n = 1; n <= 120; n++) {
    labels.push(String(n));
    if (n === 100) labels.push('100A', '100B');
    if (n === 111) labels.push('111.1', '111.2');
  }
  return {
    id: 'kabat',
    length: labels.length,
    labelAt: (i) => labels[i],
    indexOf: (l) => labels.indexOf(l),
  };
}

/** The ruler over `kabat` laid out by hand (the language declares no numbering yet), the view `[start, end)` of 700 px. */
function kabatRuler(start: number, end: number, width = 700) {
  const numbering = kabat();
  const axis = new AxisTransform('chain', [0, numbering.length], 0, width, 1, { start, end });
  const reports: string[] = [];
  const ctx = {
    box: { x: 0, y: 0, width, height: 22 },
    axis,
    section: { address: '/views/main/chain' },
    key: 'k',
    dpr: 1,
    config: renderConfig(),
    tokens: LightTokens,
    text: bundledTextMetrics(),
    settled: true,
    want: () => {},
    report: (i: { code: string }) => reports.push(i.code),
  } as unknown as LayoutContext;
  const lay = RulerType.layout(
    { id: 'kabat', numbering: { numbering, shift: 0 } },
    ctx,
    { coverage: () => undefined, keysOf: () => [], get: () => undefined },
    undefined as never,
  );
  return { lay, numbering, reports };
}

describe('the ruler: insertion labels', () => {
  it('are minor ticks, labeled only where there is room, never overlapping the majors', () => {
    // 7 positions across 700 px: room for every label
    const near = shown(kabatRuler(98, 105).lay).map((l) => l.label);
    expect(near).toEqual(expect.arrayContaining(['100', '100A', '100B']));
    // 30 positions: the minors are ticks, their labels do not fit beside 100
    const far = kabatRuler(90, 120).lay;
    const labels = shown(far).map((l) => l.label);
    expect(labels).toContain('100');
    expect(labels).not.toContain('100A');
    const minors = Array.from({ length: far.n }, (_, j) => far.major[j] === 0 && far.labels[j]).filter(Boolean);
    expect(minors).toEqual(expect.arrayContaining(['100A', '100B', '111.1', '111.2']));
    // beyond the minors' scale: no minor ticks at all
    const wide = kabatRuler(0, 126, 400).lay;
    expect(Array.from({ length: wide.n }, (_, j) => wide.major[j])).not.toContain(0);
  });

  it('every label shown, of any view, is a label of the numbering at its tick, and none overlap', () => {
    const r = rng(9);
    const numbering = kabat();
    for (let i = 0; i < 100; i++) {
      const w = 3 + r() * 100;
      const start = r() * (numbering.length - w);
      const { lay } = kabatRuler(start, start + w);
      const out = shown(lay);
      for (let k = 0; k < out.length; k++) {
        expect(numbering.labelAt(lay.axis.positionAt(out[k].x))).toBe(out[k].label);
        if (k) expect(out[k].x - out[k].w / 2 - (out[k - 1].x + out[k - 1].w / 2)).toBeGreaterThanOrEqual(4);
      }
    }
  });
});
