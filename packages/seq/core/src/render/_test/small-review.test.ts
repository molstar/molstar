/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The small layers' review: an arrowhead only where the item's whole extent is known, items inside a hidden stretch
// kept as the break's tick, the secondary classes by whole codes, a variant's stem hit, a ruler's names, a start
// triangle not under the dash, a line's label legible (_spec/README.md "Render", Layer types).

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import type { Numbering } from '@molstar/seqx-core/query';
import { describe, expect, it } from 'vitest';
import {
  AxisTransform,
  bundledTextMetrics,
  hitView,
  LightTokens,
  type LayoutContext,
  type RangesLay,
  type RulerLay,
  RulerType,
  renderConfig,
  type SecondaryLay,
  secondaryKind,
  type VariantsLay,
} from '../index.js';
import { BigRuntime } from './fixtures.js';
import { settle } from './helpers.js';
import { layerLay, readPaths, readRects } from './small-helpers.js';

const SIZE = { width: 800, height: 60 } as const;

type Rows = { start: number[]; end: number[]; id: string[]; category?: string[] };

/** A sequence of `length` letters, the rows as a layer of `type` with its params, in a section. */
function spec(
  length: number,
  rows: Rows,
  layer: Record<string, unknown>,
  options: { cyclic?: boolean; section?: Record<string, unknown>; pieces?: { start: number[]; end: number[] } } = {},
): Node {
  const b = createSeqSpec({ title: 'review' });
  const s = b.data.sequence('s', {
    alphabet: 'dna',
    letters: 'ACGT'.repeat(length / 4),
    ...(options.cyclic ? { cyclic: true } : {}),
  });
  const items = b.data.ranges('items', { target: s, rows: rows as never });
  const pieces = options.pieces ? b.data.ranges('pieces', { target: s, rows: options.pieces }) : undefined;
  b.views.view('main', (v) => {
    v.section('s', { of: s, ...(pieces ? { pieces } : {}), ...options.section })
      .track('t', { height: 36 })
      .layer('l', { data: items, ...layer } as never);
  });
  return b.build();
}

describe('an end is marked only when the item’s whole extent is known', () => {
  it('secondary: a strand across the origin of 20,000, the view holding only the row before it: no head', async () => {
    const rows: Rows = { start: [19990, 0], end: [20000, 5], id: ['x', 'x'], category: ['strand', 'strand'] };
    const s = await settle(spec(20000, rows, { type: 'secondary' }, { cyclic: true }));
    const before = layerLay<SecondaryLay>(await s.frame({ size: SIZE, viewport: { start: 19950, end: 20000 } }), 0);
    expect(before.m).toBe(1);
    expect(before.mtip[0]).toBe(0);
    // the view across the origin holds both rows: the head is on the one after it
    const w = await settle(spec(20000, rows, { type: 'secondary' }, { cyclic: true, section: { range: [19950, 20] } }));
    const both = layerLay<SecondaryLay>(await w.frame({ size: SIZE }), 0);
    const tipped = Array.from({ length: both.m }, (_, j) => [both.rows.starts[both.mitem[j]], both.mtip[j]]);
    expect(tipped.sort()).toEqual([
      [0, 1],
      [19990, 0],
    ]);
  });

  it('ranges: an item of two rows far apart, the view holding the first: no head; the second: its head', async () => {
    const rows: Rows = { start: [100, 15000], end: [200, 15100], id: ['a', 'a'] };
    const s = await settle(spec(20000, rows, { type: 'ranges', arrow: 'end' }));
    const first = layerLay<RangesLay>(await s.frame({ size: SIZE, viewport: { start: 0, end: 300 } }), 0);
    expect(Array.from({ length: first.m }, (_, j) => first.mtip[j])).toEqual([0]);
    const last = layerLay<RangesLay>(await s.frame({ size: SIZE, viewport: { start: 14900, end: 15200 } }), 0);
    expect(Array.from({ length: last.m }, (_, j) => last.mtip[j])).toEqual([2]);
  });

  it('a tiled source gives no head (its items’ extents are not known)', async () => {
    const b = createSeqSpec({ title: 'tiled' });
    const src = b.data.source('synthetic', { format: 'synthetic', uri: 'synthetic:chr' });
    const chr = b.data.sequence('chr', { alphabet: 'dna', source: src });
    const items = b.data.ranges('items', { target: chr, source: src });
    b.views.view('main', (v) => {
      v.section('chr', { of: chr }).track('t').layer('l', { type: 'ranges', data: items, arrow: 'end' });
    });
    const s = await settle(b.build(), [BigRuntime]);
    const lay = layerLay<RangesLay>(await s.frame({ size: SIZE, viewport: { start: 990, end: 1100 } }), 0);
    expect(lay.m).toBeGreaterThan(0);
    for (let j = 0; j < lay.m; j++) expect(lay.mtip[j]).toBe(0);
  });
});

describe('items inside a hidden stretch are not lost: the break’s tick', () => {
  const pieces = { start: [0, 40], end: [10, 50] };
  const hidden: Rows = { start: [20, 22], end: [21, 23], id: ['a', 'b'], category: ['strand', 'helix'] };

  for (const type of ['secondary', 'variants'] as const)
    it(`${type}: a tick for the items, hit as the break, with their count`, async () => {
      const s = await settle(spec(100, hidden, { type }, { pieces }));
      const f = await s.frame({ size: SIZE });
      const lay = layerLay<SecondaryLay | VariantsLay>(f, 0);
      expect(lay.tn).toBe(1);
      expect(lay.ticks[2]).toBe(2);
      const box = f.layout.sections[0].tracks[0].box;
      const h = hitView(f.layout, lay.ticks[0], box.y + box.height / 2);
      expect(h?.target).toBe('break');
      expect(h?.count).toBe(2);
      expect(readRects(f.lists).length).toBeGreaterThan(0);
    });
});

describe('secondary classes are whole codes and words', () => {
  it('DSSP, mmCIF and words; a prefix is not a class', () => {
    for (const k of ['H', 'G', 'I', 'helix', 'HELX_P', 'HELX_RH_AL_P']) expect(secondaryKind(k), k).toBe(0);
    for (const k of ['E', 'B', 'strand', 'STRN', 'Strand']) expect(secondaryKind(k), k).toBe(1);
    for (const k of ['T', 'S', 'turn', 'TURN_P', 'BEND', 'bend', 'break', 'in', 'helices', 'coil', 'loop'])
      expect(secondaryKind(k), k).toBe(2);
    expect(secondaryKind(undefined)).toBe(2);
  });

  it('UniProt’s words: Helix, Beta strand, Turn', () => {
    expect(['Helix', 'Beta strand', 'Turn', 'Beta'].map(secondaryKind)).toEqual([0, 1, 2, 2]);
  });
});

describe('variants: a stem is hit', () => {
  it('below its head, as the variant', async () => {
    const rows: Rows = { start: [50], end: [51], id: ['v'], category: ['missense'] };
    const s = await settle(spec(100, rows, { type: 'variants' }));
    const f = await s.frame({ size: SIZE });
    const lay = layerLay<VariantsLay>(f, 0);
    const box = f.layout.sections[0].tracks[0].box;
    const x = lay.mx[0];
    const under = lay.my[0] + 4 + 2;
    expect(under).toBeLessThan(box.y + box.height);
    expect(hitView(f.layout, x, under)?.element?.key).toBe('v');
    expect(hitView(f.layout, x + 3, under)?.element).toBeUndefined();
  });
});

describe('the ruler labels names by the label', () => {
  it('a numbering of names: every name a tick, labeled where it fits', () => {
    const names = ['alpha', 'beta', 'gamma', 'delta', 'eps', 'zeta', 'eta', 'theta'];
    const numbering: Numbering = {
      id: 'names',
      length: names.length,
      labelAt: (i) => names[i],
      indexOf: (l) => names.indexOf(l),
    };
    const axis = new AxisTransform('c', [0, names.length], 0, 200, 1);
    const ctx = {
      box: { x: 0, y: 0, width: 200, height: 22 },
      axis,
      section: { address: '/views/main/c' },
      key: 'k',
      dpr: 1,
      config: renderConfig(),
      tokens: LightTokens,
      text: bundledTextMetrics(),
      settled: true,
      want: () => {},
      report: () => {},
    } as unknown as LayoutContext;
    const lay: RulerLay = RulerType.layout(
      { id: 'names', numbering: { numbering, shift: 0 } },
      ctx,
      { coverage: () => undefined, keysOf: () => [], get: () => undefined },
      undefined as never,
    );
    // 25 px a position: a step of 2, yet every name is a tick (a name has no number to step by)
    const shown = Array.from({ length: lay.n }, (_, j) => lay.labels[j]);
    expect(shown).toEqual(names);
  });
});

describe('ranges: lines and arrows', () => {
  const ONE: Rows = { start: [40], end: [300], id: ['a'] };

  it('a start triangle: the dash starts at the line, not under the triangle', async () => {
    const s = await settle(spec(400, ONE, { type: 'ranges', line: 'dashed', arrow: 'start' }));
    const f = await s.frame({ size: SIZE });
    const lay = layerLay<RangesLay>(f, 0);
    const p = readPaths(f.lists).find((q) => q.dash > 0)!;
    const triangle = readPaths(f.lists).find((q) => q.dash === 0)!;
    const base = Math.max(...triangle.xy.filter((_, i) => i % 2 === 0));
    expect(base).toBeGreaterThan(lay.mx0[0]);
    expect(p.xy[0]).toBeGreaterThanOrEqual(base - 1e-3);
  });

  it('the label of a line stays legible: a knockout under it in the page color, none under a block', async () => {
    const labels = async (extra: Record<string, unknown>) => {
      const s = await settle(spec(400, ONE, { type: 'ranges', ...extra }));
      const f = await s.frame({ size: { width: 800, height: 60 } });
      const rects = readRects(f.lists, 'labels');
      return { f, rects };
    };
    const line = await labels({ line: 'solid' });
    expect(line.f.lists.passes.labels.texts.n).toBe(2);
    expect(line.rects.length).toBe(1);
    expect(line.rects[0].fill).toBe(LightTokens.colors['msui-bg']);
    const block = await labels({});
    expect(block.rects.length).toBe(0);
  });
});
