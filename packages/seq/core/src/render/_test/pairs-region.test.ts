/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Pairs read by region (_spec/README.md "Render", The pairs layer): zoomed out, bin pairs each an arc stroked by count
// and hit as its bin; zoomed in past the gate, the items themselves; every bin pair in view its arc (EC-208, EC-171).

import { PCG } from '@molstar/core/data/util/hash-functions';
import { categoryColumn } from '@molstar/core/repr/tree';
import { binPairs, type PairsList, type SeqData, SeqRuntime, seqToValue } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { hitView } from '../index.js';
import type { PairsLay } from '../pairs.js';
import { settle } from './helpers.js';

const Length = 1_000_000;
const Levels = [1, 1024, 16_384];

/** `n` links between 20 hotspots of a chromosome (each end jittered 2 kb), a level of bins at each coarse level. */
function linksSource(n: number): { data: SeqData; list: PairsList } {
  const rng = new PCG(9);
  const spots = Array.from({ length: 20 }, (_, i) => 20_000 + i * 48_000);
  const as = new Int32Array(n);
  const bs = new Int32Array(n);
  for (let i = 0; i < n; i++) {
    const a = spots[Math.floor(rng.float() * spots.length)];
    const b = spots[Math.floor(rng.float() * spots.length)];
    as[i] = a + Math.floor(rng.float() * 2000);
    bs[i] = b + Math.floor(rng.float() * 2000);
  }
  const list: PairsList = {
    form: 'list',
    aStart: as,
    aEnd: as.map((s) => s + 10),
    bStart: bs,
    bEnd: bs.map((s) => s + 10),
    id: categoryColumn(Array.from({ length: n }, (_, i) => `l${i}`)),
  };
  const whole = (level: number) => ({ system: 'chr', start: 0, end: Length, level, shapes: ['pairs' as const] });
  return {
    list,
    data: {
      coverage: {
        systems: [{ name: 'chr', length: Length }],
        shapes: [{ system: 'chr', shape: 'pairs' }],
        levels: Levels,
      },
      pairs: Levels.map((level) => ({
        request: whole(level),
        answer: { ok: true as const, tables: { pairs: level === 1 ? list : binPairs(list, level) } },
      })),
    },
  };
}

function spec(n: number) {
  const b = createSeqSpec({ title: 'links by region' });
  const src = b.data.source('src', { format: 'seq', value: seqToValue(linksSource(n).data) as never });
  const chr = b.data.sequence('chr', { alphabet: 'dna', source: src, name: 'chr' });
  const links = b.data.pairs('links', { a: chr, source: src } as never);
  b.views.view('main', (v) => {
    const s = v.section('chr', { of: chr });
    s.track('links').layer('links', { type: 'pairs', data: links });
  });
  return b.build();
}

const Size = { width: 720, height: 300 };
const layOf = (f: { layout: { sections: readonly { sectionLays: readonly { lay: unknown }[] }[] } }) =>
  f.layout.sections[0].sectionLays[0].lay as PairsLay;

describe('pairs by region', () => {
  it('zoomed out: bin pairs, each an arc stroked by its count and hit as its bin (Default 2)', async () => {
    const s = await settle(spec(6000), [SeqRuntime]);
    const f = await s.frame({ size: Size, viewport: { start: 0, end: Length } });
    const lay = layOf(f);
    const prep = lay.prep!;
    expect(prep.weight).toBeDefined();
    expect(prep.level).toBe(1024);
    // every bin pair in view has its arc (fewer than `maxArcs`), its count summed into the weights
    expect(lay.n).toBe(prep.data!.n);
    expect(Array.from(prep.weight!).reduce((a, b) => a + b, 0)).toBe(6000);
    const widths = Array.from(lay.width.subarray(0, lay.n));
    expect(Math.max(...widths)).toBeCloseTo(6);
    expect(Math.min(...widths)).toBeGreaterThanOrEqual(1);
    // the widest arc's middle hits a bin pair: its level said, its extent both bins
    const i = widths.indexOf(Math.max(...widths));
    const c = lay.curve.subarray(8 * i, 8 * i + 8);
    const x = (c[0] + 3 * c[2] + 3 * c[4] + c[6]) / 8;
    const y = (c[1] + 3 * c[3] + 3 * c[5] + c[7]) / 8;
    const h = hitView(f.layout, x, y);
    expect(h?.bin).toBe(1024);
    expect(String(h?.element?.key)).toMatch(/^bin:/);
  });

  it('zoomed in under the gate: the items themselves, with no bins (Default 1)', async () => {
    const s = await settle(spec(6000), [SeqRuntime]);
    // three hotspots in view: the links with an end there fit `maxArcs`
    const f = await s.frame({ size: Size, viewport: { start: 0, end: 400_000 } });
    const lay = layOf(f);
    expect(lay.prep!.weight).toBeUndefined();
    const { list } = linksSource(6000);
    let inView = 0;
    for (let i = 0; i < list.aStart.length; i++) if (list.aStart[i] < 400_000 || list.bStart[i] < 400_000) inView++;
    expect(inView).toBeLessThanOrEqual(4096);
    expect(lay.prep!.data!.n).toBe(inView);
  });
});
