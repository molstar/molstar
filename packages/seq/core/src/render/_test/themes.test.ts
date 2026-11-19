/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Themes drawn (_spec/README.md "Themes"; "Render", Layer types): drawing checked against the data (each mark's color
// is its theme's for its datum, each painted legend entry has a mark, each drawn color is in the legend), the size
// channel never moving a hit box, a shared scale's open domain unioned over its layers at settle, the export's legend
// band and scale records, a scheme that does not fit, the reader's substitution in an export; and the performance
// representatives: a layout with a field theme at most 1.1× the same layout uniform (the two timed in turn), legends
// for 100 tracks under 1 ms.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { compileChannel, legendOf, mergeLegends } from '@molstar/seqx-core/theme/index';
import { bestOf, itWithin, ratioWithin } from '@molstar/test-harness/timing';
import { describe, expect, it } from 'vitest';
import {
  DefaultThemeContext,
  type DrawLists,
  hitView,
  LightTokens,
  type RangesLay,
  recordDrawLists,
  type ValuesLay,
  viewLegends,
} from '../index.js';
import { renderSvg } from '../render-svg.js';
import { P53, p53TracksSpec } from './fixtures.js';
import { settle } from './helpers.js';

const Size = { width: 800, height: 160 };

function domainsSpec(layer: Record<string, unknown> = {}): Node {
  const b = createSeqSpec();
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const d = b.data.ranges('domains', {
    target: p53,
    rows: {
      start: new Int32Array([0, 63, 93, 318, 323, 363]),
      end: new Int32Array([61, 94, 292, 360, 355, 393]),
      category: ['activation', 'linker', 'binding', 'signal', 'oligomer', 'regulatory'],
      score: new Float32Array([0.1, 0.5, 0.9, 0.3, 0.7, 0.2]),
    },
  });
  b.views.view('main', (v) =>
    v
      .section('p', { of: p53 })
      .track('domains', { label: 'Domains' })
      .layer('d', { type: 'ranges', data: d, ...layer }),
  );
  return b.build();
}

/** The fills of the data pass's rects inside the data box (x past the label column), from the recording. */
function fills(lists: DrawLists): number[] {
  return recordDrawLists(lists, ['data'])
    .split('\n')
    .filter((l) => l.startsWith('rect ') && Number(l.split(' ')[1]) >= 100)
    .map((l) => Number.parseInt(l.split(' ')[5].slice(1), 16) >>> 0);
}

describe('drawing checked against the data', () => {
  it('each mark’s color is its theme’s for its datum; each painted legend entry has a mark; each drawn color is in the legend', async () => {
    const s = await settle(domainsSpec());
    const f = await s.frame({ size: Size, settled: true });
    const layer = f.layout.sections[0].tracks[0].layers[0];
    const lay = layer.lay as RangesLay;
    const channel = compileChannel('color', layer.layer.theme.color as never, 'token:range', {
      tokens: LightTokens,
      schemes: DefaultThemeContext.schemes,
    });
    const cats = ['activation', 'linker', 'binding', 'signal', 'oligomer', 'regulatory'];
    for (let i = 0; i < lay.items; i++) {
      const key = lay.keys[i];
      const start = lay.starts[i];
      const row = [0, 63, 93, 318, 323, 363].indexOf(start);
      expect(lay.fills[i], String(key)).toBe(
        channel.at({ value: Number.NaN, letter: '', start, end: 0, field: () => cats[row] }),
      );
    }
    const { legends } = viewLegends(f.layout);
    expect(legends).toHaveLength(1);
    const entries = legends[0].entries;
    if (entries.kind !== 'swatches') throw new Error('swatches');
    const drawn = new Set(fills(f.lists));
    for (const e of entries.items) expect(drawn.has(e.color), e.label).toBe(true);
    const legendColors = new Set(entries.items.map((e) => e.color));
    for (let i = 0; i < lay.items; i++) expect(legendColors.has(lay.fills[i])).toBe(true);
    expect(legends[0]).toMatchObject({ title: 'Domains', source: 'vocabulary', notes: [] });
  });
});

describe('channels', () => {
  it('the size channel shrinks a mark, never its lane or its hit box', async () => {
    const plain = await settle(domainsSpec());
    const sized = await settle(domainsSpec({ size: theme.field('score', { domain: [0, 1], range: [0.2, 1] }) }));
    const a = await plain.frame({ size: Size, settled: true });
    const b = await sized.frame({ size: Size, settled: true });
    for (let x = 130; x < 790; x += 7)
      for (let y = 20; y < 140; y += 3) {
        const ha = hitView(a.layout, x, y);
        const hb = hitView(b.layout, x, y);
        expect(hb?.element, `${x},${y}`).toEqual(ha?.element);
      }
    const lb = b.layout.sections[0].tracks[0].layers[0].lay as RangesLay;
    expect(lb.sizes.slice(0, lb.items).some((z) => z < 1)).toBe(true);
  });

  it('opacity scales a fill’s alpha; a scheme that does not fit draws the default and says so', async () => {
    const r = await renderSvg(domainsSpec({ opacity: 0.5, color: theme.scheme('zappo') }), { size: Size });
    expect(r.issues.map((i) => i.code)).toContain('seq.render.theme.scheme');
    expect(r.svg).toContain('fill-opacity');
  });
});

describe('a shared scale', () => {
  function shared(): Node {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'M'.repeat(20) });
    const v1 = b.data.values('a', { target: s, values: Float32Array.from({ length: 20 }, (_, i) => i) });
    const v2 = b.data.values('b', { target: s, values: Float32Array.from({ length: 20 }, (_, i) => 50 + i) });
    const sc = b.views.scale('cons', { label: 'Conservation' });
    b.views.view('main', (v) => {
      const sec = v.section('p', { of: s });
      sec.track('a').layer('a', { type: 'values', data: v1, scale: sc, color: theme.field('value', { scale: sc }) });
      sec.track('b').layer('b', { type: 'values', data: v2, scale: sc, color: theme.field('value', { scale: sc }) });
    });
    return b.build();
  }

  it('unions an open domain over the layers naming it at settle (y and color); one legend for the scale', async () => {
    const s = await settle(shared());
    const f = await s.frame({ size: Size, settled: true });
    const lays = f.layout.sections[0].tracks.map((t) => t.layers[0].lay as ValuesLay);
    expect(lays[0].domain).toEqual([0, 69]);
    expect(lays[1].domain).toEqual([0, 69]);
    // the same value has one color in both
    const { legends } = viewLegends(f.layout);
    expect(legends).toHaveLength(1);
    expect(legends[0]).toMatchObject({ scale: '/views/cons', title: 'Conservation', open: true, domain: [0, 69] });
  });

  it('an export records the scale’s domain as it drew it', async () => {
    const r = await renderSvg(shared(), { size: Size });
    expect(r.svg).toContain('data-seq-scale-cons="0,69"');
  });
});

describe('the export’s legends', () => {
  it('none by default; below or right when asked or when the view says so, outside the data box', async () => {
    const none = await renderSvg(p53TracksSpec(), { size: Size });
    const below = await renderSvg(p53TracksSpec(), { size: Size, legend: 'below' });
    const right = await renderSvg(p53TracksSpec(), { size: Size, legend: 'right' });
    expect(none.lists.height).toBe(Size.height);
    expect(below.lists.height).toBeGreaterThan(Size.height);
    expect(below.svg).toContain('Domains');
    expect(below.svg).toContain('activation');
    expect(right.lists.width).toBeGreaterThan(Size.width);
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: P53 });
    b.views.view('main', { export: { legend: 'below' } }, (v) =>
      v.section('p', { of: s }).track('t').layer('l', { type: 'letters' }),
    );
    expect((await renderSvg(b.build(), { size: Size })).lists.height).toBeGreaterThan(Size.height);
    expect((await renderSvg(b.build(), { size: Size, legend: 'none' })).lists.height).toBe(Size.height);
  });

  it('draws the reader’s substitution only when asked', async () => {
    const zappo = () => {
      const b = createSeqSpec();
      const s = b.data.sequence('s', { alphabet: 'protein', letters: 'CCCC' });
      b.views.view('main', (v) =>
        v
          .section('p', { of: s })
          .track('t')
          .layer('l', { type: 'letters', color: theme.scheme('zappo') }),
      );
      return b.build();
    };
    const plain = await renderSvg(zappo(), { size: Size });
    const safe = await renderSvg(zappo(), { size: Size, cvdSafe: true });
    expect(plain.svg).toContain('#ffff00');
    expect(safe.svg).not.toContain('#ffff00');
    expect(safe.svg).toContain('#56b4e9');
  });

  it('a substituted scheme’s view draws `classes-safe`’s group of each letter, in both themes', async () => {
    const groups: [string, number][] = [
      ['AVLIMC', 0x56b4e9ff],
      ['FWY', 0xcc79a7ff],
      ['STNQ', 0x009e73ff],
      ['KRH', 0x0072b2ff],
      ['DE', 0xd55e00ff],
      ['GP', 0xf0e442ff],
    ];
    const letters = groups.map((g) => g[0]).join('');
    const expected = [...letters].map((l) => groups.find((g) => g[0].includes(l))![1]);
    for (const name of ['classes', 'zappo', 'clustal-x', 'taylor'])
      for (const theme of ['light', 'dark'] as const) {
        const b = createSeqSpec();
        const s = b.data.sequence('s', { alphabet: 'protein', letters });
        b.views.view('main', (v) =>
          v
            .section('p', { of: s })
            .track('t')
            .layer('l', { type: 'letters', color: { kind: 'scheme', params: { name } } }),
        );
        const { lists } = await renderSvg(b.build(), { size: Size, theme, cvdSafe: true, legend: 'none' });
        const cells = recordDrawLists(lists, ['data'])
          .split('\n')
          .filter((l) => l.startsWith('rect '))
          .map((l) => l.split(' '))
          .map((p) => ({ x: Number(p[1]), fill: Number.parseInt(p[5].slice(1), 16) >>> 0 }))
          .filter((c) => expected.includes(c.fill))
          .sort((a, b) => a.x - b.x)
          .map((c) => c.fill);
        expect(cells, `${name} ${theme}`).toEqual(expected);
      }
  });
});

describe('performance representatives', () => {
  const N = 10_000;
  function many(color?: Node): Node {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(N * 2) });
    const start = Int32Array.from({ length: N }, (_, i) => i * 8);
    const r = b.data.ranges('r', {
      target: s,
      rows: {
        start,
        end: start.map((x) => x + 6),
        category: Array.from({ length: N }, (_, i) => `c${i % 8}`),
      },
    });
    b.views.view('main', (v) =>
      v
        .section('p', { of: s })
        .track('t', { height: 200 })
        .layer('l', { type: 'ranges', data: r, color: color ?? '#7a8fa6' }),
    );
    return b.build();
  }

  itWithin(
    'a layout with a field theme costs at most 1.1× the same layout uniform (10,000 ranges by a category)',
    async (within) => {
      const uniform = await settle(many());
      const themed = await settle(many(theme.field('category')));
      // a viewport of its own each frame: every frame lays out
      const at = (k: number) => ({ size: { width: 1600, height: 240 }, viewport: { start: k, end: N * 8 + k } });
      for (let k = 0; k < 6; k++) {
        await uniform.frame(at(k));
        await themed.frame(at(k));
      }
      let k = 10;
      // the two measured in turn, the themed held to 1.1× the uniform, each its best of seven
      let f: ReturnType<typeof themed.renderer.frame> | undefined;
      ratioWithin(
        within,
        1.1,
        () => uniform.renderer.frame(at(k++)),
        () => {
          f = themed.renderer.frame(at(k++));
        },
        7,
      );
      expect(f!.rebuilt).toBe(true);
      expect((f!.layout.sections[0].tracks[0].layers[0].lay as RangesLay).items).toBeGreaterThan(1000);
    },
  );

  itWithin('legends for 100 tracks are derived in under 1 ms', (within) => {
    const ch = compileChannel('color', theme.field('kind', { list: false }) as never, 'token:range', {
      tokens: LightTokens,
      schemes: DefaultThemeContext.schemes,
    });
    const run = () => {
      const ls = [];
      for (let i = 0; i < 100; i++)
        ls.push(legendOf(ch, { title: `Track ${i % 10}`, layer: `l${i}`, painted: ['a', 'b', 'c'] })!);
      return mergeLegends(ls);
    };
    run();
    expect(bestOf(within, 1, run).length).toBe(1);
  });
});
