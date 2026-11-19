/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The small layers' themes (_spec/README.md "Themes"; "Render", Layer types): secondary structure in its tokens by
// class whatever the data's spelling, keyed by the `secondary` scheme's legend, recolored by a form that reads `kind` as
// the class word; variants in `seq-variant`, or by their consequence (the builder's default for a consequence column)
// with a legend of the values drawn; the ruler in `seq-ruler`.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { LightTokens, recordDrawLists, type SecondaryLay, type VariantsLay, viewLegends } from '../index.js';
import { P53 } from './fixtures.js';
import { settle } from './helpers.js';

const Size = { width: 800, height: 120 };

function spec(
  layers: (b: ReturnType<typeof createSeqSpec>, p: string) => Record<string, Record<string, unknown>>,
): Node {
  const b = createSeqSpec();
  const p = b.data.sequence('p', { alphabet: 'protein', letters: P53 });
  const ls = layers(b, p);
  b.views.view('main', (v) => {
    const s = v.section('p', { of: p });
    for (const [id, params] of Object.entries(ls)) s.track(id, { label: id }).layer(id, params as never);
  });
  return b.build();
}

const ss = (b: ReturnType<typeof createSeqSpec>, p: string) =>
  b.data.ranges('ss', {
    target: p,
    rows: { start: [100, 140, 180, 220], end: [120, 160, 200, 240], kind: ['H', 'E', 'T', 'helix'] },
  });

async function layOf<T>(
  node: Node,
  id: string,
): Promise<{ lay: T; layout: Awaited<ReturnType<Awaited<ReturnType<typeof settle>>['frame']>>['layout'] }> {
  const s = await settle(node);
  const f = await s.frame({ size: Size, settled: true });
  const l = f.layout.sections[0].tracks.find((t) => t.track.label === id)!.layers[0];
  return { lay: l.lay as T, layout: f.layout };
}

describe('secondary structure', () => {
  it('draws its classes in the seq-ss tokens, keyed by the secondary scheme’s legend', async () => {
    const { lay, layout } = await layOf<SecondaryLay>(
      spec((b, p) => ({ ss: { type: 'secondary', data: ss(b, p) } })),
      'ss',
    );
    const c = LightTokens.colors;
    expect([...lay.fills.slice(0, 4)]).toEqual([
      c['seq-ss-helix'],
      c['seq-ss-strand'],
      c['seq-ss-turn'],
      c['seq-ss-helix'],
    ]);
    const legend = viewLegends(layout).legends.find((l) => l.layers.includes('/views/main/p/ss/ss'));
    expect(legend?.entries.kind === 'swatches' ? legend.entries.items.map((e) => e.label) : []).toEqual([
      'helix',
      'strand',
      'turn',
    ]);
  });

  it('a form reads `kind` as the class word; a uniform color paints every class', async () => {
    const red = '#ff0000';
    const themed = await layOf<SecondaryLay>(
      spec((b, p) => ({
        ss: {
          type: 'secondary',
          data: ss(b, p),
          color: theme.field('kind', { palette: theme.categorical({ colors: { strand: red } }) }),
        },
      })),
      'ss',
    );
    const strand = themed.lay.fills[1];
    expect(strand >>> 8).toBe(0xff0000);
    const uniform = await layOf<SecondaryLay>(
      spec((b, p) => ({ ss: { type: 'secondary', data: ss(b, p), color: red } })),
      'ss',
    );
    expect([...uniform.lay.fills.slice(0, 4)].every((f) => f >>> 8 === 0xff0000)).toBe(true);
  });
});

describe('variants', () => {
  const rows = {
    start: [10, 50, 90, 130],
    end: [11, 51, 91, 131],
    consequence: ['missense', 'nonsense', 'missense', 'synonymous'],
  };

  it('by consequence when the builder sees one (a legend of the listed values); seq-variant when there is none', async () => {
    const { lay, layout } = await layOf<VariantsLay>(
      spec((b, p) => ({ v: { type: 'variants', data: b.data.ranges('v', { target: p, rows }) } })),
      'v',
    );
    expect(lay.fills[0]).toBe(lay.fills[2]);
    expect(lay.fills[0]).not.toBe(lay.fills[1]);
    const legend = viewLegends(layout).legends.find((l) => l.field === 'consequence');
    expect(legend?.entries.kind === 'swatches' ? legend.entries.items.map((e) => e.label).sort() : []).toEqual([
      'missense',
      'nonsense',
      'synonymous',
    ]);
    const plain = await layOf<VariantsLay>(
      spec((b, p) => ({
        v: { type: 'variants', data: b.data.ranges('v', { target: p, rows: { start: rows.start, end: rows.end } }) },
      })),
      'v',
    );
    expect([...plain.lay.fills.slice(0, 4)].every((f) => f === LightTokens.colors['seq-variant'])).toBe(true);
  });
});

describe('the ruler', () => {
  it('draws its ticks and numbers in seq-ruler', async () => {
    const s = await settle(spec(() => ({ r: { type: 'ruler' } })));
    const f = await s.frame({ size: Size, settled: true });
    const hex = (LightTokens.colors['seq-ruler'] >>> 0).toString(16).padStart(8, '0');
    expect(recordDrawLists(f.lists, ['data'])).toContain(hex);
  });
});
