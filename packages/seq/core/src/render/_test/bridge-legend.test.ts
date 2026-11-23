/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A bridge's legend (_spec/README.md "Render", Bridges; "Themes", Legends): its ribbons coloured by the built-in
// `orientation` (same, inverted: both keyed whichever is drawn), a data column of that name winning over it, and the
// width note only where a drawn ribbon's width is its aligned length (never a line).

import type { Node } from '@molstar/core/repr/tree';
import { theme } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { viewLegends } from '../legends.js';
import { bridgeSpec, type BridgeSpecOptions } from './bridge-fixtures.js';
import { settle } from './helpers.js';

const Size = { width: 900, height: 600 };
const Pairs: BridgeSpecOptions['pairs'] = [
  { a: [10, 50], b: [20, 60] },
  { a: [100, 160], b: [120, 180], strand: '-' },
];

async function legendsOf(o: Partial<BridgeSpecOptions>, edit?: (spec: Node) => void) {
  const spec = bridgeSpec({ pairs: Pairs, layer: { color: theme.field('orientation', { list: false }) }, ...o });
  edit?.(spec);
  const s = await settle(spec);
  const f = await s.renderer.frame({ size: Size, settled: true });
  const lay = f.layout.bridges!.list[0].lays[0].lay as { n: number; color: Uint32Array };
  return Object.assign(
    viewLegends(f.layout).legends.filter((l) => l.layers.includes('/views/v/br/ribbons')),
    { colours: new Set(Array.from(lay.color.subarray(0, lay.n))) },
  );
}

const labels = (l: { entries: { kind: string; items?: readonly { label: string }[] } }) =>
  l.entries.kind === 'swatches' ? l.entries.items!.map((i) => i.label) : [];

describe('a bridge’s legend', () => {
  it('keys its ribbons’ orientation, same and inverted, titled by the bridge, noting a ribbon’s width', async () => {
    const legends = await legendsOf({});
    expect(labels(legends[0])).toEqual(['same', 'inverted']);
    expect(legends[0].notes).toContain('width');
    expect(legends.colours.size).toBe(2);
  });

  it('a data column named `orientation` wins over the computed field', async () => {
    const legends = await legendsOf({}, (spec) => {
      const p = spec.children!.find((c) => c.kind === 'data')!.children!.find((c) => c.id === 'p')!;
      (p.params!.rows as Record<string, unknown>).orientation = ['forward', 'forward'];
    });
    // one value in the data: one colour for both, the inverted one too; never the computed `same` and `inverted`
    expect(legends.colours.size).toBe(1);
    expect(labels(legends[0])).not.toContain('inverted');
  });

  it('says nothing of the width where a ribbon is drawn as a line (`shape: line`)', async () => {
    const [l] = await legendsOf({ layer: { color: theme.field('orientation', { list: false }), shape: 'line' } });
    expect(labels(l)).toEqual(['same', 'inverted']);
    expect(l.notes).not.toContain('width');
  });
});
