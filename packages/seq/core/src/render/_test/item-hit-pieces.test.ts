/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An item's hit over a cyclic axis carries the item's own pieces (_spec/README.md "Render", Lanes and Hits): the row
// across the origin and its other piece, whatever is in view; a repeated key elsewhere (exons of one gene) is not one.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { hitView, type RangesLay } from '../index.js';
import { settle } from './helpers.js';
import { layerLay } from './small-helpers.js';

const SIZE = { width: 400, height: 60 } as const;

/** A ring of 100: `a` across the origin as one row past the end, `b` as two rows, `c` twice in the middle. */
function spec(): Node {
  const b = createSeqSpec({ title: 'pieces' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(25), cyclic: true });
  const rows = {
    start: [80, 0, 20, 40],
    end: [100, 5, 30, 50],
    id: ['b', 'b', 'c', 'c'],
  };
  const items = b.data.ranges('items', { target: s, rows: rows as never });
  b.views.view('main', (v) => {
    v.section('s', { of: s }).track('items').layer('items', { type: 'ranges', data: items });
  });
  return b.build();
}

async function frameOf(viewport?: { start: number; end: number }) {
  const f = await (await settle(spec())).frame({ size: SIZE, settled: true, ...(viewport ? { viewport } : {}) });
  return { f, lay: layerLay<RangesLay>(f, 0) };
}

/** The hit at the middle of the mark of the row `row`'s item piece with this key and start. */
function hitOf(lay: RangesLay, layout: Parameters<typeof hitView>[0], key: string, start: number) {
  for (let j = 0; j < lay.hn; j++) {
    if (String(lay.keys[lay.hitem[j]]) !== key || lay.starts[lay.hitem[j]] !== start) continue;
    return hitView(
      layout,
      (lay.hx0[j] + lay.hx1[j]) / 2,
      lay.laneY + lay.hlane[j] * lay.laneHeight + lay.laneHeight / 2,
    );
  }
  return undefined;
}

const rangesOf = (h: ReturnType<typeof hitView>) => [...(h?.extent?.parts[0]?.ranges ?? [])].sort((x, y) => x - y);

describe('an item hit over a cyclic axis', () => {
  it('an item of two rows across the origin: both pieces, the hit on either', async () => {
    const { f, lay } = await frameOf();
    expect(rangesOf(hitOf(lay, f.layout, 'b', 80))).toEqual([0, 5, 80, 100]);
    expect(rangesOf(hitOf(lay, f.layout, 'b', 0))).toEqual([0, 5, 80, 100]);
  });

  it('zoomed in on one piece: the hit still carries the whole item', async () => {
    const { f, lay } = await frameOf({ start: 82, end: 98 });
    const h = hitOf(lay, f.layout, 'b', 80);
    expect(h).toBeDefined();
    expect(rangesOf(h)).toHaveLength(4);
  });

  it('a key repeated away from the origin is not one item’s pieces', async () => {
    const { f, lay } = await frameOf();
    expect(rangesOf(hitOf(lay, f.layout, 'c', 20))).toEqual([20, 30]);
    expect(rangesOf(hitOf(lay, f.layout, 'c', 40))).toEqual([40, 50]);
  });
});
