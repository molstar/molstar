/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A pairs node's numeric extra columns (_spec/README.md "Resolve", Pairs and grids; "Render", the pairs layer): kept as
// numbers beside the text ones, read by the layer's color channel as numbers (a continuous scale over an occupancy).

import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { readPairs } from '@molstar/seqx-core/resolve/held';
import { describe, expect, it } from 'vitest';
import { recordDrawLists } from '../index.js';
import { settle } from './helpers.js';
import { letters, pairRows } from './pairs-fixtures.js';

const Pairs: [number, number][] = [
  [2, 40],
  [5, 30],
  [8, 22],
  [12, 18],
];

describe('a pairs node’s numeric columns', () => {
  it('are kept as numbers, a list or a typed column, NaN missing; text columns stay text', () => {
    const d = readPairs({
      kind: 'pairs',
      params: {
        rows: {
          // the ends as the description normalizes them (typed columns)
          ...Object.fromEntries(
            Object.entries(pairRows(Pairs)).map(([k, v]) => [k, k === 'id' ? v : Int32Array.from(v as number[])]),
          ),
          occupancy: [0.1, 0.2, Number.NaN, 0.4],
          weight: Float32Array.of(1, 2, 3, 4),
          share: ['a', 'b', 'c', 'd'],
        },
      },
    });
    expect(d.extra).toEqual({
      occupancy: [0.1, 0.2, null, 0.4],
      weight: [1, 2, 3, 4],
      share: ['a', 'b', 'c', 'd'],
    });
  });

  it('color the arcs on a continuous scale: each occupancy its own color', async () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'rna', letters: letters(48) });
    const p = b.data.pairs('p', { a: s, rows: { ...pairRows(Pairs), occupancy: [0.05, 0.3, 0.6, 0.95] } });
    b.views.view('v', (v) =>
      v
        .section('s', { of: s })
        .track('pairs')
        .layer('arcs', {
          type: 'pairs',
          data: p,
          color: theme.field('occupancy', { domain: [0, 1], palette: theme.continuous({ colors: 'viridis' }) }),
        }),
    );
    const f = await (await settle(b.build())).frame({ size: { width: 800, height: 200 } });
    expect(f.layout.issues).toEqual([]);
    const strokes = new Set(
      recordDrawLists(f.lists, ['data'])
        .split('\n')
        .filter((l) => l.startsWith('path'))
        .map((l) => /#[0-9a-f]{8}/.exec(l)?.[0]),
    );
    expect(strokes.size).toBe(Pairs.length);
  });
});
