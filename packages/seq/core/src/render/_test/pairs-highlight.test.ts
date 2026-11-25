/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A highlight channel's element keys on the pairs layer (_spec/README.md "Render", The pairs layer, Highlight): the
// pairs a channel lists drawn over in its color in the `highlight` pass, the others not; arcs over a track, chords on a
// ring, dots on a grid, and the SVG export alike; the last channel listing a pair wins.

import { EmptyLocation, keySet } from '@molstar/seqx-core/query';
import { describe, expect, it } from 'vitest';
import { type EmphasisHighlight, recordDrawLists } from '../index.js';
import { renderSvg } from '../render-svg.js';
import { arcsSpec, gridSpec } from './pairs-fixtures.js';

const Red = 0xff0000ff;
const Green = 0x00ff00ff;
const Pairs: [number, number][] = [
  [0, 20],
  [2, 18],
  [5, 9],
  [11, 15],
];

const channel = (keys: string[], color: number, owner = '/data/p'): EmphasisHighlight => ({
  location: EmptyLocation,
  color,
  elements: keySet(keys.map((key) => ({ owner, key }))),
});

/** The highlight pass's strokes by color: `#rrggbb` → how many marks. */
function strokes(lists: Parameters<typeof recordDrawLists>[0]): Record<string, number> {
  const out: Record<string, number> = {};
  for (const l of recordDrawLists(lists, ['highlight']).split('\n')) {
    const m = / stroke (#[0-9a-f]{6})/.exec(l);
    if (m) out[m[1]] = (out[m[1]] ?? 0) + 1;
  }
  return out;
}

const Size = { width: 600, height: 400 };

describe('the pairs layer draws a highlight channel’s keys', () => {
  it('arcs: the listed pairs in the channel’s color, the last channel listing a pair winning; the SVG too', async () => {
    const highlights = [channel(['0-20', '5-9'], Red), channel(['5-9'], Green)];
    const r = await renderSvg(arcsSpec({ n: 24, pairs: Pairs }), { size: Size, emphasis: { highlights } });
    expect(strokes(r.lists)).toEqual({ '#ff0000': 1, '#00ff00': 1 });
    expect(r.svg).toContain('#00ff00');
    // keys of another owner, or none: nothing
    const other = await renderSvg(arcsSpec({ n: 24, pairs: Pairs }), {
      size: Size,
      emphasis: { highlights: [channel(['0-20'], Red, '/data/q')] },
    });
    expect(strokes(other.lists)).toEqual({});
  });

  it('chords on a ring and dots on a grid', async () => {
    const ring = await renderSvg(arcsSpec({ n: 24, pairs: Pairs, layout: 'rings', side: 'below' }), {
      size: Size,
      emphasis: { highlights: [channel(['2-18', '11-15'], Red)] },
    });
    expect(strokes(ring.lists)).toEqual({ '#ff0000': 2 });
    const grid = await renderSvg(gridSpec({ n: 24, pairs: Pairs }), {
      size: { width: 600, height: 900 },
      emphasis: { highlights: [channel(['11-15'], Red)] },
    });
    expect(strokes(grid.lists)['#ff0000']).toBe(1);
  });
});
