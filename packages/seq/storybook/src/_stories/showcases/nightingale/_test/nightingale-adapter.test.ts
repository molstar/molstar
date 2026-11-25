/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqContext } from '@molstar/seqx-core/context/default';
import { DataStore, SeqRuntime } from '@molstar/seqx-core/data';
import { rangesLocation } from '@molstar/seqx-core/query';
import { openModel } from '@molstar/seqx-core/render';
import { P53Views } from '@molstar/seqx-examples';
import { describe, expect, it } from 'vitest';
import {
  type RangeRow,
  readSection,
  themeColors,
  toFeatures,
  toHighlight,
  toLineData,
} from '../nightingale-adapter.js';

// The adapter's mapping on a fixture (no DOM, no Nightingale): what is read from the resolved view and the store, and the
// shapes Nightingale's elements are fed.

const Ranges: readonly RangeRow[] = [
  { id: 'A', label: 'First', start: 0, end: 10, fill: '#112233' },
  { id: 'B', label: 'Second', start: 20, end: 21, fill: '#445566' },
];

describe('Nightingale adapter, the mapping', () => {
  it('ranges become features: 0-based half-open to 1-based inclusive, in their colors', () => {
    const f = toFeatures(Ranges);
    expect(f.map((x) => x.locations[0].fragments[0])).toEqual([
      { start: 1, end: 10 },
      { start: 21, end: 21 },
    ]);
    expect(f.map((x) => x.color)).toEqual(['#112233', '#445566']);
    expect(f[0].tooltipContent).toBe('First, 1–10');
  });

  it('values become one line: positions 1-based, a missing one left out, the range from zero', () => {
    const [line] = toLineData(
      [
        { position: 0, value: 0.5 },
        { position: 1, value: Number.NaN },
        { position: 2, value: 2 },
      ],
      'conservation',
      '#000000',
      'none',
    );
    expect(line.values).toEqual([
      { position: 1, value: 0.5 },
      { position: 3, value: 2 },
    ]);
    expect(line.range).toEqual([0, 2]);
  });

  it('a location is the highlight attribute: its ranges on the system, 1-based inclusive, others ignored', () => {
    const loc = rangesLocation('/data/p53', [4, 5, 9, 12]);
    expect(toHighlight(loc, '/data/p53')).toBe('5:5,10:12');
    expect(toHighlight(loc, '/data/other')).toBe('');
    expect(toHighlight(undefined, '/data/p53')).toBe('');
  });

  it('the theme tokens differ by theme', () => {
    expect(themeColors('light').text).not.toBe(themeColors('dark').text);
    expect(themeColors('dark').text).toMatch(/^#[0-9a-f]{6}$/);
  });
});

describe('Nightingale adapter, reading the resolved view', () => {
  it('reads p53 tracks: its letters, its six domains by their category colors, its values', async () => {
    const ctx = createSeqContext({ registry: [{ kind: 'format-runtime', runtime: SeqRuntime }] });
    const store = new DataStore();
    const opened = await openModel(P53Views.find((v) => v.id === 'tracks')!.build(), store, {
      runtimes: ctx.runtimes,
      config: ctx.resolveConfig,
      assets: ctx.assets,
    });
    const read = await readSection(opened.model, opened.view!, store);
    ctx.retire();
    expect(read.pending).toBe(0);
    expect(read.length).toBe(393);
    expect(read.letters).toHaveLength(393);
    expect(read.letters.startsWith('MEEPQSDPSV')).toBe(true);
    expect(read.ranges.map((r) => r.id)).toEqual(['TAD', 'PRR', 'DBD', 'NLS', 'TET', 'CTD']);
    expect(read.ranges.find((r) => r.id === 'DBD')).toMatchObject({ label: 'DNA-binding', start: 93, end: 292 });
    expect(new Set(read.ranges.map((r) => r.fill)).size).toBeGreaterThan(1);
    expect(read.values).toHaveLength(393);
    expect(read.values.filter((v) => Number.isNaN(v.value))).toHaveLength(2);
    expect(read.system).toMatch(/p53$/);
  });
});
