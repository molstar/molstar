/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The order of a section's tracks is the order of the bands (seqspec:viewer#the-text-view): a ranges track declared
// before the letters track draws above the letters' line on every line of the text view (under `wrappedTracks: lanes`
// and under `full`, with a height or without) and above the letters in the tracks view. The Seq*/Text view story's
// Annotations panel relies on it for its "Above" set.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import { describe, expect, it } from 'vitest';
import { DefaultPresets, KernelRenderers, type LayerRenderer, type RenderProps } from '../src/index.js';
import { opened } from './layout-fixtures.js';

const BUILT_IN = {
  presets: false as const,
  renderers: [...KernelRenderers] as LayerRenderer[],
  registry: [{ presets: [...DefaultPresets] }],
};

/** A protein with two ranges tracks, `above` declared before the letters and `below` after; `height` on both. */
function doc(preset: string, height?: number): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  const value = 'MSEEKLKELAEKAIRLAKEGNPEEVLRALQ'.repeat(4);
  const a = b.data.sequence({ ref: 'a', alphabet: 'protein', value });
  const up = b.data.annotation({
    ref: 'up',
    geometry: 'ranges',
    target: a,
    value: { start: [3], stop: [30], label: ['Above'] },
  } as never);
  const down = b.data.annotation({
    ref: 'down',
    geometry: 'ranges',
    target: a,
    value: { start: [40], stop: [70], label: ['Below'] },
  } as never);
  const sec = b.view({ ref: 'v', label: 'V', layout: preset === 'tracks' ? 'tracks' : 'text' }).section({
    ref: 'v.main',
    style: { preset, wrapped_tracks: 'lanes' },
  } as never);
  const h = height === undefined ? {} : { height };
  sec.track({ ref: 'above', label: 'Above', ...h }).ranges({ data: up, labels: true } as never);
  sec.track({ ref: 'letters', label: 'Sequence' }).letters();
  sec.track({ ref: 'below', label: 'Below', ...h }).ranges({ data: down, labels: true } as never);
  return b.getDocument();
}

async function tops(preset: string, wrappedTracks: RenderProps['wrappedTracks'], height?: number) {
  const t = await opened(doc(preset, height), {
    ...BUILT_IN,
    size: { 'v.main': [600, 2000] },
    render: { wrappedGroupLabels: 'none', wrappedGutter: 0, wrappedTracks },
  });
  const e = t.entries('v.main');
  const lines = t.layout('v.main').lines?.count ?? 1;
  return { e, lines, of: (ref: string) => e.filter((x) => x[1] === ref) };
}

describe('a ranges track ordered before the letters draws above them', () => {
  for (const wrappedTracks of ['lanes', 'full'] as const) {
    for (const height of [undefined, 1]) {
      it(`text view, wrapped_tracks ${wrappedTracks}${height ? `, track height ${height}` : ', no height'}: on every line`, async () => {
        const { of, lines } = await tops('wrapped', wrappedTracks, height);
        const [above, letters, below] = [of('above'), of('letters'), of('below')];
        expect(lines).toBeGreaterThan(1);
        expect(letters).toHaveLength(lines);
        // Every line has a band of the above track directly over its letters, none under them; the below track follows.
        for (let n = 0; n < lines; n++) {
          const [from, to] = [n > 0 ? letters[n - 1][3] + letters[n - 1][4] : 0, letters[n][3]];
          const bands = above.filter((x) => x[3] >= from && x[3] < to);
          expect(bands.length, `line ${n}: a band above the letters`).toBe(1);
          expect(bands[0][3] + bands[0][4]).toBeLessThanOrEqual(to);
        }
        expect(above).toHaveLength(lines);
        expect(below.length).toBeGreaterThan(0);
        for (const x of below) expect(letters.some((l) => l[3] + l[4] <= x[3])).toBe(true);
      });
    }
  }

  it('tracks view: the ranges track declared first is the row above the letters', async () => {
    const { of } = await tops('tracks', 'lanes');
    const [above, letters, below] = [of('above')[0], of('letters')[0], of('below')[0]];
    expect(above[3] + above[4]).toBeLessThanOrEqual(letters[3]);
    expect(letters[3] + letters[4]).toBeLessThanOrEqual(below[3]);
  });
});
