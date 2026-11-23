/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Legends state their domain (_spec/README.md "Legends"): a continuous legend over an open domain says so ("in view")
// and its ends are the values the view settled on, derived again when a pan settles; a closed one says its own ends; a
// palette's domain is its own, never the y domain beside it (EC-205); a categorical legend says which values it lists:
// every listed one (the vocabulary), or the painted ones with "values in view" (EC-204), a value no palette lists one
// hashed colour on every track (EC-203).

import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { afterEach, describe, expect, it } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { legendRows } from '../legends.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Section = '/views/main/p';
const Letters = 'MEEPQSDPSVEPPLSQETFSDLWKLL';

function spec() {
  const b = createSeqSpec();
  const p = b.data.sequence('p', { alphabet: 'protein', letters: Letters });
  const values = b.data.values('v', { target: p, values: Float32Array.from({ length: 26 }, (_, i) => i) });
  const d = b.data.ranges('d', {
    target: p,
    rows: {
      start: new Int32Array([0, 6, 14]),
      end: new Int32Array([5, 12, 20]),
      id: ['A', 'B', 'C'],
      category: ['binding', 'linker', 'binding'],
    },
  });
  // the same values first met in the other order: a colour by first occurrence would differ, a hash cannot
  const d2 = b.data.ranges('d2', {
    target: p,
    rows: {
      start: new Int32Array([0, 6, 14]),
      end: new Int32Array([5, 12, 20]),
      id: ['A', 'B', 'C'],
      category: ['linker', 'binding', 'linker'],
    },
  });
  b.views.view('main', (v) => {
    const s = v.section('p', { of: p });
    s.track('open', { label: 'Open' }).layer('v', { type: 'values', data: values, color: theme.field('value') });
    // a palette's domain apart from the bars' y domain (EC-205)
    s.track('closed', { label: 'Closed' }).layer('v', {
      type: 'values',
      data: values,
      domain: [0, 100],
      color: theme.field('value', { domain: [0, 10], palette: theme.continuous({ colors: 'viridis' }) }),
    } as never);
    // the vocabulary: a palette listing its values (one not drawn here)
    s.track('listed', { label: 'Listed' }).layer('d', {
      type: 'ranges',
      data: d,
      color: theme.field('category', {
        palette: theme.categorical({ values: ['binding', 'linker', 'motif'], colors: 'categories' }),
      }),
    });
    // painted: no list, hashed colours, on two tracks of two data
    for (const [t, data] of [
      ['Painted', d],
      ['Again', d2],
    ] as const)
      s.track(t.toLowerCase(), { label: t }).layer('d', {
        type: 'ranges',
        data,
        color: theme.field('category', { list: false }),
      });
  });
  return b.build();
}

let frames: Frames | undefined;
let viewer: SeqViewer | undefined;
afterEach(() => {
  viewer?.retire();
  frames?.restore();
  viewer = undefined;
  frames = undefined;
});

async function make() {
  frames = installFrames();
  const v = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [[viewerItem('viewport.maxPxPerPosition'), 1000]],
  });
  viewer = v;
  await v.load(spec());
  v.resize(600);
  await v.settle();
  frames.flush();
  return v;
}

const row = (v: SeqViewer, title: string) => legendRows(v).find((r) => r.title === title)!;

describe('legends state their domain', () => {
  it('an open domain says "in view", its ends the values the view settled on; a pan settled derives them again', async () => {
    const v = await make();
    const open = row(v, 'Open');
    expect(open.note).toBe('in view');
    expect(open.entries[0].label).toBe('0 to 25');
    expect(open.entries[0].name).toBe('Open: 0 to 25, in view');
    v.view.setRange(Section, 10, 20);
    await v.settle();
    frames!.flush();
    const [lo, hi] = v.legends.legends.value.find((l) => l.title === 'Open')!.domain!;
    // the values in view, inside the whole data's
    expect(lo).toBeGreaterThanOrEqual(9);
    expect(hi).toBeLessThanOrEqual(21);
    expect(hi - lo).toBeLessThan(25);
    expect(row(v, 'Open').entries[0].label).toBe(`${lo} to ${hi}`);
  });

  it('a closed domain states its own ends, never "in view"; a palette’s domain apart from the y domain (EC-205)', async () => {
    const v = await make();
    const closed = row(v, 'Closed');
    expect(closed.note).toBeUndefined();
    expect(closed.entries[0].label).toBe('0 to 10');
    expect(closed.entries[0].name).toBe('Closed: 0 to 10');
  });

  it('a categorical legend says which values it lists: the vocabulary whole, or the painted ones as "values in view" (EC-204)', async () => {
    const v = await make();
    const listed = v.legends.legends.value.find((l) => l.title === 'Listed')!;
    expect(listed.source).toBe('vocabulary');
    expect(row(v, 'Listed').entries.map((e) => e.label)).toEqual(['binding', 'linker', 'motif']);
    expect(row(v, 'Listed').note).toBeUndefined();
    const painted = v.legends.legends.value.find((l) => l.title === 'Painted')!;
    expect(painted.source).toBe('painted');
    expect(row(v, 'Painted').note).toBe('values in view');
  });

  it('a value no palette lists takes one hashed colour, the same on every track (EC-203)', async () => {
    const v = await make();
    const colours = (title: string) =>
      Object.fromEntries(
        legendRows(v)
          .filter((r) => r.layers.some((l) => l.address.includes(`/${title}/`)))
          .flatMap((r) => r.entries.map((e) => [e.label, e.color])),
      );
    const painted = colours('painted');
    expect(Object.keys(painted).sort()).toEqual(['binding', 'linker']);
    // each value the same colour on both tracks, whichever it met first
    expect(colours('again')).toEqual(painted);
  });
});
