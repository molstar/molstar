/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { afterEach, describe, expect, it } from 'vitest';
import { dataFilter, itemsAddress } from '../behaviors/theming.js';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

// A feature-type track (seqx-core _spec "Data", Records): its layer's data holds a leading form (the track), and the
// legend's select, hover and hide act on the items of the track, whatever form leads.

const rows = {
  start: new Int32Array([0, 6, 14, 20]),
  end: new Int32Array([5, 12, 20, 25]),
  id: ['A', 'B', 'C', 'D'],
  category: ['binding', 'linker', 'binding', 'linker'],
  label: ['k1', 'k1', 'k2', 'k2'],
};

const Item = '(item :label k1)';
const forms = { item: Item, and: `(intersect :of [${Item}])` };

function spec(lead: string): Node {
  const b = createSeqSpec();
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLL' });
  const d = b.data.ranges('domains', { target: p, rows });
  b.views.view('main', (v) => {
    v.section('p', { of: p })
      .track('domains', { label: 'Domains' })
      .layer('d', { type: 'ranges', data: `${d}/${lead}`, color: { kind: 'field', params: { name: 'category' } } });
  });
  return b.build();
}

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];
afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

async function make(lead: string) {
  frames = installFrames();
  const viewer = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [[viewerItem('viewport.maxPxPerPosition'), 1000]],
  });
  viewers.push(viewer);
  await viewer.load(spec(lead));
  viewer.resize(600);
  await viewer.settle();
  frames.flush();
  return viewer;
}

const Layer = '/views/main/p/domains/d';
const data = (v: SeqViewer) => {
  let n = v.spec.root.value;
  for (const id of Layer.split('/').filter(Boolean)) n = n.children!.find((c) => c.id === id)!;
  return n.params!.data as string;
};

describe('the legend on a track with a leading form', () => {
  it('builds the items address from the annotation and the base: one query that resolves', () => {
    expect(itemsAddress('/data/f/(item :label k1)', 'category', 'binding')).toBe(
      '/data/f/(item :category binding :label k1)',
    );
    expect(itemsAddress('/data/f/(item :label k1)/(except :of (item :id x))', 'category', ['a', 'b c'])).toBe(
      '/data/f/(item :category [a "b c"] :label k1)',
    );
    expect(itemsAddress('/data/f/(intersect :of [(item :label k1)])', 'category', 'binding')).toBe(
      '/data/f/(intersect :of [(intersect :of [(item :label k1)]) (item :category binding)])',
    );
    expect(itemsAddress('/data/f/(range :from 1 :to 3)/(item :a b)', 'category', 'x')).toBeUndefined();
  });

  it.each(Object.entries(forms))('selects and hovers an entry’s items of the track (%s)', async (_, lead) => {
    const v = await make(lead);
    await v.execute('legend.select', { address: Layer, entry: 'binding' });
    // A is binding and of kind k1; C is binding of k2: not the track's
    expect(v.selection.value.value.elements.keys.map((k) => k.key).sort()).toEqual(['A']);
    await v.execute('legend.select', { address: Layer, entry: 'linker', toggle: true });
    expect(v.selection.value.value.elements.keys.map((k) => k.key).sort()).toEqual(['A', 'B']);
    await v.execute('legend.hover', { address: Layer, entry: 'linker' });
    expect(v.highlight.names.value).toContain('legend');
  });

  it.each(Object.entries(forms))('hides a value on it: the base kept (%s)', async (_, lead) => {
    const v = await make(lead);
    await v.execute('layer.filter', { address: Layer, entry: 'linker' });
    expect(data(v)).toBe(`/data/domains/${lead}/(except :of (item :category linker))`);
    expect(dataFilter(data(v), 'category')).toEqual({ base: `/data/domains/${lead}`, hidden: ['linker'] });
    await v.execute('layer.filter', { address: Layer });
    expect(data(v)).toBe(`/data/domains/${lead}`);
  });
});
