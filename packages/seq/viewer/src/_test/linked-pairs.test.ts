/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `linkedPairs` behavior (_spec/README.md "Behaviors", Linked pairs): a hovered or selected range element lights
// the pairs with an end in its extent through its highlight channel, drawn by the pairs layer in the highlight pass;
// cleared when the hover or selection goes; a pair or a position hovered is the layer's own rule; its params (`on`,
// `layers`); a viewer composed without it lights none; the hover's query made once an element.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { recordDrawLists } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it } from 'vitest';
import { LinkedPairsChannel, type LinkedPairsParams, LinkedPairsProvider } from '../behaviors/linked-pairs.js';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

/** An RNA of 24 with four pairs, two features over it, the pairs as arcs. */
function spec(): Node {
  const b = createSeqSpec({ title: 'linked pairs' });
  const s = b.data.sequence('s', { alphabet: 'rna', letters: 'ACGUACGUACGUACGUACGUACGU' });
  b.data.ranges('f', { target: s, rows: { start: [0, 9], end: [4, 12], id: ['stem', 'loop'] } });
  const p = b.data.pairs('p', {
    a: s,
    rows: { a_start: [0, 2, 5, 11], a_end: [1, 3, 6, 12], b_start: [20, 18, 9, 15], b_end: [21, 19, 10, 16] },
  });
  b.views.view('v', (v) => {
    const sec = v.section('s', { of: s });
    sec.track('features').layer('f', { type: 'ranges', data: '/data/f' });
    sec.track('pairs').layer('arcs', { type: 'pairs', data: p });
    sec.track('letters').layer('l', { type: 'letters' });
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

async function make(linked: LinkedPairsParams | false = {}): Promise<SeqViewer> {
  frames = installFrames();
  const base = DefaultSeqViewerSpec();
  const behaviors = base.behaviors.filter((b) => b.provider !== LinkedPairsProvider);
  viewer = createSeqViewer({
    ...base,
    behaviors: linked === false ? behaviors : [...behaviors, { provider: LinkedPairsProvider, params: linked }],
    ui: undefined,
    config: [[viewerItem('navigator.form'), 'none']],
  });
  await viewer.load(spec());
  viewer.resize(600);
  await viewer.settle();
  frames.flush();
  return viewer;
}

const Owner = {};
/** Hovers a feature (as its layer's hit says it). */
const hoverFeature = (v: SeqViewer, key: string) =>
  v.hover.set(
    { layer: '/views/v/s/features/f', system: '/data/s', element: { owner: '/data/f', key }, position: 1 },
    Owner,
  );
const lit = (v: SeqViewer) => v.highlight.get(LinkedPairsChannel)?.elements.keys.map((k) => k.key);
/** The highlight pass's curves: the arcs lit. */
const curves = (v: SeqViewer) => {
  v.draw();
  return recordDrawLists(v.frame!.lists, ['highlight'])
    .split('\n')
    .filter((l) => l.startsWith('path') && l.includes(' C ')).length;
};

describe('linkedPairs', () => {
  it('a hovered range lights the pairs with an end in it, drawn in the highlight pass; cleared when the hover goes', async () => {
    const v = await make();
    expect(curves(v)).toBe(0);
    hoverFeature(v, 'stem');
    // [0, 4): pairs 0 (0–20) and 1 (2–18)
    expect(lit(v)).toEqual([0, 1]);
    expect(curves(v)).toBe(2);
    // [9, 12): pair 2 by its b end (9), pair 3 by its a end (11)
    hoverFeature(v, 'loop');
    expect(lit(v)).toEqual([2, 3]);
    v.hover.clear(Owner);
    expect(lit(v)).toBeUndefined();
    expect(curves(v)).toBe(0);
    // a position hovered is the layer's own rule (its partners), never the channel
    v.hover.set({ layer: '/views/v/s/letters/l', system: '/data/s', position: 0 }, Owner);
    expect(lit(v)).toBeUndefined();
  });

  it('a selected range lights them on the channel, the hover’s added; cleared with the selection', async () => {
    const v = await make();
    v.addresses.select('/data/f/(item :id loop)');
    expect(lit(v)).toEqual([2, 3]);
    // the selection's pass holds the feature alone: the pairs are related to it, not in it
    expect(v.selection.value.value.elements.keys.map((k) => k.key)).toEqual(['loop']);
    hoverFeature(v, 'stem');
    expect(lit(v)).toEqual([0, 1, 2, 3]);
    v.hover.clear(Owner);
    v.selection.clear();
    expect(lit(v)).toBeUndefined();
  });

  it('its params: `on` the selection alone; `layers` by address', async () => {
    const sel = await make({ on: ['selection'] });
    hoverFeature(sel, 'stem');
    expect(lit(sel)).toBeUndefined();
    sel.addresses.select('/data/f/(item :id stem)');
    expect(lit(sel)).toEqual([0, 1]);
    viewer?.retire();
    frames?.restore();
    const none = await make({ layers: ['/views/v/s/other/arcs'] });
    hoverFeature(none, 'stem');
    expect(lit(none)).toBeUndefined();
    viewer?.retire();
    frames?.restore();
    const named = await make({ layers: ['/views/v/s/pairs/arcs'] });
    hoverFeature(named, 'stem');
    expect(lit(named)).toEqual([0, 1]);
  });

  it('a viewer composed without it lights none (the opt-out)', async () => {
    const v = await make(false);
    hoverFeature(v, 'stem');
    expect(v.highlight.names.value).toEqual([]);
    expect(curves(v)).toBe(0);
  });

  it('the hover’s query is made once an element: the same keys, the same channel value', async () => {
    const v = await make();
    hoverFeature(v, 'stem');
    const first = v.highlight.get(LinkedPairsChannel);
    const version = v.highlight.version.value;
    v.hover.set(
      { layer: '/views/v/s/features/f', system: '/data/s', element: { owner: '/data/f', key: 'stem' }, position: 3 },
      Owner,
    );
    expect(v.highlight.version.value).toBe(version);
    hoverFeature(v, 'loop');
    hoverFeature(v, 'stem');
    expect(v.highlight.get(LinkedPairsChannel)).toBe(first);
  });
});
