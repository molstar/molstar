/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The view across layouts (_spec/README.md "ViewModel", Layouts), headless: a step of the history over a layout patch
// restored by position (back from a text lands on the window's positions, and back to a text on its anchor), groups
// only within one layout, a text's moves by its own arithmetic before and after a frame, and the letters copied as a
// reversed section shows them (complemented, in its reading order) or as stored.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import type { TextPlacement } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Letters = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDD'.repeat(30);
const Main = '/views/main/main';
const Twin = '/views/main/twin';

function spec(twin = true): Node {
  const b = createSeqSpec({ title: 'layouts' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: Letters });
  const d = b.data.sequence('d', { alphabet: 'dna', letters: 'AACCGGTTAC' });
  b.views.view('main', (v) => {
    v.section('main', { of: p }).track('letters').layer('l', { type: 'letters' });
    // on the same system, the same layout at first: one scroll group
    if (twin) v.section('twin', { of: p }).track('letters').layer('l', { type: 'letters' });
    v.section('dna', { of: d, direction: 'reversed' }).track('letters').layer('l', { type: 'letters' });
  });
  return b.build();
}

let frames: Frames | undefined;
let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  frames?.restore();
  frames = undefined;
  vi.useRealTimers();
});

async function make(twin = true) {
  vi.useFakeTimers({ toFake: ['Date', 'performance'] });
  frames = installFrames();
  const v = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [[viewerItem('view.jumpMs'), 0]],
  });
  viewer = v;
  await v.load(spec(twin));
  v.resize(600);
  await v.settle();
  frames.flush();
  return v;
}

const asText = (v: SeqViewer, layout: string) =>
  v.spec.applyPatch({ ops: [{ op: 'set', at: ['views', 'main', 'main'], key: 'layout', value: layout }] });

const textOf = (v: SeqViewer) => v.frame!.layout.sections[0].placement as TextPlacement;

describe('the view across layouts', () => {
  it('groups only within one layout: a text and a window on one system move apart', async () => {
    const v = await make();
    expect(v.view.groupOf(Main)).toEqual([Main, Twin]);
    asText(v, 'text');
    frames!.flush();
    expect(v.view.groupOf(Main)).toEqual([Main]);
    v.view.setRange(Twin, 100, 120);
    frames!.flush();
    expect(v.view.get(Main)!.start).toBe(0);
  });

  it('a step of the history over a layout patch is restored by position, never dropped', async () => {
    // alone on its system: a step is restored by its section's own group, not a twin's
    const v = await make(false);
    v.view.setRange(Main, 100, 150);
    frames!.flush();
    // a jump from the window [100, 150): the step back holds it
    expect(v.view.jump(() => v.view.setRange(Main, 300, 340))).toBe('moved');
    asText(v, 'text');
    frames!.flush();
    expect(v.view.back()).toBe(true);
    frames!.flush();
    // the text shows the window's positions: its anchor at the window's first (its line at the top, a fraction in)
    const p = textOf(v);
    expect(v.view.get(Main)!.start).toBeCloseTo(100, 6);
    expect(p.window.start).toBeLessThanOrEqual(100);
    expect(p.window.end).toBeGreaterThan(100);
    // and forward again, then back to a window: the text's lines in view as its range
    expect(v.view.forward()).toBe(true);
    frames!.flush();
    v.view.panPx(Main, 0, 40);
    frames!.flush();
    const anchor = v.view.get(Main)!.start;
    expect(v.view.jump(() => v.view.fit(Main))).toBe('moved');
    asText(v, 'tracks');
    frames!.flush();
    expect(v.view.back()).toBe(true);
    frames!.flush();
    const w = v.view.get(Main)!;
    expect(w.start).toBeGreaterThan(anchor - 1);
    expect(w.start).toBeLessThan(anchor + 1);
    expect(w.end).toBeGreaterThan(w.start);
  });

  it('a text’s letter size is kept across a patch of its params, and dropped by a fit', async () => {
    const v = await make();
    asText(v, 'text');
    frames!.flush();
    expect(v.view.zoom(Main, 1.5, 0)).toBe(true);
    frames!.flush();
    const scale = v.view.get(Main)!.scale!;
    expect(scale).toBeGreaterThan(12);
    expect(textOf(v).geometry!.cell).toBe(scale);
    v.spec.applyPatch({ ops: [{ op: 'set', at: ['views', 'main', 'main'], key: 'group', value: 5 }] });
    frames!.flush();
    expect(v.view.get(Main)!.scale).toBe(scale);
    v.view.fit(Main);
    frames!.flush();
    expect(v.view.get(Main)!.scale).toBeUndefined();
  });

  it('copies a reversed DNA section’s letters complemented, in its reading order; stored, as stored', async () => {
    const v = await make();
    await v.execute('selection.add', { range: { system: '/data/d', start: 2, end: 6 } });
    // stored CCGG at 2–5; read right to left on the reversed section, each complemented: CCGG
    expect(await v.execute('selection.copy-letters', {})).toBe('CCGG');
    await v.execute('selection.clear', {});
    await v.execute('selection.add', { range: { system: '/data/d', start: 0, end: 3 } });
    // AAC read backwards (CAA), complemented: GTT
    expect(await v.execute('selection.copy-letters', {})).toBe('GTT');
    expect(await v.execute('selection.copy-letters', { stored: true })).toBe('AAC');
  });
});
