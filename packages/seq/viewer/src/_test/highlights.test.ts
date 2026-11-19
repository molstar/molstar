// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Highlight channels through the viewer: a channel's style is set with it (`address.highlight`'s `style`, or the
// model's value) and kept onto the sections it is carried to; overlapping channels draw as stripes; the Legend part
// lists the active channels with their colors and clears each by name.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { recordDrawLists } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { mountLegend } from '../mount-legend.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

/** 10 residues on a 40-nucleotide gene (one exon, 5 to 35). */
function spec(): Node {
  const b = createSeqSpec({ title: 'highlights' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSV' });
  const g = b.data.sequence('g', { alphabet: 'dna', letters: 'ACGT'.repeat(10) });
  b.data.map('m', { from: p, to: g, method: 'codon', exons: { start: [5], end: [35] }, strand: '+' });
  b.views.view('main', (v) => {
    const protein = v.section('protein', { of: p });
    protein.track('letters').layer('l', { type: 'letters' });
    v.section('gene', { of: g, follow: protein }).track('letters').layer('l', { type: 'letters' });
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

async function make(): Promise<SeqViewer> {
  frames = installFrames();
  viewer = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [
      [viewerItem('viewport.maxPxPerPosition'), 1000],
      [viewerItem('navigator.form'), 'none'],
    ],
  });
  await viewer.load(spec());
  viewer.resize(600);
  await viewer.settle();
  frames.flush();
  return viewer;
}

const Range = (a: number, b: number) => `/data/p/(range :from ${a} :to ${b})`;

/** The highlight pass's rects in a section: `[y, h, paint]`. */
function rectsIn(v: SeqViewer, section: number) {
  const s = v.frame!.layout.sections[section];
  return recordDrawLists(v.frame!.lists, ['highlight'])
    .split('\n')
    .filter((l) => l.startsWith('rect '))
    .map((l) => l.split(' '))
    .filter((p) => Number(p[2]) >= s.y - 4 && Number(p[2]) < s.y + s.height)
    .map((p) => ({ y: Number(p[2]), h: Number(p[4]), paint: p.slice(5).join(' ') }));
}

describe('a channel’s style', () => {
  it('is set with the channel, and carried onto the linked section with it', async () => {
    const v = await make();
    await v.execute('address.highlight', { address: Range(2, 5), channel: 'a', style: 'band' });
    frames!.flush();
    expect(v.highlight.get('a')?.style).toBe('band');
    for (const section of [0, 1]) {
      const rects = rectsIn(v, section);
      expect(rects.length).toBeGreaterThan(0);
      // a bar along the foot, in the channel's color at full opacity: no tint
      expect(rects.every((r) => r.h === 3 && !r.paint.includes('4d'))).toBe(true);
    }
  });

  it('defaults to the tint, and a changed style is a changed channel', async () => {
    const v = await make();
    v.addresses.highlightAddress(Range(2, 5), 'a');
    frames!.flush();
    expect(rectsIn(v, 0).every((r) => r.paint.endsWith('4d'))).toBe(true);
    const version = v.highlight.version.value;
    v.addresses.highlightAddress(Range(2, 5), 'a', 'outline');
    frames!.flush();
    expect(v.highlight.version.value).toBeGreaterThan(version);
    expect(rectsIn(v, 0).every((r) => r.paint.includes('stroke'))).toBe(true);
    expect(rectsIn(v, 1).every((r) => r.paint.includes('stroke'))).toBe(true);
  });

  it('overlapping channels stripe the band on both sections', async () => {
    const v = await make();
    v.addresses.highlightAddress(Range(2, 6), 'a');
    v.addresses.highlightAddress(Range(4, 8), 'b');
    frames!.flush();
    for (const section of [0, 1]) {
      const heights = new Set(rectsIn(v, section).map((r) => r.h));
      // the whole band where one lights, half of it where both do
      expect(heights.size).toBe(2);
    }
  });
});

describe('the Legend part lists the highlight channels', () => {
  it('each with its color and style, each cleared by name', async () => {
    const v = await make();
    v.addresses.highlightAddress(Range(2, 5), 'a');
    v.addresses.highlightAddress(Range(4, 8), 'b', 'band');
    const el = document.createElement('div');
    document.body.append(el);
    const dispose = mountLegend(el, v);
    const rows = [...el.querySelectorAll<HTMLElement>('[data-highlight]')];
    expect(rows.map((r) => r.dataset.highlight)).toEqual(['a', 'b']);
    expect(rows.map((r) => r.querySelector<HTMLElement>('[data-style]')!.dataset.style)).toEqual(['tint', 'band']);
    // the swatch is the channel's color
    const hex = (c: number) => `#${(c >>> 0).toString(16).padStart(8, '0')}`;
    expect(rows[0].querySelector<HTMLElement>('[data-style]')!.dataset.color).toBe(hex(v.highlightColor('a')));
    const clear = (name: string) =>
      [...el.querySelectorAll('button')].find((b) => b.getAttribute('aria-label') === `Clear the ${name} highlight`)!;
    clear('a').click();
    await v.settle();
    expect(v.highlight.names.value).toEqual(['b']);
    expect([...el.querySelectorAll<HTMLElement>('[data-highlight]')].map((r) => r.dataset.highlight)).toEqual(['b']);
    clear('b').click();
    expect(el.querySelectorAll('[data-highlight]').length).toBe(0);
    dispose();
  });
});
