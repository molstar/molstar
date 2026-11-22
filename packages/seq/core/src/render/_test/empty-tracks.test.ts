/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A track with no items never draws a full-height blank well (_spec/README.md "Render", Tracks layout, Empty tracks):
// by default a one-line row naming it with a quiet "no items" (`track.empty: 'note'`), the row alone (`collapse`), or
// its well as any track's (`keep`); a track still waiting for its source is pending, never empty; a track with a
// ruler, letters or any layer of a fixed height is never empty.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { GffRuntime } from '../../runtimes/gff.js';
import { type BlockPlacement, DefaultRenderConfig, recordDrawLists } from '../index.js';
import { lanesHeight } from '../sizing.js';
import { settle } from './helpers.js';

const Size = { width: 800, height: 600 };

/** A section of a full track, an empty one (no rows) and a ruler; `empty` the empty track's own. */
function spec(empty?: 'note' | 'collapse' | 'keep'): Node {
  const b = createSeqSpec({ title: 'empty tracks' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(100) });
  const full = b.data.ranges('full', { target: s, rows: { start: [1, 20], end: [10, 40], id: ['a', 'b'] } });
  const none = b.data.ranges('none', { target: s, rows: { start: [], end: [], id: [] } });
  b.views.view('main', (v) => {
    const sec = v.section('s', { of: s });
    sec.track('full', { label: 'Full' }).layer('l', { type: 'ranges', data: full });
    sec
      .track('none', { label: 'Promoters', ...(empty ? { empty } : {}) } as never)
      .layer('l', { type: 'ranges', data: none });
    sec.track('ruler', { label: '' }).layer('r', { type: 'ruler' });
  });
  return b.build();
}

const heights = async (n: Node, config = {}) => {
  const s = await settle(n, [], config);
  const f = await s.frame({ size: Size, settled: true });
  return { f, h: f.layout.sections[0].tracks.map((t) => t.box.height) };
};

describe('a track with no items is one line, not a blank well', () => {
  it('collapses to its row with a quiet note by default; its neighbours keep theirs', async () => {
    const { f, h } = await heights(spec());
    expect(h[1]).toBe(DefaultRenderConfig.track.emptyPx);
    expect(h[1]).toBeLessThan(lanesHeight(1, DefaultRenderConfig));
    expect(h[0]).toBe(lanesHeight(1, DefaultRenderConfig));
    expect(f.layout.sections[0].tracks.map((t) => t.empty)).toEqual([undefined, 'note', undefined]);
    expect(recordDrawLists(f.lists)).toContain('text "no items"');
  });

  it('`collapse` draws the row alone; `keep` the well as any track’s; the config’s default is the host’s', async () => {
    const c = await heights(spec('collapse'));
    expect(c.h[1]).toBe(DefaultRenderConfig.track.emptyPx);
    expect(recordDrawLists(c.f.lists)).not.toContain('"no items"');
    const k = await heights(spec('keep'));
    expect(k.h[1]).toBe(lanesHeight(1, DefaultRenderConfig));
    expect(recordDrawLists(k.f.lists)).not.toContain('"no items"');
    const host = await heights(spec(), { track: { empty: 'keep' } });
    expect(host.h[1]).toBe(lanesHeight(1, DefaultRenderConfig));
    // the track's own wins over the host's
    const own = await heights(spec('collapse'), { track: { empty: 'keep' } });
    expect(own.h[1]).toBe(DefaultRenderConfig.track.emptyPx);
  });

  it('a track whose source has not answered is pending, never empty', async () => {
    const b = createSeqSpec({ title: 'pending' });
    const file = b.data.source('gff', {
      format: 'gff3',
      value: '##gff-version 3\n##sequence-region c 1 400\nc\t.\tmRNA\t1\t90\t.\t+\t.\tID=x;Name=T\n' as never,
    });
    const chr = b.data.sequence('c', { alphabet: 'dna', source: file, name: 'c' });
    const tx = b.data.ranges('tx', { target: chr, source: file });
    b.views.view('main', (v) => {
      v.section('c', { of: chr }).track('tx', {}).layer('tx', { type: 'ranges', data: tx });
    });
    const s = await settle(b.build(), [GffRuntime]);
    const pending = await s.frame({ size: Size, settled: false });
    expect(pending.layout.sections[0].tracks[0].empty).toBeUndefined();
    expect(pending.layout.sections[0].tracks[0].box.height).toBe(lanesHeight(1, DefaultRenderConfig));
  });
});

/** A stack of two blocks, a track of items on the first block's system only, a track of both, and the letters. */
function stack(): Node {
  const b = createSeqSpec({ title: 'stack' });
  const a = b.data.sequence('a', { alphabet: 'dna', letters: 'ACGT'.repeat(10) });
  const c = b.data.sequence('c', { alphabet: 'dna', letters: 'ACGT'.repeat(10) });
  const ra = b.data.ranges('ra', { target: a, rows: { start: [1, 20], end: [10, 30], id: ['a', 'b'] } });
  const rc = b.data.ranges('rc', { target: c, rows: { start: [2], end: [8], id: ['x'] } });
  b.views.view('main', (v) => {
    const sec = v.section('s', { arrange: 'stack' } as never);
    sec.block('a', { of: a, label: 'A' });
    sec.block('c', { of: c, label: 'C' });
    sec.track('only-a', { label: 'Only A' }).layer('l', { type: 'ranges', data: ra });
    const t = sec.track('both', { label: 'Both' });
    t.layer('a', { type: 'ranges', data: ra });
    t.layer('c', { type: 'ranges', data: rc });
    sec.track('letters', { label: 'Sequence' }).layer('l', { type: 'letters' });
  });
  return b.build();
}

describe('a stacked block leaves out its tracks with no items for it (`blocks.emptyTracks`)', () => {
  const bandsOf = async (config = {}) => {
    const s = await settle(stack(), [], config);
    const f = await s.frame({ size: Size, settled: true });
    const p = f.layout.sections[0].placement as BlockPlacement;
    const st = p.strips;
    return {
      p,
      f,
      h: Array.from({ length: st.n }, (_, i) => Array.from({ length: st.bands }, (_, j) => st.bandH[i * st.bands + j])),
    };
  };

  it('the second block draws the tracks it has items in and its letters; the first keeps every track', async () => {
    const { h, p } = await bandsOf();
    expect(h[0].every((x) => x > 0)).toBe(true);
    expect(h[1][0]).toBe(0);
    expect(h[1][1]).toBeGreaterThan(0);
    expect(h[1][2]).toBeGreaterThan(0);
    // its row is shorter by the track left out, its label row kept
    expect(p.strips.height[1]).toBeLessThan(p.strips.height[0]);
    expect(p.labels[1]).toBe('C');
  });

  it('`keep` draws every track in every row', async () => {
    const { h } = await bandsOf({ blocks: { emptyTracks: 'keep' } });
    expect(h.every((row) => row.every((x) => x > 0))).toBe(true);
  });
});
