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
import { DefaultRenderConfig, recordDrawLists } from '../index.js';
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
