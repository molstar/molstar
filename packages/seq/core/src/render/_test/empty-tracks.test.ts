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
import { createSeqSpec, SeqLanguage } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { GffRuntime } from '../../runtimes/gff.js';
import { rangesLocation } from '@molstar/seqx-core/query';
import { type BlockPlacement, DefaultRenderConfig, hitView, recordDrawLists } from '../index.js';
import { renderSvg } from '../render-svg.js';
import { lanesHeight } from '../sizing.js';
import { settle } from './helpers.js';

const Size = { width: 800, height: 600 };

/** A section of a full track, an empty one (no rows) and a ruler; `height` the empty track's own. */
function spec(height?: number): Node {
  const b = createSeqSpec({ title: 'empty tracks' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(100) });
  const full = b.data.ranges('full', { target: s, rows: { start: [1, 20], end: [10, 40], id: ['a', 'b'] } });
  const none = b.data.ranges('none', { target: s, rows: { start: [], end: [], id: [] } });
  b.views.view('main', (v) => {
    const sec = v.section('s', { of: s });
    sec.track('full', { label: 'Full' }).layer('l', { type: 'ranges', data: full });
    sec
      .track('none', { label: 'Promoters', ...(height !== undefined ? { height } : {}) })
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

describe('the language has no `empty` on a track (the host’s config alone)', () => {
  it('the track kind’s params carry none', () => {
    const params = (SeqLanguage.kinds.track as { params?: Record<string, unknown> }).params ?? {};
    expect(Object.keys(params)).not.toContain('empty');
  });
});

describe('a track with no items is one line, not a blank well', () => {
  it('collapses to its row with a quiet note by default; its neighbours keep theirs', async () => {
    const { f, h } = await heights(spec());
    expect(h[1]).toBe(DefaultRenderConfig.track.emptyPx);
    expect(h[1]).toBeLessThan(lanesHeight(1, DefaultRenderConfig));
    expect(h[0]).toBe(lanesHeight(1, DefaultRenderConfig));
    expect(f.layout.sections[0].empties).toEqual([undefined, 'note', undefined]);
    expect(recordDrawLists(f.lists)).toContain('text "no items"');
  });

  it('`collapse` draws the row alone; `keep` the well as any track’s (the host’s config: no spec param)', async () => {
    const c = await heights(spec(), { track: { empty: 'collapse' } });
    expect(c.h[1]).toBe(DefaultRenderConfig.track.emptyPx);
    expect(recordDrawLists(c.f.lists)).not.toContain('"no items"');
    const k = await heights(spec(), { track: { empty: 'keep' } });
    expect(k.h[1]).toBe(lanesHeight(1, DefaultRenderConfig));
    expect(recordDrawLists(k.f.lists)).not.toContain('"no items"');
  });

  it('an explicit `height` and no items is empty too: the note row, unless `keep`', async () => {
    const own = await heights(spec(60));
    expect(own.h[1]).toBe(DefaultRenderConfig.track.emptyPx);
    expect(own.f.layout.sections[0].empties[1]).toBe('note');
    expect((await heights(spec(60), { track: { empty: 'keep' } })).h[1]).toBe(60);
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
    expect(pending.layout.sections[0].empties[0]).toBeUndefined();
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

describe('a source track settled with no items is empty, not a loading lane', () => {
  const none = () => geneSpecOf(['##gff-version 3', '##sequence-region c 1 400'].join('\n'));

  it('pending keeps the placeholder lane (EC-169); settled with no item it is the note row', async () => {
    const s = await settle(none(), [GffRuntime]);
    const pending = await s.frame({ size: Size, settled: false });
    expect(pending.layout.sections[0].empties[0]).toBeUndefined();
    expect(pending.layout.sections[0].tracks[0].box.height).toBe(lanesHeight(1, DefaultRenderConfig));
    const f = await s.frame({ size: Size, settled: true });
    const t = f.layout.sections[0].tracks[0];
    expect(f.layout.sections[0].empties[0]).toBe('note');
    expect(t.box.height).toBe(DefaultRenderConfig.track.emptyPx);
  });

  it('a source with items is never empty', async () => {
    const s = await settle(
      geneSpecOf(
        ['##gff-version 3', '##sequence-region c 1 400', 'c\t.\tmRNA\t1\t90\t.\t+\t.\tID=x;Name=T'].join('\n'),
      ),
      [GffRuntime],
    );
    const f = await s.frame({ size: Size, settled: true });
    expect(f.layout.sections[0].empties[0]).toBeUndefined();
  });
});

/** One section, one gene-model track over a GFF3 file's text. */
function geneSpecOf(text: string): Node {
  const b = createSeqSpec({ title: 'source' });
  const file = b.data.source('gff', { format: 'gff3', value: text as never });
  const chr = b.data.sequence('c', { alphabet: 'dna', source: file, name: 'c' });
  const tx = b.data.ranges('tx', { target: chr, source: file });
  b.views.view('main', (v) => {
    v.section('c', { of: chr })
      .track('tx', {})
      .layer('tx', { type: 'ranges', data: tx, line: 'gene' } as never);
  });
  return b.build();
}

describe('over a collapsed track: the emphasis and the export', () => {
  it('the hover, the selection and the cursor draw around the note row, and hit nothing in it', async () => {
    const s = await settle(spec());
    const frame = (emphasis: object) => s.renderer.frame({ size: Size, settled: true, emphasis } as never);
    await s.frame({ size: Size, settled: true });
    const sel = { selection: rangesLocation('/data/s', [0, 30]) };
    const f = frame({ ...sel, cursor: { system: '/data/s', position: 5 } });
    const [full, none] = f.layout.sections[0].tracks;
    // the empty row has no layer to hit: a point in it is no element
    const h = hitView(f.layout, none.box.x + 40, none.box.y + none.box.height / 2);
    expect(h?.element).toBeUndefined();
    // the selection tints the full track's row, none of the empty one's marks
    const rects = recordDrawLists(f.lists, ['selection'])
      .split('\n')
      .filter((l) => l.startsWith('rect '));
    expect(rects.length).toBeGreaterThan(0);
    // every tint inside the section (the one-line row is a band like any other, never a well of its old height)
    const sec = f.layout.sections[0];
    expect(
      rects.every((l) => Number(l.split(' ')[2]) >= sec.y - 4 && Number(l.split(' ')[2]) < sec.y + sec.height),
    ).toBe(true);
    expect(full.box.height).toBeGreaterThan(none.box.height);
    // a hover on the full track's item still hits it
    const item = hitView(f.layout, full.box.x + 12, full.box.y + full.box.height / 2);
    expect(item?.element).toBeDefined();
  });

  it('the SVG export draws the note row (and the canvas lists the same text)', async () => {
    const r = await renderSvg(spec(), { size: Size });
    expect(r.svg.split('>no items<').length - 1).toBe(1);
    expect(r.svg).toContain('>Promoters<');
    const kept = await renderSvg(spec(), { size: Size, config: { track: { empty: 'keep' } } } as never);
    expect(kept.svg).not.toContain('>no items<');
  });
});

/** A stack of three blocks over two tracks: each track has items on one block's system, the third block on neither. */
function stack3(): Node {
  const b = createSeqSpec({ title: 'stack3' });
  const a = b.data.sequence('a', { alphabet: 'dna', letters: 'ACGT'.repeat(10) });
  const c = b.data.sequence('c', { alphabet: 'dna', letters: 'ACGT'.repeat(10) });
  const e = b.data.sequence('e', { alphabet: 'dna', letters: 'ACGT'.repeat(10) });
  const ra = b.data.ranges('ra', { target: a, rows: { start: [1], end: [10], id: ['a'] } });
  const rc = b.data.ranges('rc', { target: c, rows: { start: [2], end: [8], id: ['x'] } });
  b.views.view('main', (v) => {
    const sec = v.section('s', { arrange: 'stack' } as never);
    sec.block('c', { of: c, label: 'C' });
    sec.block('a', { of: a, label: 'A' });
    sec.block('e', { of: e, label: 'E' });
    sec.track('t-a', { label: 'On A' }).layer('l', { type: 'ranges', data: ra });
    sec.track('t-c', { label: 'On C' }).layer('l', { type: 'ranges', data: rc });
  });
  return b.build();
}

describe('a stack applies the rule to every block, the first too', () => {
  it('an empty first block leaves its tracks out (a track sits in the first row that keeps it); a block with every band empty is its heading row alone', async () => {
    const s = await settle(stack3());
    const f = await s.frame({ size: Size, settled: true });
    const sec = f.layout.sections[0];
    const p = sec.placement as BlockPlacement;
    const st = p.strips;
    const band = (i: number, j: number) => st.bandH[i * st.bands + j];
    // block C (first): its `On A` track left out, `On C` kept
    expect([band(0, 0) > 0, band(0, 1) > 0]).toEqual([false, true]);
    // block A: `On A` kept, `On C` left out
    expect([band(1, 0) > 0, band(1, 1) > 0]).toEqual([true, false]);
    // block E: neither, its row the label row
    expect([band(2, 0), band(2, 1)]).toEqual([0, 0]);
    expect(st.height[2]).toBe(p.labelRow);
    // `On A` is laid where block A keeps it: its well and label there, not over block C's row
    expect(sec.tracks[0].box.y).toBe(st.bandY[1 * st.bands] + st.dy[1]);
    expect(sec.tracks[1].box.y).toBe(st.bandY[1] + st.dy[0]);
    // the first block's row is shorter than a full one
    expect(st.height[0]).toBeLessThan(st.height[1] + band(1, 0) + 1 + p.labelRow);
  });

  it('keeps one set of flags across frames (nothing made a frame), and none outside a stack', async () => {
    const s = await settle(stack3());
    const a = (await s.frame({ size: Size, settled: true })).layout.sections[0].placement as BlockPlacement;
    const flags = a.gone;
    expect(flags.length).toBeGreaterThan(0);
    const b = (await s.frame({ size: Size, settled: true, viewport: { start: 0, end: 20 } })).layout.sections[0]
      .placement as BlockPlacement;
    expect(b.gone).toBe(flags);
    const row = await settle(rowOfTwo());
    const rp = (await row.frame({ size: Size, settled: true })).layout.sections[0].placement as BlockPlacement;
    expect(rp.gone.length).toBe(0);
  });
});

/** The same blocks in a row (no `arrange`). */
function rowOfTwo(): Node {
  const b = createSeqSpec({ title: 'row' });
  const a = b.data.sequence('a', { alphabet: 'dna', letters: 'ACGT'.repeat(10) });
  const c = b.data.sequence('c', { alphabet: 'dna', letters: 'ACGT'.repeat(10) });
  b.views.view('main', (v) => {
    const sec = v.section('s', {});
    sec.block('a', { of: a, label: 'A' });
    sec.block('c', { of: c, label: 'C' });
    sec.track('letters', { label: 'Sequence' }).layer('l', { type: 'letters' });
  });
  return b.build();
}
