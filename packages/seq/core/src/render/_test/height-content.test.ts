/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `height: 'content'` (_spec/README.md "Render", Track heights and Section scroll): a track as tall as its lanes with
// no `track.maxPx` cap, a section as tall as its content with no `section.maxPx` cap; absent `height` unchanged.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { DefaultRenderConfig, type RangesLay } from '../index.js';
import { RowHeight } from '../rows-model.js';
import { renderSvg } from '../render-svg.js';
import { lanesHeight } from '../sizing.js';
import { settle } from './helpers.js';
import { geneTracks } from './section-scroll-fixtures.js';

const Size = { width: 800, height: 3000 };

/** One ranges track of `n` stacked lanes, its `height` as given. */
function lanesSpec(n: number, height?: number | 'content'): Node {
  const b = createSeqSpec({ title: 'content height' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(50) });
  b.views.view('main', (v) => {
    const d = b.data.ranges('r', {
      target: s,
      rows: {
        start: Array.from({ length: n }, (_, i) => i),
        end: Array.from({ length: n }, () => 150),
        id: Array.from({ length: n }, (_, i) => `r${i}`),
      },
    });
    v.section('s', { of: s })
      .track('t', { label: 'T', ...(height !== undefined ? { height } : {}) } as never)
      .layer('l', { type: 'ranges', data: d });
  });
  return b.build();
}

async function layoutOf(node: Node, config = {}) {
  const s = await settle(node, [], config);
  const f = await s.frame({ size: Size, settled: true });
  return f.layout.sections[0];
}

describe('a track at height: content', () => {
  it('is as tall as its 40 lanes, no inner scroll; absent, it stops at the cap', async () => {
    const full = lanesHeight(40, DefaultRenderConfig);
    expect(full).toBeGreaterThan(DefaultRenderConfig.track.maxPx);
    const content = await layoutOf(lanesSpec(40, 'content'));
    expect(content.tracks[0].box.height).toBe(full);
    const capped = await layoutOf(lanesSpec(40));
    expect(capped.tracks[0].box.height).toBe(DefaultRenderConfig.track.maxPx);
  });

  it('under the cap it is its content, as without it', async () => {
    const a = await layoutOf(lanesSpec(3, 'content'));
    const b = await layoutOf(lanesSpec(3));
    expect(a.tracks[0].box.height).toBe(b.tracks[0].box.height);
  });

  it('the canvas layout and the SVG export draw the same lanes', async () => {
    const r = await renderSvg(lanesSpec(40, 'content'), { size: { width: 800, height: 'auto' } } as never);
    expect(r.svg).toContain('>r39<');
  });
});

describe('a section at height: content', () => {
  it('past 600 px it is its content, never scrolling inside; absent, it scrolls at the cap', async () => {
    const capped = await layoutOf(geneTracks(30));
    expect(capped.placement.contentHeight).toBeGreaterThan(600);
    expect(capped.placement.box.height).toBe(600);
    expect(capped.placement.clip).toBe(true);
    const content = await layoutOf(geneTracks(30, 'content'));
    expect(content.placement.box.height).toBe(content.placement.contentHeight);
    expect(content.placement.clip).toBe(false);
    expect(content.sectionRail).toBeUndefined();
  });

  it('a text section at content shows every line, past section.maxPx; absent, it stops at the cap', async () => {
    const textSpec = (height?: 'content') => {
      const b = createSeqSpec({ title: 'text' });
      const p = b.data.sequence('p', { alphabet: 'protein', letters: 'ACDEFGHIKLMNPQRSTVWY'.repeat(300) });
      b.views.view('main', { layout: 'text' }, (v) => {
        v.section('p', { of: p, ...(height ? { height } : {}) } as never)
          .track('letters')
          .layer('letters', { type: 'letters' });
      });
      return b.build();
    };
    const capped = await layoutOf(textSpec());
    expect(capped.placement.box.height).toBe(DefaultRenderConfig.section.maxPx);
    const content = await layoutOf(textSpec('content'));
    expect(content.placement.box.height).toBeGreaterThan(DefaultRenderConfig.section.maxPx);
    expect(content.placement.box.height).toBe(content.placement.contentHeight);
  });
});

describe('a rows band at height: content', () => {
  const alignment = (height?: 'content') => {
    const b = createSeqSpec({ title: 'rows' });
    const a = b.data.alignment('a', {
      rows: {
        label: Array.from({ length: 60 }, (_, k) => `R${k}`),
        letters: Array.from({ length: 60 }, () => 'ACGTACGTACGTACGTACGT'),
      },
      alphabet: 'dna',
    } as never);
    b.views.view('main', (v) => {
      v.section('s', { of: a } as never)
        .track('rows', height ? { height } : {})
        .layer('r', { type: 'rows', data: a } as never);
    });
    return b.build();
  };
  const rowsBox = async (spec: Node, size: { width: number; height: number }) => {
    const s = await settle(spec);
    const f = await s.frame({ size, settled: true });
    return f.layout.sections[0].tracks.at(-1)!.box.height;
  };

  it('is every row, past rows.maxPx, in a tall frame and in a short one; absent, it stops at the cap', async () => {
    const all = 60 * RowHeight;
    expect(all).toBeGreaterThan(DefaultRenderConfig.rows.maxPx);
    expect(await rowsBox(alignment(), { width: 800, height: 1400 })).toBeLessThanOrEqual(
      DefaultRenderConfig.rows.maxPx,
    );
    expect(await rowsBox(alignment('content'), { width: 800, height: 1400 })).toBeGreaterThanOrEqual(all);
    expect(await rowsBox(alignment('content'), { width: 800, height: 300 })).toBeGreaterThanOrEqual(all);
  });
});

describe('a rows band at height: content with its rows not counted', () => {
  it('keeps the band’s own height, never the room (no rows.maxPx rule applies either)', async () => {
    const b = createSeqSpec({ title: 'pending' });
    // rows not counted: the alignment names no inline letters (as a source still pending)
    const a = b.data.alignment('a', { rows: { label: ['R0'] }, alphabet: 'dna' } as never);
    b.views.view('main', (v) => {
      v.section('s', { of: a } as never)
        .track('rows', { height: 'content' } as never)
        .layer('r', { type: 'rows', data: a } as never);
    });
    const s = await settle(b.build());
    const f = await s.frame({ size: { width: 800, height: 1400 }, settled: true });
    expect(f.layout.sections[0].tracks.at(-1)!.box.height).toBe(DefaultRenderConfig.heights.rows);
  });
});

describe('a ring’s track at height: content', () => {
  it('keeps rings.trackMaxPx: a ring has no page to grow', async () => {
    const b = createSeqSpec({ title: 'ring' });
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(100) });
    const r = b.data.ranges('r', {
      target: s,
      rows: {
        start: Array.from({ length: 20 }, (_, i) => i),
        end: Array.from({ length: 20 }, () => 150),
        id: Array.from({ length: 20 }, (_, i) => `r${i}`),
      },
    });
    b.views.view('main', (v) => {
      v.section('s', { of: s, layout: 'rings' } as never)
        .track('t', { height: 'content' } as never)
        .layer('l', { type: 'ranges', data: r });
    });
    const st = await settle(b.build());
    const f = await st.frame({ size: { width: 700, height: 700 }, settled: true });
    const lay = f.layout.sections[0].tracks[0].layers[0].lay as RangesLay;
    expect(lay.box.height).toBeLessThanOrEqual(DefaultRenderConfig.rings.trackMaxPx);
  });
});
