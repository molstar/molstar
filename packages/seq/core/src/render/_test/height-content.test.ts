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
import { DefaultRenderConfig } from '../index.js';
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
  const sectionOf = (n: number, height?: 'content') => {
    const node = geneTracks(n);
    return height === undefined
      ? node
      : (JSON.parse(JSON.stringify(node).replace('"of":', `"height":"content","of":`)) as Node);
  };

  it('past 600 px it is its content, never scrolling inside; absent, it scrolls at the cap', async () => {
    const capped = await layoutOf(sectionOf(30));
    expect(capped.placement.contentHeight).toBeGreaterThan(600);
    expect(capped.placement.box.height).toBe(600);
    expect(capped.placement.clip).toBe(true);
    const content = await layoutOf(sectionOf(30, 'content'));
    expect(content.placement.box.height).toBe(content.placement.contentHeight);
    expect(content.placement.clip).toBe(false);
    expect(content.sectionRail).toBeUndefined();
  });
});
