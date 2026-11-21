/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { letterCap, letterSize } from '../config.js';
import { DefaultRenderConfig, fontSize, type Frame, recordDrawLists } from '../index.js';
import { p53TracksSpec } from './fixtures.js';
import { settle } from './helpers.js';

// A letter is never larger than the chrome's label (the foundations' type: optically matched, capped by
// `letters.maxPx`, 0.9 of the label's size for the mono face when absent): at the closest zoom, where the cell would
// take any size, every path that draws letters (a track, a ring, an alignment's rows, a text figure) stays at the cap,
// in the canvas's lists and so in the SVG's; a host's `maxPx` lifts it.

const Label = fontSize(DefaultRenderConfig.labelFont);
const Cap = letterSize(DefaultRenderConfig, 1000, 1000);
const Size = { width: 640, height: 700 };

/** The sizes of the glyph runs a frame's data pass draws. */
function glyphSizes(f: Frame): number[] {
  return recordDrawLists(f.lists, ['data'])
    .split('\n')
    .map((l) => /^glyphs "(\d+(?:\.\d+)?)px mono"/.exec(l))
    .filter((m): m is RegExpExecArray => !!m)
    .map((m) => Number(m[1]));
}

async function closest(spec: Node, config = {}, viewport?: { start: number; end: number }) {
  const s = await settle(spec, [], config as never);
  return glyphSizes(await s.frame({ size: Size, ...(viewport ? { viewport } : {}) }));
}

function section(layout?: 'rings' | 'text'): Node {
  const b = createSeqSpec({ title: 'cap' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGTTGCAAC'.repeat(12), cyclic: layout === 'rings' });
  b.views.view('v', layout === 'text' ? { layout: 'text' } : {}, (v) => {
    v.section('sec', { of: s, ...(layout === 'rings' ? { layout: 'rings' } : {}) })
      .track('letters')
      .layer('letters', { type: 'letters' });
  });
  return b.build();
}

function rows(): Node {
  const b = createSeqSpec();
  const aln = b.data.alignment('aln', {
    alphabet: 'protein',
    rows: { label: ['A/1-6', 'B/1-6'], letters: ['MVLSPA', 'MVHLTP'] },
  });
  b.views.view('v', (v) =>
    v.section('s', { of: aln }, (s) =>
      s.track('rows', { label: '' }, (t) => t.layer('r', { type: 'rows', data: aln } as never)),
    ),
  );
  return b.build();
}

describe('letters are capped at the label size', () => {
  it('the cap is 0.9 of the label’s size by default, and the scale’s largest size within it', () => {
    expect(letterCap(DefaultRenderConfig)).toBeCloseTo(Label * 0.9, 6);
    expect(Cap).toBeLessThanOrEqual(Label);
    expect(Cap).toBe(10);
    expect(letterCap({ ...DefaultRenderConfig, letters: { ...DefaultRenderConfig.letters, maxPx: 13 } })).toBe(13);
  });

  const paths: [string, () => Node, { start: number; end: number } | undefined][] = [
    ['a track', () => p53TracksSpec(), { start: 100, end: 108 }],
    ['a ring', () => section('rings'), undefined],
    ['a text figure', () => section('text'), undefined],
    ['an alignment’s rows', rows, undefined],
  ];
  for (const [name, spec, viewport] of paths) {
    it(`${name}: no letter exceeds the label size at the closest zoom`, async () => {
      const sizes = await closest(spec(), {}, viewport);
      expect(sizes.length).toBeGreaterThan(0);
      for (const size of sizes) expect(size).toBeLessThanOrEqual(Label);
    });
    it(`${name}: a host’s maxPx lifts the cap`, async () => {
      const sizes = await closest(spec(), { letters: { maxPx: 14 } }, viewport);
      expect(Math.max(...sizes)).toBeGreaterThan(Label);
    });
  }
});
