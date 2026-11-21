/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { letterCap, letterSize } from '../config.js';
import { createSeqContext } from '../../context/default.js';
import { DefaultRenderConfig, fontSize, type Frame, recordDrawLists } from '../index.js';
import { p53TracksSpec } from './fixtures.js';
import { settle } from './helpers.js';

// A letter's cap height equals a label's (the foundations' type): its size is capped by `letters.maxPx`, else by the
// label's size times `letters.capRatio`, the label face's cap height over the mono face's, measured by the host's text
// metrics port (1 where it measures none). At the closest zoom, where the cell would take any size, every path that
// draws letters (a track, a ring, an alignment's rows) stays at the cap itself, in the canvas's lists and so in the
// SVG's; a host's `maxPx` lifts it. A text figure keeps the scale: it is read as text.

const Label = fontSize(DefaultRenderConfig.labelFont);
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

/** Faces with cap heights per em: the label face's 0.72, the mono face's 0.70 (what a canvas measures of `H`). */
const Faces = { sans: 0.72, mono: 0.7 };
const measured = (font: string) => (/mono/.test(font) ? Faces.mono : Faces.sans) * 100;

describe('letters are capped at the label size', () => {
  it('without measured faces the cap is the label’s size, the cap itself a size (not the scale’s)', () => {
    expect(letterCap(DefaultRenderConfig)).toBe(Label);
    expect(letterSize(DefaultRenderConfig, 1000, 1000)).toBe(Label);
    // a text figure keeps the scale
    expect(letterSize(DefaultRenderConfig, 1000, 1000, false)).toBe(14);
    expect(letterCap({ ...DefaultRenderConfig, letters: { ...DefaultRenderConfig.letters, maxPx: 13 } })).toBe(13);
    // a narrow cell takes the scale's largest that fits, never past the cap
    expect(letterSize(DefaultRenderConfig, 8, 1000)).toBe(10);
  });

  it('with measured faces the cap heights are equal: the letter’s within half a pixel of the label’s', () => {
    const context = createSeqContext({ ports: { textMetrics: () => ({ measure: () => 0, capHeight: measured }) } });
    const config = context.renderConfig;
    expect(config.letters.capRatio).toBeCloseTo(Faces.sans / Faces.mono, 6);
    const letterCapHeight = letterCap(config) * Faces.mono;
    const labelCapHeight = Label * Faces.sans;
    expect(Math.abs(letterCapHeight - labelCapHeight)).toBeLessThan(0.5);
    // computed once a config
    expect(letterCap(config)).toBe(letterCap(config));
  });

  const paths: [string, () => Node, { start: number; end: number } | undefined][] = [
    ['a track', () => p53TracksSpec(), { start: 100, end: 108 }],
    ['a ring', () => section('rings'), undefined],
    ['an alignment’s rows', rows, undefined],
  ];
  for (const [name, spec, viewport] of paths) {
    it(`${name}: at the closest zoom a letter is the cap, no larger, in the lists`, async () => {
      const sizes = await closest(spec(), {}, viewport);
      expect(sizes.length).toBeGreaterThan(0);
      for (const size of sizes) expect(size).toBeLessThanOrEqual(Label);
      expect(Math.max(...sizes)).toBe(Label);
    });
    it(`${name}: a host’s maxPx lifts the cap`, async () => {
      const sizes = await closest(spec(), { letters: { maxPx: 14 } }, viewport);
      expect(Math.max(...sizes)).toBeGreaterThan(Label);
    });
  }
});
