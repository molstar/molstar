/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// EC-138 (.seq/sources/seqstar/edge-cases.md): the hover outline (`seq-hover`) holds 3:1 (WCAG's for graphics) against
// what it outlines, in both themes: the surfaces the outline is drawn on (the page, the track, the sequence's cells,
// the missing and pending fills), the residue and nucleotide class fills a letter cell has, and a range's fill.

import { contrastRgba } from '@molstar/seqx-core/theme/index';
import { describe, expect, it } from 'vitest';
import { DarkTokens, LightTokens, type SeqTokenName, SeqTokenNames } from '../index.js';

const Outlined: readonly SeqTokenName[] = [
  'msui-bg',
  'seq-track',
  'seq-sequence',
  'seq-missing',
  'seq-pending',
  'seq-range',
  ...SeqTokenNames.filter((n) => n.startsWith('seq-aa-') || n.startsWith('seq-nt-')),
];

describe('EC-138: the hover outline against what it outlines', () => {
  for (const [name, tokens] of [
    ['light', LightTokens],
    ['dark', DarkTokens],
  ] as const) {
    it(`${name}: seq-hover holds 3:1 against each surface and fill it rings`, () => {
      expect(Outlined.length).toBeGreaterThan(12);
      const hover = tokens.colors['seq-hover'];
      const low = Outlined.map((n) => ({ n, ratio: contrastRgba(hover, tokens.colors[n]) })).filter((e) => e.ratio < 3);
      expect(low).toEqual([]);
    });
  }
});
