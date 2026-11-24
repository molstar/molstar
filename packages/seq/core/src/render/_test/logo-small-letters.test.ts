/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { type Frame, type LogoLay, MinLogoGlyphPx, recordDrawLists } from '../index.js';
import { settle } from './helpers.js';
import { layerLay } from './small-helpers.js';

// The logo's small letters (_spec/README.md "The logo", the defaults audit): a stacked letter is drawn as a letter only
// from about 9 px of glyph height; below it the stack's bottom is a flat block of the letter's color, its column wide
// and its height tall, never a pile of glyphs.

const Size = { width: 400, height: 400 };

/** Ten rows of one residue, five of a second, three of a third and two of a fourth: a low stack that is tall in bits. */
function spec(trackHeight: number) {
  const col = 'A'.repeat(10) + 'C'.repeat(5) + 'D'.repeat(3) + 'E'.repeat(2);
  const rows = { label: [] as string[], letters: [] as string[] };
  for (let i = 0; i < col.length; i++) {
    rows.label.push(`r${i}`);
    rows.letters.push(col[i].repeat(12));
  }
  const b = createSeqSpec();
  const aln = b.data.alignment('aln', { alphabet: 'protein', rows });
  b.views.view('v', (v) =>
    v.section('s', { of: aln }, (s) =>
      s.track('logo', { label: 'logo', height: trackHeight }, (t) => t.layer('l', { type: 'logo', data: aln })),
    ),
  );
  return b.build();
}

async function drawn(height: number) {
  const s = await settle(spec(height));
  const f: Frame = await s.frame({ size: Size });
  const lay = layerLay<LogoLay>(f, 0);
  const lines = recordDrawLists(f.lists, ['data']).split('\n');
  return { lay, glyphRuns: lines.filter((l) => l.startsWith('glyphs')).length, lines };
}

describe('the logo’s small letters', () => {
  it('a stack too short for glyphs is flat blocks: no glyph run, a block of each letter’s color and height', async () => {
    const { lay, glyphRuns } = await drawn(24);
    expect(lay.n).toBeGreaterThan(0);
    expect(lay.flatN).toBe(lay.n);
    expect(glyphRuns).toBe(0);
    for (let i = 0; i < lay.n; i++) {
      const h = lay.vs[i] * lay.fontSize * 0.72;
      expect(h).toBeLessThan(MinLogoGlyphPx);
    }
    expect(lay.cellW).toBeGreaterThan(0);
  });

  it('a tall stack keeps its letters as glyphs, and only the small ones at its bottom are blocks', async () => {
    const { lay, glyphRuns } = await drawn(400);
    expect(glyphRuns).toBeGreaterThan(0);
    let tall = 0;
    for (let i = 0; i < lay.n; i++) {
      const h = lay.vs[i] * lay.fontSize * 0.72;
      expect(lay.flat[i] === 1).toBe(h < MinLogoGlyphPx);
      if (!lay.flat[i]) tall++;
    }
    expect(tall).toBeGreaterThan(0);
  });
});
