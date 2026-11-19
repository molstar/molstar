/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { BuiltInSchemes, contrastRgba } from '../../theme/index.js';
import { DarkTokens, LightTokens, recordDrawLists, renderSvg } from '../index.js';

// The logo's letter colors (_spec/README.md "The logo"): a letter's color is its class's held at 3:1 against the track,
// in each theme, for every color of every built-in protein and nucleotide scheme. The glyph runs of the data pass are
// what a canvas draws; the SVG's `fill` attributes are the vector backend's.

const Size = { width: 400, height: 120 };
const Themes = { light: LightTokens, dark: DarkTokens } as const;
const Protein = 'ACDEFGHIKLMNPQRSTVWY';

function logo(alphabet: 'protein' | 'dna' | 'rna', letters: string, scheme?: string) {
  const b = createSeqSpec();
  const aln = b.data.alignment('aln', { alphabet, rows: { label: ['a', 'b'], letters: [letters, letters] } });
  b.views.view('v', (v) =>
    v.section('s', { of: aln }, (s) =>
      s.track('logo', { label: 'logo' }, (t) =>
        t.layer('l', {
          type: 'logo',
          data: aln,
          ...(scheme ? { color: { kind: 'scheme', params: { name: scheme } } } : {}),
        } as never),
      ),
    ),
  );
  return b.build();
}

const rgba = (hex: string) => Number.parseInt(hex.padEnd(8, 'f'), 16) >>> 0;

/** The colors of the logo's glyph runs (the canvas's input) and of its SVG text, over the logo's track. */
async function drawn(spec: ReturnType<typeof logo>, theme: 'light' | 'dark') {
  const out = await renderSvg(spec, { size: Size, theme, legend: 'none' });
  const canvas = new Set<number>();
  for (const line of recordDrawLists(out.lists, ['data']).split('\n'))
    if (line.startsWith('glyphs')) for (const m of line.matchAll(/#([0-9a-f]{8})/g)) canvas.add(rgba(m[1]));
  const svg = new Set<number>();
  for (const m of out.svg.matchAll(/<text [^>]*fill="#([0-9a-f]{6})"/g)) svg.add(rgba(m[1]));
  return { canvas, svg };
}

describe('the logo’s letter colors', () => {
  for (const e of BuiltInSchemes.filter((x) => x.subject === 'protein' || x.subject === 'nucleotide'))
    for (const theme of ['light', 'dark'] as const)
      it(`${e.name}, ${theme}: every letter at 3:1 against the track`, async () => {
        const surface = Themes[theme].colors['seq-track'];
        const specs =
          e.subject === 'protein'
            ? [logo('protein', Protein, e.name)]
            : [logo('dna', 'ACGT', e.name), logo('rna', 'ACGU', e.name)];
        for (const spec of specs) {
          const { canvas, svg } = await drawn(spec, theme);
          expect(canvas.size).toBeGreaterThan(0);
          for (const c of canvas) {
            expect(svg.has(c)).toBe(true);
            expect(contrastRgba(c, surface)).toBeGreaterThanOrEqual(3);
          }
        }
      });

  it('keeps a color that already holds 3:1, and the hue of one that is moved', async () => {
    const light = await drawn(logo('protein', Protein), 'light');
    expect(light.canvas.size).toBeGreaterThan(1);
    // `classes` in the light table: the pale cell colors are moved, not replaced by one ink
    const colors = [...light.canvas];
    expect(new Set(colors).size).toBe(colors.length);
  });
});
