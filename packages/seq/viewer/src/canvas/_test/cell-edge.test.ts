/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A letter cell's edge (core _spec/README.md "Render", Tokens: `seq-cell-edge`) drawn against the data, on canvas and
// in SVG: in dark, every cell drawn with its letter (the letters layer, a rows layer's letter cells) has a hairline of
// one device pixel in the token; in light, none.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { builtInTokens, colorHex, type DrawLists, recordDrawLists, renderSvg } from '@molstar/seqx-core/render';
import { describe, expect, it } from 'vitest';
import { CanvasBackend } from '../backend.js';
import { FakeCanvas } from './recording.js';

const Dpr = 2;
const Size = { width: 600, height: 200 };
const Letters = 'MEEPQSDPSV';

function lettersSpec() {
  const b = createSeqSpec();
  const s = b.data.sequence('s', { alphabet: 'protein', letters: Letters });
  b.views.view('main', (v) => v.section('p', { of: s }).track('t').layer('l', { type: 'letters' }));
  return b.build();
}

function rowsSpec() {
  const b = createSeqSpec();
  const aln = b.data.alignment('aln', {
    alphabet: 'protein',
    rows: { label: ['a', 'b', 'c'], letters: ['MVLSPA', 'MVHLTP', 'MGLS--'] },
  });
  b.views.view('v', (v) =>
    v.section('s', { of: aln }, (s) =>
      s.track('rows', { label: '', height: 90 }, (t) => t.layer('r', { type: 'rows', data: aln } as never)),
    ),
  );
  return b.build();
}

/** The data pass's rects with a letter's fill, as `[fill, stroke, width]` (the recording's spelling). */
function cells(lists: DrawLists): string[][] {
  return recordDrawLists(lists, ['data'])
    .split('\n')
    .filter((l) => l.startsWith('rect '))
    .map((l) => l.slice(5).split(' '));
}

/** The canvas calls drawing the data pass. */
function canvasCalls(lists: DrawLists): string[] {
  const base = new FakeCanvas();
  const backend = new CanvasBackend(base, new FakeCanvas(), new FakeCanvas(), { createCanvas: () => new FakeCanvas() });
  backend.resize(lists.width, lists.height, Dpr);
  base.context.take();
  backend.draw(lists, true);
  return base.context.take();
}

describe('a letter cell’s edge', () => {
  for (const [name, spec, count] of [
    ['letters', lettersSpec, Letters.length],
    // three rows of six, two gaps: 16 letter cells
    ['rows', rowsSpec, 16],
  ] as const)
    it(`${name}: a hairline in the token in dark, none in light, on canvas and in SVG`, async () => {
      const edge = colorHex(builtInTokens('dark').colors['seq-cell-edge']);
      const dark = await renderSvg(spec(), { size: Size, theme: 'dark', dpr: Dpr, legend: 'none' });
      const light = await renderSvg(spec(), { size: Size, theme: 'light', dpr: Dpr, legend: 'none' });

      // the draw lists: every cell drawn with its letter has the edge, one device pixel wide
      const edged = cells(dark.lists).filter((r) => r.join(' ').includes(edge));
      expect(edged.length, name).toBe(count);
      for (const r of edged) expect(r[r.indexOf(edge) + 1], r.join(' ')).toBe('0.5');
      expect(cells(light.lists).some((r) => r.join(' ').includes(edge))).toBe(false);

      // SVG: a stroke in the token, 1 / dpr wide, a cell each; none in light
      const rgb = edge.slice(0, 7);
      const strokes = dark.svg.match(new RegExp(`stroke="${rgb}" stroke-width="0.5"`, 'g')) ?? [];
      expect(strokes.length).toBe(edged.length);
      expect(light.svg).not.toContain(`stroke="${rgb}"`);

      // canvas: a stroked rect a cell, in the token at 1 / dpr
      const calls = canvasCalls(dark.lists);
      expect(calls.filter((c) => c.startsWith('strokeRect(')).length).toBe(edged.length);
      expect(calls).toContain(`lineWidth=${1 / Dpr}`);
      const lc = canvasCalls(light.lists);
      expect(lc.filter((c) => c.startsWith('strokeRect(')).length).toBe(0);
    });
});
