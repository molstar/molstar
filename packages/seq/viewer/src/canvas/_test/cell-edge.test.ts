/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A letter cell's edge (core _spec/README.md "Render", Tokens: `seq-cell-edge`) drawn against the data, on canvas and
// in SVG: in dark, every cell drawn with its letter (the letters layer, a rows layer's letter cells) has a hairline of
// one device pixel in the token; in light, none.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { builtInTokens, colorHex, type DrawLists, recordDrawLists } from '@molstar/seqx-core/render';
import { renderSvg } from '@molstar/seqx-core/render/render-svg';
import { describe, expect, it } from 'vitest';
import { CanvasBackend } from '../backend.js';
import { FakeCanvas } from './recording.js';

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

/** The data pass's rects as numbers and colors: `[x, y, w, h]`, the fill, the stroke and its width. */
function rects(lists: DrawLists) {
  return recordDrawLists(lists, ['data'])
    .split('\n')
    .filter((l) => l.startsWith('rect '))
    .map((l) => {
      const p = l.slice(5).split(' ');
      const k = p.indexOf('stroke');
      return {
        box: p.slice(0, 4).map(Number),
        fill: p[4],
        stroke: k > 0 ? p[k + 1] : undefined,
        width: k > 0 ? Number(p[k + 2]) : 0,
      };
    });
}

/** The canvas calls drawing the data pass. */
function canvasCalls(lists: DrawLists, dpr: number): string[] {
  const base = new FakeCanvas();
  const backend = new CanvasBackend(base, new FakeCanvas(), new FakeCanvas(), { createCanvas: () => new FakeCanvas() });
  backend.resize(lists.width, lists.height, dpr);
  base.context.take();
  backend.draw(lists, true);
  return base.context.take();
}

const args = (call: string) =>
  call
    .slice(call.indexOf('(') + 1, -1)
    .split(',')
    .map(Number);

/** `inner` is `outer` inset by `d` on every side. */
function expectInset(inner: readonly number[], outer: readonly number[], d: number, what: string) {
  expect(inner[0], what).toBeCloseTo(outer[0] + d, 1);
  expect(inner[1], what).toBeCloseTo(outer[1] + d, 1);
  expect(inner[2], what).toBeCloseTo(outer[2] - 2 * d, 1);
  expect(inner[3], what).toBeCloseTo(outer[3] - 2 * d, 1);
}

const SvgRect = /<rect x="([^"]+)" y="([^"]+)" width="([^"]+)" height="([^"]+)" ([^/]*)\/>/g;

describe('a letter cell’s edge', () => {
  for (const [name, spec, count] of [
    ['letters', lettersSpec, Letters.length],
    // three rows of six, two gaps: 16 letter cells
    ['rows', rowsSpec, 16],
  ] as const)
    for (const dpr of [1, 1.5, 2])
      it(`${name} at dpr ${dpr}: a hairline inside each cell in dark, none in light; canvas and SVG`, async () => {
        const edge = colorHex(builtInTokens('dark').colors['seq-cell-edge']);
        const px = 1 / dpr;
        const dark = await renderSvg(spec(), { size: Size, theme: 'dark', dpr, legend: 'none' });
        const light = await renderSvg(spec(), { size: Size, theme: 'light', dpr, legend: 'none' });

        // the draw lists: each edge right after its cell's fill, unfilled, a device pixel wide, on the cell inset by
        // half a device pixel (inside the cell: no neighbour's fill paints over it, nothing spills out)
        const all = rects(dark.lists);
        const edges = all.map((r, i) => [r, i] as const).filter(([r]) => r.stroke === edge);
        expect(edges.length, name).toBe(count);
        for (const [r, i] of edges) {
          const cell = all[i - 1];
          expect(cell.stroke).toBeUndefined();
          expect(r.fill).toBe('none');
          expect(r.width).toBeCloseTo(px, 1);
          expectInset(r.box, cell.box, px / 2, `${name} ${dpr}`);
        }
        expect(rects(light.lists).some((r) => r.stroke === edge)).toBe(false);

        // SVG: an unfilled rect in the token a cell, inset the same; none in light
        const rgb = edge.slice(0, 7);
        const svgRects = [...dark.svg.matchAll(SvgRect)];
        const svgEdges = svgRects.map((m, i) => [m, i] as const).filter(([m]) => m[5].includes(`stroke="${rgb}"`));
        expect(svgEdges.length).toBe(count);
        for (const [m, i] of svgEdges) {
          expect(m[5]).toContain('fill="none"');
          const prev = svgRects[i - 1];
          expectInset(m.slice(1, 5).map(Number), prev.slice(1, 5).map(Number), px / 2, `svg ${name} ${dpr}`);
        }
        expect(light.svg).not.toContain(`stroke="${rgb}"`);

        // canvas: each cell filled, then its edge stroked a device pixel wide on the inset rect
        const calls = canvasCalls(dark.lists, dpr);
        const strokes = calls.map((c, i) => [c, i] as const).filter(([c]) => c.startsWith('strokeRect('));
        expect(strokes.length).toBe(count);
        expect(calls.some((c) => c.startsWith('lineWidth=') && Math.abs(Number(c.slice(10)) - px) < 1e-3)).toBe(true);
        for (const [c, i] of strokes) {
          let f = i - 1;
          while (!calls[f].startsWith('fillRect(')) f--;
          expectInset(args(c), args(calls[f]), px / 2, `canvas ${name} ${dpr}`);
        }
        expect(canvasCalls(light.lists, dpr).filter((c) => c.startsWith('strokeRect(')).length).toBe(0);
      });
});
