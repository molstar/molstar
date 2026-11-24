/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { drawListsToSvg, type Frame, recordDrawLists } from '../index.js';
import { settle } from './helpers.js';

// A letter is centred in its whole cell, not in the part of it the view shows: panning moves each letter by the
// pan, and the spacing between letters stays one advance a cell, at the edges too. The fills, the hover and the
// selection stay cut to the view; a letter that overhangs is clipped to the band's box.

const Size = { width: 640, height: 200 };

async function frameAt(start: number, end: number): Promise<Frame> {
  const b = createSeqSpec({ title: 'whole cell' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGTTGCAAC'.repeat(12) });
  b.views.view('v', {}, (v) => {
    v.section('sec', { of: s }).track('letters').layer('letters', { type: 'letters' });
  });
  const st = await settle(b.build(), []);
  return st.frame({ size: Size, viewport: { start, end } });
}

/** The x of each letter of the data pass's glyph runs, and the data box the clips name. */
function glyphXs(f: Frame): number[] {
  const out: number[] = [];
  for (const line of recordDrawLists(f.lists, ['data']).split('\n')) {
    const m = /^glyphs "[^"]+" [\d.]+ per (.*)$/.exec(line);
    if (!m) continue;
    for (const part of m[1].split(', ')) out.push(Number(part.split(' ')[1]));
  }
  return out;
}

function dataBox(f: Frame): { left: number; right: number } {
  // the letters' own clip: the one pushed right before their glyph run
  const m = /^clip ([\d.-]+) [\d.-]+ ([\d.-]+) [\d.-]+\nglyphs /m.exec(recordDrawLists(f.lists, ['data']))!;
  return { left: Number(m[1]), right: Number(m[1]) + Number(m[2]) };
}

describe('a letter is centred in its whole cell', () => {
  it('at a pan where the first and last cells are partly out, the spacing is one advance a cell, edges included', async () => {
    const f = await frameAt(10.4, 20.6);
    const xs = glyphXs(f);
    expect(xs.length).toBeGreaterThanOrEqual(10);
    const advance = dataBox(f).right - dataBox(f).left;
    const cell = advance / 10.2;
    for (let i = 1; i < xs.length; i++) expect(xs[i] - xs[i - 1]).toBeCloseTo(cell, 0);
    // the first cell (position 10) starts 0.4 cells left of the box: its letter is at its whole cell's middle
    expect(xs[0]).toBeCloseTo(dataBox(f).left + (0.5 - 0.4) * cell, 0);
  });

  it('every letter moves by the pan, the edge cells too', async () => {
    const a = glyphXs(await frameAt(10.4, 20.6));
    const b = glyphXs(await frameAt(10.6, 20.8));
    expect(b.length).toBe(a.length);
    const cell = a[1] - a[0];
    for (let i = 0; i < a.length; i++) expect(Math.abs(a[i] - b[i] - 0.2 * cell)).toBeLessThanOrEqual(1);
  });

  it('a cell wider than the view centres its letter off-screen, not in view', async () => {
    const f = await frameAt(10.05, 10.15);
    const { left, right } = dataBox(f);
    const xs = glyphXs(f);
    expect(xs.length).toBe(1);
    // the cell spans the whole view and more: its middle is a few hundred pixels right of the box
    expect(xs[0]).toBeGreaterThan(right);
    expect(xs[0]).toBeGreaterThan(left);
  });

  it('the SVG is the canvas draw list: the same glyph x, the letters clipped to the box', async () => {
    const f = await frameAt(10.4, 20.6);
    const svg = drawListsToSvg(f.lists);
    const xs = glyphXs(f);
    const block = /<text y="[^"]+" text-anchor="middle"[^>]*>((?:<tspan[^>]*>[^<]*<\/tspan>)+)<\/text>/.exec(svg);
    expect(block).not.toBeNull();
    const fromSvg = [...block![1].matchAll(/<tspan x="([\d.-]+)"/g)].map((t) => Number(t[1]));
    expect(fromSvg.length).toBe(xs.length);
    for (let i = 0; i < xs.length; i++) expect(fromSvg[i]).toBeCloseTo(xs[i], 1);
    // the edge cell's letter at its whole cell's middle (0.4 of the cell is out on the left)
    const { left, right } = dataBox(f);
    expect(fromSvg[0]).toBeCloseTo(left + (0.5 - 0.4) * ((right - left) / 10.2), 0);
    // the letters' text sits in a clip group whose box is the row's data box (the wells' clip, the first of the pass)
    const group = /<g clip-path="url\(#(c\d+)\)"><text y="[^"]+" text-anchor="middle"[^>]*><tspan/.exec(svg);
    expect(group).not.toBeNull();
    const box = new RegExp(`<clipPath id="${group![1]}"><rect x="([\\d.-]+)" y="[\\d.-]+" width="([\\d.-]+)"`).exec(
      svg,
    );
    expect(box).not.toBeNull();
    const well = /^clip ([\d.-]+) [\d.-]+ ([\d.-]+) /m.exec(recordDrawLists(f.lists, ['data']))!;
    expect(Number(box![1])).toBeCloseTo(Number(well[1]), 1);
    expect(Number(box![2])).toBeCloseTo(Number(well[2]), 1);
  });
});
