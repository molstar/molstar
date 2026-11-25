/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The wells and every emphasis as squircles in the canvas as in the SVG (the core's _spec/README.md "Render", Shapes):
// a frame's lists drawn by both backends, each squircle the core's curve (a cubic a round corner), each clip a well's.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { rangesLocation } from '@molstar/seqx-core/query';
import { drawListsToSvg, type Emphasis, type PassName, recordDrawLists } from '@molstar/seqx-core/render';
import { renderSvg } from '@molstar/seqx-core/render/render-svg';
import { describe, expect, it } from 'vitest';
import { CanvasBackend } from '../backend.js';
import { FakeCanvas } from './recording.js';

const Protein = '/data/p';

function spec() {
  const b = createSeqSpec({ title: 'squircles' });
  const p = b.data.sequence('p', {
    alphabet: 'protein',
    letters: 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP',
  });
  b.views.view('v', (v) => {
    const s = v.section('s', { of: p });
    s.track('letters', { label: 'p', padding: { top: 6 } }).layer('letters', { type: 'letters' });
  });
  return b.build();
}

/** A frame's lists drawn by the canvas (its base and overlay calls) and as SVG, and the pass's record. */
async function both(emphasis: Emphasis, pass: PassName) {
  const r = await renderSvg(spec(), { size: { width: 600, height: 80 }, viewport: { start: 0, end: 20 }, emphasis });
  const base = new FakeCanvas();
  const overlay = new FakeCanvas();
  const backend = new CanvasBackend(base, overlay, new FakeCanvas(), { createCanvas: () => new FakeCanvas() });
  backend.resize(600, 80, 1);
  base.context.take();
  overlay.context.take();
  backend.draw(r.lists, true);
  return {
    record: recordDrawLists(r.lists, [pass]),
    base: base.context.take(),
    over: overlay.context.take(),
    svg: drawListsToSvg(r.lists),
  };
}

describe('squircles in the canvas and the SVG', () => {
  it('the well and its clip: the same curve in both', async () => {
    const { record, base, svg } = await both({}, 'data');
    expect(record).toMatch(/clip 96 8 496 26 squircle 8\nrect 96 14 /);
    // the well 26 px (the data's 20 and a top padding of 6), the data below its padding; the clip traced as four
    // cubics before the canvas's clip, and as the SVG's clip path
    const at = base.indexOf('clip()');
    expect(at).toBeGreaterThan(0);
    const before = base.slice(0, at);
    expect(before.slice(before.lastIndexOf('beginPath()')).filter((c) => c.startsWith('bezierCurveTo')).length).toBe(4);
    expect(svg).toMatch(/<clipPath id="[^"]*"><path d="M 104 8 /);
  });

  it('a cell’s hover and a projected hover: squircles in both', async () => {
    for (const emphasis of [
      { hover: { layer: '/views/v/s/letters/letters', section: 's', system: Protein, position: 4 } },
      { projected: { hover: rangesLocation(Protein, [3, 7]), hoverSection: '/views/w/x' } },
    ] satisfies Emphasis[]) {
      const { record, over, svg } = await both(emphasis, 'hover');
      const rects = record.split('\n').filter((l) => l.startsWith('rect '));
      expect(rects.length).toBe(1);
      expect(rects[0]).toMatch(/ squircle [\d.]+$/);
      // the canvas: the outline traced as four cubics and stroked
      expect(over.filter((c) => c.startsWith('bezierCurveTo')).length).toBe(4);
      expect(over.filter((c) => c === 'stroke()').length).toBe(1);
      // the SVG: the outline a path of four cubics
      const strokes = [...svg.matchAll(/<path d="([^"]*)" fill="none" stroke="#[0-9a-f]{6}"[^>]*\/>/g)];
      expect(strokes.length).toBeGreaterThan(0);
      expect(strokes.at(-1)![1].match(/C/g)?.length).toBe(4);
    }
  });
});
