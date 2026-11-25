/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Letters on a ring (_spec/README.md "Render", Section layouts: the polar warp): each letter's center on its cell's
// middle angle and its band's middle radius, upright or along the ring, in the canvas's lists and so in the SVG's.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { fontSize, type Frame, type RingPlacement } from '../index.js';
import { settle } from './helpers.js';

const Size = { width: 640, height: 700 };
const Length = 120;

function spec(letters?: 'upright' | 'tangent') {
  const b = createSeqSpec({ title: 'ring letters' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGTTGCAAC'.repeat(Length / 10), cyclic: true });
  b.views.view('v', (v) => {
    v.section('ring', { of: s, layout: 'rings' }).track('letters').layer('letters', { type: 'letters' });
  });
  return { spec: b.build(), config: letters ? { rings: { letters } } : {} };
}

/** Each placed letter's center (its baseline's point turned back up by the center's height) on the ring. */
function centers(f: Frame) {
  const ring = f.layout.sections.find((x) => x.layout.name === 'rings')!;
  const p = ring.placement as RingPlacement;
  const g = f.lists.passes.data.glyphs;
  const out: { offset: number; radius: number }[] = [];
  for (let i = 0; i < g.n; i++) {
    const k = 0.35 * fontSize(f.lists.string(g.font[i]));
    for (let j = g.start[i]; j < g.start[i] + g.count[i]; j++) {
      const t = g.angle[j];
      const x = g.x[j] + Math.sin(t) * g.vs[j] * k;
      const y = g.ly[j] - Math.cos(t) * g.vs[j] * k;
      out.push({ offset: p.offsetAt(x, y), radius: Math.hypot(x - p.cx, y - p.cy) });
    }
  }
  return out;
}

describe('letters on a ring', () => {
  for (const letters of ['upright', 'tangent'] as const)
    it(`${letters}: each letter's center on its cell's middle angle and every center at one radius`, async () => {
      const { spec: node, config } = spec(letters);
      const s = await settle(node, [], config as never);
      const f = await s.frame({ size: Size });
      const c = centers(f);
      expect(c.length).toBe(Length);
      // within a twentieth of a cell: the cells' edges are snapped to device pixels in the strip, half a pixel
      for (const { offset } of c) expect(Math.abs((offset % 1) - 0.5)).toBeLessThan(0.05);
      const r = c.map((x) => x.radius);
      expect(Math.max(...r) - Math.min(...r)).toBeLessThan(0.05);
    });
});
