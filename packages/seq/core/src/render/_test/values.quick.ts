/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The values layer's warm frames (_spec/README.md "Render", Layer types, `values`): a frame of each encoding, bars, line,
// area and points, by value and uniform, at a per-position step (dots) and at a device-column step, and the hover on a
// dot, allocate nothing once warm. By hand, never in CI: `pnpm -C packages/seq/core test:perf`.

import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { drawView, type FrameInput, type Hit, hitView } from '../index.js';
import { settle } from './helpers.js';
import { reusedHover } from './scenarios.js';

const Size = { width: 800, height: 80 };

function spec(n: number, layer: Record<string, unknown>) {
  const b = createSeqSpec({ title: 'values' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'A'.repeat(n) });
  const values = Float32Array.from({ length: n }, (_, i) => (i % 17 === 5 ? Number.NaN : 5 + 4 * Math.sin(i / 7)));
  const v = b.data.values('v', { target: s, values });
  b.views.view('main', (view) =>
    view
      .section('s', { of: s })
      .track('v')
      .layer('values', { type: 'values', data: v, ...layer }),
  );
  return b.build();
}

const Cases: [string, number, Record<string, unknown>][] = [];
for (const encoding of ['bars', 'line', 'area', 'points'])
  for (const [n, step] of [
    [60, 'dots'],
    [4000, 'columns'],
  ] as const)
    for (const color of ['uniform', 'by value'])
      Cases.push([
        `${encoding}, ${color}, ${step}`,
        n,
        { encoding, ...(color === 'by value' ? { color: theme.field('value', { domain: [0, 10] }) } : {}) },
      ]);

// two reference lines (a dashed and a labelled, a solid one) beside each encoding: nothing allocated once warm
for (const encoding of ['bars', 'line', 'area', 'points'])
  for (const [n, step] of [
    [60, 'dots'],
    [4000, 'columns'],
  ] as const)
    Cases.push([
      `${encoding}, two reference lines, ${step}`,
      n,
      {
        encoding,
        reference: [
          { at: 7.3, label: 'genome-wide' },
          { at: 5, style: 'solid', color: '#d62728' },
        ],
      },
    ]);

describe('a warm frame of the values layer allocates nothing', () => {
  for (const [name, n, layer] of Cases)
    it(name, async ({ task }) => {
      const s = await settle(spec(n, layer));
      const input: FrameInput = { size: Size, dpr: 2 };
      for (let i = 0; i < 4; i++) await s.frame(input);
      const { r } = quickCase(task, `values ${name}`, () => void s.renderer.frame(input));
      console.info(`values ${name}: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
      expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
      // its draw alone, every mark written
      const f = s.renderer.frame(input);
      const d = quickCase(task, `values ${name} draw`, () => drawView(f.layout, f.lists), {});
      console.info(`values ${name} draw: ${d.r.bytesPerStep.toFixed(0)} bytes a frame (${d.r.note})`);
      expect(d.r.bytesPerStep, d.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    });

  it('the hover on a dot of points', async ({ task }) => {
    const s = await settle(spec(60, { encoding: 'points' }));
    const input: FrameInput = { size: Size, dpr: 2 };
    for (let i = 0; i < 4; i++) await s.frame(input);
    const f = s.renderer.frame(input);
    const hits: Hit[] = [];
    for (let x = 120; x < 780; x += 11) {
      const h = hitView(f.layout, x, 30);
      if (h) hits.push(h);
    }
    expect(hits.length).toBeGreaterThan(20);
    const hover = reusedHover(hits);
    const at: FrameInput = { size: Size, dpr: 2, emphasis: hover.emphasis };
    const { r } = quickCase(task, 'values points hover', (i) => {
      hover.set(i);
      s.renderer.frame(at);
    });
    console.info(`values points hover: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});
