/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { DataStore } from '@molstar/seqx-core/data';
import { AlignmentFileRuntime } from '@molstar/seqx-core/data/alignment-readers';
import { attachSources, bundledTextMetrics, LightTokens, Renderer, recordDrawLists } from '@molstar/seqx-core/render';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import * as globins from '../globins/family.js';
import { nodeFixtureAssets } from '../node-assets.js';

// The Globins figure's labels (_spec/README.md "Render", Layer types, ruler; Text layout, Numbers): a position carries
// one label where a group's end number and the ruler's tick label fall together, and no two labels overlap.

async function drawn(width: number, spec: Node = globins.figure()) {
  const model = resolveSpec(spec, { runtimes: [AlignmentFileRuntime] });
  const view = resolveViews(model)[0];
  const store = new DataStore();
  await attachSources(model, store, { assets: nodeFixtureAssets().reader });
  const renderer = new Renderer({ model, view, data: store, tokens: LightTokens, text: bundledTextMetrics() });
  const input = { size: { width, height: 1200 }, dpr: 1, settled: true, figure: true };
  let f = renderer.frame(input);
  for (let i = 0; i < 6; i++) {
    store.want(f.wants);
    await store.whenSettled();
    f = renderer.frame(input);
  }
  return recordDrawLists(f.lists, ['labels'])
    .split('\n')
    .filter((l) => l.startsWith('text '));
}

interface Label {
  readonly text: string;
  readonly a: number;
  readonly b: number;
  readonly y: number;
}

/** The numbers drawn (a label of digits only) with their x extents and baselines. */
function numbers(lines: readonly string[]): Label[] {
  const out: Label[] = [];
  for (const l of lines) {
    const m = /^text "(\d+)" (\S+) (\S+) (start|middle|end) "[^"]*" \S+ width (\S+)/.exec(l);
    if (!m) continue;
    const x = Number(m[2]);
    const w = Number(m[5]);
    const a = m[4] === 'start' ? x : m[4] === 'end' ? x - w : x - w / 2;
    out.push({ text: m[1], a, b: a + w, y: Number(m[3]) });
  }
  return out;
}

describe('the figure’s labels', () => {
  for (const width of [420, 640, 900, 1300])
    it(`at ${width} px: a position carries one label and no two labels overlap`, async () => {
      const labels = numbers(await drawn(width));
      expect(labels.length).toBeGreaterThan(8);
      const clash: string[] = [];
      for (let i = 0; i < labels.length; i++)
        for (let j = i + 1; j < labels.length; j++) {
          const p = labels[i];
          const q = labels[j];
          // labels of two lines are never within a label's height of each other
          if (Math.abs(p.y - q.y) < 10 && p.a < q.b && q.a < p.b) clash.push(`${p.text}@${p.a} ${q.text}@${q.a}`);
        }
      expect(clash).toEqual([]);
    });

  it('the line’s start number is on the ruler’s label baseline: one row of numbers across the line', async () => {
    const lines = await drawn(900);
    const ys = (re: RegExp) => lines.filter((l) => re.test(l)).map((l) => Number(l.split(' ')[3]));
    const ruler = ys(/"mono"|10px mono/);
    const start = ys(/^text "1" .* end "12px sans"/);
    expect(ruler.length).toBeGreaterThan(0);
    expect(start).toHaveLength(1);
    // the first line's labels
    expect(Math.min(...ruler)).toBe(start[0]);
  });

  it('inline group labels: the end number in the gap, the ruler’s label of that position dropped', async () => {
    const b = createSeqSpec({ title: 'inline' });
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'M'.repeat(100) });
    b.views.view('main', (v) => {
      const sec = v.section('fig', { of: s, layout: 'text', group: 10, group_labels: 'inline' } as never);
      sec.track('ruler', { label: '' }).layer('ruler', { type: 'ruler' });
      sec.track('seq', { label: '' }).layer('l', { type: 'letters', data: s });
    });
    const lines = await drawn(900, b.build());
    expect(lines.filter((l) => l.startsWith('text "10" '))).toHaveLength(1);
  });

  it('a ruler of another numbering than the group labels print keeps its labels', async () => {
    const b = createSeqSpec({ title: 'other' });
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'M'.repeat(100) });
    b.views.view('main', (v) => {
      const sec = v.section('fig', { of: s, layout: 'text', group: 5 } as never);
      // the first layer naming a numbering is the one the text prints: 0-based; the second ruler counts from 1
      sec.track('r0', { label: '' }).layer('r0', { type: 'ruler', numbering: 'index' } as never);
      sec.track('r1', { label: '' }).layer('r1', { type: 'ruler' });
      sec.track('seq', { label: '' }).layer('l', { type: 'letters', data: s });
    });
    const lines = await drawn(1500, b.build());
    // "10" of the 0-based ruler (its unit 10) and of the 1-based one (the group's end, printed there as 9)
    expect(lines.filter((l) => l.startsWith('text "10" ') && l.includes('10px mono'))).toHaveLength(2);
    expect(lines.filter((l) => l.startsWith('text "9" ') && l.includes('10px sans'))).toHaveLength(1);
  });

  it('a group’s end number is drawn once on the first line, in its own style', async () => {
    const tens = (await drawn(900)).filter((l) => l.startsWith('text "10" '));
    expect(tens).toHaveLength(1);
    expect(tens[0]).toContain('10px sans');
  });
});
