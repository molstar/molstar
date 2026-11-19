/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DataStore } from '@molstar/seqx-core/data';
import { attachSources, bundledTextMetrics, LightTokens, Renderer, recordDrawLists } from '@molstar/seqx-core/render';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import * as globins from '../globins/family.js';

// The Globins figure's labels (_spec/README.md "Render", Layer types, ruler; Text layout, Numbers): a position carries
// one label where a group's end number and the ruler's tick label fall together, and no two labels overlap.

async function drawn(width: number) {
  const model = resolveSpec(globins.figure());
  const view = resolveViews(model)[0];
  const store = new DataStore();
  await attachSources(model, store);
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

  it('a group’s end number is drawn once on the first line, in its own style', async () => {
    const tens = (await drawn(900)).filter((l) => l.startsWith('text "10" '));
    expect(tens).toHaveLength(1);
    expect(tens[0]).toContain('10px sans');
  });
});
