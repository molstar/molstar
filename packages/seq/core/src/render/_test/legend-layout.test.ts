/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Legend } from '@molstar/seqx-core/theme/index';
import { describe, expect, it } from 'vitest';
import type { PassWriter } from '../draw-lists.js';
import { LightTokens } from '../index.js';
import { layoutLegends } from '../legends.js';
import { createTextMetrics } from '../text.js';

// The export's legend band (_spec/README.md "Render", Export): a legend that wraps puts each text on its own line,
// title, entries and "+n" never over one another; an export lists 12 entries and never opens.

interface Box {
  s: string;
  x0: number;
  x1: number;
  y: number;
}

function lay(width: number, n = 12, more = 28) {
  const legend = {
    key: 'k',
    title: 'Chromosomes',
    channel: 'color',
    entries: {
      kind: 'swatches',
      items: Array.from({ length: n }, (_, i) => ({ label: `chr${i + 1}`, color: 0xff0000ff, value: `chr${i + 1}` })),
    },
    open: false,
    source: 'painted',
    more,
    notes: ['unlisted'],
    layers: ['/l'],
  } as unknown as Legend;
  const texts: Box[] = [];
  const out = {
    text: (s: string, x: number, y: number, _a: unknown, _f: string, _c: number, w: number) =>
      texts.push({ s, x0: x, x1: x + w, y }),
    fill() {},
    image() {},
  } as unknown as PassWriter;
  const metrics = createTextMetrics(() => (t: string) => 6 * t.length);
  const height = layoutLegends([legend], { x: 0, y: 0, width, height: 0 }, LightTokens, metrics, out);
  return { texts, height };
}

describe('the export legend band', () => {
  it('lays the title, the entries and "+n" out without their text boxes meeting, wrapping onto lines', () => {
    const { texts, height } = lay(300);
    expect(height).toBeGreaterThan(18);
    for (let i = 0; i < texts.length; i++)
      for (let j = i + 1; j < texts.length; j++) {
        const a = texts[i];
        const b = texts[j];
        if (a.y === b.y) expect(a.x1 <= b.x0 || b.x1 <= a.x0, `${a.s} / ${b.s}`).toBe(true);
      }
    expect(new Set(texts.map((t) => t.y)).size).toBeGreaterThan(1);
    expect(texts.at(-1)!.s).toBe('+28 more');
  });

  it('keeps its 12 entries and the count of the rest', () => {
    const { texts } = lay(300);
    expect(texts.filter((t) => t.s.startsWith('chr'))).toHaveLength(12);
  });
});
