/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A merged block's copies each with their own secondary structure (3J3Q: 1,356 chains, 18,965 rows, 253 distinct):
// an opaque shape equal to one drawn at its start is skipped, so the frame draws what one copy of each variant draws,
// pixel for pixel; a translucent one is drawn every time (its copies darken each other).

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { describe, expect, it } from 'vitest';
import type { Theme } from '../src/themes.js';
import { drawn, layerOf, record } from './draw-fixtures.js';

/** Two variants of a chain's elements: the second's helix one residue longer. */
const VARIANTS = [
  { start: [2, 12, 24], stop: [9, 20, 30], category: ['strand', 'helix', 'turn'] },
  { start: [2, 12, 24], stop: [9, 21, 30], category: ['strand', 'helix', 'turn'] },
];

/**
 * `n` copies merged into one block, copy k with variant `variant(k)`; the secondary track's paths and rectangles. `alpha`:
 * its color theme's alpha over every feature (what a `rules` theme with an alpha gives).
 */
async function copies(n: number, variant: (k: number) => number, alpha?: number) {
  const b = createSeqViewSpec();
  const s = b.data.sequence({ ref: 's', value: 'MKVLATGHEK'.repeat(4) });
  const units = Array.from({ length: n }, (_, k) => b.data.unit({ ref: `u${k}`, sequence: s }));
  const value = { target: [] as string[], start: [] as number[], stop: [] as number[], category: [] as string[] };
  units.forEach((u, k) => {
    const x = VARIANTS[variant(k)];
    value.target.push(...x.start.map(() => u));
    value.start.push(...x.start);
    value.stop.push(...x.stop);
    value.category.push(...x.category);
  });
  const ss = b.data.annotation({ ref: 'ss', geometry: 'ranges', value });
  const v = b.view({ ref: 'v', label: 'V', height: 'box' });
  for (const u of units) v.block({ unit: u });
  v.track({ label: 'S' }).secondary({ data: ss });
  const t = await drawn(b.getDocument(), { width: 600, height: 200 });
  if (alpha !== undefined) {
    const at = t.ctx.graph.assembled.value?.view?.prepared[layerOf(t.ctx, 'S').k] as { themes: (Theme | null)[] };
    const color = at.themes[0] as Theme;
    at.themes[0] = { ...color, alpha: () => alpha };
  }
  const { rec } = record(t.ctx, 'S');
  const p = rec.pathItems;
  const paths: string[] = [];
  for (let i = 0; i < p.n; i++) {
    const c1 = i + 1 < p.n ? p.cstart[i + 1] : p.nc;
    paths.push(`${p.color[i]}:${Array.from(p.coords.subarray(p.cstart[i], c1))}`);
  }
  const r = rec.rectItems;
  const rects = Array.from({ length: r.n }, (_, i) => `${r.color[i]}:${r.x[i]},${r.y[i]},${r.w[i]},${r.h[i]}`);
  return { paths, rects };
}

describe('secondary: copies alike drawn once', () => {
  it('200 copies alike draw what one does', async () => {
    const one = await copies(1, () => 0);
    const many = await copies(200, () => 0);
    // A strand, a helix and a turn.
    expect(one.paths.length).toBe(3);
    expect(many).toEqual(one);
  });

  it('copies of two variants draw each variant once, as two copies do', async () => {
    const two = await copies(2, (k) => k);
    const many = await copies(200, (k) => ((k * 7) % 3 === 0 ? 1 : 0));
    // The strand and the turn alike in both, the helices a residue apart (two).
    expect(two.paths.length).toBe(4);
    expect(new Set(many.paths)).toEqual(new Set(two.paths));
    expect(many.paths.length).toBe(two.paths.length);
  });

  it('translucent copies are each drawn', async () => {
    const many = await copies(20, () => 0, 128);
    expect(many.paths.length).toBe(20 * 3);
  });
});
