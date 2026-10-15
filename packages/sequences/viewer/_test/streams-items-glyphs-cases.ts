/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The streamed glyphs' and placeholders' cases, for their heap guard and its quick twin (streams-items-glyphs.heap.ts,
// streams-items-glyphs.quick.ts): at 5 kb with the regions held, a pan's pass and composer frame (variants as glyphs
// region by region, stacked or as points, genes as items, letters as cells); and at 100 kb with every region pending, a
// pan's frame drawing the placeholder bands. Each call loads the viewer afresh (`vi.resetModules`).

import { vi } from 'vitest';

export const GLYPH_CASES = ['glyphs', 'points', 'pending'] as const;
export type GlyphCase = (typeof GLYPH_CASES)[number];

/** The case's frame: `glyphs` (5 kb, stacked), `points` (5 kb, unstacked) or `pending` (100 kb, every region pending). */
export async function glyphsCase(c: GlyphCase) {
  vi.resetModules();
  const { itemsFrame } = await import('./stream-items-fixtures.js');
  if (c === 'glyphs') return itemsFrame(100_040_000, 100_045_000);
  if (c === 'points') return itemsFrame(100_040_000, 100_045_000, false, 'points');
  return itemsFrame(100_000_000, 100_100_000, true);
}

/** The placeholder bands the frame drew. */
export async function pendingBands(g: Awaited<ReturnType<typeof glyphsCase>>): Promise<number> {
  const { SkinToken } = await import('../src/draw/backend.js');
  let bands = 0;
  for (let k = 0; k < g.backend.rectItems.n; k++) if (g.backend.rectItems.color[k] === SkinToken.pending) bands++;
  return bands;
}
