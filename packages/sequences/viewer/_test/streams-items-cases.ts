/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The streamed items' cases, for their heap guard and its quick twin (streams-items.heap.ts, streams-items.quick.ts): a
// pan's pass and composer frame once the regions are held (`itemsFrame`, stream-items-fixtures.ts) at 100 kb (genes as
// items region by region, variants as density), at 1 Mb (both as density from `count` bins) and at 200 Mb (the letters
// gate's hatch and hint). A step is a pan, its pass and the frame. Each call loads the viewer afresh (`vi.resetModules`).

import { vi } from 'vitest';

/** The views: label, start and stop on chr2. */
export const ITEMS_VIEWS = [
  ['100 kb', 100_000_000, 100_100_000],
  ['1 Mb', 100_000_000, 101_000_000],
  ['200 Mb', 20_000_000, 220_000_000],
] as const;

/** The genome composed at [start, stop) with its regions held: `frame` is the step. */
export async function itemsCase(start: number, stop: number) {
  vi.resetModules();
  const { itemsFrame } = await import('./stream-items-fixtures.js');
  return itemsFrame(start, stop);
}
